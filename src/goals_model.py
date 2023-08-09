import math
import numpy as np
import scipy as sp
import openpyxl as xlsx
import src.goals_const as CONST
import src.goals_utils as Utils
import lib.Release.GoalsARM as Goals # At last check, calculation was ~100-fold slower with the Debug version of the library

## TODO:
## Create an Excel reader that just loads the raw inputs from Excel into member
## variables. Then the Model here would be responsible for doing any
## transformations on those variables before passing them to the calculation
## engine. The C++ transfer layer and calculation engine ideally should not do
## any input transformations.

class Model:
    """! Goals model class. This is wraps an external Goals ARM core projection object
    so that calling applications should not need to care about the Python-C++ API
    """

    def __init__(self):
        self._dtype = np.float64
        self._order = "C"
        self._initialized = False # True if projection inputs have been initialized, False otherwise
        self._projected   = -1    # The latest year that the projection has been calculated through (-1 if not done)
    
    def is_initialized(self):
        """! Check if the projection has been initialized"""
        return self._initialized
    
    def last_valid_year(self):
        """! Return the latest year for which the projection has been calculated, or -1 if uncalculated"""
        return self._projected
    
    def init_from_xlsx(self, xlsx_name):
        """! Create and initialize a Goals ARM model instance from inputs stored in Excel
        @param xlsx_name An Excel workbook with Goals ARM inputs
        @return An initialized Goals ARM model instance
        """

        wb = xlsx.load_workbook(filename=xlsx_name, read_only=True)
        cfg_opts = Utils.xlsx_load_config(wb[CONST.XLSX_TAB_CONFIG])
        self.epi_pars = Utils.xlsx_load_epi(wb[CONST.XLSX_TAB_EPI])

        # Conver % epi parameters to proportions
        self.epi_pars[CONST.EPI_INITIAL_PREV   ] *= 0.01
        self.epi_pars[CONST.EPI_TRANSMIT_F2M   ] *= 0.01
        self.epi_pars[CONST.EPI_EFFECT_VMMC    ] *= 0.01
        self.epi_pars[CONST.EPI_EFFECT_CONDOM  ] *= 0.01
        self.epi_pars[CONST.EPI_ART_MORT_WEIGHT] *= 0.01

        self.year_first = cfg_opts[CONST.CFG_FIRST_YEAR]
        self.year_final = cfg_opts[CONST.CFG_FINAL_YEAR]

        num_years = self.year_final - self.year_first + 1
        shp_adult_neg = (num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP)
        shp_adult_hiv = (num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP, CONST.N_HIV_ADULT, CONST.N_DTX)
        shp_child_neg = (num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD)
        shp_child_hiv = (num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD, CONST.N_HIV_CHILD, CONST.N_DTX)

        self.pop_adult_neg = np.zeros(shp_adult_neg, dtype=self._dtype, order=self._order)
        self.pop_adult_hiv = np.zeros(shp_adult_hiv, dtype=self._dtype, order=self._order)
        self.pop_child_neg = np.zeros(shp_child_neg, dtype=self._dtype, order=self._order)
        self.pop_child_hiv = np.zeros(shp_child_hiv, dtype=self._dtype, order=self._order)

        self.deaths_adult_neg = np.zeros(shp_adult_neg, dtype=self._dtype, order=self._order)
        self.deaths_adult_hiv = np.zeros(shp_adult_hiv, dtype=self._dtype, order=self._order)
        self.deaths_child_neg = np.zeros(shp_child_neg, dtype=self._dtype, order=self._order)
        self.deaths_child_hiv = np.zeros(shp_child_hiv, dtype=self._dtype, order=self._order)

        self.births = np.zeros((num_years, CONST.N_SEX), dtype=self._dtype, order=self._order)
        self.births_exposed = np.zeros((num_years), dtype=self._dtype, order=self._order)

        self.new_infections = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE, CONST.N_POP), dtype=self._dtype, order=self._order)

        self._proj = Goals.Projection(self.year_first, self.year_final)
        self._proj.initialize(cfg_opts[CONST.CFG_UPD_NAME])
        self._proj.share_output_population(self.pop_adult_neg, self.pop_adult_hiv, self.pop_child_neg, self.pop_child_hiv)
        self._proj.share_output_births(self.births)
        self._proj.share_output_deaths(self.deaths_adult_neg, self.deaths_adult_hiv, self.deaths_child_neg, self.deaths_child_hiv)
        self._proj.share_output_new_infections(self.new_infections)
        self._proj.share_output_births_exposed(self.births_exposed)

        med_age_debut, med_age_union, avg_dur_union, kp_size, kp_stay, kp_turnover = Utils.xlsx_load_popsize(wb[CONST.XLSX_TAB_POPSIZE])
        self._initialize_population_sizes(med_age_debut, med_age_union, avg_dur_union, kp_size, kp_stay, kp_turnover)

        if not cfg_opts[CONST.CFG_USE_UPD_PASFRS]:
            pasfrs = Utils.xlsx_load_pasfrs(wb[CONST.XLSX_TAB_PASFRS])
            self._proj.init_pasfrs_from_5yr(pasfrs)

        if not cfg_opts[CONST.CFG_USE_UPD_MIGR]:
            migr_net, migr_dist_m, migr_dist_f = Utils.xlsx_load_migr(wb[CONST.XLSX_TAB_MIGR])
            self._proj.init_migr_from_5yr(migr_net, migr_dist_f, migr_dist_m)

        self._proj.init_effect_vmmc(self.epi_pars[CONST.EPI_EFFECT_VMMC])
        self._proj.init_effect_condom(self.epi_pars[CONST.EPI_EFFECT_CONDOM])
        if cfg_opts[CONST.CFG_USE_DIRECT_INCI]:
            inci, sirr, airr_m, airr_f, rirr_m, rirr_f = Utils.xlsx_load_inci(wb[CONST.XLSX_TAB_INCI])
            self._proj.use_direct_incidence(True)
            self._proj.init_direct_incidence(0.01 * inci, sirr, airr_f, airr_m, rirr_f, rirr_m)
        else:
            self.partner_time_trend, self.partner_age_params, self.partner_pop_ratios = Utils.xlsx_load_partner_rates(wb[CONST.XLSX_TAB_PARTNER])
            age_prefs, pop_prefs, self.p_married = Utils.xlsx_load_partner_prefs(wb[CONST.XLSX_TAB_PARTNER])
            mix_raw = Utils.xlsx_load_mixing_levels(wb[CONST.XLSX_TAB_MIXNG_MATRIX])
            self.sex_acts, condom_freq, self.pwid_force, needle_sharing = Utils.xlsx_load_contact_params(wb[CONST.XLSX_TAB_CONTACT])
            self.partner_rate = self.calc_partner_rates(self.partner_time_trend, self.partner_age_params, self.partner_pop_ratios)
            self.age_mixing = self.calc_partner_prefs(age_prefs)
            self.pop_assort = self.calc_pop_assort(pop_prefs)
            self.mix_levels = self.calc_mix_levels(mix_raw)
            self.condom_freq = 0.01 * condom_freq
            self.needle_sharing = 0.01 * needle_sharing
            self.p_married = 0.01 * np.array([self.p_married[CONST.SEX_FEMALE, CONST.POP_PWID - CONST.POP_KEY_MIN],
                                              self.p_married[CONST.SEX_MALE,   CONST.POP_PWID - CONST.POP_KEY_MIN],
                                              self.p_married[CONST.SEX_FEMALE, CONST.POP_FSW  - CONST.POP_KEY_MIN],
                                              self.p_married[CONST.SEX_MALE,   CONST.POP_CSW  - CONST.POP_KEY_MIN],
                                              self.p_married[CONST.SEX_MALE,   CONST.POP_MSM  - CONST.POP_KEY_MIN],
                                              self.p_married[CONST.SEX_FEMALE, CONST.POP_TGW  - CONST.POP_KEY_MIN]])
            self._proj.share_input_partner_rate(self.partner_rate)
            self._proj.share_input_age_mixing(self.age_mixing)
            self._proj.share_input_pop_assort(self.pop_assort)
            self._proj.share_input_pwid_risk(self.pwid_force, self.needle_sharing)
            self._proj.use_direct_incidence(False)
            self._proj.init_epidemic_seed(self.epi_pars[CONST.EPI_INITIAL_YEAR] - self.year_first, self.epi_pars[CONST.EPI_INITIAL_PREV])
            self._proj.init_transmission(
                self.epi_pars[CONST.EPI_TRANSMIT_F2M],
                self.epi_pars[CONST.EPI_TRANSMIT_M2F],
                self.epi_pars[CONST.EPI_TRANSMIT_M2M],
                self.epi_pars[CONST.EPI_TRANSMIT_PRIMARY],
                self.epi_pars[CONST.EPI_TRANSMIT_CHRONIC],
                self.epi_pars[CONST.EPI_TRANSMIT_SYMPTOM],
                self.epi_pars[CONST.EPI_TRANSMIT_ART_VS],
                self.epi_pars[CONST.EPI_TRANSMIT_ART_VF])
            self._proj.init_keypop_married(self.p_married)
            self._proj.init_mixing_matrix(self.mix_levels)
            self._proj.init_sex_acts(self.sex_acts)
            self._proj.init_condom_freq(self.condom_freq)

        if cfg_opts[CONST.CFG_USE_DIRECT_CLHIV]:
            direct_clhiv = Utils.xlsx_load_direct_clhiv(wb[CONST.XLSX_TAB_DIRECT_CLHIV])
            self._proj.init_clhiv_agein(direct_clhiv)

        self.hiv_frr = Utils.xlsx_load_hiv_fert(wb[CONST.XLSX_TAB_HIV_FERT])
        dist, prog, mort, art1, art2, art3 = Utils.xlsx_load_adult_prog(wb[CONST.XLSX_TAB_ADULT_PROG])
        art_elig, art_num, art_pct, art_stop, art_mrr, art_vs = Utils.xlsx_load_adult_art(wb[CONST.XLSX_TAB_ADULT_ART])
        uptake_mc = Utils.xlsx_load_mc_uptake(wb[CONST.XLSX_TAB_MALE_CIRC])

        self.anc_par = Utils.xlsx_load_fit_pars(wb[CONST.XLSX_TAB_FITTING_PARS])

        self._proj.init_hiv_fertility(self.hiv_frr['age'] * self.hiv_frr['laf'], self.hiv_frr['cd4'], self.hiv_frr['art'] * self.hiv_frr['laf'])
        self._proj.init_adult_prog_from_10yr(0.01 * dist, prog, mort)
        self._proj.init_adult_art_mort_from_10yr(art1, art2, art3, art_mrr)
        self._proj.init_adult_art_eligibility(art_elig)
        self._proj.init_adult_art_curr(art_num, 0.01 * art_pct)
        self._proj.init_adult_art_allocation(self.epi_pars[CONST.EPI_ART_MORT_WEIGHT])
        self._proj.init_adult_art_interruption(-np.log(1.0 - 0.01 * art_stop)) # convert %/year to an event rate
        self._proj.init_adult_art_suppressed(0.01 * art_vs)

        self._proj.init_male_circumcision_uptake(uptake_mc)
        self._initialized = True

        wb.close()

    def project(self, year_stop):
        """! Calculate the projection from the first year to the requested final year. The
        projection must be initialized (e.g., via init_from_xlsx) and the year_final must
        not exceed 
        """
        self._proj.project(year_stop)
        self._projected = year_stop

    def invalidate(self, year):
        """! Invalidate projections from a given year onward. Call this after project(year_stop) if
        you need to recalculate indicators for years before year_stop, otherwise projection will
        resume from year_stop. """
        self._proj.invalidate(year)
        self._projected = min(year, self._projected)
        
    def _initialize_population_sizes(self, med_age_debut, med_age_union, avg_dur_union, kp_size, kp_stay, kp_turnover):
        """! Convenience function for initializing model population sizes
        """
        self._proj.init_median_age_debut(med_age_debut[CONST.SEX_FEMALE], med_age_debut[CONST.SEX_MALE])
        self._proj.init_median_age_union(med_age_union[CONST.SEX_FEMALE], med_age_union[CONST.SEX_MALE])
        self._proj.init_mean_duration_union(avg_dur_union)
        self._proj.init_keypop_size_params(0.01 * kp_size, kp_stay, kp_turnover)

    def calc_partner_rates(self, time_trend, age_params, pop_ratios):
        """! Calculate partnership rates by year, sex, age, and behavioral risk group
        @param time_trend lifetime partnership-years by sex and year
        @param age_params beta distribution mean and size parameters that specify partner rates by age
        @param pop_params rate ratios by behavioral risk group, excluding the sexually naive group
        """
        
        ## Calculate age-specific rate ratios from age_params
        age_ratios = np.zeros((CONST.N_SEX, CONST.N_AGE_ADULT), dtype=self._dtype, order=self._order)
        raw_ages = np.array(range(CONST.AGE_ADULT_MIN, CONST.AGE_ADULT_MAX + 1))
        std_ages = (raw_ages - CONST.AGE_ADULT_MIN) / (CONST.AGE_ADULT_MAX - CONST.AGE_ADULT_MIN)
        for s in [CONST.SEX_FEMALE, CONST.SEX_MALE]: # not idiomatic...
            raw_mean = age_params[0,s]
            std_mean = (raw_mean - CONST.AGE_ADULT_MIN) / (CONST.AGE_ADULT_MAX - CONST.AGE_ADULT_MIN)
            dist = sp.stats.beta(age_params[1,s] * std_mean, age_params[1,s] * (1.0 - std_mean))

            ## This intentionally excludes CONST.AGE_ADULT_MAX so that its age_ratio is 0
            age_ratios[s,0:(CONST.N_AGE_ADULT - 1)] = np.diff(dist.cdf(std_ages))

        ## Calculate the full partner rate matrix
        num_yrs = self.year_final - self.year_first + 1
        yr_bgn = self.year_first - CONST.XLSX_FIRST_YEAR
        yr_end = self.year_final - CONST.XLSX_FIRST_YEAR + 1
        partner_rate = np.zeros((num_yrs, CONST.N_SEX, CONST.N_AGE_ADULT, CONST.N_POP), dtype=self._dtype, order=self._order)

        ## Reorganize pop_ratios to include sexually inactive people and to map gender
        ## identity to assigned sex at birth
        pop_ratios_aug = np.zeros((CONST.N_POP, CONST.N_SEX), dtype=self._dtype, order=self._order)
        pop_ratios_aug[CONST.POP_NEVER:(CONST.POP_FSW+1), CONST.SEX_FEMALE] = pop_ratios[0:5, CONST.SEX_FEMALE]
        pop_ratios_aug[CONST.POP_NEVER:(CONST.POP_MSM+1), CONST.SEX_MALE  ] = pop_ratios[0:6, CONST.SEX_MALE  ] 
        pop_ratios_aug[CONST.POP_TGW, CONST.SEX_MALE] = pop_ratios[6, CONST.SEX_FEMALE]

        # ## Loop variant (readable)
        # for s in range(CONST.N_SEX):
        #     for a in range(CONST.N_AGE_ADULT):
        #         for r in range(CONST.N_POP):
        #             partner_rate[:,s,a,r] = time_trend[s,yr_bgn:yr_end] * age_ratios[s,a] * pop_ratios_aug[r,s]

        ## Vectorized variant (much faster)
        partner_rate = np.zeros((num_yrs, CONST.N_SEX, CONST.N_AGE_ADULT, CONST.N_POP), dtype=self._dtype, order=self._order)
        for s in range(CONST.N_SEX):
            partner_rate[:,s,:,:] = np.outer(time_trend[s,yr_bgn:yr_end], np.outer(age_ratios[s,:], pop_ratios_aug[:,s])).reshape((num_yrs, CONST.N_AGE_ADULT, CONST.N_POP))

        return partner_rate
    
    def calc_partner_prefs(self, age_prefs):
        ## age differences mean and variance
        oppo_diff_avg, oppo_diff_var = age_prefs[0], age_prefs[1] # male-female partnerships
        male_diff_var = age_prefs[2]                              # male-male partnerships
        mix = np.zeros((CONST.N_SEX, CONST.N_AGE_ADULT, CONST.N_SEX, CONST.N_AGE_ADULT), dtype=self._dtype, order=self._order)

	    ## Use Newton-Raphson to approximate shape and scale parameters of the Fisk
	    ## distribution given the specified mean and variance.
        num_iter = 5
        shift = -10 # assuming negligible partnerships with females who are 10 years older than their male partners
        m, v = oppo_diff_avg - shift, oppo_diff_var
        target = m * m / (m * m + v)
        x = 0.5 * math.pi - target
        for k in range(num_iter):
            cot_x = 1.0 / math.tan(x)
            csc_x = 1.0 / math.sin(x)
            fx = x * cot_x - target
            dx = cot_x - x * csc_x * csc_x
            x = x - fx / dx

        shape = math.pi / x
        scale = m * math.sin(x) / x
        oppo_dist = sp.stats.fisk(shape, shift, scale)
        same_dist = sp.stats.norm(0.0, math.sqrt(male_diff_var))
        
        ## This could be vectorized for speed. Beware Brian Kernighan's debugging insight.
        ## Calculate unnormalized mixing preferences 
        oppo_raw = np.zeros((CONST.N_AGE_ADULT, CONST.N_AGE_ADULT), dtype=self._dtype, order=self._order)
        same_raw = np.zeros((CONST.N_AGE_ADULT, CONST.N_AGE_ADULT), dtype=self._dtype, order=self._order)
        for a in range(CONST.AGE_ADULT_MIN, CONST.AGE_ADULT_MAX): # intentionally omits the 80+ age group
            b = a - CONST.AGE_ADULT_MIN
            min_age_diff = CONST.AGE_ADULT_MIN - a
            max_age_diff = CONST.AGE_ADULT_MAX - a
            oppo_raw[b,:-1] = np.diff(oppo_dist.cdf(range(min_age_diff, max_age_diff + 1)))
            same_raw[b,:-1] = np.diff(same_dist.cdf(range(min_age_diff, max_age_diff + 1)))

        ## Fill in normalized mixing matrix. This could be vectorized further
        for b in range(CONST.N_AGE_ADULT - 1):
            mix[CONST.SEX_FEMALE, b, CONST.SEX_MALE,   :] = oppo_raw[b,:] / oppo_raw[b,:].sum()
            mix[CONST.SEX_MALE,   b, CONST.SEX_FEMALE, :] = oppo_raw[:,b] / oppo_raw[:,b].sum()
            mix[CONST.SEX_MALE,   b, CONST.SEX_MALE,   :] = same_raw[b,:] / same_raw[b,:].sum()
        
        return mix
    
    def calc_pop_assort(self, pop_prefs):
        """"! Convert raw assortativity inputs from Excel into usable inputs by converting
        them from percentages to proportions, adding a row for POP_NOSEX, and mapping TGW
        from female identity to male sex at birth
        @param pop_prefs array of assortativity parameters by sex and population, excluding POP_NOSEX
        """
        assort = np.zeros((CONST.N_SEX, CONST.N_POP), dtype=self._dtype, order=self._order)
        assort[CONST.SEX_FEMALE, CONST.POP_NEVER:(CONST.POP_FSW+1)] = 0.01 * pop_prefs[0:5, CONST.SEX_FEMALE]
        assort[CONST.SEX_MALE, CONST.POP_NEVER:(CONST.POP_MSM+1)] = 0.01 * pop_prefs[0:6, CONST.SEX_MALE]
        assort[CONST.SEX_MALE, CONST.POP_TGW] = 0.01 * pop_prefs[6, CONST.SEX_FEMALE]

        return assort
    
    def calc_mix_levels(self, mix_levels):
        """! Reshape the raw mixing levels read Excel into a more usable layout"""

        ## 1. Augment with rows ignored by Excel (never had sex, female versions of male-only populations)
        mix = np.zeros((CONST.N_SEX * CONST.N_POP, CONST.N_SEX * CONST.N_POP), dtype=np.int32, order=self._order)
        mix[1:7, 1:7 ] = mix_levels[0:6, 0:6 ] # female-female block
        mix[1:7, 9:15] = mix_levels[0:6, 6:12] # female-male block
        mix[9:15,1:7 ] = mix_levels[6:12,0:6 ] # male-female block
        mix[9:15,9:15] = mix_levels[6:12,6:12] # male-male block

        ## Swap TGW from females (gender) to males (assigned sex at birth)
        mix[[15,6],:] = mix[[6,15],:] # swaps TGW rows into place
        mix[:,[15,6]] = mix[:,[6,15]] # swaps TGW cols into place

        return mix.reshape((CONST.N_SEX, CONST.N_POP, CONST.N_SEX, CONST.N_POP))
