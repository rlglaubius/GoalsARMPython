import numpy as np
import openpyxl as xlsx
import src.goals_const as CONST
import src.goals_utils as Utils
import lib.Debug.GoalsARM as Goals

## TODO:
## Truncate inputs from XLSX to the requested years before passing them to GoalsARMCore. Use x[t0:t1,:] syntax
## [x] Config
## [x] MigrInputs
## [x] PopSizeInputs
## [x] DirectIncidenceInputs
## [ ] PartnershipInputs
## [ ] MixingMatrix
## [ ] ContactInputs
## [x] EpiInputs
## [x] HIVDiseaseInputs
## [x] HIVFertilityInputs
## [x] ARTAdultInputs
## [x] MCInputs
## [x] DirectCLHIV

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
        pop_pars = Utils.xlsx_load_popsize(wb[CONST.XLSX_TAB_POPSIZE])
        epi_pars = Utils.xlsx_load_epi(wb[CONST.XLSX_TAB_EPI])

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

        self._proj = Goals.Projection(self.year_first, self.year_final)
        self._proj.initialize(cfg_opts[CONST.CFG_UPD_NAME])
        self._proj.setup_storage_population(self.pop_adult_neg, self.pop_adult_hiv, self.pop_child_neg, self.pop_child_hiv)
        self._proj.setup_storage_deaths(self.deaths_adult_neg, self.deaths_adult_hiv, self.deaths_child_neg, self.deaths_child_hiv)

        self._initialize_population_sizes(pop_pars)

        if not cfg_opts[CONST.CFG_USE_UPD_PASFRS]:
            pasfrs = Utils.xlsx_load_pasfrs(wb[CONST.XLSX_TAB_PASFRS])
            self._proj.init_pasfrs_from_5yr(pasfrs)

        if not cfg_opts[CONST.CFG_USE_UPD_MIGR]:
            migr_net, migr_dist_m, migr_dist_f = Utils.xlsx_load_migr(wb[CONST.XLSX_TAB_MIGR])
            self._proj.init_migr_from_5yr(migr_net, migr_dist_f, migr_dist_m)

        if cfg_opts[CONST.CFG_USE_DIRECT_INCI]:
            inci, sirr, airr_m, airr_f, rirr_m, rirr_f = Utils.xlsx_load_inci(wb[CONST.XLSX_TAB_INCI])
            self._proj.use_direct_incidence(True)
            self._proj.init_direct_incidence(inci, sirr, airr_f, airr_m, rirr_f, rirr_m)
        else:
            self._proj.use_direct_incidence(False)
            self._proj.init_epidemic_seed(epi_pars[CONST.EPI_INITIAL_YEAR], epi_pars[CONST.EPI_INITIAL_PREV])
            self._proj.init_transmission(
                epi_pars[CONST.EPI_TRANSMIT_F2M],
                epi_pars[CONST.EPI_TRANSMIT_M2F],
                epi_pars[CONST.EPI_TRANSMIT_M2M],
                epi_pars[CONST.EPI_TRANSMIT_PRIMARY],
                epi_pars[CONST.EPI_TRANSMIT_CHRONIC],
                epi_pars[CONST.EPI_TRANSMIT_SYMPTOM],
                epi_pars[CONST.EPI_TRANSMIT_ART_VS],
                epi_pars[CONST.EPI_TRANSMIT_ART_VF])

        if cfg_opts[CONST.CFG_USE_DIRECT_CLHIV]:
            direct_clhiv = Utils.xlsx_load_direct_clhiv(wb[CONST.XLSX_TAB_DIRECT_CLHIV])
            self._proj.init_clhiv_agein(direct_clhiv)

        frr_age_no_art, frr_cd4_no_art, frr_age_on_art = Utils.xlsx_load_hiv_fert(wb[CONST.XLSX_TAB_HIV_FERT])
        dist, prog, mort, art1, art2, art3 = Utils.xlsx_load_adult_prog(wb[CONST.XLSX_TAB_ADULT_PROG])
        art_elig, art_num, art_pct, art_drop, art_mrr, art_vs = Utils.xlsx_load_adult_art(wb[CONST.XLSX_TAB_ADULT_ART])
        uptake_mc = Utils.xlsx_load_mc_uptake(wb[CONST.XLSX_TAB_MALE_CIRC])

        self._proj.init_hiv_fertility(frr_age_no_art, frr_cd4_no_art, frr_age_on_art)
        self._proj.init_adult_prog_from_10yr(dist, prog, mort)
        self._proj.init_adult_art_mort_from_10yr(art1, art2, art3, art_mrr)
        self._proj.init_adult_art_eligibility(art_elig)
        self._proj.init_adult_art_curr(art_num, art_pct)
        self._proj.init_adult_art_allocation(epi_pars[CONST.EPI_ART_MORT_WEIGHT])
        self._proj.init_adult_art_dropout(art_drop)
        self._proj.init_adult_art_suppressed(art_vs)

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
        
    def _initialize_population_sizes(self, pop_pars):
        """! Convenience function for initializing model population sizes
        """
        FEMALE, MALE = 1, 0
        self._proj.init_median_age_debut(pop_pars[CONST.POP_FIRST_SEX  ][FEMALE], pop_pars[CONST.POP_FIRST_SEX  ][MALE])
        self._proj.init_median_age_union(pop_pars[CONST.POP_FIRST_UNION][FEMALE], pop_pars[CONST.POP_FIRST_UNION][MALE])
        self._proj.init_mean_duration_union(pop_pars[CONST.POP_FIRST_UNION][MALE])
        self._proj.init_mean_duration_pwid(pop_pars[CONST.POP_DUR_PWID][FEMALE], pop_pars[CONST.POP_DUR_PWID][MALE])
        self._proj.init_mean_duration_fsw(pop_pars[CONST.POP_DUR_KEYPOP][FEMALE])
        self._proj.init_mean_duration_msm(pop_pars[CONST.POP_DUR_KEYPOP][MALE])
        self._proj.init_size_pwid(pop_pars[CONST.POP_SIZE_PWID][FEMALE], pop_pars[CONST.POP_SIZE_PWID][MALE])
        self._proj.init_size_fsw(pop_pars[CONST.POP_SIZE_KEYPOP][FEMALE])
        self._proj.init_size_msm(pop_pars[CONST.POP_SIZE_KEYPOP][MALE])
        self._proj.init_size_trans(pop_pars[CONST.POP_SIZE_TRANS][FEMALE], pop_pars[CONST.POP_SIZE_TRANS][MALE])
        self._proj.init_age_pwid(
            pop_pars[CONST.POP_PWID_LOC][FEMALE],
            pop_pars[CONST.POP_PWID_SHP][FEMALE],
            pop_pars[CONST.POP_PWID_LOC][MALE],
            pop_pars[CONST.POP_PWID_SHP][MALE])
        self._proj.init_age_fsw(pop_pars[CONST.POP_KEYPOP_LOC][FEMALE], pop_pars[CONST.POP_KEYPOP_SHP][FEMALE])
        self._proj.init_age_msm(pop_pars[CONST.POP_KEYPOP_LOC][MALE], pop_pars[CONST.POP_KEYPOP_SHP][MALE])
