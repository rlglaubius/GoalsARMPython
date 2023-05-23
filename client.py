import numpy as np
import openpyxl as xlsx
import os
import sys
import lib.Debug.GoalsARM as Goals

import src.goals_const as CONST

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

# Return the contents of a range in an excel tab as a numpy array
# The numpy array has "C" ordering
def xlsx_load_range(tab, cell_first, cell_final, dtype=np.float64):
    return np.array([[cell.value for cell in row] for row in tab[cell_first:cell_final]], dtype=dtype, order="C")

def xlsx_load_config(tab_config):
    """! Load configuration options
    @param tab_config an openpyxl workbook tab
    @return a dict mapping configuration tags to parameter values (e.g., x["first.year"] = 1970)
    """
    vals = [cell[0].value for cell in tab_config['B2:B8']]
    keys = [cell[0].value for cell in tab_config['D2:D8']]
    return dict(zip(keys,vals))

def xlsx_load_popsize(tab_popsize):
    """! Load population size inputs
    @param tab_popsize an openpyxl workbook tab
    @return a dict mapping population size tags to parameter values. dict values are tuples: (male_value, female_value) 
    """
    vals = [tuple([cell.value for cell in row]) for row in tab_popsize['B2:C16']]
    keys = [cell[0].value for cell in tab_popsize['E2:E16']]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if keys != None} # Prune empty rows

def xlsx_load_epi(tab_epi):
    """! Load various epidemiological inputs
    @param tab_epi an openpyxl workbook tab
    @return a dict mapping parameter tags to values
    """
    vals = [cell[0].value for cell in tab_epi['B3:B19']]
    keys = [cell[0].value for cell in tab_epi['D3:D19']]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if keys != None} # Prune empty rows

def xlsx_load_pasfrs(tab_fert):
    p = xlsx_load_range(tab_fert, 'B2', 'CD8')
    return 0.01 * p.transpose()

def xlsx_load_migr(tab_migr):
    migr   = xlsx_load_range(tab_migr, 'B2',  'CD3')
    dist_m = xlsx_load_range(tab_migr, 'B6',  'CD22')
    dist_f = xlsx_load_range(tab_migr, 'B25', 'CD41')
    return migr.transpose(), 0.01 * dist_m.transpose(), 0.01 * dist_f.transpose()

def xlsx_load_inci(tab_inci):
    inci   = xlsx_load_range(tab_inci, 'B2',  'CD2')
    sirr   = xlsx_load_range(tab_inci, 'B5',  'CD5')
    airr_m = xlsx_load_range(tab_inci, 'B9',  'CD25').transpose()
    airr_f = xlsx_load_range(tab_inci, 'B27', 'CD43').transpose()
    rirr_m = xlsx_load_range(tab_inci, 'B47', 'CD53').transpose()
    rirr_f = xlsx_load_range(tab_inci, 'B55', 'CD61').transpose()
    return inci[0], sirr[0], airr_m, airr_f, rirr_m, rirr_f

def xlsx_load_hiv_fert(tab_frr):
    age_no_art = xlsx_load_range(tab_frr, 'B2',  'CD8')
    cd4_no_art = xlsx_load_range(tab_frr, 'B11', 'B17')[:,0] # X[:,0] converts to 1d array
    age_on_art = xlsx_load_range(tab_frr, 'B20', 'B26')[:,0]
    adj = tab_frr['B29'].value # local adjustment factor
    return age_no_art.transpose() * adj, cd4_no_art, age_on_art * adj

def xlsx_load_adult_prog(tab_prog):
    dist = xlsx_load_range(tab_prog, 'B4',  'I9')  # CD4 distribution after primary infection
    prog = xlsx_load_range(tab_prog, 'B14', 'I19') # HIV progression rates with untreated HIV
    mort = xlsx_load_range(tab_prog, 'B24', 'I30') # HIV-related mortality rates with untreated HIV
    art1 = xlsx_load_range(tab_prog, 'B35', 'I41') # HIV-related mortality rates if on ART for [0,6) months
    art2 = xlsx_load_range(tab_prog, 'B46', 'I52') # HIV-related mortality rates if on ART for [6,12) months
    art3 = xlsx_load_range(tab_prog, 'B57', 'I63') # HIV-related mortality rates if on ART for [12,\infty) months
    return dist, prog, mort, art1, art2, art3

def xlsx_load_adult_art(tab_art):
    art_elig = xlsx_load_range(tab_art, 'B2',  'CD2', dtype=np.int32)[0] # CD4-based eligibility threshold
    art_num  = xlsx_load_range(tab_art, 'B4',  'CD5').transpose()        # PLHIV on ART, #
    art_pct  = xlsx_load_range(tab_art, 'B7',  'CD8').transpose()        # PLHIV on ART, %
    art_drop = xlsx_load_range(tab_art, 'B10', 'CD11').transpose()       # Annual dropout, %
    art_mrr  = xlsx_load_range(tab_art, 'B13', 'CD14').transpose()       # Mortality time trend rate ratio over time
    art_vs   = xlsx_load_range(tab_art, 'B16', 'CD23').transpose()       # Viral suppression on ART, %
    return art_elig, art_num, art_pct, art_drop, art_mrr, art_vs

def xlsx_load_mc_uptake(tab_mc):
    return xlsx_load_range(tab_mc, 'B3', 'CD19').transpose()

def xlsx_load_direct_clhiv(tab_clhiv):
    return xlsx_load_range(tab_clhiv, 'D3', 'CF86').transpose()

def initialize_population_sizes(model, pop_pars):
    FEMALE, MALE = 1, 0
    model.init_median_age_debut(pop_pars[CONST.POP_FIRST_SEX  ][FEMALE], pop_pars[CONST.POP_FIRST_SEX  ][MALE])
    model.init_median_age_union(pop_pars[CONST.POP_FIRST_UNION][FEMALE], pop_pars[CONST.POP_FIRST_UNION][MALE])
    model.init_mean_duration_union(pop_pars[CONST.POP_FIRST_UNION][MALE])
    model.init_mean_duration_pwid(pop_pars[CONST.POP_DUR_PWID][FEMALE], pop_pars[CONST.POP_DUR_PWID][MALE])
    model.init_mean_duration_fsw(pop_pars[CONST.POP_DUR_KEYPOP][FEMALE])
    model.init_mean_duration_msm(pop_pars[CONST.POP_DUR_KEYPOP][MALE])
    model.init_size_pwid(pop_pars[CONST.POP_SIZE_PWID][FEMALE], pop_pars[CONST.POP_SIZE_PWID][MALE])
    model.init_size_fsw(pop_pars[CONST.POP_SIZE_KEYPOP][FEMALE])
    model.init_size_msm(pop_pars[CONST.POP_SIZE_KEYPOP][MALE])
    model.init_size_trans(pop_pars[CONST.POP_SIZE_TRANS][FEMALE], pop_pars[CONST.POP_SIZE_TRANS][MALE])
    model.init_age_pwid(
        pop_pars[CONST.POP_PWID_LOC][FEMALE],
        pop_pars[CONST.POP_PWID_SHP][FEMALE],
        pop_pars[CONST.POP_PWID_LOC][MALE],
        pop_pars[CONST.POP_PWID_SHP][MALE])
    model.init_age_fsw(pop_pars[CONST.POP_KEYPOP_LOC][FEMALE], pop_pars[CONST.POP_KEYPOP_SHP][FEMALE])
    model.init_age_msm(pop_pars[CONST.POP_KEYPOP_LOC][MALE], pop_pars[CONST.POP_KEYPOP_SHP][MALE])

def init_from_xlsx(xlsx_name):
    """! Create and initialize a Goals ARM model instance from inputs stored in Excel
    @param xlsx_name An Excel workbook with Goals ARM inputs
    @return An initialized Goals ARM model instance
    """

    wb = xlsx.load_workbook(filename=xlsx_name, read_only=True)
    cfg_opts = xlsx_load_config(wb[CONST.XLSX_TAB_CONFIG])
    pop_pars = xlsx_load_popsize(wb[CONST.XLSX_TAB_POPSIZE])
    epi_pars = xlsx_load_epi(wb[CONST.XLSX_TAB_EPI])

    first_year = cfg_opts[CONST.CFG_FIRST_YEAR]
    final_year = cfg_opts[CONST.CFG_FINAL_YEAR]

    num_years = final_year - first_year + 1
    pop_adult_neg = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP), dtype=np.float64, order="C")
    pop_adult_hiv = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP, CONST.N_HIV_ADULT, CONST.N_DTX), dtype=np.float64, order="C")
    pop_child_neg = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD), dtype=np.float64, order="C")
    pop_child_hiv = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD, CONST.N_HIV_CHILD, CONST.N_DTX), dtype=np.float64, order="C")

    deaths_adult_neg = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP), dtype=np.float64, order="C")
    deaths_adult_hiv = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_ADULT, CONST.N_POP, CONST.N_HIV_ADULT, CONST.N_DTX), dtype=np.float64, order="C")
    deaths_child_neg = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD), dtype=np.float64, order="C")
    deaths_child_hiv = np.zeros((num_years, CONST.N_SEX_MC, CONST.N_AGE_CHILD, CONST.N_HIV_CHILD, CONST.N_DTX), dtype=np.float64, order="C")

    model = Goals.Projection(first_year, final_year)
    model.initialize(cfg_opts[CONST.CFG_UPD_NAME])
    model.setup_storage_population(pop_adult_neg, pop_adult_hiv, pop_child_neg, pop_child_hiv)
    model.setup_storage_deaths(deaths_adult_neg, deaths_adult_hiv, deaths_child_neg, deaths_child_hiv)

    initialize_population_sizes(model, pop_pars)

    if not cfg_opts[CONST.CFG_USE_UPD_PASFRS]:
        pasfrs = xlsx_load_pasfrs(wb[CONST.XLSX_TAB_PASFRS])
        model.init_pasfrs_from_5yr(pasfrs)

    if not cfg_opts[CONST.CFG_USE_UPD_MIGR]:
        migr_net, migr_dist_m, migr_dist_f = xlsx_load_migr(wb[CONST.XLSX_TAB_MIGR])
        model.init_migr_from_5yr(migr_net, migr_dist_f, migr_dist_m)

    if cfg_opts[CONST.CFG_USE_DIRECT_INCI]:
        inci, sirr, airr_m, airr_f, rirr_m, rirr_f = xlsx_load_inci(wb[CONST.XLSX_TAB_INCI])
        model.use_direct_incidence(True)
        model.init_direct_incidence(inci, sirr, airr_f, airr_m, rirr_f, rirr_m)
    else:
        model.use_direct_incidence(False)
        model.init_epidemic_seed(epi_pars[CONST.EPI_INITIAL_YEAR], epi_pars[CONST.EPI_INITIAL_PREV])
        model.init_transmission(
            epi_pars[CONST.EPI_TRANSMIT_F2M],
            epi_pars[CONST.EPI_TRANSMIT_M2F],
            epi_pars[CONST.EPI_TRANSMIT_M2M],
            epi_pars[CONST.EPI_TRANSMIT_PRIMARY],
            epi_pars[CONST.EPI_TRANSMIT_CHRONIC],
            epi_pars[CONST.EPI_TRANSMIT_SYMPTOM],
            epi_pars[CONST.EPI_TRANSMIT_ART_VS],
            epi_pars[CONST.EPI_TRANSMIT_ART_VF])

    if cfg_opts[CONST.CFG_USE_DIRECT_CLHIV]:
        direct_clhiv = xlsx_load_direct_clhiv(wb[CONST.XLSX_TAB_DIRECT_CLHIV])
        model.init_clhiv_agein(direct_clhiv)

    frr_age_no_art, frr_cd4_no_art, frr_age_on_art = xlsx_load_hiv_fert(wb[CONST.XLSX_TAB_HIV_FERT])
    dist, prog, mort, art1, art2, art3 = xlsx_load_adult_prog(wb[CONST.XLSX_TAB_ADULT_PROG])
    art_elig, art_num, art_pct, art_drop, art_mrr, art_vs = xlsx_load_adult_art(wb[CONST.XLSX_TAB_ADULT_ART])
    uptake_mc = xlsx_load_mc_uptake(wb[CONST.XLSX_TAB_MALE_CIRC])

    model.init_hiv_fertility(frr_age_no_art, frr_cd4_no_art, frr_age_on_art)
    model.init_adult_prog_from_10yr(dist, prog, mort)
    model.init_adult_art_mort_from_10yr(art1, art2, art3, art_mrr)
    model.init_adult_art_eligibility(art_elig)
    model.init_adult_art_curr(art_num, art_pct)
    model.init_adult_art_allocation(epi_pars[CONST.EPI_ART_MORT_WEIGHT])
    model.init_adult_art_dropout(art_drop)
    model.init_adult_art_suppressed(art_vs)

    model.init_male_circumcision_uptake(uptake_mc)

    wb.close()
    return model

def main(xlsx_name):
    """! Main program entry point"""
    model = init_from_xlsx(xlsx_name)
    model.project(1980)
    pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "inputs\\example_unversioned_nohiv.xlsx"
    main(xlsx_name)
