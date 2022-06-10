import numpy as np
import openpyxl as xlsx
import os
import sys
import lib.Debug.GoalsARM as Goals

import src.goals_const as CONST

## TODO:
## [ ] Config [partially done; direct CLHIV input]
## [x] MigrInputs
## [x] PopSizeInputs
## [x] DirectIncidenceInputs
## [ ] PartnershipInputs
## [ ] MixingMatrix
## [ ] ContactInputs
## [ ] EpiInputs        <- next 3
## [x] HIVDiseaseInputs
## [ ] HIVFertilityInputs
## [ ] ARTAdultInputs   <- next 1
## [ ] MCInputs         <- next 2
## [ ] DirectCLHIV      <- next 4

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
    art_drop = xlsx_load_range(tab_art, 'B10', 'CD11')                   # Annual dropout, %
    art_mrr  = xlsx_load_range(tab_art, 'B13', 'CD14').transpose()       # Mortality time trend rate ratio over time
    art_vs   = xlsx_load_range(tab_art, 'B16', 'CD23')                   # Viral suppression on ART, %
    return art_elig, art_num, art_pct, art_drop, art_mrr, art_vs

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
    model.init_age_pwid(pop_pars[CONST.POP_PWID_LOC][FEMALE], pop_pars[CONST.POP_PWID_SHP][FEMALE],
                        pop_pars[CONST.POP_PWID_LOC][MALE  ], pop_pars[CONST.POP_PWID_SHP][MALE])
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

    first_year = cfg_opts[CONST.CFG_FIRST_YEAR]
    final_year = cfg_opts[CONST.CFG_FINAL_YEAR]

    model = Goals.Projection(first_year, final_year)
    model.initialize(cfg_opts[CONST.CFG_UPD_NAME])

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

    dist, prog, mort, art1, art2, art3 = xlsx_load_adult_prog(wb[CONST.XLSX_TAB_ADULT_PROG])
    art_elig, art_num, art_pct, art_drop, art_mrr, art_vs = xlsx_load_adult_art(wb[CONST.XLSX_TAB_ADULT_ART])

    model.init_adult_prog_from_10yr(dist, prog, mort)
    model.init_adult_art_mort_from_10yr(art1, art2, art3, art_mrr)
    model.init_adult_art_eligibility(art_elig)
    model.init_adult_art_curr(art_num, art_pct)
    model.init_adult_art_dropout(art_drop)

    ## TODO: TO USE
    ## art_vs

    wb.close()
    return model

def main(xlsx_name):
    """! Main program entry point"""
    model = init_from_xlsx(xlsx_name)

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "C:\\Proj\\Repositories\\GoalsARM\\tests\\example_unversioned.xlsx"
    main(xlsx_name)
