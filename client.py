import numpy as np
import openpyxl as xlsx
import os
import sys
import lib.Debug.GoalsARM as Goals

import src.goals_const as CONST

## TODO:
## Config [partially done; TODO: net migration, direct incidence, direct CLHIV input]
## MigrInputs [requires implementing Beers disaggregation here or C-side]
## PopSizeInputs
## DirectIncidenceInputs
## PartnershipInputs
## MixingMatrix
## ContactInputs
## EpiInputs
## HIVDiseaseInputs
## HIVFertilityInputs
## ARTAdultInputs
## MCInputs
## DirectCLHIV

# Return the contents of a range in an excel tab as a numpy array
# The numpy array has float64 data type and "C" ordering
def xlsx_load_range(tab, cell_first, cell_final):
    return np.array([[cell.value for cell in row] for row in tab[cell_first:cell_final]], dtype=np.float64, order="C")

# Return a dict() mapping configuration tags to their values
# (e.g., x["first.year"] = 1970)
def xlsx_load_config(tab_config):
    vals = [cell[0].value for cell in tab_config['B2':'B8']]
    keys = [cell[0].value for cell in tab_config['D2':'D8']]
    return dict(zip(keys,vals))

def xlsx_load_pasfrs(tab_fert):
    p = xlsx_load_range(tab_fert, 'B2', 'CD8')
    return 0.01 * p.transpose()

def xlsx_load_migr(tab_migr):
    migr   = xlsx_load_range(tab_migr, 'B2',  'CD3')
    dist_m = xlsx_load_range(tab_migr, 'B6',  'CD22')
    dist_f = xlsx_load_range(tab_migr, 'B25', 'CD41')
    return migr, dist_m, dist_f

def xlsx_load_inci(tab_inci):
    inci   = xlsx_load_range(tab_inci, 'B2',  'CD2')
    sirr   = xlsx_load_range(tab_inci, 'B5',  'CD5')
    airr_m = xlsx_load_range(tab_inci, 'B9',  'CD25')
    airr_f = xlsx_load_range(tab_inci, 'B27', 'CD43')
    rirr_m = xlsx_load_range(tab_inci, 'B47', 'CD53')
    rirr_f = xlsx_load_range(tab_inci, 'B55', 'CD61')
    return inci[0], sirr[0], airr_m, airr_f, rirr_m, rirr_f

def init_from_xlsx(xlsx_name):
    wb = xlsx.load_workbook(filename=xlsx_name, read_only=True)
    config = xlsx_load_config(wb[CONST.XLSX_TAB_CONFIG])

    first_year = config[CONST.CFG_FIRST_YEAR]
    final_year = config[CONST.CFG_FINAL_YEAR]

    model = Goals.Projection(first_year, final_year)
    model.initialize(config[CONST.CFG_UPD_NAME])

    if not config[CONST.CFG_USE_UPD_PASFRS]:
        pasfrs = xlsx_load_pasfrs(wb[CONST.XLSX_TAB_PASFRS])
        model.init_pasfrs_from_5yr(pasfrs)

    if not config[CONST.CFG_USE_UPD_MIGR]:
        migr_net, migr_dist_m, migr_dist_f = xlsx_load_migr(wb[CONST.XLSX_TAB_MIGR])
        Warning("Net migration input from Excel is not supported yet")

    if config[CONST.CFG_USE_DIRECT_INCI]:
        inci, sirr, airr_m, airr_f, rirr_m, rirr_f = xlsx_load_inci(wb[CONST.XLSX_TAB_INCI])
        Warning("Direct incidence input from Excel is not supported yet")

    wb.close()
    return model

def main(xlsx_name):
    model = init_from_xlsx(xlsx_name)

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "C:\\Proj\\Repositories\\GoalsARM\\tests\\example.xlsx"
    main(xlsx_name)
