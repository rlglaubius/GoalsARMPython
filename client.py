import numpy as np
import openpyxl as xlsx
import os
import sys
import lib.Debug.GoalsARM as Goals

import src.goals_const as CONST

# Return a dict() mapping configuration tags to their values
# (e.g., x["first.year"] = 1970)
def xlsx_load_config(tab_config):
    vals = [cell[0].value for cell in tab_config['B2':'B8']]
    keys = [cell[0].value for cell in tab_config['D2':'D8']]
    return dict(zip(keys,vals))

def xlsx_load_pasfrs(tab_fert):
    p = np.array([[cell.value for cell in row] for row in tab_fert['B2':'CD8']], dtype=np.float64, order="C")
    return 0.01 * p.transpose()

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

    wb.close()
    return model

def main(xlsx_name):
    model = init_from_xlsx(xlsx_name)

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "C:\\Proj\\Repositories\\GoalsARM\\tests\\example.xlsx"
    main(xlsx_name)
