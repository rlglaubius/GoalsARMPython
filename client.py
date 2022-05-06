import numpy as np
import openpyxl as xlsx
import os
import sys
import lib.Debug.GoalsARM as Goals

import src.goals_const as CONST

def xlsx_load_config(cfgtab):
    vals = [cell[0].value for cell in cfgtab['B2':'B8']]
    keys = [cell[0].value for cell in cfgtab['D2':'D8']]
    return dict(zip(keys,vals))

def init_from_xlsx(xlsx_name):
    wb = xlsx.load_workbook(filename=xlsx_name, read_only=True)
    config = xlsx_load_config(wb[CONST.XLSX_TAB_CONFIG])

    first_year = config[CONST.CFG_FIRST_YEAR]
    final_year = config[CONST.CFG_FINAL_YEAR]
    model = Goals.Projection(first_year, final_year)

    wb.close()
    return model

def main(xlsx_name):
    model = init_from_xlsx(xlsx_name)

if __name__ == "__main__":
    sys.stderr.write("Process %d" % (os.getpid()))
    xlsx_name = "C:\\Proj\\Repositories\\GoalsARM\\tests\\example.xlsx"
    main(xlsx_name)
