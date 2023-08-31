import numpy as np
import os
import pandas as pd
import plotnine
import sys
import src.goals_const as CONST
from src.goals_model import Model
from src.goals_results import Results

def main(xlsx_name, data_path):
    """! Main program entry point
    @param xlsx_name Excel file with Goals ARM inputs
    @param data_path Path to write output CSV files
    """
    model = Model()
    model.init_from_xlsx(xlsx_name)
    model.project(2030)

    results = Results(model)

    bin = list(range(CONST.AGE_MIN, CONST.AGE_MAX + 1)) + [np.inf]
    pop_total = results.pop_total(bin)
    pop_risk  = results.pop_risk([0,15,25,50,np.inf])

    ## TODO: replace sex indices with names
    ## TODO: stacked proportions instead of lines? 
    # pop_risk = pop_risk.to_frame()
    # pop_risk.reset_index(inplace=True)
    # pop_risk[CONST.STR_POP].astype("category")
    # pop_risk[CONST.STR_POP].rename_categories([])

    p = (plotnine.ggplot(pop_risk)
         + plotnine.aes(x=CONST.STR_YEAR, y=CONST.STR_VALUE, fill=CONST.STR_POP, color=CONST.STR_POP)
         + plotnine.geom_col()
         + plotnine.facet_grid('Age~Sex', scales="free_y")
         + plotnine.theme_bw())
    p.save("temp.tiff", dpi=600, units="in", width=6.5, height=4*2.25+0.25, pil_kwargs={"compression" : "tiff_lzw"})

    ## TODO: visualization stuff should be in a different branch (changes to goals_results.py and addition of visualize.py)
    ## TODO: use cProfile to check how expensive Results operations are

    pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    if len(sys.argv) == 1:
        xlsx_name = "inputs\\example-inputs.xlsx"
        data_path = "."
        main(xlsx_name, data_path)
    elif len(sys.argv) < 3:
        sys.stderr.write("USAGE: %s <input_param>.xlsx <output_path>" % (sys.argv[0]))
    else:
        xlsx_name = sys.argv[1]
        data_path = sys.argv[2]
        main(xlsx_name, data_path)
