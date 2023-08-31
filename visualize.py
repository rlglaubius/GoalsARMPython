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

    p = (plotnine.ggplot(pop_risk)
         + plotnine.aes(x=CONST.STR_YEAR, y=CONST.STR_VALUE, fill=CONST.STR_POP, color=CONST.STR_POP)
         + plotnine.geom_col()
         + plotnine.scale_fill_brewer(type="qual", palette="Dark2")
         + plotnine.scale_color_brewer(type="qual", palette="Dark2")
         + plotnine.facet_grid('Age~Gender', scales="free_y")
         + plotnine.theme_bw())
    p.save("temp.tiff", dpi=600, units="in", width=6.5, height=4*2.25+0.25, pil_kwargs={"compression" : "tiff_lzw"})

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
