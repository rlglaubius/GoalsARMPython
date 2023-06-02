import os
import pandas as pd
import sys
import src.goals_const as CONST
import time
from src.goals_model import Model
from src.goals_results import Results

def main(xlsx_name):
    """! Main program entry point"""
    t0 = time.time()
    model = Model()
    t1 = time.time()
    model.init_from_xlsx(xlsx_name)
    t2 = time.time()
    model.project(2030)
    t3 = time.time()

    sys.stdout.write("Construct\t%0.2fs\nInitialize\t%0.2fs\nProject\t%0.2fs\n" % (t1-t0, t2-t1, t3-t2))

    results = Results(model)
    pop = results.bigpop()
    
    pop_names = ['Year', 'Sex', 'Age']
    pop_index = pd.MultiIndex.from_product([range(s) for s in pop.shape], names=pop_names)
    pop_frame = pd.DataFrame({'Value' : pop.flatten()}, index=pop_index)['Value']
    pop_frame.to_csv("bigpop.csv")

    pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "inputs\\example-unversioned-nohiv-dbg.xlsx"
    main(xlsx_name)
