import os
import sys
from src.goals_model import Model
from src.goals_results import Results

def main(xlsx_name):
    """! Main program entry point"""
    model = Model()
    model.init_from_xlsx(xlsx_name)
    model.project(1980)

    results = Results(model)
    pop = results.bigpop()
    
    pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    xlsx_name = "inputs\\example-unversioned-nohiv.xlsx"
    main(xlsx_name)
