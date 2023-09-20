import os
import pandas as pd
import sys
import time
from src.goals_model import Model

## Convert a numpy ndarray to a long data frame
## @param array a numpy ndarray
## @param names a list of names, one per dimension of ndarray
## @return a long data frame with one column per dimension of ndarray
def array2frame(array, names):
    array_index = pd.MultiIndex.from_product([range(s) for s in array.shape], names=names)
    array_frame = pd.DataFrame({'Value' : array.flatten()}, index=array_index)['Value']
    return array_frame

def main(xlsx_name, data_path):
    """! Main program entry point
    @param xlsx_name Excel file with Goals ARM inputs
    @param data_path Path to write output CSV files
    """
    t0 = time.time()
    model = Model()
    t1 = time.time()
    model.init_from_xlsx(xlsx_name)
    t2 = time.time()
    model.project(2030)
    t3 = time.time()

    birth_frame = array2frame(model.births, ['Year', 'Sex'])
    pop_child_neg = array2frame(model.pop_child_neg, ['Year', 'Sex', 'Age'])
    pop_child_hiv = array2frame(model.pop_child_hiv, ['Year', 'Sex', 'Age', 'CD4', 'ART'])
    pop_adult_neg = array2frame(model.pop_adult_neg, ['Year', 'Sex', 'Age', 'Risk'])
    pop_adult_hiv = array2frame(model.pop_adult_hiv, ['Year', 'Sex', 'Age', 'Risk', 'CD4', 'ART'])
    new_hiv = array2frame(model.new_infections, ['Year', 'Sex', 'Age', 'Risk'])
    t4 = time.time()

    birth_frame.to_csv(data_path + "/births.csv")
    pop_child_neg.to_csv(data_path + "/child-neg.csv")
    pop_child_hiv.to_csv(data_path + "/child-hiv.csv")
    pop_adult_neg.to_csv(data_path + "/adult-neg.csv")
    pop_adult_hiv.to_csv(data_path + "/adult-hiv.csv")
    new_hiv.to_csv(data_path + "/new-hiv.csv")
    t5 = time.time()

    sys.stdout.write("Construct\t%0.2fs\nInitialize\t%0.2fs\nProject\t\t%0.2fs\nAnalysis\t%0.2fs\nCSV write\t%0.2fs\n" % (t1-t0, t2-t1, t3-t2, t4-t3, t5-t4))

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
