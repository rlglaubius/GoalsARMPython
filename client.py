import numpy as np
import os
import sys
import lib.Debug.GoalsARM as Goals

def main():
    first_year = 1970
    final_year = 2030
    model = Goals.Projection(first_year, final_year)

if __name__ == "__main__":
    sys.stderr.write("Process %d" % (os.getpid()))
    main()
    