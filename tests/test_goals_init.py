import numpy as np
import unittest
import src.goals_const as CONST
import src.goals_proj.x64.Release.goals_proj as GoalsProj

class Test_TestGoalsInit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.dtype = np.float64
        self.order = "C"

        self.year_first = 1970
        self.year_final = 1971
        self.num_years = self.year_final - self.year_first + 1
        self.proj = GoalsProj.Projection(self.year_first, self.year_final)
        
    def test_init_pasfrs_from_5yr(self):
        pasfrs5y = 0.01 * np.array([[12.02, 21.08, 20.45, 17.70, 14.93,  9.99, 3.83],
                                    [12.05, 21.10, 20.38, 17.65, 14.96, 10.03, 3.83]], dtype=self.dtype, order=self.order)
        try:
            self.proj.init_pasfrs_from_5yr(pasfrs5y)
        except RuntimeError:
            self.fail("Unexpected runtime error during init_pasfrs_from_5yr(births)")

if __name__ == "__main__":
    unittest.main()
