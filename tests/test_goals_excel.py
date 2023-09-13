import numpy as np
import unittest
import src.goals_const as CONST
from src.goals_model import Model

## Unit tests for Goals input initialization from Excel

class Test_TestGoalsExcel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.goals = Model()
        
    def test_external_clhiv(self):
        self.goals.init_from_xlsx("tests/test-external-clhiv.xlsx")
        self.goals.project(2049)

        ## Reference values
        ref_by_year = np.array([   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                                   0,    0,    0,    0,    0,    0,    0,    0,    8,   16,
                                  31,   50,   76,  118,  176,  268,  381,  549,  762, 1052,
                                1417, 1863, 2392, 2994, 3644, 4328, 4952, 5518, 6003, 6362,
                                6645, 6840, 6958, 7060, 7170, 7199, 7132, 7050, 6930, 6811,
                                6650, 6385, 5806, 5528, 5202, 4348, 4170, 4008, 3655, 3131,
                                2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583,
                                2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583, 2583])
        ref_by_sex = np.array([107550, 105748])
        ref_by_age = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 213298])
        ref_by_cd4 = np.array([19761, 21823, 36314, 46410, 44504, 44486])
        ref_by_art = np.array([0, 0, 94155, 9741, 0, 109402])

        ## Output values
        out_by_year = self.goals.pop_child_hiv.sum((1,2,3,4))
        out_by_sex = [self.goals.pop_child_hiv[:,0, :,:].sum(),
                      self.goals.pop_child_hiv[:,1:,:,:].sum()] # aggregate circumcised and uncircumcised males
        out_by_age = self.goals.pop_child_hiv.sum((0,1,3,4))
        out_by_cd4 = self.goals.pop_child_hiv.sum((0,1,2,4))
        out_by_art = self.goals.pop_child_hiv.sum((0,1,2,3))

        [self.assertAlmostEqual(x,y) for x, y in zip(ref_by_year, out_by_year)]
        [self.assertAlmostEqual(x,y) for x, y in zip(ref_by_sex,  out_by_sex )]
        [self.assertAlmostEqual(x,y) for x, y in zip(ref_by_age,  out_by_age )]
        [self.assertAlmostEqual(x,y) for x, y in zip(ref_by_cd4,  out_by_cd4 )]
        [self.assertAlmostEqual(x,y) for x, y in zip(ref_by_art,  out_by_art )]

if __name__ == "__main__":
    unittest.main()
