import numpy as np
import unittest
import src.goals_const as CONST
import src.goals_proj.x64.Release.goals_proj as GoalsProj

class Test_TestMemorySharing(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.dtype = np.float64
        self.order = "C"

        self.year_first = 1970
        self.year_final = 2030
        self.num_years = self.year_final - self.year_first + 1
        self.proj = GoalsProj.Projection(self.year_first, self.year_final)
    
    
    def test_births(self):
        births = np.zeros((self.num_years, CONST.N_SEX), dtype=self.dtype, order=self.order)
        try:
            self.proj.share_output_births(births)
        except RuntimeError:
            self.fail("Unexpected runtime error during share_output_births(births)")

    def test_births_ndim(self):
        births = np.zeros((self.num_years), dtype=self.dtype, order=self.order)
        self.assertRaises(RuntimeError, self.proj.share_output_births, births)

    def test_births_shape(self):
        births = np.zeros((CONST.N_SEX, self.num_years), dtype=self.dtype, order=self.order)
        self.assertRaises(RuntimeError, self.proj.share_output_births, births)

    def test_births_layout(self):
        births = np.zeros((self.num_years, CONST.N_SEX), dtype=self.dtype, order="F")
        self.assertRaises(RuntimeError, self.proj.share_output_births, births)
    
if __name__ == "__main__":
    unittest.main()
