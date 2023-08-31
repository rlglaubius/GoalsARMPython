import numpy as np
import pandas as pd
import src.goals_const as CONST
import src.goals_model as Goals

def array2frame(array, names):
    """! Convert a numpy ndarray to a long data frame
    @param array a numpy ndarray
    @param names a list of names, one per dimension of ndarray
    @return a long data frame with one column per dimension of ndarray
    """
    array_index = pd.MultiIndex.from_product([range(s) for s in array.shape], names=names)
    array_frame = pd.DataFrame({'Value' : array.flatten()}, index=array_index)['Value']
    return array_frame.to_frame()

class Results:
    def __init__(self, model):
        """! Initialize the object with a Goals model
        @param model The Goals model
        """
        self._model = model
        self._dtype = model._dtype
        self._order = model._order
        self._years = self._model.year_final - self._model.year_first + 1

        self.strmap_pop = {CONST.POP_NOSEX : CONST.STR_POP_NOSEX,
                           CONST.POP_NEVER : CONST.STR_POP_NEVER,
                           CONST.POP_UNION : CONST.STR_POP_UNION,
                           CONST.POP_SPLIT : CONST.STR_POP_SPLIT,
                           CONST.POP_PWID  : CONST.STR_POP_PWID,
                           CONST.POP_BOTH  : CONST.STR_POP_BOTH,
                           CONST.POP_MSM   : CONST.STR_POP_MSM,
                           CONST.POP_TGW   : CONST.STR_POP_TGW}
        
        self.strmap_gender = {CONST.GENDER_FEM : CONST.STR_GENDER_FEM,
                              CONST.GENDER_MSC : CONST.STR_GENDER_MSC}

        # Process the population for standardized reporting. This may be
        # trading a lot of time and memory for (hopefully) simpler code here
        self._pop = self._pop_merge()
    
    def _pop_merge(self):
        # Define lists of column names
        str_child_neg = [CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE]
        str_child_hiv = [CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE, CONST.STR_HIV, CONST.STR_DTX]
        str_adult_neg = [CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE, CONST.STR_POP]
        str_adult_hiv = [CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE, CONST.STR_POP, CONST.STR_HIV, CONST.STR_DTX]

        # Convert population arrays to data frames
        pop_child_neg = array2frame(self._model.pop_child_neg, str_child_neg)
        pop_child_hiv = array2frame(self._model.pop_child_hiv, str_child_hiv)
        pop_adult_neg = array2frame(self._model.pop_adult_neg, str_adult_neg)
        pop_adult_hiv = array2frame(self._model.pop_adult_hiv, str_adult_hiv)

        # Convert indices to columns
        pop_child_neg.reset_index(inplace=True)
        pop_child_hiv.reset_index(inplace=True)
        pop_adult_neg.reset_index(inplace=True)
        pop_adult_hiv.reset_index(inplace=True)

        # Standardize population structure columns. -1 indicates HIV-negative. Sexual
        # activity before age 15 is assumed negligible.
        pop_child_neg.insert(3, CONST.STR_POP, CONST.POP_NOSEX)
        pop_child_neg.insert(4, CONST.STR_HIV, -1)
        pop_child_neg.insert(5, CONST.STR_DTX, -1)
        pop_child_hiv.insert(3, CONST.STR_POP, CONST.POP_NOSEX)
        pop_adult_neg.insert(4, CONST.STR_HIV, -1)
        pop_adult_neg.insert(5, CONST.STR_DTX, -1)

        # Convert adult age indices to age values
        pop_adult_neg[CONST.STR_AGE] = pop_adult_neg[CONST.STR_AGE] + CONST.AGE_ADULT_MIN
        pop_adult_hiv[CONST.STR_AGE] = pop_adult_hiv[CONST.STR_AGE] + CONST.AGE_ADULT_MIN

        # Concatenate and convert year indices to year values
        pop = pd.concat([pop_child_neg, pop_child_hiv, pop_adult_neg, pop_adult_hiv])
        pop[CONST.STR_YEAR] = pop[CONST.STR_YEAR] + self._model.year_first

        return pop
    
    def _pop_remap_age(self, age_breaks = [0, np.inf]):
        """! Return a copy of the population with ages replaced by age groups"""
        pop = self._pop.copy()
        pop[CONST.STR_AGE] = pd.cut(pop[CONST.STR_AGE], age_breaks, include_lowest=True, right=False)
        return pop
    
    def _pop_remap_age_gender(self, age_breaks = [0, np.inf]):
        """! Return a copy of the population with ages replaced by age groups
        and sex + male circumcision status replaced by gender """
        pop = self._pop_remap_age(age_breaks)
        pop[CONST.STR_GENDER] = ((pop[CONST.STR_POP] == CONST.POP_TGW) | (pop[CONST.STR_SEX] == CONST.SEX_FEMALE))
        pop[CONST.STR_GENDER] = pop[CONST.STR_GENDER].map({True  : CONST.GENDER_FEM,
                                                           False : CONST.GENDER_MSC})
        pop.drop(CONST.STR_SEX, axis=1)
        return pop

    def _pop_remap_age_sex(self, age_breaks = [0, np.inf]):
        """! Return a copy of the population with ages replaced by age groups
        and sex + male circumcision status replaced by sex alone """
        pop = self._pop_remap_age(age_breaks)
        pop[CONST.STR_SEX] = pop[CONST.STR_SEX].map({CONST.SEX_FEMALE : CONST.SEX_FEMALE,
                                                     CONST.SEX_MALE_U : CONST.SEX_MALE,
                                                     CONST.SEX_MALE_C : CONST.SEX_MALE})
        return pop
        
    def _series2df(self, series):
        """! Convert a pandas series to a pandas dataframe. This converts
        indices into columns to allow R style plotting and analysis"""
        return series.to_frame().reset_index()
    
    def pop_total(self, age_breaks=[0, np.inf]):
        """! Calculate the total population by year, sex, age
        @param age_breaks cut points for age aggregation
        """
        pop = self._pop_remap_age_gender(age_breaks)
        rval = pop.groupby([CONST.STR_YEAR, CONST.STR_GENDER, CONST.STR_AGE])[CONST.STR_VALUE].sum()
        return self._series2df(rval)
    
    def pop_risk(self, age_breaks=[0, np.inf]):
        """! Calculate the total population by year, sex, age, risk
        @param age_breaks cut points for age aggregation
        """
        aggr = self._pop_remap_age_gender(age_breaks)
        rval = aggr.groupby([CONST.STR_YEAR, CONST.STR_GENDER, CONST.STR_AGE, CONST.STR_POP])[CONST.STR_VALUE].sum()
        rval = self._series2df(rval)

        rval[CONST.STR_POP] = rval[CONST.STR_POP].astype("category")
        rval[CONST.STR_POP] = rval[CONST.STR_POP].cat.rename_categories(self.strmap_pop)

        rval[CONST.STR_GENDER] = rval[CONST.STR_GENDER].astype("category")
        rval[CONST.STR_GENDER] = rval[CONST.STR_GENDER].cat.rename_categories(self.strmap_gender)

        return rval
