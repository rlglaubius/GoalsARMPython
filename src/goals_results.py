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

        self.pop_strmap = {CONST.POP_NOSEX : CONST.STR_POP_NOSEX,
                           CONST.POP_NEVER : CONST.STR_POP_NEVER,
                           CONST.POP_UNION : CONST.STR_POP_UNION,
                           CONST.POP_SPLIT : CONST.STR_POP_SPLIT,
                           CONST.POP_PWID  : CONST.STR_POP_PWID,
                           CONST.POP_BOTH  : CONST.STR_POP_BOTH,
                           CONST.POP_MSM   : CONST.STR_POP_MSM,
                           CONST.POP_TGW   : CONST.STR_POP_TGW}

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
    
    def _pop_aggregate(self, age_breaks=[0, np.inf]):
        """! Return a copy of the population with age and circumcision status aggregated away"""
        dict_sex = {CONST.SEX_FEMALE : CONST.SEX_FEMALE,
                    CONST.SEX_MALE_U : CONST.SEX_MALE,
                    CONST.SEX_MALE_C : CONST.SEX_MALE}

        pop = self._pop.copy()
        pop[CONST.STR_AGE] = pd.cut(pop[CONST.STR_AGE], age_breaks, include_lowest=True, right=False)
        pop[CONST.STR_SEX] = pop[CONST.STR_SEX].map(dict_sex)
        return pop
    
    def _series2df(self, series):
        """! Convert a pandas series to a pandas dataframe. This converts
        indices into columns to allow R style plotting and analysis"""
        return series.to_frame().reset_index()

    ## TODO: phase out bigpop. pop_total does the same job more flexibly with less code
    def bigpop(self):
        """! Calculate the total population by year, sex, age"""
        rval = np.zeros((self._years, CONST.N_SEX, CONST.N_AGE), dtype=self._dtype, order=self._order)

        ## Add up children. We must sum males across circumcision states
        rval[:, CONST.SEX_MALE, CONST.AGE_CHILD_MIN:(CONST.AGE_CHILD_MAX+1)] \
            = self._model.pop_child_neg[:,CONST.SEX_MALE_U:,:].sum((1)) \
            + self._model.pop_child_hiv[:,CONST.SEX_MALE_U:,:,:,:].sum((1,3,4))
        rval[:, CONST.SEX_FEMALE, CONST.AGE_CHILD_MIN:(CONST.AGE_CHILD_MAX+1)] \
            = self._model.pop_child_neg[:,CONST.SEX_FEMALE,:] \
            + self._model.pop_child_hiv[:,CONST.SEX_FEMALE,:,:,:].sum((2,3))
        
        ## Add up adults. We must sum males across circumcision states
        rval[:, CONST.SEX_MALE, CONST.AGE_ADULT_MIN:] \
            = self._model.pop_adult_neg[:,CONST.SEX_MALE_U:,:,:].sum((1,3)) \
            + self._model.pop_adult_hiv[:,CONST.SEX_MALE_U:,:,:,:,:].sum((1,3,4,5))
        rval[:, CONST.SEX_FEMALE, CONST.AGE_ADULT_MIN:] \
            = self._model.pop_adult_neg[:,CONST.SEX_FEMALE,:,:].sum((2)) \
            + self._model.pop_adult_hiv[:,CONST.SEX_FEMALE,:,:,:,:].sum((2,3,4))
        
        return(rval)
    
    def pop_total(self, age_breaks=[0, np.inf]):
        """! Calculate the total population by year, sex, age"""
        pop = self._pop_aggregate(age_breaks)
        rval = pop.groupby([CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE])[CONST.STR_VALUE].sum()
        return self._series2df(rval)
    
    def pop_risk(self, age_breaks=[0, np.inf]):
        """! Calculate the total population by year, sex, age, risk"""
        aggr = self._pop_aggregate(age_breaks)
        rval = aggr.groupby([CONST.STR_YEAR, CONST.STR_SEX, CONST.STR_AGE, CONST.STR_POP])[CONST.STR_VALUE].sum()
        rval = self._series2df(rval)
        rval[CONST.STR_POP] = rval[CONST.STR_POP].astype("category")
        rval[CONST.STR_POP] = rval[CONST.STR_POP].cat.rename_categories(self.pop_strmap)
        return rval
