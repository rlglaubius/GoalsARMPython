import numpy as np
import src.goals_const as CONST
import src.goals_model as Goals

class Results:
    def __init__(self, model):
        """! Initialize the object with a Goals model
        @param model The Goals model
        """
        self._model = model
        self._dtype = model._dtype
        self._order = model._order
    
    def bigpop(self):
        """! Calculate the total population by year, sex, age"""
        rval = np.zeros((self._model.year_final - self._model.year_first + 1, CONST.N_SEX, CONST.N_AGE),
                        dtype=self._dtype, order=self._order)

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
