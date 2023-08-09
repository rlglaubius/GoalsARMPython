import numpy as np
import openpyxl as xlsx
import os
import pandas as pd
import plotnine
import scipy.optimize as optimize
import scipy.stats as stats
import sys
import src.goals_model as Goals
import src.goals_const as CONST
import time
from percussion import ancprev, hivprev

## TODO: change to use a dict of parameters that may be fitted. Add a column to 
## the input workbook "FittingParams" tab that allows parameter prior specification
## and a True/False indicator for whether to vary each input. Extend the workbook reader in goals_model
## and goals_util as needed to support this selection. This is important since
## we may only want to vary (e.g.) behavioral parameters for a key population if
## there are data to fit to for that population.

## Install percussion, the likelihood model package:
## python -m pip install git+https://rlglaubius:${token}@github.com/rlglaubius/percussion.git

def fill_hivprev_template(hivsim, template):
    # Poor practice first cut: loop. Since the template rows are invariant from
    # one simulation to the next, we may be able to precompute and store the
    # indices needed to calculate the model's prevalence estimates
    nobs = len(template.index)
    pop = template['Population']
    sex = template['Gender']
    yidx = template['Year'] - hivsim.year_first     # convert from years to indices relative to the first year of projection
    amin = template['AgeMin'] - CONST.AGE_ADULT_MIN # convert from ages to years since entry into the adult model
    amax = template['AgeMax'] - CONST.AGE_ADULT_MIN + 1
    for row in range(nobs):
        match pop[row]:
            case 'All': pop_min, pop_max = CONST.POP_MIN, CONST.POP_MAX + 1
            case 'FSW': pop_min, pop_max = CONST.POP_FSW, CONST.POP_FSW + 1
            case 'MSM': pop_min, pop_max = CONST.POP_MSM, CONST.POP_MSM + 1
            case _: sys.stderr.write("Error: Unrecognized population %s" % (pop[row]))

        match sex[row]:
            case 'All':    sex_min, sex_max = CONST.SEX_MC_MIN, CONST.SEX_MC_MAX + 1
            case 'Female': sex_min, sex_max = CONST.SEX_FEMALE, CONST.SEX_FEMALE + 1
            case 'Male':   sex_min, sex_max = CONST.SEX_MALE_U, CONST.SEX_MALE_C + 1
            case _: sys.stderr.write("Error: Unrecognized gender %s" % (sex[row]))

        pop_hiv = hivsim.pop_adult_hiv[yidx[row], sex_min:sex_max, amin[row]:amax[row], pop_min:pop_max, :, :].sum()
        pop_neg = hivsim.pop_adult_neg[yidx[row], sex_min:sex_max, amin[row]:amax[row], pop_min:pop_max].sum()
        template.at[row,'Prevalence'] = pop_hiv / (pop_hiv + pop_neg)

def plot_fit_anc(hivsim, ancdat, tiffname):
    anc_data = ancdat.anc_data.copy()
    anc_data['Source'] = ['Census' if site=='Census' else 'ANC-%s' % (kind) for site, kind in zip(anc_data['Site'], anc_data['Type'])]

    mod_data = pd.DataFrame({'Year'     : list(range(hivsim.year_first, hivsim.year_final + 1)),
                             'Prevalence' : hivsim.births_exposed / hivsim.births.sum((1)),
                             'Site'       : 'Goals',
                             'Source'     : 'Goals'})

    p = (plotnine.ggplot(anc_data[anc_data['Site'] != 'Census'])
         + plotnine.aes(x='Year', y='Prevalence', color='Source', group='Site')
         + plotnine.geom_line()
         + plotnine.geom_point()
         + plotnine.geom_line(data=anc_data[anc_data['Site'] == 'Census'])
         + plotnine.geom_point(data=anc_data[anc_data['Site'] == 'Census'])
         + plotnine.geom_line(data=mod_data)
         + plotnine.theme_bw())
    p.save(filename=tiffname, dpi=600, units="in", width=6.5, height=5.0, pil_kwargs={"compression" : "tiff_lzw"})

def plot_fit_hiv(hivsim, hivdat, tiffname):
    hiv_data = hivdat.hiv_data.copy()
    hiv_data['Age'] = ['%s-%s' % (amin, amax) for (amin, amax) in zip(hiv_data['AgeMin'], hiv_data['AgeMax'])]
    hiv_data.rename(columns={'Value' : 'Prevalence'}, inplace=True)

    # Capture HIV prevalence every year for the populations with data
    pop_frame = hiv_data.groupby(['Population', 'Gender', 'AgeMin', 'AgeMax']).size().reset_index(name='Prevalence')
    yrs_frame = pd.DataFrame({'Year' : range(hivsim.year_first, hivsim.year_final + 1)})
    mod_data = yrs_frame.join(pop_frame, how='cross')
    fill_hivprev_template(hivsim, mod_data)
    mod_data['Age'] = ['%s-%s' % (amin, amax) for (amin, amax) in zip(mod_data['AgeMin'], mod_data['AgeMax'])]

    p = (plotnine.ggplot(hiv_data[hiv_data['AgeMax'] > 14])
         + plotnine.aes(x='Year', y='Prevalence', color='Gender')
         + plotnine.geom_point()
         + plotnine.geom_line(data=mod_data[mod_data['AgeMax'] > 14])
         + plotnine.facet_grid('Population~Age', scales='free_y')
         + plotnine.theme_bw()
         + plotnine.theme(axis_text_x = plotnine.element_text(angle=90)))
    p.save(filename=tiffname, dpi=600, units="in", width=16, height=9, pil_kwargs={"compression" : "tiff_lzw"})

class GoalsFitter:
    def __init__(self, par_xlsx, anc_csv, hiv_csv):
        self.init_hivsim(par_xlsx)
        self.init_data_anc(anc_csv)
        self.init_data_hiv(hiv_csv)
        self.eval_count = 0

    def init_hivsim(self, par_xlsx):
        self.hivsim = Goals.Model()
        self.hivsim.init_from_xlsx(par_xlsx)
        self.year_first = self.hivsim.year_first
        self.year_final = self.hivsim.year_final

    def init_data_anc(self, anc_csv):
        self._ancdat = ancprev.ancprev(self.year_first)
        self._ancdat.read_csv(anc_csv)

    def init_data_hiv(self, hiv_csv):
        self._hivdat = hivprev.hivprev(self.year_first)
        self._hivdat.read_csv(hiv_csv)
        self._hivest = self._hivdat.projection_template()

    def prior(self, params):
        """! Prior density on log scale """
        return (stats.norm.logpdf(params[0], 0.15, 1.00) +          # ANC-SS bias term
                stats.norm.logpdf(params[1], 0.00, 1.00) +          # ANC-RT calibration term
                stats.expon.logpdf(params[2], scale=1.0 / 0.015) +  # ANC variance inflation term, site
                stats.expon.logpdf(params[3], scale=1.0 / 0.015))   # ANC variance inflation term, census

    def likelihood(self, params):
        """! Log-likelihood """
        # self.project(params)
        # fill_hivprev_template(self.hivsim, self._hivest)
        self.eval_count += 1
        self._ancdat.set_parameters(params[0], params[1], params[2], params[3]) # TODO: once testing with real parameters, cut this line (redundant with self.project)
        self._ancest = self.hivsim.births_exposed / self.hivsim.births.sum((1))
        lhood_hiv = self._hivdat.likelihood(self._hivest)
        lhood_anc = self._ancdat.likelihood(self._ancest)
        return lhood_hiv + lhood_anc

    def posterior(self, params):
        """"! Posterior density on log scale """
        return self.prior(params) + self.likelihood(params)
    
    def project(self, params):
        """! Set fitting parameter values into the model then run a projection """
        self._ancdat.set_parameters(params[0], params[1], params[2], params[3])
        self.hivsim.project(self.year_final)

    def calibrate(self, ancss_bias, ancrt_bias, var_infl_site, var_infl_census, method='Nelder-Mead'):
        """! Calibrate the model to HIV prevalence data
        ancss_bias      -- Initial ANC-SS bias parameter value
        ancrt_bias      -- Initial ANC-SS bias parameter value
        var_infl_site   -- Variance inflation to account for non-sampling error in ANC data at site level
        var_infl_census -- Variance inflation to account for non-sampling error in ANC data at census level
        """
        self.eval_count = 0

        par_init = np.array([ancss_bias, ancrt_bias, var_infl_site, var_infl_census])

        ## TODO: As a special case when calibrating ANC likelihood parameters
        ## only, we only need to calculate the projection once. Once we start
        ## estimating more parameters we will need to reproject repeatedly in
        ## the likelihood calculation instead
        self.project(par_init)
        fill_hivprev_template(self.hivsim, self._hivest)

        bounds = optimize.Bounds(lb = [-np.inf, -np.inf, 0.0, 0.0],
                                 ub = [+np.inf, +np.inf, +np.inf, +np.inf])
        par = optimize.minimize(lambda p : -self.posterior(p),
                                par_init,
                                method = method,
                                bounds = bounds,
                                callback = lambda res : print(res))
        return par

def main(par_file, anc_file, hiv_file, out_file):
    Fitter = GoalsFitter(par_file, anc_file, hiv_file)
    print(Fitter.eval_count)

    # Get initial conditions from the fitter's model instance
    ancss_bias    = Fitter.hivsim.anc_par['ancss.bias']
    ancrt_bias    = Fitter.hivsim.anc_par['ancrt.bias']
    varinf_site   = Fitter.hivsim.anc_par['var.infl.site']
    varinf_census = Fitter.hivsim.anc_par['var.infl.census']

    # TODO: Adjust lifetime partners levels for males and females

    par = Fitter.calibrate(ancss_bias, ancrt_bias, varinf_site, varinf_census, method='L-BFGS-B')

    print(par.x)
    print(Fitter.eval_count)
    ## TODO: use the fitted model to update the goodness-of-fit plots
    pass


    # TODO: outfile - xlsx workbook with fitted parameter estimates

    # # Initialize GoalsARM
    # hivsim = Goals.Model()
    # hivsim.init_from_xlsx(par_file)

    # # Initialize likelihood evaluation objects
    # ancdat = ancprev.ancprev(hivsim.year_first)
    # hivdat = hivprev.hivprev(hivsim.year_first)
    # ancdat.read_csv(anc_file)
    # hivdat.read_csv(hiv_file)

    # hiv_proj = hivdat.projection_template()    

    # # Run a model simulation
    # hivsim.project(2030)

    # plot_fit_anc(hivsim, ancdat, "ancfit.tiff")
    # plot_fit_hiv(hivsim, hivdat, "hivfit.tiff")

    # ## Set ANC likleihood parameters
    # ancdat.bias_ancss = hivsim.anc_par['ancss.bias']
    # ancdat.bias_ancrt = hivsim.anc_par['ancrt.bias']
    # ancdat.var_inflate_site = hivsim.anc_par['var.infl.site']
    # ancdat.var_inflate_census = hivsim.anc_par['var.infl.census']

    # ## Evaluate the ANC likelihood
    # anc_proj = hivsim.births_exposed / hivsim.births.sum((1))
    # lnlhood_anc = ancdat.likelihood(anc_proj)

    # # Evaluate the HIV prevalence likelihood
    # fill_hivprev_template(hivsim, hiv_proj)
    # lnlhood_hiv = hivdat.likelihood(hiv_proj)

    # print(lnlhood_anc)
    # print(lnlhood_hiv)

    # pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
    time_start = time.time()
    if len(sys.argv) == 1:
        par_file = "inputs/mwi-2023-inputs.xlsx"
        anc_file = "inputs/mwi-2023-anc-prev.csv"
        hiv_file = "inputs/mwi-2023-hiv-prev.csv"
        out_file = ""
        main(par_file, anc_file, hiv_file, out_file)
    elif len(sys.argv) < 4:
        sys.stderr.write("USAGE: %s <input_param>.xlsx <anc_data>.csv <hiv_data>.csv <output_param>.xlsx" % (sys.argv[0]))
    else:
        par_file = sys.argv[1]
        anc_file = sys.argv[2]
        hiv_file = sys.argv[3]
        out_file = sys.argv[4]
        main(par_file, anc_file, hiv_file, out_file)
    print("Completed in %s seconds" % (time.time() - time_start))
