import openpyxl as xlsx
import os
import pandas as pd
import plotnine
import sys
import src.goals_model as Goals
import src.goals_const as CONST
from percussion import ancprev, hivprev


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
            case 'Male':   sex_min, sex_max = CONST.SEX_MALE, CONST.SEX_MALE + 1
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

def plot_fit_hiv(hivsim, hivdat):
    pass

def main(par_file, anc_file, hiv_file, out_file):
    # Initialize GoalsARM
    hivsim = Goals.Model()
    hivsim.init_from_xlsx(par_file)

    # Initialize likelihood evaluation objects
    ancdat = ancprev.ancprev(hivsim.year_first)
    hivdat = hivprev.hivprev(hivsim.year_first)
    ancdat.read_csv(anc_file)
    hivdat.read_csv(hiv_file)

    hiv_proj = hivdat.projection_template()    

    # Run a model simulation
    hivsim.project(2030)

    plot_fit_anc(hivsim, ancdat, "ancfit.tiff")

    ## Set ANC likleihood parameters
    ancdat.bias_ancss = hivsim.anc_par['ancss.bias']
    ancdat.bias_ancrt = hivsim.anc_par['ancrt.bias']
    ancdat.var_inflate_site = hivsim.anc_par['var.infl.site']
    ancdat.var_inflate_census = hivsim.anc_par['var.infl.census']

    ## Evaluate the ANC likelihood
    anc_proj = hivsim.births_exposed / hivsim.births.sum((1))
    lnlhood_anc = ancdat.likelihood(anc_proj)

    # Evaluate the HIV prevalence likelihood
    fill_hivprev_template(hivsim, hiv_proj)
    lnlhood_hiv = hivdat.likelihood(hiv_proj)

    print(lnlhood_anc)
    print(lnlhood_hiv)

    pass

if __name__ == "__main__":
    sys.stderr.write("Process %d\n" % (os.getpid()))
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
