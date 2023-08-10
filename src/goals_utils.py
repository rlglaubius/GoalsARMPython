import numpy as np
import src.goals_const as CONST

def xlsx_load_range(tab, cell_first, cell_final, dtype=np.float64, order="C"):
    """! Return the contents of a range in an Excel tab as a numpy array
    @param tab an openpyxl workbook tab
    @param cell_first A string specifying the first cell of the range (e.g., "A1")
    @param cell_final A string specifying the final cell of the range (e.g., "B2")
    @return a numpy 2-d array
    """
    return np.array([[cell.value for cell in row] for row in tab[cell_first:cell_final]], dtype=dtype, order=order)

def xlsx_load_config(tab_config):
    """! Load configuration options
    @param tab_config an openpyxl workbook tab
    @return a dict mapping configuration tags to parameter values (e.g., x["first.year"] = 1970)
    """
    vals = [cell[0].value for cell in tab_config['B2:B8']]
    keys = [cell[0].value for cell in tab_config['D2:D8']]
    return dict(zip(keys,vals))

def xlsx_load_popsize(tab_popsize):
    med_age_debut = xlsx_load_range(tab_popsize, 'B3', 'C3')
    med_age_union = xlsx_load_range(tab_popsize, 'B4', 'C4')
    avg_dur_union = tab_popsize['B7'].value
    kp_size = xlsx_load_range(tab_popsize, 'B10', 'G10')
    kp_stay = xlsx_load_range(tab_popsize, 'B11', 'G11', dtype=np.int32)
    kp_turnover = xlsx_load_range(tab_popsize, 'B12', 'G14')
    return med_age_debut[0], med_age_union[0], avg_dur_union, kp_size[0], kp_stay[0], kp_turnover

def xlsx_load_epi(tab_epi):
    """! Load various epidemiological inputs
    @param tab_epi an openpyxl workbook tab
    @return a dict mapping parameter tags to values
    """
    vals = [cell[0].value for cell in tab_epi['B3:B19']]
    keys = [cell[0].value for cell in tab_epi['D3:D19']]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if keys != None} # Prune empty rows

def xlsx_load_pasfrs(tab_fert):
    p = xlsx_load_range(tab_fert, 'B2', 'CD8')
    return 0.01 * p.transpose()

def xlsx_load_migr(tab_migr):
    migr   = xlsx_load_range(tab_migr, 'B2',  'CD3')
    dist_m = xlsx_load_range(tab_migr, 'B6',  'CD22')
    dist_f = xlsx_load_range(tab_migr, 'B25', 'CD41')
    return migr.transpose(), 0.01 * dist_m.transpose(), 0.01 * dist_f.transpose()

def xlsx_load_inci(tab_inci):
    inci   = xlsx_load_range(tab_inci, 'B2',  'CD2')
    sirr   = xlsx_load_range(tab_inci, 'B5',  'CD5')
    airr_m = xlsx_load_range(tab_inci, 'B9',  'CD25').transpose()
    airr_f = xlsx_load_range(tab_inci, 'B27', 'CD43').transpose()
    rirr_raw_m = xlsx_load_range(tab_inci, 'B47', 'CD53').transpose()
    rirr_raw_f = xlsx_load_range(tab_inci, 'B55', 'CD61').transpose()

    ## Reorganize incidence rate ratios by population for Goals input
    n_years = rirr_raw_m.shape[0]
    rirr_m = np.zeros((n_years, CONST.N_POP), dtype=np.double, order="C")
    rirr_f = np.zeros((n_years, CONST.N_POP), dtype=np.double, order="C")
    rirr_f[:,CONST.POP_NOSEX:(CONST.POP_FSW+1)] = rirr_raw_f[:,CONST.POP_NOSEX:(CONST.POP_FSW+1)]
    rirr_m[:,CONST.POP_NOSEX:(CONST.POP_MSM+1)] = rirr_raw_m[:,CONST.POP_NOSEX:(CONST.POP_MSM+1)]
    rirr_m[:,CONST.POP_TGW] = rirr_raw_f[:,6]

    return inci[0], sirr[0], airr_m, airr_f, rirr_m, rirr_f

def xlsx_load_hiv_fert(tab_frr):
    frr_age = xlsx_load_range(tab_frr, 'B2',  'CD8')
    frr_cd4 = xlsx_load_range(tab_frr, 'B11', 'B17')[:,0] # X[:,0] converts to 1d array
    frr_art = xlsx_load_range(tab_frr, 'B20', 'B26')[:,0]
    frr_laf = tab_frr['B29'].value # local adjustment factor
    return {'age' : frr_age.transpose(),
            'cd4' : frr_cd4,
            'art' : frr_art,
            'laf' : frr_laf}

def xlsx_load_adult_prog(tab_prog):
    dist = xlsx_load_range(tab_prog, 'B4',  'I9')  # CD4 distribution after primary infection
    prog = xlsx_load_range(tab_prog, 'B14', 'I19') # HIV progression rates with untreated HIV
    mort = xlsx_load_range(tab_prog, 'B24', 'I30') # HIV-related mortality rates with untreated HIV
    art1 = xlsx_load_range(tab_prog, 'B35', 'I41') # HIV-related mortality rates if on ART for [0,6) months
    art2 = xlsx_load_range(tab_prog, 'B46', 'I52') # HIV-related mortality rates if on ART for [6,12) months
    art3 = xlsx_load_range(tab_prog, 'B57', 'I63') # HIV-related mortality rates if on ART for [12,\infty) months
    return dist, prog, mort, art1, art2, art3

def xlsx_load_adult_art(tab_art):
    art_elig = xlsx_load_range(tab_art, 'B2',  'CD2', dtype=np.int32)[0] # CD4-based eligibility threshold
    art_num  = xlsx_load_range(tab_art, 'B4',  'CD5').transpose()        # PLHIV on ART, #
    art_pct  = xlsx_load_range(tab_art, 'B7',  'CD8').transpose()        # PLHIV on ART, %
    art_stop = xlsx_load_range(tab_art, 'B10', 'CD11').transpose()       # Annual ART interruption, %
    art_mrr  = xlsx_load_range(tab_art, 'B13', 'CD14').transpose()       # Mortality time trend rate ratio over time
    art_vs   = xlsx_load_range(tab_art, 'B16', 'CD23').transpose()       # Viral suppression on ART, %
    return art_elig, art_num, art_pct, art_stop, art_mrr, art_vs

def xlsx_load_mc_uptake(tab_mc):
    return xlsx_load_range(tab_mc, 'B3', 'CD19').transpose()

def xlsx_load_partner_rates(tab_partners):
    """! Load parameters that specify numbers of partners per year """
    time_trend = xlsx_load_range(tab_partners, 'B2',  'CD3') # sex x year matrix
    age_params = xlsx_load_range(tab_partners, 'B6',  'C7')  # par x sex matrix, par = (mean, scale)
    pop_ratios = xlsx_load_range(tab_partners, 'B10', 'C16') # pop x sex matrix, excludes sexually-inactive
    return time_trend, age_params, pop_ratios

def xlsx_load_partner_prefs(tab_partners):
    """! Load parameters that specify age-based mixing and married/union status """
    age_prefs = xlsx_load_range(tab_partners, 'B19', 'B21')
    pop_prefs = xlsx_load_range(tab_partners, 'B24', 'C30')
    p_married = xlsx_load_range(tab_partners, 'B33', 'C36')
    return age_prefs, pop_prefs, p_married.transpose()

def xlsx_load_mixing_levels(tab_mixing):
    """"! Load the behavioral risk group mixing matrix structure"""
    mix_levels = xlsx_load_range(tab_mixing, 'C3', 'N14')
    return mix_levels

def xlsx_load_contact_params(tab_contact):
    """"! Load coital frequency and condom use inputs"""
    sex_acts    = xlsx_load_range(tab_contact, 'B2',  'B5')
    condom_freq = xlsx_load_range(tab_contact, 'B8',  'CD11')
    pwid_force_infection = xlsx_load_range(tab_contact, 'B14', 'CD15')
    needle_sharing = xlsx_load_range(tab_contact, 'B18', 'CD18')
    return sex_acts.reshape((4)), condom_freq.transpose(), pwid_force_infection.transpose(), needle_sharing[0]

def xlsx_load_direct_clhiv(tab_clhiv):
    return xlsx_load_range(tab_clhiv, 'D3', 'CF86').transpose()

def xlsx_load_likelihood_pars(tab_lhood):
    """! Load parameters used in likelihood calculations but not in simulation
    @param tab_lhood an openpyxl workbook tab
    @return a dict mapping parameter tags to values
    """
    vals = [cell[0].value for cell in tab_lhood['B2:B5']]
    keys = [cell[0].value for cell in tab_lhood['D2:D5']]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if keys != None} # Prune empty rows

def xlsx_load_fitting_pars(tab_fit):
    """! Load initial values and fitting metadata for parameters varied during model fitting
    @param tab_fit an openpyxl workbook tab
    """
    last_row = 15
    keys = [cell[0].value for cell in tab_fit['H2:H%d' % (last_row)]]
    vals = [tuple(cell.value for cell in row) for row in tab_fit['B2:F%d' % (last_row)]]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if key != None} # Prune empty rows
