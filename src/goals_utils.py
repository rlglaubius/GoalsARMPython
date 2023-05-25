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
    """! Load population size inputs
    @param tab_popsize an openpyxl workbook tab
    @return a dict mapping population size tags to parameter values. dict values are tuples: (male_value, female_value) 
    """
    vals = [tuple([cell.value for cell in row]) for row in tab_popsize['B2:C16']]
    keys = [cell[0].value for cell in tab_popsize['E2:E16']]
    rval = dict(zip(keys, vals))
    return {key : rval[key] for key in keys if keys != None} # Prune empty rows

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
    rirr_m = xlsx_load_range(tab_inci, 'B47', 'CD53').transpose()
    rirr_f = xlsx_load_range(tab_inci, 'B55', 'CD61').transpose()
    return inci[0], sirr[0], airr_m, airr_f, rirr_m, rirr_f

def xlsx_load_hiv_fert(tab_frr):
    age_no_art = xlsx_load_range(tab_frr, 'B2',  'CD8')
    cd4_no_art = xlsx_load_range(tab_frr, 'B11', 'B17')[:,0] # X[:,0] converts to 1d array
    age_on_art = xlsx_load_range(tab_frr, 'B20', 'B26')[:,0]
    adj = tab_frr['B29'].value # local adjustment factor
    return age_no_art.transpose() * adj, cd4_no_art, age_on_art * adj

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
    art_drop = xlsx_load_range(tab_art, 'B10', 'CD11').transpose()       # Annual dropout, %
    art_mrr  = xlsx_load_range(tab_art, 'B13', 'CD14').transpose()       # Mortality time trend rate ratio over time
    art_vs   = xlsx_load_range(tab_art, 'B16', 'CD23').transpose()       # Viral suppression on ART, %
    return art_elig, art_num, art_pct, art_drop, art_mrr, art_vs

def xlsx_load_mc_uptake(tab_mc):
    return xlsx_load_range(tab_mc, 'B3', 'CD19').transpose()

def xlsx_load_partner_rates(tab_partners):
    """! Load parameters that specify numbers of partners per year """
    time_trend = xlsx_load_range(tab_partners, 'B2',  'CD3') # sex x year matrix
    age_params = xlsx_load_range(tab_partners, 'B6',  'C7')  # par x sex matrix, par = (mean, scale)
    pop_ratios = xlsx_load_range(tab_partners, 'B10', 'C15') # pop x sex matrix, excludes sexually-inactive
    # reorder males and females to Goals ARM ordering
    return time_trend[[1,0],:], age_params[:,[1,0]], pop_ratios[:,[1,0]]

def xlsx_load_partner_prefs(tab_partners):
    """! Load parameters that specify age-based mixing and married/union status """
    age_prefs = xlsx_load_range(tab_partners, 'B18', 'B20')
    pop_prefs = xlsx_load_range(tab_partners, 'B23', 'C28')
    p_married = xlsx_load_range(tab_partners, 'B31', 'C33')
    return age_prefs, 0.01 * pop_prefs[:,[1,0]], 0.01 * p_married[:,[1,0]]

def xlsx_load_direct_clhiv(tab_clhiv):
    return xlsx_load_range(tab_clhiv, 'D3', 'CF86').transpose()
