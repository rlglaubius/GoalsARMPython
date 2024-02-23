## +===+ Excel workbook constants +============================================+
## Names of input workbook tabs
XLSX_TAB_CONFIG       = "Config"
XLSX_TAB_PASFRS       = "FertilityInputs"
XLSX_TAB_MIGR         = "MigrInputs"
XLSX_TAB_INCI         = "DirectIncidenceInputs"
XLSX_TAB_PARTNER      = "PartnershipInputs"
XLSX_TAB_MIXNG_MATRIX = "MixingMatrix"
XLSX_TAB_CONTACT      = "ContactInputs"
XLSX_TAB_POPSIZE      = "PopSizeInputs"
XLSX_TAB_EPI          = "EpiInputs"
XLSX_TAB_STIPREV      = "STIPrevInputs"
XLSX_TAB_HIV_FERT     = "HIVFertilityInputs"
XLSX_TAB_ADULT_PROG   = "HIVDiseaseInputs"
XLSX_TAB_ADULT_ART    = "ARTAdultInputs"
XLSX_TAB_MALE_CIRC    = "MCInputs"
XLSX_TAB_DIRECT_CLHIV = "DirectCLHIV"
XLSX_TAB_LIKELIHOOD   = "LikelihoodInputs"
XLSX_TAB_FITTING      = "FittingInputs"

## The Excel file specifies inputs for 1970-2050, the projection
## uses a subset of these
XLSX_FIRST_YEAR = 1970
XLSX_FINAL_YEAR = 2050

## Configuration tab tags
CFG_FIRST_YEAR       = "first.year"
CFG_FINAL_YEAR       = "final.year"
CFG_UPD_NAME         = "upd.file"
CFG_USE_UPD_PASFRS   = "use.upd.pasfrs"
CFG_USE_UPD_MIGR     = "use.upd.migr"
CFG_USE_DIRECT_INCI  = "use.direct.inci"
CFG_USE_DIRECT_CLHIV = "use.direct.clhiv"

## EpiInputs tab tags
EPI_TRANSMIT_F2M     = "transmit.f2m"
EPI_TRANSMIT_M2F     = "transmit.m2f"
EPI_TRANSMIT_M2M     = "transmit.m2m"
EPI_TRANSMIT_CHRONIC = "transmit.chronic"
EPI_TRANSMIT_PRIMARY = "transmit.primary"
EPI_TRANSMIT_SYMPTOM = "transmit.symptom"
EPI_TRANSMIT_ART_VS  = "transmit.art.vs"
EPI_TRANSMIT_ART_VF  = "transmit.art.vf"
EPI_TRANSMIT_STI_NEG = "effect.sti.neg"
EPI_TRANSMIT_STI_POS = "effect.sti.pos"
EPI_EFFECT_VMMC      = "effect.vmmc"
EPI_EFFECT_CONDOM    = "effect.condom"
EPI_ART_MORT_WEIGHT  = "art.mort.weight"
EPI_INITIAL_YEAR     = "seed.time"
EPI_INITIAL_PREV     = "seed.prev"

## LikelihoodInputs tab tags
LHOOD_ANCSS_BIAS     = "ancss.bias"
LHOOD_ANCRT_BIAS     = "ancrt.bias"
LHOOD_VARINFL_SITE   = "var.infl.site"
LHOOD_VARINFL_CENSUS = "var.infl.census"

## FittingInputs tab tags
FIT_INITIAL_PREV        = EPI_INITIAL_PREV
FIT_TRANSMIT_F2M        = EPI_TRANSMIT_F2M
FIT_TRANSMIT_M2F        = EPI_TRANSMIT_M2F
FIT_TRANSMIT_STI_NEG    = EPI_TRANSMIT_STI_NEG
FIT_TRANSMIT_STI_POS    = EPI_TRANSMIT_STI_POS
FIT_FORCE_PWID          = "force.pwid"
FIT_LT_PARTNER_F        = "lt.partner.f"
FIT_LT_PARTNER_M        = "lt.partner.m"
FIT_PARTNER_AGE_MEAN_F  = "partner.age.mean.f"
FIT_PARTNER_AGE_MEAN_M  = "partner.age.mean.m"
FIT_PARTNER_AGE_SCALE_F = "partner.age.scale.f"
FIT_PARTNER_AGE_SCALE_M = "partner.age.scale.m"
FIT_PARTNER_POP_FSW     = "partner.pop.fsw"
FIT_PARTNER_POP_CLIENT  = "partner.pop.client"
FIT_PARTNER_POP_MSM     = "partner.pop.msm"
FIT_PARTNER_POP_TGW     = "partner.pop.tgw"
FIT_ASSORT_GEN          = "assort.gen"
FIT_ASSORT_FSW          = "assort.fsw"
FIT_ASSORT_MSM          = "assort.msm"
FIT_ASSORT_TGW          = "assort.tgw"
FIT_HIV_FRR_LAF         = "hiv.frr.laf"
FIT_ANCSS_BIAS          = LHOOD_ANCSS_BIAS
FIT_ANCRT_BIAS          = LHOOD_ANCRT_BIAS
FIT_VARINFL_SITE        = LHOOD_VARINFL_SITE
FIT_VARINFL_CENSUS      = LHOOD_VARINFL_CENSUS

## +===+ Model constants +=====================================================+
## Model constants are aligned with values in GoalsARM_Core DPConst.H

## +-+ Sex and male circumcision status +--------------------------------------+
SEX_FEMALE = 0 # female
SEX_MALE   = 1
SEX_MALE_U = 1 # male, uncircumcised
SEX_MALE_C = 2 # male, circumcised
SEX_MIN    = 0
SEX_MAX    = 1
SEX_MC_MIN = 0
SEX_MC_MAX = 2
N_SEX    = SEX_MAX - SEX_MIN + 1
N_SEX_MC = SEX_MC_MAX - SEX_MC_MIN + 1

## Map from sex and circumcision status (SEX_FEMALE, SEX_MALE_U, SEX_MALE_C) to sex alone
sex = [SEX_FEMALE, SEX_MALE, SEX_MALE]

## +-+ Age constants +---------------------------------------------------------+
AGE_MIN = 0
AGE_MAX = 80
AGE_CHILD_MIN = AGE_MIN
AGE_CHILD_MAX = 14
AGE_ADULT_MIN = AGE_CHILD_MAX + 1
AGE_ADULT_MAX = AGE_MAX
AGE_BIRTH_MIN = AGE_ADULT_MIN # minimum reproductive age
AGE_BIRTH_MAX = 49            # maximum reproductive age

N_AGE = AGE_MAX - AGE_MIN + 1
N_AGE_CHILD = AGE_CHILD_MAX - AGE_CHILD_MIN + 1 # number of child ages
N_AGE_ADULT = AGE_ADULT_MAX - AGE_ADULT_MIN + 1 # number of adult ages
N_AGE_BIRTH = AGE_BIRTH_MAX - AGE_BIRTH_MIN + 1 # number of reproductive ages

## +-+ Behavioral risk constants +---------------------------------------------+
POP_NOSEX = 0 # male or female, never had sex
POP_NEVER = 1 # male or female, never married
POP_UNION = 2 # male or female, married or in stable union
POP_SPLIT = 3 # male or female, previously married
POP_PWID = 4  # male or female, people who inject drugs
POP_FSW = 5   # female only, female sex workers
POP_CSW = 5   # male only, male clients of female sex workers
POP_BOTH = 5  # male or female, female sex workers or their clients
POP_MSM = 6   # male only, men who have sex with men
POP_TGW = 7   # male only, transgender women

POP_MIN = 0
POP_MAX = 7
N_POP = POP_MAX - POP_MIN + 1

POP_KEY_MIN = POP_PWID
POP_KEY_MAX = POP_MAX
N_POP_KEY = POP_KEY_MAX - POP_KEY_MIN + 1

## +-+ HIV infection stage constants +-----------------------------------------+
## Child categories (CD4 percent), ages 0-4
HIV_PRC_GT_30 = 0
HIV_PRC_26_30 = 1
HIV_PRC_21_25 = 2
HIV_PRC_16_20 = 3
HIV_PRC_11_15 = 4
HIV_PRC_05_10 = 5
HIV_PRC_LT_05 = 6

## Child categories, ages 5-14
HIV_NUM_GEQ_1000 = 0
HIV_NUM_750_1000 = 1
HIV_NUM_500_750  = 2
HIV_NUM_350_500  = 3
HIV_NUM_200_350  = 4
HIV_NUM_LT_200   = 5
HIV_NUM_INVALID  = 6 # Not used

## Adult (ages 15+) categories, defined by CD4 counts
HIV_PRIMARY = 0
HIV_GEQ_500 = 1
HIV_350_500 = 2
HIV_200_350 = 3
HIV_100_200 = 4
HIV_050_100 = 5
HIV_000_050 = 6

HIV_CHILD_PRC_MIN = 0 # Children aged 0-4
HIV_CHILD_PRC_MAX = 6 # Children aged 0-4
HIV_CHILD_NUM_MIN = 0 # Children aged 5-14
HIV_CHILD_NUM_MAX = 5 # Children aged 5-14

HIV_CHILD_MIN = 0 # Children 0-14 (will include HIV_NUM_INVALID for ages 5-14)
HIV_CHILD_MAX = 6

HIV_ADULT_MIN = 0
HIV_ADULT_MAX = 6

HIV_MIN = 0
HIV_MAX = 6

N_HIV_CHILD_PRC = HIV_CHILD_PRC_MIN - HIV_CHILD_PRC_MAX + 1
N_HIV_CHILD_NUM = HIV_CHILD_NUM_MAX - HIV_CHILD_NUM_MIN + 1
N_HIV_CHILD     = HIV_CHILD_MAX - HIV_CHILD_MIN + 1
N_HIV_ADULT     = HIV_ADULT_MAX - HIV_ADULT_MIN + 1
N_HIV           = HIV_MAX - HIV_MIN + 1

## +-+ HIV diagnosis (Dx) and treatment (Tx) constants +---------------------+
		
DTX_UNAWARE = 0 # HIV-positive but status-unaware
DTX_AWARE   = 1 # HIV-positive and status-aware but not on ART
DTX_PREV_TX = 2 # HIV-positive and previously on ART
DTX_ART1    = 3 # Months [0,6) on ART
DTX_ART2    = 4 # Months [6,12) on ART
DTX_ART3    = 5 # Months [12,\infty) on ART

DTX_MIN = 0
DTX_MAX = 5

# Bounds for looping just for off-ART states
DTX_OFF_MIN = 0
DTX_OFF_MAX = 2

# Bounds for looping just over on-ART states
DTX_ART_MIN = 3
DTX_ART_MAX = 5

N_DTX = DTX_MAX - DTX_MIN + 1
N_ART = DTX_ART_MAX - DTX_ART_MIN + 1

## +===+ Fitting constants +===================================================+
DIST_LOGNORMAL = 'Lognormal'
DIST_NORMAL    = 'Normal'
DIST_GAMMA     = 'Gamma'
DIST_BETA      = 'Beta'
