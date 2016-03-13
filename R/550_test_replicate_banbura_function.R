# test for replicate banbura function


source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")
source("500_replicate_banbura_function.R")

# need to run only once
source("200_load_after_eviews.R")


df <- readr::read_csv("../data/df_2015_final.csv")

temp_data <- replicate_banbura(df)



# to step inside replicate_banbura() in case of fail:

parallel = "off" # c("off", "unix", "windows")
ncpu = 30
h_max = 12
fast_forecast = TRUE
keep = 5000
verbose = FALSE
testing_mode = FALSE
carriero_hack = FALSE
num_AR_lags = NULL
set_delta_by = "KPSS"
c_0 = 0.5
c_1 = 1
set_A = c("ind_prod", "cpi", "ib_rate")
set_B = c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price")
set_C = c("employment", "ind_prod", "cpi", 
          "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
          "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
          "agriculture", "retail", "gov_balance", "export", "import")
v_prior = "m+2"
T_common = 120 
p_max = 12



