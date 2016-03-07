# this file replicates table 1/2 in Demeshev/Malakhovskaya

# this file should be in /R subfolder
# Session - Set wd - To source file location

dir.create("../estimation/tables_rmsfe/")
# will throw warning if the folder exists
data_for_tables_file <- "../estimation/tables_rmsfe/tables_rmsfe_raw.Rds"
# we overwrite file after each variable


source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")
source("500_replicate_banbura_function.R")

# need to run only once
source("200_load_after_eviews.R")


# if ind_prod/cpi/ib_rate - nothing to change
# if var \in set_6, but not var \in set_3 - add var to set_3
# if var \in set_23, but not var \in set_6 - add var to set_3 and set_6

base_set_A <- c("ind_prod", "cpi", "ib_rate")
base_set_B <- c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price")
base_set_C <- c("employment", "ind_prod", "cpi", 
    "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
    "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
    "agriculture", "retail", "gov_balance", "export", "import")



all_vars <- base_set_C

all_results <- NULL

for (analysed_variable in all_vars) {
  set_A <- base_set_A
  set_B <- base_set_B
  set_C <- base_set_C

  if (!analysed_variable %in% base_set_A) {
    set_A <- c(base_set_A, analysed_variable)
  }
  if (!analysed_variable %in% base_set_B) {
    set_B <- c(base_set_B, analysed_variable)
  }
  
  message("Calculating RMSFE for variable: ", analysed_variable)
  
  
  set.seed(36)
  temp_data <- replicate_banbura(set_A = set_A, set_B = set_B, set_C = set_C)
  temp_data$analysed <- analysed_variable
  
  all_results <- dplyr::bind_rows(all_results, temp_data)
  saveRDS(all_results, file = data_for_tables_file)
}




