# this file replicates table 1/2 in Demeshev/Malakhovskaya

# this file should be in /R subfolder
# Session - Set wd - To source file location


working_folder <- "../estimation/tables_rmsfe/"
data_for_tables_file_short <- "tables_rmsfe_raw_2.Rds"
forecast_filename_prefix <- "forecasts_" # name of analysed variable will be added automatically
model_info_filename_prefix <- "model_info_" # name of analysed variable will be added automatically


dir.create(working_folder)
# will throw warning if the folder exists
data_for_tables_file <- paste0(working_folder, data_for_tables_file_short)
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

df <- readr::read_csv("../data/df_2015_final.csv")
df <- dplyr::select(df, -time_y)

base_set_A <- c("ind_prod", "cpi", "ib_rate")
base_set_B <- c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price")
base_set_C <- c("employment", "ind_prod", "cpi", 
    "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
    "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
    "agriculture", "retail", "gov_balance", "export", "import")



all_vars <- base_set_C

all_rmsfe_results <- NULL

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
  
  part_A <- data_frame(var_set = "set_A", variable = set_A)
  part_B <- data_frame(var_set = "set_B", variable = set_B)
  part_C <- data_frame(var_set = "set_C", variable = set_C)
  var_set_info <- bind_rows(part_A, part_B, part_C)
  
  
  message("---------------------------------------------------")
  message("Calculating RMSFE for variable: ", analysed_variable)
  message("---------------------------------------------------")
  
  
  
  set.seed(36)
  banbura_data <- replicate_banbura(df, var_set_info = var_set_info, 
                                 set_delta_by = "AR1", num_AR_lags = 1)
  
  # saving useful info:
  rmsfe_one_variable <- banbura_data$rmsfe_wide
  rmsfe_one_variable$analysed <- analysed_variable
  all_rmsfe_results <- dplyr::bind_rows(all_rmsfe_results, rmsfe_one_variable)
  saveRDS(all_rmsfe_results, file = data_for_tables_file)
  
  forecast_filename <- paste0(working_folder, forecast_filename_prefix, analysed_variable, ".Rds")
  saveRDS(banbura_data$forecasts, file = forecast_filename)
  
  model_info_filename <- paste0(working_folder, model_info_filename_prefix, analysed_variable, ".Rds")
  saveRDS(banbura_data$model_info, file = model_info_filename)
  
}




