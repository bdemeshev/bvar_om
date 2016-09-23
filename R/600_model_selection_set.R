# model selection set

library("dplyr")
library("ggplot2")
library("MCS")
library("reshape2")
library("parallel")

statistic <- "TR" # we've done with "Tmax"
best_model_prefix <- "best_models_TR_"

working_folder <- "../estimation/tables_rmsfe/"
forecast_filename_prefix <- "forecasts_"
model_info_filename_prefix <- "model_info_"

all_vars <- c("employment", "ind_prod", "cpi", 
              "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
              "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
              "agriculture", "retail", "gov_balance", "export", "import")

#analysed_variable <- all_vars[2]
for (analysed_variable in all_vars) {
  
  analysed_variable
  
  filename <- paste0(working_folder, forecast_filename_prefix, analysed_variable, ".Rds")
  
  
  d <- readRDS(filename)
  
  str(d)
  
  # d$actual_obs$variable
  
  T_min <- 133
  T_max <- max(d$actual_obs$t)
  
  actual <- d$actual_obs %>% filter(t %in% T_min:T_max, 
                                    variable == analysed_variable)
  bvar <- d$bvar_out_forecasts %>% filter(t %in% T_min:T_max, 
                                          variable == analysed_variable)
  var <- d$rwwn_var_out_forecasts %>% filter(t %in% T_min:T_max, 
                                             variable == analysed_variable)
  # 
  # t_count_var <- var %>% group_by(t) %>% summarize(count = n())
  # t_count_var %>% head(15)
  # t_count_bvar <- bvar %>% group_by(t) %>% summarize(count = n())
  # t_count_bvar %>% head(15)
  
  # total number of models (for t >= 144)
  # bvar: 864
  # var + rwwn: 312
  # we have h=1, ..., 12
  # for fixed h the number of models is 
  # bvar: 72
  # var: 26
  
  filename <- paste0(working_folder, model_info_filename_prefix, analysed_variable, ".Rds")
  info <- readRDS(filename)
  # 
  info$bvar_out_list -> blist
  blist <- mutate(blist, model_eq = interaction("b", as.numeric(factor(var_set)), as.numeric(factor(fit_set)), n_lag))
  
  
  blist %>% select(model_id = id, fit_set, var_set, n_lag,model_eq) -> blist
  bvar_new <- left_join(bvar, blist, by = "model_id") %>% select(-type)
  bvar_new %>% filter(fit_set == "ind+cpi") -> bvar_new
  
  
  info$rwwn_var_out_list -> vlist
  vlist <- mutate(vlist, model_eq = interaction("v", as.numeric(factor(type)), as.numeric(factor(var_set)), n_lag))
  vlist %>% select(model_id = id, type, var_set, n_lag, model_eq) -> vlist
  var_new <- left_join(var, vlist, by = "model_id") 
  
  for (hh in unique(bvar_new$h)) {
    message("hh = ", hh)
    bvar_new %>% filter(h == hh) -> bvar_h
    bvar_wide <- dcast(bvar_h, t ~ model_eq, value.var = "forecast")
    
    var_new %>% filter(h == hh) -> var_h
    var_wide <- dcast(var_h, t ~ model_eq, value.var = "forecast")
    
    all_wide <- left_join(bvar_wide, var_wide, by = "t") %>% select(-t) 
    actual_partial <- tail(actual$actual, nrow(all_wide))
    all_wide_actual <- matrix(actual_partial, nrow = nrow(all_wide), ncol = ncol(all_wide))
    loss <- (all_wide - all_wide_actual) ^ 2
    
    n_cores <- detectCores()
    cluster <- makeCluster(n_cores)
    best_models <- MCSprocedure(loss, cl = cluster, statistic = statistic)
    stopCluster(cluster)
    filename <- paste0(working_folder, "loss_", analysed_variable, "_", hh, ".Rds")
    saveRDS(loss, file = filename)
    filename <- paste0(working_folder, best_model_prefix, analysed_variable, "_", hh, ".Rds")
    saveRDS(best_models, file = filename)
    diff <- best_models@Info$elapsed.time
    message("time = ", diff, " ", attr(diff, "units"))
  }
}