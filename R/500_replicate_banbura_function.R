# 500_replicate_banbura_function.R

# input: '../data/df_2015_final.csv' output:

# needs the following functions:
# source("400_model_funs.R")
# source("400_model_lists.R")
# source("500_banbura_funs.R")

create_var_set_info <- function() {
  add_A <- dplyr::data_frame(var_set = "set_A", variable = c("ind_prod", "cpi", "ib_rate"))
  add_B <- dplyr::data_frame(var_set = "set_B", variable = c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))
  add_C <- dplyr::data_frame(var_set = "set_C", variable = c("employment", "ind_prod", "cpi", 
                                                             "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                                                             "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                                                             "agriculture", "retail", "gov_balance", "export", "import"))
  
  var_set_info <- dplyr::bind_rows(add_A, add_B, add_C)
  
  return(var_set_info)
}


create_fit_set_info <- function() {
  fit_set_2vars <- dplyr::data_frame(variable = c("ind_prod", "cpi"), fit_set = "ind+cpi")
  fit_set_3vars <- dplyr::data_frame(variable = c("ind_prod", "cpi", "ib_rate"), fit_set = "ind+cpi+rate")
  
  fit_set_info <- dplyr::bind_rows(fit_set_2vars, fit_set_3vars)
  return(fit_set_info)
}

#' Returns the use of variable sets. The smallest variable set will be market as
#' `target`, BVAR will be estimated for all var-sets, VAR only for these with 
#' 10 variables or less.
#' @param var_set_info data.frame with `var_set` and `variable` columns 
create_var_set_use <- function(var_set_info) {
  var_set_use <- dplyr::group_by(var_set_info, var_set) %>%
    dplyr::summarise(n_vars = n())
  var_set_use$estimate_bvar <- TRUE
  var_set_use <- dplyr::mutate(var_set_use, estimate_var = (n_vars <= 10))
  var_set_use <- dplyr::arrange(var_set_use, n_vars)
  var_set_use$target <- FALSE
  var_set_use$target[1] <- TRUE
  return(var_set_use)
}
  
#' @param parallel either off/unix/windows
#' @param ncpu number of cpu for parallel computations, ignored if parallel=='off'
#' @param h_max maximum forecast horizont for VAR and BVAR 
#' @param fast_forecast TRUE/FALSE. TRUE means that posterior means of coefficients are used for forecast
#' @param keep number of simulations from posterior (used only if fast_forecast is FALSE). 
#' keep is automatically set to zero if fast_forecast is TRUE
#' @param verbose turn on/off messages from functions 
#' @param testing_mode TRUE/FALSE. If TRUE then less lambdas are estimated, see 400_model_lists.R
#' @param carriero_hack FALSE/TRUE. If TRUE then we use wrong formulas from carriero code 
#' for dummy cNIW without square root for sigma^2
#' @param v_prior simple formula for prior nu. By default is "m+2", may be "T_dummy" or a constant
#' @param c_0 value for delta for variables classified as stationary (default is 0.5)
#' @param c_1 value for delta for variables classified as non-stationary (default is 1)
#' @param set_delta_by how we classify stationary and non-stationary variables
#' maybe: "KPSS", "ADF", "global AR1", "AR1" or a number between 0 and 1.
#' This influences the delta hyper-parameter in priors.
#' "AR1" --- AR(1) is estimated using only observations for estimated model
#' "global AR1" --- AR(1) is estimated using all available observations (does not depend on the model)
#' If it is numeric then delta is equal for stationary and non-stationary series
#' @param num_AR_lags number or NULL (by default). Number of lags in AR() model used to estimate sigma^2 
#' If NULL then p (number of lags in VAR/BVAR) will be used
#' @param T_common (by default 120) number of observations for in-sample forecast
#' @param p_max (by default 12) maximum number of lags
#' @param 
#' @param var_set_info data.frame with `var_set` and `variable` columns
#' @param fit_set_info data.frame with `fit_set` and `variable` columns
#' @param var_set_use data.frame with `var_set`, `n_vars`, `estimate_bvar`, `estimate_var`, `target` columns.
#' The `estimate_bvar`, `estimate_var`, `target` columns are TRUE/FALSE. Unique TRUE in `target` column is required.
#' This unique `target` var_set serves as a base to calculate lambda.
replicate_banbura <- function(df, parallel = c("off", "unix", "windows"), ncpu = 30,
                              h_max = 12,
                              fast_forecast = TRUE, keep = 5000,
                              verbose = FALSE, testing_mode = FALSE,
                              carriero_hack = FALSE, num_AR_lags = NULL,
                              set_delta_by = "KPSS",
                              c_0 = 0.5, c_1 = 1,
                              var_set_info = create_var_set_info(),
                              fit_set_info = create_fit_set_info(),
                              var_set_use = create_var_set_use(var_set_info),
                              v_prior = "m+2",
                              T_common = 120, p_max = 12) {
  
  ########################################
  ########################## begin set-up part #######################
  
  parallel <- match.arg(parallel)
  
  if ("t" %in% colnames(df)) {
    message("Data set contains the variable named `t`, it will be overwritten.")
  }

  if (p_max != 12) {
    message("Function has not been tested for p_max != 12. Maybe 12 is still hardcoded somewhere")
  }
  
  # subset t:
  # df <- tail(df, -p_max) ### ????????????????

  
  T_available <- nrow(df)  # number of observations
  
  if (fast_forecast) {
    keep <- 0
  }

  
  
  #############################################################
  ##### end set-up part #########################
  
  ####### step 0 (before banbura procedure) melting actual observations
  df <- dplyr::mutate(df, t = row_number())
  actual_obs <- reshape2::melt(df, id.vars = "t") %>% rename(actual = value) %>% mutate(variable = as.character(variable))
  
  
  ##### banbura step 1 calculate msfe-0. Estimate RWWN (random walk OR white noise
  ##### model)
  
  # set delta (prior hyperparameter)
  
  if (set_delta_by %in% c("ADF", "KPSS")) {
    deltas <- delta_i_prior(df, test = set_delta_by, remove_vars = "t", 
                            c_0 = c_0, c_1 = c_1)
  }
  
  if (set_delta_by == "global AR1") {
    deltas <- delta_i_from_ar1(df, remove_vars = "t")
  }
  
  if ((set_delta_by == "AR1") | is.finite(set_delta_by)) {
    # if 'AR1' then will be calculated automatically by bvar_conj_setup
    deltas <- data_frame(variable = setdiff(colnames(df), "t"), delta = set_delta_by)
  }
  
  
  
  # Important!!!  value of delta determines prior 
  # RW/WN determines model for comparison
  
  # always compare with Random Walk:
  deltas <- mutate(deltas, rw_wn = "rw", variable = as.character(variable))

  
  # estimate all RW and WN models
  rwwn_list <- create_rwwn_list(p_max = p_max, T_common = T_common)
  message("Estimate RW and WN")
  rwwn_list <- estimate_models(rwwn_list, parallel = parallel, verbose = verbose,
                               var_set_info = var_set_info, df = df, 
                               deltas = deltas, num_AR_lags = num_AR_lags,
                               carriero_hack = carriero_hack,
                               v_prior = v_prior, keep = keep)
  
  # forecast all RW and WN models
  rwwn_forecast_list <- data_frame(model_id = c(1, 2), h = NA, type = "in-sample")
  message("Forecast RW and WN")
  rwwn_forecasts <- forecast_models(rwwn_forecast_list, rwwn_list, verbose = verbose,
                                    var_set_info = var_set_info, df = df,
                                    fast_forecast = fast_forecast)
  
  # calculate all msfe-0 two ways to calculate msfe (we use b) a) use all available
  # predictions for each model b) use only common available predictions for all
  # models
  
  # rwwn_forecasts %>% group_by(model_id) %>% summarise(Tf_start = min(t), Tf_end = max(t))
  
  # plot_forecast(rwwn_forecasts, var_name = "m2", mod_id = 2, actual_obs = actual_obs)
  
  
  # build table with corresponding msfe-0 (RW or WN)
  msfe0_all <- get_msfe(rwwn_forecasts, actual_obs, 
                        models = dplyr::select(rwwn_list, id, type), msfe_name = "msfe")
  
  
  # add deltas info to msfe0_all
  
  msfe0 <- left_join(deltas, msfe0_all, by = c(variable = "variable", rw_wn = "type"))
  
  ##### banbura step 2
  
  # estimate VAR
  VAR_variable_sets <- var_set_use$var_set[var_set_use$estimate_var]
  var_list <- create_var_list(T_common = T_common, p_max = p_max, var_sets = VAR_variable_sets)
  message("Estimating VAR")
  var_list <- estimate_models(var_list, parallel = parallel, verbose = verbose,
                              var_set_info = var_set_info, df = df, 
                              deltas = deltas, num_AR_lags = num_AR_lags,
                              carriero_hack = carriero_hack,
                              v_prior = v_prior, keep = keep)
  
  # forecast VAR
  var_forecast_list <- data_frame(model_id = var_list$id, h = NA, type = "in-sample")
  message("Forecasting VAR")
  var_forecasts <- forecast_models(var_forecast_list, var_list, verbose = verbose,
                                   var_set_info = var_set_info, df = df,
                                   fast_forecast = fast_forecast)
  
  
  # calculate msfe-inf
  msfe_Inf_info <- get_msfe(var_forecasts, actual_obs, 
                            models = select(var_list, id, type, var_set, n_lag), 
                            msfe_name = "msfe_Inf")
  
  

  # join msfe0
  msfe_0_Inf <- left_join(msfe_Inf_info, msfe0 %>% select(msfe, rw_wn, variable) %>% 
                            rename(msfe0 = msfe), by = "variable")
  
  msfe_0_Inf <- msfe_0_Inf %>% select(-model_id, -type) %>% mutate(msfe_ratio = msfe_Inf/msfe0)
  msfe_0_Inf
  
  
  # create empty table
  fit_inf_table <- NULL
  
  # cycle all fit_sets:
  for (current_fit_set in unique(fit_set_info$fit_set)) {
    fit_variables <- filter(fit_set_info, fit_set == current_fit_set)$variable
    block <- filter(msfe_0_Inf, variable %in% fit_variables) %>% 
      group_by(var_set, n_lag) %>% summarise(fit_inf = mean(msfe_ratio)) %>% 
      mutate(fit_set = current_fit_set)
    fit_inf_table <- bind_rows(fit_inf_table, block)
  }
  
  fit_inf_table
  
  ##### banbura step 3 goal: calculate fit-lam
  
  # create model list to find optimal lambda
  BVAR_variable_sets <- var_set_use$var_set[var_set_use$estimate_bvar]
  bvar_list_ <- create_bvar_banbura_list(T_common = T_common, p_max = p_max, 
                                         var_sets = BVAR_variable_sets)
  # write_csv(bvar_list, path = '../estimation/bvar_list.csv')
  
  # estimate models from list bvar_list <- read_csv('../estimation/bvar_list.csv')
  message("Estimating BVAR")
  bvar_list <- estimate_models(bvar_list_, parallel = parallel, ncpu = ncpu, verbose = verbose,
                               var_set_info = var_set_info, df = df, 
                               deltas = deltas, num_AR_lags = num_AR_lags,
                               carriero_hack = carriero_hack,
                               v_prior = v_prior, keep = keep)  # status and filename are updated
  # write_csv(bvar_list, path = '../estimation/bvar_list.csv')
  
  # forecast BVAR
  bvar_forecast_list <- data_frame(model_id = bvar_list$id, h = NA, type = "in-sample")
  message("Forecasting BVAR in-sample")
  bvar_forecasts <- forecast_models(bvar_forecast_list, bvar_list, verbose = verbose,
                                    var_set_info = var_set_info, df = df,
                                    fast_forecast = fast_forecast)
  message("Forecasting BVAR in-sample ok")
  
  
  # calculate msfe-lam
  msfe_lam_info <- get_msfe(bvar_forecasts, actual_obs, 
                            models = select(bvar_list, id, type, var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, n_lag), 
                            msfe_name = "msfe_lam")
  
  
  msfe_lam_info %>% head()
  
  # join msfe0
  msfe_0_lam <- left_join(msfe_lam_info, msfe0 %>% select(msfe, rw_wn, variable) %>% 
                            rename(msfe0 = msfe), by = "variable")
  
  msfe_0_lam <- msfe_0_lam %>% select(-model_id, -type) %>% mutate(msfe_ratio = msfe_lam/msfe0)
  msfe_0_lam
  
  # calculate fit-lam create empty table
  fit_lam_table <- NULL
  
  # cycle all fit_sets:
  for (current_fit_set in unique(fit_set_info$fit_set)) {
    fit_variables <- filter(fit_set_info, fit_set == current_fit_set)$variable
    block <- filter(msfe_0_lam, variable %in% fit_variables) %>% 
      group_by(var_set, n_lag, l_1, l_const, l_io, l_power, l_sc) %>% 
      summarise(fit_lam = mean(msfe_ratio)) %>% 
      mutate(fit_set = current_fit_set)
    fit_lam_table <- bind_rows(fit_lam_table, block)
  }
  
  fit_lam_table
  
  ##### banbura step 4 find optimal lambda ungroup() is needed! otherwise cannot remove
  ##### var_set (active group on var_set)
  target_var_set <- dplyr::filter(var_set_use, target == TRUE)[["var_set"]]
  fit_goal <- dplyr::filter(fit_inf_table, var_set == target_var_set) %>% ungroup() %>% 
    dplyr::select(-var_set)
  fit_goal
  
  fit_lam_table <- left_join(fit_lam_table, fit_goal, by = c("n_lag", "fit_set"))
  fit_lam_table <- mutate(fit_lam_table, delta_fit = abs(fit_lam - fit_inf))
  
  # for each var_set, n_lag, fit_set the best lambdas are calculated
  optimal_by <- c("var_set", "n_lag", "fit_set")
  best_lambda <- fit_lam_table %>% group_by_(.dots = optimal_by) %>% 
    mutate(fit_rank = dense_rank(delta_fit)) %>% 
    filter(fit_rank == 1) %>% ungroup()
  best_lambda %>% arrange(var_set, n_lag, fit_set)
  
  
  
  # check whether best lambda is unique
  check_uniqueness <- best_lambda %>% group_by_(.dots = optimal_by) %>% 
    summarise(num_of_best_lambdas = n())
  # num of best lambdas should be always one
  if (max(check_uniqueness$num_of_best_lambdas) > 1) {
    warning("**** ACHTUNG ****: non unique optimal lambdas")
    check_uniqueness %>% filter(num_of_best_lambdas > 1)
    message("First in list optimal lambda will be chosen")
    best_lambda <- best_lambda %>% group_by_(.dots = optimal_by) %>% 
      mutate(rownum = row_number()) %>% 
      filter(rownum == 1) %>% select(-rownum)
    
  }
  
  
  ##### banbura step 5: calculate omsfe for bvar forecast and evaluate using optimal
  ##### lambda
  
  # create best models lists with correct time spec we need to keep fit_set
  # variable for further comparison
  bvar_out_list <- create_bvar_out_list(best_lambda, T_available = T_available,
                                        T_common = T_common, p_max = p_max)
  # best_lambda table should have: var_set, n_lag, l_1, l_const, l_io, l_power,
  # l_sc, fit_set
  
  message("Estimate rolling BVAR")
  bvar_out_list <- estimate_models(bvar_out_list, parallel = parallel, verbose = verbose,
                                   var_set_info = var_set_info, df = df, 
                                   deltas = deltas, num_AR_lags = num_AR_lags,
                                   carriero_hack = carriero_hack,
                                   v_prior = v_prior, keep = keep)
  # write_csv(bvar_out_list, path = '../estimation/bvar_out_list.csv')
  
  # forecast BVAR check structure
  bvar_out_list %>% group_by(var_set, fit_set, n_lag) %>% summarise(n = n())
  
  bvar_out_forecast_list <- bvar_out_list %>% rowwise() %>% 
    mutate(model_id = id, h = min(h_max,T_available - T_start - T_in + 1)) %>% 
    select(model_id, h) %>% mutate(type = "out-of-sample")
  # without rowwise min will be global and always equal to 1
  
  # a lot of forecasts WARNING: maybe too many observations!!!
  message("Forecasting rolling BVAR, out-of-sample")
  bvar_out_forecasts <- forecast_models(bvar_out_forecast_list, bvar_out_list, verbose = verbose,
                                        var_set_info = var_set_info, df = df,
                                        fast_forecast = fast_forecast)
  message("Forecasting rolling BVAR out-of-sample ok")
  
  omsfe_bvar_table <- get_msfe(bvar_out_forecasts, actual_obs, 
                               models = select(bvar_out_list, id, var_set, n_lag, fit_set), 
                               plus_group_vars = "fit_set", msfe_name = "omsfe", 
                               msfe_type = "out-of-sample")
  omsfe_bvar_table %>% head()
  
  ##### banbura step 6: calculate omsfe for RW/WN/VAR look at lists:
  var_list
  rwwn_list
  
  
  rwwn_var_unique_list <- bind_rows(var_list, rwwn_list)  # %>% mutate_each('as.numeric', n_lag, T_in, T_start)
  
  # every model should be rolled
  rwwn_var_out_list <- rolling_model_replicate(rwwn_var_unique_list, T_available = T_available,
                                               T_common = T_common, p_max = p_max) %>% 
    mutate(file = paste0(type, "_", id, "_T_", T_start, "_", T_in, "_", var_set, "_lags_", n_lag, ".Rds"))
  message("Estimating rolling rw/wn/var models")
  rwwn_var_out_list <- estimate_models(rwwn_var_out_list, parallel = parallel, verbose = verbose,
                                       var_set_info = var_set_info, df = df, 
                                       deltas = deltas, num_AR_lags = num_AR_lags,
                                       carriero_hack = carriero_hack,
                                       v_prior = v_prior, keep = keep)  # takes some minutes
  
  rwwn_var_out_forecast_list <- rwwn_var_out_list %>% rowwise() %>% 
    mutate(model_id = id, h = min(h_max, T_available - T_start - T_in + 1)) %>% 
    select(model_id, h) %>% 
    mutate(type = "out-of-sample")
  # without rowwise min will be global and always equal to 1
  
  # forecast all rolling models
  message("Forecasting rolling rw/wn/var models, out-of-sample")
  rwwn_var_out_forecasts <- forecast_models(rwwn_var_out_forecast_list, rwwn_var_out_list, 
                                            verbose = verbose, 
                                            var_set_info = var_set_info, df = df,
                                            fast_forecast = fast_forecast)
  
  
  
  omsfe_rwwn_var_table <- get_msfe(forecasts = rwwn_var_out_forecasts, actuals = actual_obs, 
                                   models = select(rwwn_var_out_list, id, var_set, n_lag, model_type = type), 
                                   plus_group_vars = "model_type", 
                                   msfe_name = "omsfe", msfe_type = "out-of-sample")
  
  
  ##### banbura step 7: calculate relative omsfe
  
  # we need to join: omsfe_bvar_table, omsfe_rwwn_var_table, deltas
  omsfe_rwwn_var_table
  deltas
  omsfe_bvar_table
  
  # select rw or wn omsfe for each variable, for each h
  omsfe_selected_rwwn <- left_join(deltas, filter(omsfe_rwwn_var_table, 
                                                  model_type %in% c("wn", "rw")), 
                                   by = c(rw_wn = "model_type", variable = "variable")) %>% 
    select(-var_set, -n_lag, -delta)
  omsfe_selected_rwwn %>% glimpse()
  
  
  omsfe_bvar_table <- ungroup(omsfe_bvar_table) %>% mutate_each("as.numeric", n_lag)
  
  
  # replicate banbura table from page 79
  
  omsfe_var_banbura <- ungroup(omsfe_rwwn_var_table) %>% filter(model_type == "var")
  
  omsfe_bvar_banbura <- omsfe_bvar_table %>% mutate(model_type = "bvar")
  
  omsfe_rwwn_banbura <- omsfe_selected_rwwn %>% rename(model_type = rw_wn)
  
  omsfe_var_banbura
  omsfe_bvar_banbura
  omsfe_rwwn_banbura
  
  rmsfe_long <- bind_rows(omsfe_var_banbura, omsfe_bvar_banbura) %>% 
    left_join(omsfe_rwwn_banbura %>% select(-model_type, omsfe_rwwn = omsfe), 
              by = c("h", "variable")) %>% mutate(rmsfe = omsfe/omsfe_rwwn)
  
  # rmsfe_long
  
  
  rmsfe_wide_bvar <- rmsfe_long %>% filter(model_type == "bvar") %>% 
    select(-omsfe, -omsfe_rwwn) %>% 
    dcast(h + variable + n_lag + fit_set ~ var_set + model_type, value.var = "rmsfe")
  
  rmsfe_wide_var <- rmsfe_long %>% filter(model_type == "var") %>% select(-omsfe, -omsfe_rwwn) %>% 
    dcast(h + variable + n_lag ~ var_set + model_type, value.var = "rmsfe")
  
  rmsfe_wide <- left_join(rmsfe_wide_bvar, rmsfe_wide_var, by = c("h", "variable", 
                                                                  "n_lag"))
  
  all_forecasts <- list(actual_obs = actual_obs,
                        rwwn_forecasts = rwwn_forecasts, 
                        var_forecasts = var_forecasts,
                        bvar_forecasts = bvar_forecasts, # bvar in-sample
                        bvar_out_forecasts = bvar_out_forecasts,
                        rwwn_var_out_forecasts = rwwn_var_out_forecasts)
  
  all_model_info <- list(fit_set_info = fit_set_info,
                         var_set_info = var_set_info,
                         deltas = deltas,
                         rwwn_list = rwwn_list,
                         rwwn_forecast_list = rwwn_forecast_list,
                         var_list = var_list,
                         var_forecast_list = var_forecast_list,
                         bvar_list = bvar_list,
                         bvar_list_ = bvar_list_,
                         bvar_forecast_list = bvar_forecast_list,
                         fit_lam_table = fit_lam_table,
                         bvar_out_list = bvar_out_list,
                         bvar_out_forecast_list = bvar_out_forecast_list,
                         rwwn_var_out_list = rwwn_var_out_list,
                         rwwn_var_out_forecast_list = rwwn_var_out_forecast_list)
  
  banbura_res <- list(rmsfe_wide = rmsfe_wide, forecasts = all_forecasts, model_info = all_model_info)
  
  return(banbura_res)
}