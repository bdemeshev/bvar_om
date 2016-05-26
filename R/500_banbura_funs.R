# 500_banbura_funs.R

library("BMR")

# this function set delta=c_1 for non-stationary variables and
# delta=c_1 for stationary for ADF all three types are supported for
# KPSS only trend and constant c_0=0 value for stationary series c_1=1
# value for non-stationary series alpha=5%
delta_i_prior <- function(df_ts, varset = colnames(df_ts), test = c("ADF", 
                                                                    "KPSS"), type = c("trend", "constant", "neither"), c_0 = 0, c_1 = 1, 
                          remove_vars = NULL) {
  varset <- setdiff(varset, remove_vars)
  
  df_sel <- df_ts[, varset]
  nvar <- ncol(df_sel)
  stat_res <- stationarity(df_sel, print = FALSE)
  
  type <- match.arg(type)  # take first option if not specified
  test <- match.arg(test)
  
  if (type == "trend") {
    line <- 1
  }
  if (type == "constant") {
    line <- 2
  }
  if (type == "neither") {
    line <- 3
  }
  
  nvar <- ncol(df_sel)
  
  # use time_trend - ADF - 5%
  if (test == "ADF") {
    adf_cr <- stat_res$ADF$`5 Pct`[line]
    adf_obs <- stat_res$ADF[line, 1:nvar]
    delta <- ifelse(adf_obs < adf_cr, c_0, c_1)
  }
  
  if (test == "KPSS") {
    kpss_cr <- stat_res$KPSS$`5 Pct`[line]
    kpss_obs <- stat_res$KPSS[line, 1:nvar]
    delta <- ifelse(kpss_obs > kpss_cr, c_1, c_0)
  }
  
  answer <- data.frame(variable = varset, delta = as.vector(delta))
  
  return(answer)
}

# this function calculates AR(1) coefficient for each variable
delta_i_from_ar1 <- function(df_ts, varset = colnames(df_ts), remove_vars = NULL, 
                             replace_by_one = TRUE) {
  varset <- setdiff(varset, remove_vars)
  
  answer <- data.frame(variable = varset, delta = 0)
  for (j in varset) {
    y_uni <- df_ts[, j]
    AR_1 <- ar.ols(y_uni, aic = FALSE, order.max = 1)
    answer$delta[answer$variable == j] <- AR_1$ar
  }
  if (max(answer$delta) > 1) {
    warning("Some AR(1) coefficients estimates are greater than 1")
    if (replace_by_one) 
      answer$delta[answer$delta > 1] <- 1
  }
  return(answer)
}

#' @param models list of models, all field from models will be joined
get_msfe <- function(forecasts, actuals, models = NULL, plus_group_vars = NULL, 
                     msfe_name = "msfe", msfe_type = c("in-sample", "out-of-sample")) {
  
  msfe_type <- match.arg(msfe_type)
  
  # joining forecast and actual observations
  for_and_act <- left_join(forecasts, actuals, by = c("t", "variable"))
  for_and_act <- mutate(for_and_act, sq_error = (forecast - actual)^2)
  
  # join info from model list
  
  
  # calculate msfe
  
  if (msfe_type == "in-sample") {
    # for in-sample forecast errors are averaged for each model
    msfes <- for_and_act %>% group_by_(.dots = c("variable", "model_id", 
                                                 plus_group_vars)) %>% summarise(msfe = mean(sq_error))
    if (!is.null(models)) {
      msfes <- left_join(msfes, models, by = c(model_id = "id"))
    }
  }
  if (msfe_type == "out-of-sample") {
    if (!is.null(models)) {
      for_and_act <- left_join(for_and_act, models, by = c(model_id = "id"))
    }
    # for out-of-sample forecasts errors are averaged between models of the
    # same type evaluated for moving time slot
    msfes <- for_and_act %>% group_by_(.dots = c("variable", "var_set", 
                                                 "n_lag", "h", plus_group_vars)) %>% summarise(msfe = mean(sq_error))
  }
  
  
  
  colnames(msfes)[colnames(msfes) == "msfe"] <- msfe_name
  
  return(msfes)
}

#' @param bvar_list
#' @return bvar_list augmented with mdd variable
calculate_mdd <- function(bvar_list, parallel = c("off", "windows", "unix"), 
                          ncpu = 4, test = FALSE, do_log = FALSE, progress_bar = TRUE, verbose = FALSE) {
  start_time <- Sys.time()
  
  parallel <- match.arg(parallel)
  
  if (parallel == "off") {
    bvar_list$mdd <- NA
    
    if (progress_bar) {
      pb <- txtProgressBar(min = 1, max = nrow(bvar_list), style = 3)
    }
    for (i in 1:nrow(bvar_list)) {
      file_wpath <- paste0("../estimation/models/", bvar_list$file[i])
      model <- readRDS(file_wpath)
      if (verbose) {
        message(bvar_list$file[i])
      }
      bvar_list$mdd[i] <- bvar_conj_mdd(model)
      if (progress_bar) {
        setTxtProgressBar(pb, i)
      }
    }
    if (progress_bar) {
      close(pb)
    }
    
  }
  end_time <- Sys.time()
  message("Time elapsed: ", round(end_time - start_time, 1), " ", attr(end_time - 
                                                                         start_time, "units"))
  
  return(bvar_list)
}
