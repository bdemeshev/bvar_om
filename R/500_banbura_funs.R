# 500_banbura_funs.R

library("BMR")

# this function set delta=1 for non-stationary variables and delta=0 for stationary
# for ADF all three types are supported
# for KPSS only trend and constant
# c_0=0 value for stationary series
# c_1=1 value for non-stationary series
# alpha=5%
delta_i_prior <- function(df_ts, varset=colnames(df_ts), test=c("ADF","KPSS"),
                          type=c("trend", "constant", "neither"), c_0=0, c_1=1, remove_vars=NULL) {
  varset <- setdiff(varset, remove_vars)
  
  df_sel <- df_ts[,varset]
  nvar <- ncol(df_sel)
  stat_res <- stationarity(df_sel, print = FALSE)
  
  type <- match.arg(type) # take first option if not specified
  test <- match.arg(test)
  
  if (type=="trend") line <- 1
  if (type=="constant") line <- 2
  if (type=="neither") line <- 3
  
  nvar <- ncol(df_sel)
  
  # use time_trend - ADF - 5%
  if(test=="ADF") {
    adf_cr <- stat_res$ADF$`5 Pct`[line]
    adf_obs <- stat_res$ADF[line,1:nvar]  
    delta <- ifelse(adf_obs < adf_cr, c_0, c_1)
  }
  
  if(test=="KPSS") {
    kpss_cr <- stat_res$KPSS$`5 Pct`[line]
    kpss_obs <- stat_res$KPSS[line,1:nvar]  
    delta <- ifelse(kpss_obs > kpss_cr, c_1, c_0)
  }
  
  answer <- data.frame(variable=varset, delta=as.vector(delta))
  
  return(answer)
}

# this function calculates AR(1) coefficient for each variable
delta_i_from_ar1 <- function(df_ts, varset=colnames(df_ts), remove_vars=NULL, replace_by_one=TRUE) {
  varset <- setdiff(varset, remove_vars)
  
  answer <- data.frame(variable=varset, delta=0)
  for (j in varset) {
    y_uni <- df_ts[,j]
    AR_1 <- ar.ols(y_uni, aic=FALSE, order.max = 1)
    answer$delta[answer$variable==j] <- AR_1$ar
  }
  if (max(answer$delta)>1) {
    warning("Some AR(1) coefficients estimates are greater than 1")
    if (replace_by_one) answer$delta[answer$delta>1] <- 1
  }
  return(answer)
}


get_msfe <- function(forecasts, actuals, 
                      plus_group_vars = NULL, 
                      msfe_name="msfe" ) {
  
  # joining actual observations
  for_p_act <- left_join(forecasts, actuals, by=c("t","variable"))
  
  for_p_act <- mutate(for_p_act, sq_error=(forecast-actual)^2)
  
  # calculate msfe

  msfes <- for_p_act %>% group_by_(.dots=c("variable", "model_id", plus_group_vars)) %>% 
                                    summarise(msfe=mean(sq_error))
  
  colnames(msfes)[colnames(msfes)=="msfe"] <- msfe_name
  
  return(msfes)
}


