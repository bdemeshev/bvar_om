#' Dummy function to estimate random walk model
#' 
#' \code{estimate_rw} estimates nothing at all :)
#' 
#' This function just provides common syntax interface and should be never used.
#' 
#' @param y multivariate time series
#' @param ... dummy for further arguments 
#' @return message
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' junk <- estimate_rw(y_small)
estimate_rw <- function(y, ...) {
  return("RW model. Does not need estimation.")
}

#' Estimate auto arima model on multivariate time series
#' 
#' \code{estimate_arima} returns a list of estimated arima model
#' 
#' Just a wrapper for \code{auto.arima} function from \code{forecast} package.
#' 
#' @param y multivariate time series
#' @param ... further arguments passed to \code{auto.arima} function
#' @return list of estimated arima models
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' arima_model <- estimate_arima(y_small)
estimate_arima <- function(y, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)

  model <- list()

  for (i in 1:m) {
    model[[i]] <- forecast::auto.arima(y[, i], ...)
  }
  return(model)
}

#' Estimate VAR model on multivariate time series
#' 
#' \code{estimate_var} returns estimated VAR model
#' 
#' Simply a wrapper for \code{VAR} function from \code{vars} package.
#' 
#' @param y multivariate time series
#' @param p number of lags for VAR model
#' @param ... further arguments passed to \code{VAR} function
#' @return estimated VAR model
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' var_model <- estimate_var(y_small, p = 1)
estimate_var <- function(y, p = 1, ...) {
  y_matrix <- as.matrix(y)
  model <- vars::VAR(y_matrix, p = p, ...)
  return(model)
}



#' Estimate ETS model on multivariate time series
#' 
#' \code{estimate_ets} returns a list of estimated ETS models
#' 
#' Just a wrapper for \code{auto.arima} function from \code{forecast} package.
#' 
#' @param y multivariate time series
#' @param ... further arguments passed to \code{ets} function
#' @return list of estimated ETS models
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' ets_model <- estimate_ets(y_small)
estimate_ets <- function(y, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model <- list()
  
  for (i in 1:m) {
    model[[i]] <- forecast::ets(y[, i], ...)
  }
  return(model)
}


#' Estimate VAR-lasso model on multivariate time series
#' 
#' \code{estimate_var_lasso} returns an estimated VAR-lasso model
#' 
#' Just a wrapper for \code{cv.BigVAR} function from \code{BigVAR} package. This function does NOT normalise variables. 
#' Pre-normalisation is required by the idea.
#' 
#' @param y multivariate time series
#' @param h forecasting horizon
#' @param p number of lags
#' @param struct type of VAR-lasso
#' @param gran granularity vector, If gran = (a, b) then cross-validation checks b values form lambda to lambda/a.
#' @param T1 left border for cross validation samples
#' @param T2 right border for cross validation samples
#' @return estimated VAR-lasso model
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' y_small_scaled <- scale_to(y_small)
#' var_lasso_model <- estimate_var_lasso(y_small_scaled)
estimate_var_lasso <- function(y, h = 1, p = 12,
                               struct = "OwnOther",
                               gran = c(25, 10),
                               T1 = floor(nrow(as.matrix(y)) / 3),
                               T2 = floor(2 * nrow(as.matrix(y)) / 3)) {
  y_matrix <- as.matrix(y)
  
  # strange error for time series in BigVAR ?
  if (stats::is.ts(y_matrix)) {
    y_matrix <- zoo::coredata(y_matrix)
  }
  
  cat("Estimation of var lasso ", struct, " started...\n")
  model_spec <- BigVAR::constructModel(y_matrix, p = p, struct = struct, 
                                       T1 = T1, T2 = T2, 
                                       gran = gran, h = h, 
                                       verbose = TRUE, VARX = list())
  
  model <- BigVAR::cv.BigVAR(model_spec)
  cat("Estimation of var lasso ", struct, " done.\n")
  return(model)
}


#' Estimate TVP model on multivariate time series
#' 
#' \code{estimate_tvp_primiceri} returns an estimated TVP model a la Primiceri
#' 
#' Just a wrapper for \code{bvar.sv.tvp} function from \code{bvarsv} package.
#' ACHTUNG: probably we need to scale series for this model!
#' 
#' @param y multivariate time series
#' @param p number of lags
#' @return estimated TVP model
#' @param nrep number of MCMC draws excluding burn-in
#' @param nburn number of burn-in MCMC draws 
#' @param ... furhter arguments passed to \code{bvar.sv.tvp}
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' y_small_scaled <- scale_to(y_small)
#' primiceri_model <- estimate_tvp_primiceri(y_small_scaled, p = 1)
estimate_tvp_primiceri <- function(y, p = 12, nrep = 1000, nburn = 1000, ...) {
  y_matrix <- as.matrix(y)
  
  model <- bvarsv::bvar.sv.tvp(y_matrix, p = p, nrep = nrep, nburn = nburn, ...)
  return(model)
}
