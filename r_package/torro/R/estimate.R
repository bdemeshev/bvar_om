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
