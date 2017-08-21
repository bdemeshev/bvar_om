#' Get time of the next to last observation of time series
#' 
#' \code{next_obs_time} gives the time of the next to last observation of time series
#' 
#' Let's consider monthly time series. In the output 2015.0 means 2015 first month, 
#' 2015.083 means 2015 second month and so on.
#' 
#' @param y time series 
#' @return time of the last observation plus one time unit
#' @export
#' @examples 
#' data(rus_macro)
#' next_obs_time(rus_macro)
next_obs_time <- function(y) {
  return(stats::end(y)[1] + stats::deltat(y) * stats::end(y)[2])
  # !!! don't need to add one 
}

#' Transform mforecast object to simple matrix of forecasts.
#' 
#' \code{mforecast_to_matrix} transform complex mforecast object to simple matrix of forecasts.
#' 
#' \code{forecast} package uses complex mforecast format to store forecasts.
#' This function transforms mforecast objects to simple matrices.
#' 
#' @param mforecast mforecast object to be transformed
#' @return h by m matrix (h — number of steps, m — number of variables)
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' var_model <- vars::VAR(y_small)
#' y_for <- forecast::forecast(var_model, h = 2)
#' mforecast_to_matrix(y_for)
mforecast_to_matrix <- function(mforecast) {
  m <- length(mforecast$forecast)
  h <- length(mforecast$forecast[[1]]$mean)
  forecast_matrix <- matrix(0, nrow = h, ncol = m)
  for (i in 1:m) {
    forecast_matrix[, i] <- mforecast$forecast[[i]]$mean
  }
  colnames(forecast_matrix) <- names(mforecast$forecast)
  return(forecast_matrix)
}

#' Transform simple matrix of forecasts to mforecast object.
#' 
#' \code{matrix_to_mforecast} transform simple matrix of forecasts to complex mforecast object.
#' 
#' \code{forecast} package uses complex mforecast format to store forecasts.
#' This function transforms simple matrices to mforecast objects.
#'  
#' @param forecast_matrix h by m forecast matrix
#' @param y_before tibble with three columns: variable, mean, sd
#' @param method the name of the method
#' @return complex mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' y_forecast <- tail(y_small, 1) # trivial forecast 
#' mfor <- matrix_to_mforecast(y_forecast, y_small, method = "Trivial")
#' # ggplot2:::autoplot.mforecast(mfor) # TODO: does not work!!!
matrix_to_mforecast <- function(forecast_matrix, y_before, 
                                method = "Unspecified") {
  
  mforecast <- list()
  mforecast$forecast <- list()
  m <- ncol(forecast_matrix) # m = number of times series
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- list()
    
    # frequency of plain matrices is equal to one
    fors_freq <- stats::frequency(y_before[, i]) 
    fors_start <- next_obs_time(y_before[, i])    
    
    mforecast$forecast[[i]]$method <- method # method name
    mforecast$forecast[[i]]$x <- y_before[, i] # actual y before forecast period
    
    # forecasts with correct frequency and start:
    future_ts <- stats::ts(forecast_matrix[, i], start = fors_start, frequency = fors_freq)
    try_na_remove <- try(stats::na.omit(future_ts))
    if (class(try_na_remove) == "try-error") {
      mforecast$forecast[[i]]$mean <- future_ts
    } else {
      mforecast$forecast[[i]]$mean <- try_na_remove
    }
    
    mforecast$forecast[[i]]$series <- colnames(forecast_matrix)[i] # names of series
    class(mforecast$forecast[[i]]) <- "forecast"
  }
  
  names(mforecast$forecast) <- colnames(forecast_matrix)
  class(mforecast) <- "mforecast"
  return(mforecast)
}
