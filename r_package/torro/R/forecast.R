#' Forecast using VAR model
#' 
#' \code{forecast_var} forecasts time series using VAR model
#' 
#' Models are estimated using \code{VAR} function from \code{vars} package if necessary.
#' 
#' @param y multivariate time series
#' @param model estimated VAR model. If missing will be estimated automatically.
#' @param h forecast horizon
#' @param ... further arguments passed to \code{estimate_var} and than to \code{VAR} function
#' @return forecasts from VAR model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' var_forecast <- forecast_var(y_small, h = 2, p = 1)
forecast_var <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_var(y, ...)
  }
  
  mforecast <- forecast::forecast(model, h = h)
  return(mforecast)
}




#' Forecast using random walk model
#' 
#' \code{forecast_rw} forecasts time series using RW model
#' 
#' \code{rwf} function from \code{forecast} package is used
#' 
#' @param y multivariate time series
#' @param h forecast horizon
#' @param ... further arguments passed to \code{rwf} function
#' @return forecasts from RW model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' rw_forecast <- forecast_rw(y_small, h = 2)
forecast_rw <- function(y, h = 1, ...) {
  # drift = TRUE/FALSE (FALSE by default)
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- forecast::rwf(y_matrix[, i], h = h, ...)
    forecast_data$forecast[[i]]$series <- colnames(y_matrix)[i]
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  class(forecast_data) <- "mforecast"
  return(forecast_data)
}