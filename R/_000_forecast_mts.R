library(forecast)

estimate_arima <- function(y, h = 1, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model <- list()
  
  for (i in 1:m) {
    model[[i]] <- auto.arima(y[, i], ...)
  }
  return(model)
}

estimate_ets <- function(y, h = 1, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model <- list()
  
  for (i in 1:m) {
    model[[i]] <- ets(y[, i], ...)
  }
  return(model)
}

# we return something like mforecast (!)
forecast_arima <- function(y, h = 1, 
                           model = NULL, 
                           ...) {
  if (is.null(model)) {
    model <- estimate_arima(y, h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- forecast(model[[i]], h = h)
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  return(forecast_data)
}

forecast_ets <- function(y, h = 1, 
                           model = NULL, 
                           ...) {
  if (is.null(model)) {
    model <- estimate_ets(y, h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- forecast(model[[i]], h = h)
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  return(forecast_data)
}

scale_series <- function(y) {
  y_matrix <- as.matrix(y)
  y_scaled <- apply(y_matrix, 2, FUN = scale)
  return(y_scaled)
}


# package forecast defines 
# mforecast class for multivariate forecasts
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

matrix_to_mforecast <- function(forecast_matrix) {
  
  mforecast <- list()
  mforecast$forecast <- list()
  m <- ncol(forecast_matrix)
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- list()
    mforecast$forecast[[i]]$mean <- forecast_matrix[, i]
  }
  
  names(mforecast$forecast) <- colnames(forecast_matrix)
  return(mforecast)
}


# test block

library(vars)
data("Canada")
y <- scale_series(Canada)
fors <- forecast_arima(y, h = 3)
mods <- estimate_ets(y, h = 3)
fors <- forecast_ets(y, h = 3, model = mods)

m <- auto.arima(y)
