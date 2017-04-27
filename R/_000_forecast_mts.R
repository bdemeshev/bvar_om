library(forecast)


estimate_rw <- function(y, h = 1, ...) {
  return("RW model. Does not need estimation.")
}

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


estimate_var_lasso <- function(y, h = 1, p = 12,
                               struct = "OwnOther",
                               gran = c(10000, 3),
                               T1 = floor(nrow(as.matrix(y)) / 3),
                               T2 = floor(2 * nrow(as.matrix(y)) / 3)) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model_spec <- constructModel(y_matrix, p = p, struct = struct, 
                           T1 = T1, T2 = T2, 
                           gran = gran, h = h, 
                           verbose = TRUE, VARX = list())
  
  model <- cv.BigVAR(model_spec)
  return(model)
}


# we return something like mforecast (!)
forecast_rw <- function(y, h = 1, ...) {
  # drift = TRUE/FALSE (FALSE by default)
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- rwf(y_matrix[, i], h = h, ...)
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  return(forecast_data)
}


# we return something like mforecast (!)
forecast_var_lasso <- function(y, h = 1, 
                           model = NULL, 
                           ...) {
  if (is.null(model)) {
    model <- estimate_var_lasso(y, h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_matrix <- matrix(0, nrow = h, ncol = m)
  
  for (i in 1:h) {
    forecast_matrix[i, ] <- 
      as.vector(predict(model, n.ahead = i))
    # as.vector нужен так как predict для n.ahead = 1 возвращает строку
    # а для n.ahead > 1 возвращает столбец
  }
  
  colnames(forecast_matrix) <- colnames(y_matrix)
  mforecast <- matrix_to_mforecast(forecast_matrix)
  return(mforecast)
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
forecast_
fors <- forecast_rw(y, h = 3)

m <- auto.arima(y)

mforecast_to_matrix(fors)

mods <- estimate_var_lasso(y, h = 3)

fors <- forecast_var_lasso(y, h = 3, model = mods)
fors
mforecast_to_matrix(fors)
