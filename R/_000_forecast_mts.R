library(forecast)
library(vars)
library(tidyverse)

estimate_rw <- function(y, h = 1, ...) {
  return("RW model. Does not need estimation.")
}

estimate_arima <- function(y, h = 1, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model <- list()
  
  for (i in 1:m) {
    model[[i]] <- forecast::auto.arima(y[, i], ...)
  }
  return(model)
}



next_obs_time <- function(y) {
  return(end(y)[1] + deltat(y) * (end(y)[2] + 1))
}

estimate_ets <- function(y, h = 1, ...) {
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  model <- list()
  
  for (i in 1:m) {
    model[[i]] <- forecast::ets(y[, i], ...)
  }
  return(model)
}


estimate_var_lasso <- function(y, h = 1, p = 12,
                               struct = "OwnOther",
                               gran = c(10000, 3),
                               T1 = floor(nrow(as.matrix(y)) / 3),
                               T2 = floor(2 * nrow(as.matrix(y)) / 3)) {
  y_matrix <- as.matrix(y)
  
  # strange error for time series in BigVAR ?
  if (is.ts(y_matrix)) {
    y_matrix <- coredata(y_matrix)
  }
  
  model_spec <- BigVAR::constructModel(y_matrix, p = p, struct = struct, 
                           T1 = T1, T2 = T2, 
                           gran = gran, h = h, 
                           verbose = TRUE, VARX = list())
  
  model <- BigVAR::cv.BigVAR(model_spec)
  return(model)
}


estimate_tvp_primiceri <- function(y, h = 1, p = 12, nrep = 1000, nburn = 1000, ...) {
  y_matrix <- as.matrix(y)
  
  model <- bvarsv::bvar.sv.tvp(y_matrix, p = p, nrep = nrep, nburn = nburn, ...)
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
    forecast_data$forecast[[i]] <- forecast::rwf(y_matrix[, i], h = h, ...)
    forecast_data$forecast[[i]]$series <- colnames(y_matrix)[i]
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  class(forecast_data) <- "mforecast"
  return(forecast_data)
}


# we return something like mforecast (!)
forecast_var_lasso <- function(y, h = 1, model = NULL, ...) {
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
  mforecast <- matrix_to_mforecast(forecast_matrix, y_matrix, method = "BigVAR")
  mforecast$model <- model
  return(mforecast)
}


primiceri_draws_to_1d_forecast <- function(model, y_actual, series_no = 1, 
                                           h = 1, level = c(80, 95)) {

  
  
create_interval_borders <- function(y_future, level, 
                                    type = c("lower", "upper"),
                                    fors_start, fors_freq) {
  type <- match.arg(type)
  
  lower_probs <- (1 - level / 100) / 2
  if (type == "lower") {
    probs <- lower_probs
  } else {
    probs <- 1 - lower_probs
  }
  
  y_future_border <- map(y_future, ~ quantile(.x, probs = lower_probs)) %>% as_tibble() %>% t()
  y_future_border <- ts(y_future_border, start = fors_start, frequency = fors_freq)
  colnames(y_future_border) <- paste0(level, "%")
  
  return(y_future_border)
}  
  
  
  one_ts_forecast <- list()
  one_ts_forecast$method <- "TVP-BVAR_Primiceri"
  one_ts_forecast$level <- level
  one_ts_forecast$series <- colnames(y_actual)[series_no]
  one_ts_forecast$x <- y_actual[, series_no]
  
  y_future <- list()
  for (h_index in 1:h) {
    y_future[[h_index]] <- bvarsv::predictive.draws(model, v = series_no, h = h_index)$y
  }
  names(y_future) <- paste0("h", 1:h)
  y_future <- as_tibble(y_future)
  
  fors_freq <- frequency(y_actual[, series_no])
  fors_start <- next_obs_time(y_actual[, series_no])
  
  y_future_mean <- map_dbl(y_future, mean)
  y_future_mean <- ts(y_future_mean, start = fors_start, frequency = fors_freq)
  
  
  one_ts_forecast$mean <- y_future_mean
  one_ts_forecast$lower <- create_interval_borders(y_future, 
                                                   level = level, type = "lower",
                                                   fors_start = fors_start,
                                                   fors_freq = fors_freq)
  one_ts_forecast$upper <- create_interval_borders(y_future, 
                                                   level = level, type = "upper",
                                                   fors_start = fors_start,
                                                   fors_freq = fors_freq)
  
  class(one_ts_forecast) <- "forecast"
  return(one_ts_forecast)
}


# we return something like mforecast (!)
forecast_tvp_primiceri <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_tvp_primiceri(y, h = h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix) # number of series
  
  mforecast <- list()
  mforecast$model <- model
  
  mforecast$forecast <- list()
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- primiceri_draws_to_1d_forecast(model, y_actual = y_matrix,
                                                              series_no = i, h = h)
  }

  class(mforecast) <- "mforecast"
  return(mforecast)
}




# we return something like mforecast (!)
forecast_arima <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_arima(y, h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- forecast::forecast(model[[i]], h = h)
    forecast_data$forecast[[i]]$series <- colnames(y_matrix)[i]
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  
  class(forecast_data) <- "mforecast"
  return(forecast_data)
}

forecast_ets <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_ets(y, h, ...)
  }
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_data <- list()
  forecast_data$forecast <- list()
  
  for (i in 1:m) {
    forecast_data$forecast[[i]] <- forecast::forecast(model[[i]], h = h)
    forecast_data$forecast[[i]]$series <- colnames(y_matrix)[i]
  }
  
  names(forecast_data$forecast) <- colnames(y_matrix)
  class(forecast_data) <- "mforecast"
  return(forecast_data)
}

scale_series <- function(y) {
  if (is.ts(y)) {
    y_freq <- frequency(y)
    y_start <- start(y)
  }
  y_matrix <- as.matrix(y)
  y_scaled <- apply(y_matrix, 2, FUN = scale)
  if (is.ts(y)) {
    y_scaled <- ts(y_scaled, frequency = y_freq, start = y_start)
  } else {
    y_scaled <- ts(y_scaled, frequency = 1, start = 1)
  }
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

matrix_to_mforecast <- function(forecast_matrix, y_before, method = "Unspecified") {
  
  mforecast <- list()
  mforecast$forecast <- list()
  m <- ncol(forecast_matrix) # m = number of times series
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- list()
    
    mforecast$forecast[[i]]$method <- method # method name
    mforecast$forecast[[i]]$x <- y_before[, i] # actual y before forecast period

    fors_freq <- frequency(y_before[, i])
    fors_start <- tsp(y_before[, i])[2] + (1 / fors_freq)
    # forecasts with correct frequency and start:
    mforecast$forecast[[i]]$mean <- ts(forecast_matrix[, i], start = fors_start, frequency = fors_freq)

    mforecast$forecast[[i]]$series <- colnames(forecast_matrix)[i] # names of series
    class(mforecast$forecast[[i]]) <- "forecast"
  }
  
  names(mforecast$forecast) <- colnames(forecast_matrix)
  class(mforecast) <- "mforecast"
  return(mforecast)
}

# load data
rus_macro <- readr::read_csv("../data/df_2015_final.csv")
# head(rus_macro$time_y)
rus_macro <- dplyr::select(rus_macro, -time_y)
rus_macro <- ts(rus_macro, start = c(1995, 1), frequency = 12)
head(rus_macro)
str(rus_macro)
# test block


y <- scale_series(rus_macro)
y_subset <- y[, 1:2]
fors <- forecast_arima(y_subset, h = 3)
fors <- forecast_ets(y_subset, h = 3)
fors <- forecast_rw(y_subset, h = 3)
fors <- forecast_var_lasso(y_subset, h = 3)
autoplot(fors)


mods <- estimate_tvp_primiceri(y_subset, h = 3, p = 1)
fors <- forecast_tvp_primiceri(y_subset, h = 3, model = mods)

library(listviewer)
one_ts_forecast <- fors$forecast[[1]]
jsonedit(one_ts_forecast)

library(vars)
model_VAR <- VAR(y_subset, p = 1)
fors_var <- forecast(model_VAR, h = 3)
jsonedit(fors_var)


