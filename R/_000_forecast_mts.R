# this file contains functions used to estimate arima, ets, rw, var_lasso

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

# p — number of lags
estimate_var <- function(y, h = 1, p = 1, ...) {
  y_matrix <- as.matrix(y)
  model <- VAR(y_matrix, p = p, ...)
  return(model)
}

forecast_var <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_var(y, h = h, ...)
  }
  
  mforecast <- forecast(model, h = h)
  return(mforecast)
}



# 2015y 1m ~ 2015.0
# 2015y 2m ~ 2015.1
next_obs_time <- function(y) {
  return(end(y)[1] + deltat(y) * end(y)[2])
  # don't need to add one (!!!)
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
                               gran = c(25, 10),
                               T1 = floor(nrow(as.matrix(y)) / 3),
                               T2 = floor(2 * nrow(as.matrix(y)) / 3)) {
  y_matrix <- as.matrix(y)
  
  # strange error for time series in BigVAR ?
  if (is.ts(y_matrix)) {
    y_matrix <- coredata(y_matrix)
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
forecast_var_lasso <- function(y, h = 1, model = NULL, p = 1, 
                               struct = "OwnOther", gran = c(25, 10),
                               type = c("fast", "honest"), h_cv = 1, ...) {
  # honest: one should redo cross-validation for each h
  # h should be a vector
  # model should be missing or a list of the same length as h (or one :))
  # fast: do cross-validation for h_cv
  # and than predict for each h = 1, ..., h
  type <- match.arg(type)
  
  # var_lasso fails! for time series!
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_matrix <- matrix(NA, nrow = max(h), ncol = m)    
  
  if (type == "fast") {
    if (is.null(model)) {
      model <- estimate_var_lasso(y, h = h_cv, p = p, struct = struct, gran = gran, ...)
    }

    for (i in 1:h) {
      cat("  forecast for h = ", i, " started...\n")
      forecast_matrix[i, ] <- 
        as.vector(BigVAR::predict(model, n.ahead = i))
      # as.vector нужен так как predict для n.ahead = 1 возвращает строку
      # а для n.ahead > 1 возвращает столбец
      cat(" forecast for h = ", i, " done.\n")
    }
  }
  
  if (type == "honest") {
    if (!length(model) %in% c(0, 1, length(h))) {
      stop("Argument 'model' should be NULL, one model for all 'h' or a list of models for each 'h'.\n
           Length of 'model' is ", length(model), ", while lenght of 'h' is ", length(h),".")
    }

    # estimate if no models are estimated
    if (length(model) == 0) {
      if (length(h) == 1) {
        model <- estimate_var_lasso(y, h = h, p = p, struct = struct, gran = gran, ...)
      } else {
        model <- list()
        for (i in 1:length(h)) {
          model[[i]] <- estimate_var_lasso(y, h = h[i], p = p, struct = struct, gran = gran, ...)
        }
      }
    }
    
    # forecast 
    for (i in 1:length(h)) {
      cat("  forecast for h = ", h[i], " started...\n")
      
      if (length(model) > 1) {
        # separate model for each h
        model_i <- model[[i]]
      } else {
        # same model for each h
        model_i <- model
      }
      forecast_matrix[h[i], ] <- 
        as.vector(BigVAR::predict(model_i, n.ahead = h[i]))
      # as.vector нужен так как predict для n.ahead = 1 возвращает строку
      # а для n.ahead > 1 возвращает столбец
      cat("  forecast for h = ", h[i], " done.\n")
    }
    
  }
  
  
  
  colnames(forecast_matrix) <- colnames(y_matrix)
  mforecast <- matrix_to_mforecast(forecast_matrix, y, method = "BigVAR")
  mforecast$model <- model
  return(mforecast)
}


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
  
  y_future_border <- map(y_future, ~ quantile(.x, probs = probs)) %>% as_tibble() %>% t()
  y_future_border <- ts(y_future_border, start = fors_start, frequency = fors_freq)
  colnames(y_future_border) <- paste0(level, "%")
  
  return(y_future_border)
}  



primiceri_draws_to_1d_forecast <- function(model, y_actual, series_no = 1, 
                                           h = 1, level = c(80, 95)) {
  
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
  y_matrix <- y # y_matrix <- as.matrix(y)
  m <- ncol(y_matrix) # number of series
  
  mforecast <- list()
  mforecast$model <- model
  
  mforecast$forecast <- list()
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- primiceri_draws_to_1d_forecast(model, y_actual = y_matrix,                                                              series_no = i, h = h)
  }
  names(mforecast$forecast) <- colnames(y_matrix)

  class(mforecast) <- "mforecast"
  return(mforecast)
}




# we return something like mforecast (!)
forecast_arima <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_arima(y, h, ...)
  }
  y_matrix <- y # y_matrix <- as.matrix(y)
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

    fors_freq <- frequency(y_before[, i]) # will get one for plain matrices
    fors_start <- next_obs_time(y_before[, i])    
        
    mforecast$forecast[[i]]$method <- method # method name
    mforecast$forecast[[i]]$x <- y_before[, i] # actual y before forecast period

    # forecasts with correct frequency and start:
    future_ts <- ts(forecast_matrix[, i], start = fors_start, frequency = fors_freq)
    try_na_remove <- try(na.omit(future_ts))
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


load_rus_data <- function() {
  rus_macro <- readr::read_csv("../data/df_2015_final.csv")
  rus_macro <- dplyr::select(rus_macro, -time_y)
  rus_macro <- ts(rus_macro, start = c(1995, 1), frequency = 12)
  return(rus_macro)
}

# just an example 
# 
# # load data
# rus_macro <- load_rus_data()
# # test block
# 
# 
# y <- scale_series(rus_macro)
# y_subset <- y[, 1:2]
# 
# fors <- forecast_arima(y_subset, h = 3)
# fors <- forecast_ets(y_subset, h = 3)
# fors <- forecast_rw(y_subset, h = 3)
# fors <- forecast_var_lasso(y_subset, h = 3)
# fors <- forecast_var_lasso(y_subset, h = 1, p = 12, struct = "OwnOther", gran = c(25, 10), type = "honest")
# fors <- forecast_var_lasso(y_subset, h = 12, p = 12, struct = "OwnOther", gran = c(25, 10), type = "fast")
# fors <- forecast_var_lasso(y_subset, h = 12, p = 12, struct = "OwnOther", gran = c(25, 10), type = "honest")
# fors <- forecast_tvp_primiceri(y_subset, h = 3, p = 1)
# fors <- forecast_var(y_subset, h = 3, p = 1)
# autoplot(fors)
# plot(fors$model)
# 
# 
# 
# 
# library(listviewer)
# one_ts_forecast <- fors$forecast[[1]]
# autoplot(one_ts_forecast)
# jsonedit(one_ts_forecast)
# 




