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


#' Forecast using ARIMA model
#' 
#' \code{forecast_arima} forecasts multivariate time series using ARIMA model
#' 
#' \code{auto.arima} function from \code{forecast} package is used
#' 
#' @param y multivariate time series
#' @param h forecast horizon
#' @param model list of estimated ARIMA models. If missing will be estimated automatically.
#' @param ... further arguments passed to \code{auto.arima} function
#' @return forecasts from ARIMA model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' arima_forecast <- forecast_arima(y_small, h = 2)
forecast_arima <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_arima(y, ...)
  }
  y_matrix <- y 
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

#' Forecast using ETS model
#' 
#' \code{forecast_ets} forecasts time series using ETS model
#' 
#' \code{ets} function from \code{forecast} package is used
#' 
#' @param y multivariate time series
#' @param h forecast horizon
#' @param model list of estimated ETS models. If missing will be estimated automatically.
#' @param ... further arguments passed to \code{ets} function
#' @return forecasts from ETS model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' ets_forecast <- forecast_ets(y_small, h = 2)
forecast_ets <- function(y, h = 1, model = NULL, ...) {
  if (is.null(model)) {
    model <- estimate_ets(y, ...)
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


#' Forecast using VAR-lasso model
#' 
#' \code{forecast_var_lasso} forecasts time series using VAR-lasso model
#' 
#' \code{cv.BigVAR} function from \code{BigVAR} package is used
#' 
#' @param y multivariate time series
#' @param h forecast horizon
#' @param h_cv forecast horizon used for cross validation
#' @param p number of lags
#' @param struct type of VAR-lasso
#' @param gran granularity vector, If gran = (a, b) then cross-validation checks b values form lambda to lambda/a.
#' @param model estimated VAR-lasso model. If missing will be estimated automatically.
#' @param ... further arguments passed to \code{cv.BigVAR} function
#' @param type fast or honest cross-validation. 
#' honest: cross-validation is repeated for each h, so h should be a vector.
#' model should be missing or a list of the same length as h or the same model for all h.
#' fast: do cross-validation for h_cv and than predict for each horizon from 1 to h.
#' @return forecasts from VAR-lasso model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' var_lasso_forecast <- forecast_var_lasso(y_small, h = 2)
forecast_var_lasso <- function(y, h = 1, model = NULL, p = 1, 
                               struct = "OwnOther", gran = c(25, 10),
                               type = c("fast", "honest"), h_cv = 1, ...) {
  type <- match.arg(type)
  
  # var_lasso fails! for time series!
  y_matrix <- as.matrix(y)
  m <- ncol(y_matrix)
  
  forecast_matrix <- matrix(NA, nrow = max(h), ncol = m)    
  
  if (type == "fast") {
    if (is.null(model)) {
      model <- estimate_var_lasso(y, h = h_cv, p = p, 
                                  struct = struct, gran = gran, ...)
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
           Length of 'model' is ", length(model), 
           ", while lenght of 'h' is ", length(h), ".")
    }
    
    # estimate if no models are estimated
    if (length(model) == 0) {
      if (length(h) == 1) {
        model <- estimate_var_lasso(y, h = h, p = p, 
                                    struct = struct, gran = gran, ...)
      } else {
        model <- list()
        for (i in 1:length(h)) {
          model[[i]] <- estimate_var_lasso(y, h = h[i], p = p, 
                                           struct = struct, gran = gran, ...)
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





#' Create lower or upper interval border from future MCMC draws
#' 
#' \code{create_interval_borders} creates lower or upper interval borders from future MCMC draws
#' 
#' Just finds the appropriate quantile of future draws.
#' 
#' @param y_future data frame with n_draws rows. 
#' Each column represent future MCMC draws for specific variable for specific forecast horizon.
#' @param level vector of confidence levels
#' @param type lower or upper for left and right borders of intervals
#' @return matrix of interval borders. One column corresponds for one confidence level.
#' @export
#' @examples 
#' a <- tibble::tibble(x = rnorm(100, mean = 2, sd = 1), 
#'   y = rnorm(100, mean = -2, sd = 4))
#' create_interval_borders(a, type = "lower")
create_interval_borders <- function(y_future, level = c(80, 95), 
                                    type = c("lower", "upper")) {
  type <- match.arg(type)
  
  lower_probs <- (1 - level / 100) / 2
  if (type == "lower") {
    probs <- lower_probs
  } else {
    probs <- 1 - lower_probs
  }
  
  y_future_border <- purrr::map(y_future, ~ quantile(.x, probs = probs))
  y_future_border <- t(tibble::as_tibble(y_future_border))
  colnames(y_future_border) <- paste0(level, "%")
  
  return(y_future_border)
}  


#' Convert future MCMC draws from TVP a la Primiceri to interval borders
#' 
#' \code{primiceri_draws_to_1d_forecast} creates lower and upper interval borders, mean  
#' from future MCMC draws a la Primiceri for one specific variable.
#' 
#' For every h future MCMC draws are obtained, than standard forecast object is created.
#' 
#' @param y_actual original multivariate time series. One column is extracted and is stored in forecast object. 
#' @param model estimated TVP a la primiceri model
#' @param series_no number of time series to predict
#' @param h forecasting horizon
#' @param level vector of confidence levels
#' @return one dimensional complex forecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' y_small_scaled <- scale_to(y_small)
#' primiceri_model <- estimate_tvp_primiceri(y_small_scaled, p = 1)
#' draws <- primiceri_draws_to_1d_forecast(primiceri_model, 
#'   y_actual = y_small_scaled, h = 1, series_no = 1)
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
  y_future <- tibble::as_tibble(y_future)
  
  fors_freq <- stats::frequency(y_actual[, series_no])
  fors_start <- next_obs_time(y_actual[, series_no])
  
  y_future_mean <- purrr::map_dbl(y_future, mean)
  y_future_mean <- stats::ts(y_future_mean, start = fors_start, frequency = fors_freq)
  one_ts_forecast$mean <- y_future_mean
  
  lower_border <- create_interval_borders(y_future, 
                                          level = level, type = "lower")
  lower_border <- stats::ts(lower_border, start = fors_start, frequency = fors_freq)
  one_ts_forecast$lower <- lower_border
  
  upper_border <- create_interval_borders(y_future, 
                                          level = level, type = "upper")
  upper_border <- stats::ts(upper_border, start = fors_start, frequency = fors_freq)
  one_ts_forecast$upper <- upper_border
  
  class(one_ts_forecast) <- "forecast"
  return(one_ts_forecast)
}


#' Forecast using TVP a la Primiceri model
#' 
#' \code{forecast_tvp_primiceri} forecasts time series using TVP model
#' 
#' One dimensional forecasts are obtained from future MCMC draws. 
#' Mean and interval forecasts are calculated from the draws.
#' ACHTUNG: probably we need to scale series for this model!
#' 
#' @param y multivariate time series
#' @param h forecast horizon
#' @param level vector of confidence levels
#' @param model estimated TVP models. If missing will be estimated automatically.
#' @param ... further arguments passed to \code{bvar.sv.tvp} function
#' @return forecasts from TVP model as mforecast object
#' @export
#' @examples 
#' data(rus_macro)
#' y_small <- rus_macro[, c("cpi", "employment", "m2")]
#' y_small_scaled <- scale_to(y_small)
#' tvp_forecast <- forecast_tvp_primiceri(y_small_scaled, h = 2, p = 1)
forecast_tvp_primiceri <- function(y, h = 1, model = NULL, level = c(80, 95), ...) {
  if (is.null(model)) {
    model <- estimate_tvp_primiceri(y, ...)
  }
  y_matrix <- y 
  m <- ncol(y_matrix) # number of series
  
  mforecast <- list()
  mforecast$model <- model
  
  mforecast$forecast <- list()
  
  for (i in 1:m) {
    mforecast$forecast[[i]] <- primiceri_draws_to_1d_forecast(model, 
                            y_actual = y_matrix, series_no = i, h = h, level = level)
  }
  names(mforecast$forecast) <- colnames(y_matrix)
  
  class(mforecast) <- "mforecast"
  return(mforecast)
}


