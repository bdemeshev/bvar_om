#' Means and standard deviations of time series
#' 
#' \code{get_scales} returns tibble with means and standard deviations of times series
#' 
#' The function calculates means and standard deviations of multivariate time series or data.frame object. 
#' Missing values are omitted. The result is output as tibble.
#' 
#' @param y multivariate time series, matrix or data.frame
#' @return tibble with three columns: variable, mean, sd
#' @export
#' @examples 
#' get_scales(cars)
get_scales <- function(y) {
  variable <- value <- NULL # black magic to remove NOTE in R CMD check
  
  y_tibble <- tibble::as.tibble(y)
  y_long <- reshape2::melt(y_tibble, id.vars = NULL, na.rm = TRUE)
  y_long <- dplyr::mutate(y_long, variable = as.character(variable))
  y_sum <-  dplyr::summarise(dplyr::group_by(y_long, variable),
              mu = mean(value), sd = stats::sd(value))
  return(y_sum)
}

#' Scale multivariate time series to specific means and standard deviations
#' 
#' \code{scale_to} scales multivariate time series to specific means and standard deviations
#' 
#' The function scales multivariate time series to specific means and standard deviations.
#' If mu_sd tibble is not specified then times series are scaled to zero mean and unit standard deviation.
#' 
#' @param y multivariate time series to be scales
#' @param mu_sd tibble with three columns: variable, mean, sd
#' @return scaled mutlivatiate time series
#' @export
#' @examples 
#' scale_to(cars)
#' data(rus_macro)
#' scale_to(rus_macro)
scale_to <- function(y, mu_sd = NULL) {
  # black magic to remove NOTE in R CMD check
  row_number <- variable <- value <- sd <- mu <- .id <-  NULL 
  
  y_tibble <- dplyr::mutate(tibble::as.tibble(y), .id = row_number())

  if (is.null(mu_sd)) {
    mu_sd <- tibble::tibble(variable = colnames(y_tibble), mu = 0, sd = 1)
  }

  y_long <- reshape2::melt(y_tibble, id.vars = ".id")
  y_long <- dplyr::mutate(y_long, variable = as.character(variable))
  y_long <- dplyr::left_join(y_long, mu_sd, by = "variable")
  y_long <- dplyr::mutate(dplyr::group_by(y_long, variable),
                          new_value = scale(value) * sd + mu)
  y_scaled <- reshape2::dcast(y_long, .id ~ variable, value.var = "new_value")
  y_scaled <-  dplyr::select(y_scaled, -.id)

  if (stats::is.ts(y)) {
    y_scaled <- stats::ts(y_scaled,
                  frequency = stats::frequency(y), start = stats::start(y))
  } 
  return(y_scaled)
}
