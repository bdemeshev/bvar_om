#' Coerce gran_1 and gran_2 scalars in list into vector gran 
#' 
#' \code{granulatiry_to_vector} coerces gran_1 and gran_2 scalars in list into vector gran in the same list.
#' 
#' Coerces gran_1 and gran_2 scalars in list into vector gran in the same list.
#'  
#' @param df list with gran_1 and gran_2 scalar numbers
#' @return list with gran vector 
#' @export
#' @examples 
#' a <- list(gran_1 = 25, gran_2 = 10)
#' granularity_to_vector(a)
granularity_to_vector <- function(df) {
  df$gran <- c(df$gran_1, df$gran_2)
  df$gran_1 <- NULL
  df$gran_2 <- NULL
  return(df)
}




#' Set maximal horizon for h-agnostic models 
#' 
#' \code{set_agnostic_max_h} sets maximal horizon for h-agnostic models.
#' 
#' There are two types of models. Horizon agnostic that do not require horizon h for estimation (ARIMA, RW, ETS...). 
#' Horizon gnostic models that require h both for estimation and forecasting (cross validated VAR-lasso, ...). 
#' This function removes redundant estimation of h-agnostic models. If one requires estimation of h-agnostic
#' model for h = 1, h = 2, h = 3 this function will leave only estimation of h = 3 in the request.
#' List of models to estimate should contain at least
#' \itemize{
#' \item h forecasting horizon
#' \item h_required logical TRUE/FALSE, FALSE for h-agnostic models
#' }
#' All other variables except \code{ignored_vars} are used to identify whether models are different.
#'    
#' @param fits_long list of models to estimate
#' @param ignored_vars names of variables that are ignored 
#' when \code{set_agnostic_max_h} groups rows correponding to same model. 
#' @return shorter list of models to estimate 
#' @export
#' @examples 
#' fits_long <- tibble::tibble(h_required = c(TRUE, TRUE, FALSE, FALSE, FALSE),
#'    h = c(1, 2, 1, 2, 3), params = c(1, 1, 2, 2, 5))
#' set_agnostic_max_h(fits_long)
set_agnostic_max_h <- function(fits_long, ignored_vars = "model_filename") {
  h_required <- h <- NULL # black magic to remove NOTE in R CMD check
  
  fits_h_gnostic <- dplyr::filter(fits_long, h_required == TRUE)
  fits_h_agnostic <- dplyr::filter(fits_long, h_required == FALSE)
  
  group_by_vars <- setdiff(colnames(fits_long), c("h", ignored_vars))
  
  fits_h_agnostic <- dplyr::group_by_at(fits_h_agnostic, .vars = group_by_vars) 
  fits_h_agnostic <- dplyr::ungroup(dplyr::filter(fits_h_agnostic, h == max(h)))

  fits_long <- dplyr::bind_rows(fits_h_gnostic, fits_h_agnostic)
  return(fits_long)
}


#' Transform shifts data frame to samples data frame
#' 
#' \code{shifts_to_samples} transforms shifts data frame to samples data frame.
#' 
#' Transforms shifts data frame to samples data frame. Shifts data frame should contain
#' \itemize{
#' \item shift_name name of the shift
#' \item win_expanding logical, TRUE for growing window, FALSE for moving window
#' \item shift_T_start first T of the first window
#' \item win_start_length length of the first window
#' \item n_shifts number of window shifts
#' }
#' Samples data frame contains the same columns as shifts data frame plus three: T_start, T_end, sample_name. 
#'  
#' @param shifts shifts data frame
#' @return samples data frame. 
#' @export
#' @examples 
#' shifts <- tibble::tribble(~shift_name, ~shift_T_start, ~win_expanding, 
#'     ~win_start_length, ~n_shifts, "moving_120", 13, FALSE, 120, 2) 
#' shifts_to_samples(shifts)
shifts_to_samples <- function(shifts) {
  # black magic to remove NOTE in R CMD check:
  row_number <- shift_name <- shift_T_start <- NULL 
  shift_no <- T_end <- win_expanding <- NULL 
  win_start_length <- NULL
  
  line_no <- rep(1:nrow(shifts), times = shifts$n_shifts)
  samples <- shifts[line_no, ]
  
  samples <-  dplyr::group_by(samples, shift_name) 
  samples <- dplyr::mutate(samples, shift_no = row_number())
  
  samples <- dplyr::mutate(dplyr::ungroup(samples), 
           T_end = shift_T_start + win_start_length + shift_no - 2,
           T_start = T_end - win_start_length - win_expanding * (shift_no - 1) + 1,
           sample_name = paste0(shift_name, "_", shift_no))
  return(samples)
}
