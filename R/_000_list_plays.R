# main cycle of estimation

source("_000_forecast_mts.R")
fits_folder <- "../estimation/bvar_alternatives/fits/"
bvar_alt_folder <- "../estimation/bvar_alternatives/"
save_fits_long_every <- 10 # after every 10 moves update fits_long

rus_macro <- load_rus_data()

# create all shifts (like moving/expanding...)
shifts <- tribble(~shift_name, ~shift_T_start, ~win_expanding, ~win_start_length, ~n_shifts,
                  "moving_120", 1, FALSE, 120, 10,
                  "expanding_120", 1, TRUE, 120, 10)

# describe non lasso models (automatic)
models <- tribble(~model_type, ~comment,
                  "arima", "auto arima applied series wise",
                  "ets", "auto ets applied series wise",
                  "rw", "random walk applied series wise")
fake_arg <- tibble(pars_id = "automatic")
models$model_args <- rep(list(fake_arg), 3)
models$h_required <- rep(FALSE, 3)

# lasso params:
var_lasso_pars <- crossing(p = 1:12, struct = c("SparseLag", "OwnOther", "SparseOO"))
# three top players according to 
# https://arxiv.org/pdf/1508.07497.pdf page 23
# "SparseLag" (Lag Sparse Group VARX-L)
# "OwnOther" (Own/Other Group VARX-L)
# "SparseOO" (Own/Other Sparse Group VARX-L)

# we have grid of 10 points from lambda_max/25 to lambda_max
# 25 and 10 suggested by BigVAR user's guide, 
# http://www.wbnicholson.com/BigVAR.pdf, page 9
# ideally we would like to do this:
# var_lasso_pars$granularity <- rep(list(c(25, 10)), nrow(var_lasso_pars))
# but group_by does not work with ugly variables 
# so we introduce gran_1 and gran_2 and unite them after group_by
var_lasso_pars$gran_1 <- 25
var_lasso_pars$gran_2 <- 10


var_lasso_models <- tribble(~model_type, ~model_args, ~h_required, ~comment,
                            "var_lasso", var_lasso_pars, TRUE, "var lasso with cross validation for each h")
models <- bind_rows(var_lasso_models, models)

# create all variable sets
create_var_set_info <- function() {
  add_A <- dplyr::tibble(var_set = "set_A", variable = c("ind_prod", "cpi", "ib_rate"))
  add_B <- dplyr::tibble(var_set = "set_B", variable = c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))
  add_C <- dplyr::tibble(var_set = "set_C", variable = c("employment", "ind_prod", "cpi", 
                                                             "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                                                             "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                                                             "agriculture", "retail", "gov_balance", "export", "import"))
  
  var_set_info <- dplyr::bind_rows(add_A, add_B, add_C)
  var_set_info <- dplyr::mutate(var_set_info, pre_transform = "none", post_transform = "none")
  
  return(var_set_info)
}

var_sets <- create_var_set_info()


# create all horizons
horizons <- tibble(h = c(1, 3, 6, 12))



# GO!
shifts_to_samples <- function(shifts) {
  line_no <- rep(1:nrow(shifts), times = shifts$n_shifts)
  samples <- shifts[line_no, ]
  samples <- samples %>% group_by(shift_name) %>% mutate(shift_no = row_number())
  samples <- samples %>% ungroup() %>% 
    mutate(T_end = shift_T_start + win_start_length + shift_no - 2,
           T_start = T_end - win_start_length - win_expanding * (shift_no - 1) + 1,
           sample_name = paste0(shift_name, "_", shift_no))
  return(samples)
}

samples <- shifts_to_samples(shifts)

models_long <- unnest(models) # do we need it?


fits <- crossing(samples, models, var_set = var_sets$var_set, horizons)
fits_long <- unnest(fits) 

is_not_h <- function(col_names) {
  return(setdiff(col_names, c("h", "fit_id", "model_filename")))
}


# the function sets maximal h for h-agnostic models 
# (models that do not require h during estimation phase)
set_agnostic_max_h <- function(fits_long) {
  fits_h_gnostic <- filter(fits_long, h_required == TRUE)
  fits_h_agnostic <- filter(fits_long, h_required == FALSE)
  group_by_vars <- is_not_h(colnames(fits_long))
  fits_h_agnostic <- fits_h_agnostic %>% 
    group_by_at(.vars = group_by_vars) %>% 
    filter(h == max(h)) %>% ungroup()
  fits_long <- bind_rows(fits_h_gnostic, fits_h_agnostic)
  return(fits_long)
}

fits_long <- set_agnostic_max_h(fits_long)


# add id just before cycle
fits_long <- fits_long %>% mutate(fit_id = paste0("fit_", row_number()),
                model_filename = paste0("fit_", row_number(), ".Rds"))



non_parameter_colnames <- setdiff(c(colnames(fits), "model_filename", "pars_id", "fit_id"), c("model_args"))

fit_no <- 42 # for testing. remove later!!!!!

granulatiry_to_vector <- function(df) {
  df$gran <- list(c(df$gran_1, df$gran_2))
  df <- dplyr::select(df, -gran_1, -gran_2)
  return(df)
}


forecast_one_fit <- function(fits_long, fit_no) {
  fit_row <- fits_long[fit_no, ]
  
  parameter_colnames <- setdiff(colnames(fit_row), non_parameter_colnames)
  parameters <- fit_row[, parameter_colnames]
  # remove NA parameters that are not required for actual fit:
  parameters <- parameters[, as.vector(!is.na(parameters[1, ]))]
  parameters <- granulatiry_to_vector(parameters)
  
  forecast_fun_name <- paste0("forecast", "_", fit_row$model_type)
  forecast_fun <- eval(parse(text = forecast_fun_name))
  
  vars <- var_sets %>% filter(var_set == fit_row$var_set) %>% .$variable
  y <- rus_macro[fit_row$T_start:fit_row$T_end, vars]
  parameters$y <- list(y)
  
  attempt <- try(do.call(forecast_fun, parameters))
  
  if ("try-error" %in% class(attempt)) {
    result <- as.character(attempt) 
  } else {
    # write file
    readr::write_rds(attempt, path = paste0(fits_folder, fit_row$model_filename))
    result <- "OK"
  }
  return(result)
}

readr::write_rds(fits_long, path = paste0(bvar_alt_folder, "fits_long.Rds"))
forecast_one_fit(fits_long, 1)


# forecasting
for (fit_no in 1:nrow(fits_long)) {
  fits_long$result <- forecast_one_fit(fits_long, fit_no)
  if (fit_no %% save_fits_long_every == 0) {
    readr::write_rds(fits_long, path = paste0(bvar_alt_folder, "fits_long.Rds"))
  }
}

readr::write_rds(fits_long, path = paste0(bvar_alt_folder, "fits_long.Rds"))


