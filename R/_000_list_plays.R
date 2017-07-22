# main cycle of estimation

source("_000_forecast_mts.R")

# here are separated files for each fit:
fits_folder <- "../estimation/bvar_alternatives/fits/"
# here lies the list of all fits:
bvar_alt_folder <- "../estimation/bvar_alternatives/"
# here lie all generated files
estimation_folder <- "../estimation/"

# create folders if they do not exist
if (!dir.exists(estimation_folder)) {
  dir.create(estimation_folder)
}

if (!dir.exists(bvar_alt_folder)) {
  dir.create(bvar_alt_folder)
}

if (!dir.exists(fits_folder)) {
  dir.create(fits_folder)
}

rus_macro <- load_rus_data()

# create all shifts (like moving/expanding...)
shifts <- tribble(~shift_name, ~shift_T_start, ~win_expanding, ~win_start_length, ~n_shifts,
                  "moving_120", 1, FALSE, 120, 10,
                  "expanding_120", 1, TRUE, 120, 10)

# describe non lasso models (automatic)
fake_args <- tibble(pars_id = "automatic")
models <- tribble(~model_type, ~comment, ~h_required, ~model_args, 
                  "arima", "auto arima applied series wise", FALSE, fake_args,
                  "ets", "auto ets applied series wise", FALSE, fake_args,
                  "rw", "random walk applied series wise", FALSE, fake_args)


# lasso params:
var_lasso_args <- crossing(p = 1:12, 
                           struct = c("OwnOther"), 
                           gran_1 = 25, gran_2 = 10, 
                           pars_id = paste0(struct, "_", p))
# lasso type:
# three top players according to 
# https://arxiv.org/pdf/1508.07497.pdf page 23
# "SparseLag" (Lag Sparse Group VARX-L)
# "OwnOther" (Own/Other Group VARX-L)
# "SparseOO" (Own/Other Sparse Group VARX-L)
# but:
# "SparseLag" hangs indefinetely
# "SparseOO" probably hangs indefinetely

# granularity:
# we have grid of 10 points from lambda_max/25 to lambda_max
# 25 and 10 suggested by BigVAR user's guide, 
# http://www.wbnicholson.com/BigVAR.pdf, page 9
# ideally we would like to do this:
# var_lasso_args$granularity <- rep(list(c(25, 10)), nrow(var_lasso_pars))
# but group_by does not work with ugly variables 
# so we introduce gran_1 and gran_2 and unite them after group_by



var_lasso_models <- tribble(~model_type, ~model_args, ~h_required, ~comment,
                            "var_lasso", var_lasso_args, TRUE, "var lasso with cross validation for each h")
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


models_long <- unnest(models) # not used for computations but maybe useful for user
readr::write_rds(models_long, path = paste0(bvar_alt_folder, "models_long.Rds"))


fits <- crossing(samples, models, var_set = var_sets$var_set, horizons)
fits_long <- unnest(fits) 



# the function sets maximal h for h-agnostic models 
# (models that do not require h during estimation phase)
set_agnostic_max_h <- function(fits_long) {
  fits_h_gnostic <- filter(fits_long, h_required == TRUE)
  fits_h_agnostic <- filter(fits_long, h_required == FALSE)

  group_by_vars <- setdiff(colnames(fits_long), c("h", "model_filename"))
  
  fits_h_agnostic <- fits_h_agnostic %>% 
    group_by_at(.vars = group_by_vars) %>% 
    filter(h == max(h)) %>% ungroup()
  
  fits_long <- bind_rows(fits_h_gnostic, fits_h_agnostic)
  return(fits_long)
}

fits_long <- set_agnostic_max_h(fits_long)


# add model_filename (used as id) just before cycle
fits_long <- fits_long %>% mutate(result = "Non-estimated",
                model_filename = paste0("fit_", row_number(), ".Rds"))



non_parameter_colnames <- setdiff(c(colnames(fits), "model_filename", "pars_id", "result"), c("model_args"))


granulatiry_to_vector <- function(df) {
  df$gran <- c(df$gran_1, df$gran_2)
  df$gran_1 <- NULL
  df$gran_2 <- NULL
  return(df)
}


forecast_one_fit <- function(fits_long, fit_no) {
  fit_row <- fits_long[fit_no, ]
  
  parameter_colnames <- setdiff(colnames(fit_row), non_parameter_colnames)
  parameters <- fit_row[, parameter_colnames]
  # remove NA parameters that are not required for actual fit:
  parameters <- as.list(parameters[, as.vector(!is.na(parameters[1, ]))])
  # concatenate two granularity parameters in one vector:
  if (fit_row$model_type == "var_lasso")  {
    parameters <- granulatiry_to_vector(parameters)
  }
  if (fit_row$h_required) {
    parameters$h <- fit_row$h
  }
  
  forecast_fun_name <- paste0("forecast", "_", fit_row$model_type)
  forecast_fun <- eval(parse(text = forecast_fun_name))
  
  vars <- var_sets %>% filter(var_set == fit_row$var_set) %>% .$variable
  y <- rus_macro[fit_row$T_start:fit_row$T_end, vars]
  parameters$y <- y
  
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
# testing:
# forecast_one_fit(fits_long, 1)


# forecasting
for (fit_no in 1:nrow(fits_long)) {
  cat("Processing fit no ", fit_no, " out of ", nrow(fits_long), ", ", fits_long$model_type[fit_no], "...\n")
  
  fits_long$result <- forecast_one_fit(fits_long, fit_no)

  readr::write_rds(fits_long, path = paste0(bvar_alt_folder, "fits_long.Rds"))
  
  cat("Processing fit no ", fit_no, " out of ", nrow(fits_long), ", ", fits_long$model_type[fit_no], "done.\n")
}


readr::write_rds(fits_long, path = paste0(bvar_alt_folder, "fits_long.Rds"))


