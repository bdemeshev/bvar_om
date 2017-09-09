# main cycle of estimation

# TODO: test left border for cv in lasso!


# torro may be installed via devtools::install_github("bdemeshev/torro")


library(torro)
library(tidyverse)
library(forecast)


basefolder <- "~/Documents/bvar_om/estimation/bvar_alternatives/"
n_cores <- 4 
# 1 for non-parallel version
# max is equal to parallel::detectCores()
# but maybe non-optimal

# look at torro built-in data frames
# almost every data frame has "toy" little brother: shifts and shifts_toy etc
glimpse(rus_macro)
glimpse(shifts)
glimpse(horizons)
glimpse(var_sets)

glimpse(arguments_arima)
glimpse(arguments_var_lasso)
glimpse(arguments_var_lasso_toy)


glimpse(models)


fits_long <- create_fits_long(shifts, models, var_sets, horizons)
fits_long_toy <- create_fits_long(shifts_toy, models_toy, var_sets_toy, horizons_toy)


readr::write_rds(fits_long, path = paste0(basefolder, "/fits_long.Rds"))



# forecasting

fits_long_subset <- fits_long[1:100, ] # estimate on machine 1 for example

fits_long_toy_2 <- forecast_all(fits_long_toy, var_sets, basefolder = basefolder,
                            n_cores = n_cores)



get_status_from_fit_files(basefolder, fits_long_toy$model_filename)




