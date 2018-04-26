library(torro) # estimation of multivariate models
library(tidyverse)
library(forecast)
library(skimr)
# library(stringr)

basefolder <- "~/Documents/bvar_om/estimation/refits/"

fits_long <- create_fits_long(shifts, models, var_sets, horizons)
# fits_long_toy <- create_fits_long(shifts_toy, models_toy, var_sets_toy, horizons_toy)

actual_files <- list.files(paste0(basefolder, "fits/"))
due_files <- fits_long %>% filter(!model_type == "var_lasso") %>% .$model_filename
setdiff(due_files, actual_files)

fits_long_lasso <- mutate(fits_long, 
    corner_optim = get_optimum_at_border_flags(model_filename, basefolder)) %>% 
  filter(model_type == "var_lasso")

glimpse(fits_long_lasso)
table(fits_long_lasso$corner_optim, fits_long_lasso$p)
table(fits_long_lasso$corner_optim, fits_long_lasso$var_set)

# help(package = "torro")


# save all forecasts to 1 file
all_forecasts <- get_forecasts_from_fit_files(basefolder = basefolder)
write_rds(all_forecasts, paste0(basefolder, "all_forecasts.Rds"))

# go on - quick entry point provided all_forecasts.Rds exists
fits_long <- create_fits_long(shifts, models, var_sets, horizons)
all_forecasts <- read_rds(path = paste0(basefolder, "all_forecasts.Rds"))

all_f2 <- augment_forecasts(all_forecasts, fits_long, rus_macro)
# row_number means forecasting horizon (sorry, bad name)

rmse_table <- calculate_rmse(all_f2, fits_long, 
                             T_rmse_min = 133, 
                             T_rmse_max = 243)


# 133 # start 2006
# 243 # last obs in rus_macro, in 2015


df <- rmse_table %>% filter(model_type == "arima", 
                            var_set == "set_C", 
                            variable == "ind_prod") 
qplot(data = df, x = row_number, y = rmse)


df <- rmse_table %>% filter(model_type == "ets", 
                            var_set == "set_C", 
                            variable == "ind_prod") 
qplot(data = df, x = row_number, y = rmse)

df <- rmse_table %>% filter(model_type == "rw", 
                            var_set == "set_C", 
                            variable == "ind_prod") 
qplot(data = df, x = row_number, y = rmse)

df <- rmse_table %>% filter(model_type == "var_lasso",
                            pars_id == "OwnOther_1",
                            var_set == "set_B", 
                            variable == "ind_prod") 
qplot(data = df, x = row_number, y = mae)




rel_rmse_table <- transform_to_relative_measure(rmse_table)

rel_rmse_subset <- rel_rmse_table %>% filter(
  ((model_type %in% c("arima", "ets")) & (var_set == "set_C")) |
     ((model_type %in% c("var_lasso")) & (var_set == "set_C") & (pars_id == "OwnOther_12")),
  row_number %in% c(1, 3, 6, 9, 12)
  ) 

rel_comparison2 <- reshape2::dcast(rel_rmse_subset, 
                                   variable ~ row_number + model_type, 
                                   value.var = "mse")

