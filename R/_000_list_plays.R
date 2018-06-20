# main cycle of estimation

# TODO: test left border for cv in lasso!

# torro may be installed via 
# devtools::install_github("bdemeshev/torro")

# this R script is available at
# https://github.com/bdemeshev/bvar_om/blob/master/R/_000_list_plays.R

library(torro) # estimation of multivariate models
library(tidyverse)
library(forecast)
library(reshape2)
# library(stringr)

missing <- c(3773L, 3774L, 3775L)

basefolder <- "~/Documents/bvar_om/estimation/refits/"
basefolder <- "D:\\Research\\Projects\\ProjectXIV\\Code_2018\\estimation\\refits"

n_cores <- 1
# 1 for non-parallel version
# max is equal to parallel::detectCores()
# but maybe non-optimal

# look at torro built-in data frames
# almost every data frame has "toy" little brother: shifts and shifts_toy etc
# glimpse(rus_macro)
# glimpse(shifts)
# glimpse(horizons)
# glimpse(var_sets)
# 
# glimpse(arguments_arima)
# glimpse(arguments_var_lasso)
# glimpse(arguments_var_lasso_toy)
# 
# 
# glimpse(models)


# create list containing all models fit requested
fits_long <- create_fits_long(shifts, models, var_sets, horizons)
# fits_long_toy <- create_fits_long(shifts_toy, models_toy, var_sets_toy, horizons_toy)

# forecasting

fits_long_subset <- fits_long[missing, ] # estimate on machine 1 for example

fits_long_2 <- forecast_all(fits_long_subset, var_sets, basefolder = basefolder,
                            n_cores = n_cores)


# save all forecasts to 1 file
all_forecasts <- get_forecasts_from_fit_files(basefolder = basefolder)
write_rds(all_forecasts, paste0(basefolder, "all_forecasts_new.Rds"))

# go on
fits_long <- create_fits_long(shifts, models, var_sets, horizons)
all_forecasts <- read_rds(path = paste0(basefolder, "all_forecasts_new.Rds"))

all_f2 <- augment_forecasts(all_forecasts, fits_long, rus_macro)
# row_number means forecasting horizon (sorry, bad name)

rmse_table <- calculate_rmse(all_f2, fits_long, 
                             T_rmse_min = 133, 
                             T_rmse_max = 243)
# row_number means forecasting horizon (sorry, bad name)

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
                            var_set == "set_A", 
                            variable == "ind_prod") 
qplot(data = df, x = row_number, y = mae)




rel_rmse_table <- transform_to_relative_measure(rmse_table)

rel_rmse_subset <- rel_rmse_table %>% filter(
                          model_type %in% c("arima", "ets"),
                          row_number %in% c(1, 3, 6, 9, 12),
                          var_set == "set_C") 

rel_comparison2 <- reshape2::dcast(rel_rmse_subset, 
                                   variable ~ row_number + model_type, 
                                   value.var = "mse")


write_rds(rel_rmse_table, 'rel_rms_table.Rds')

# univariate leader  ------------------------------------------------------
univariate <- rel_rmse_table %>% rename(h = row_number) %>%
  filter(model_type %in% c('arima', 'ets', 'rw')) %>% select(-rmse, -mae) %>% filter(var_set == 'set_C') %>%
  select(-var_set, -pars_id) %>% filter(h %in% c(1, 3, 6, 9, 12))

uni_leader_long <- univariate %>% group_by(h, variable) %>% top_n(-1, mse) %>% ungroup() 


uni_leader_long %>% dcast(variable ~ h, value.var = 'model_type')
uni_leader_long %>% dcast(variable ~ h, value.var = 'mse')


ggplot(data = uni_leader_long, aes(x = factor(h), y = variable)) +
  geom_tile(aes(fill = model_type)) +
  geom_text(aes(label = round(mse, 2))) +
  theme_minimal() + labs(x = 'Forecasting horizon', 
                         title = 'Univariate models MSE relative to random walk MSE',
                         y = 'Variable', 
                         fill = 'Model')




# bad performance of var-lasso
# var-lasso vs RW

# var-lasso border problem


# old leader vs new univariate models




