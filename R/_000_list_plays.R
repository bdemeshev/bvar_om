# main cycle of estimation

# TODO: test left border for cv in lasso!


# torro may be installed via devtools::install_github("bdemeshev/torro")


library(torro)
library(tidyverse)
library(forecast)


missing <- c(3773L, 3774L, 3775L, 3776L, 3777L, 3778L, 3779L, 3780L, 3781L, 
             3782L, 3783L, 3784L, 3785L, 3786L, 3787L, 3788L, 3789L, 3790L, 
             3866L, 3867L, 3868L, 3869L, 3870L, 3871L, 3872L, 3873L, 3874L, 
             3875L, 3876L, 3877L, 3878L, 3879L, 3880L, 3881L, 3882L, 3883L, 
             3884L, 3885L, 3886L, 3887L, 3888L, 3889L, 3890L, 3891L, 3892L, 
             3893L, 3894L, 3895L, 3896L, 3897L, 3898L, 3899L, 3900L, 3901L, 
             3902L, 3903L, 3904L, 3905L, 3906L, 3907L, 3908L, 3909L, 3910L, 
             3911L, 3912L, 3913L, 3914L, 3915L, 3916L, 3917L, 3918L, 3919L, 
             3920L, 3921L, 3922L, 3923L, 3924L, 3925L, 3926L, 3927L, 3928L, 
             3929L, 3930L, 3931L, 3932L, 3933L, 3934L, 3935L, 3936L, 3937L, 
             3938L, 3939L)



basefolder <- "~/Documents/bvar_om/estimation/bvar_alternatives/"
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


fits_long <- create_fits_long(shifts, models, var_sets, horizons)
# fits_long_toy <- create_fits_long(shifts_toy, models_toy, var_sets_toy, horizons_toy)


# readr::write_rds(fits_long, path = paste0(basefolder, "/fits_long.Rds"))



# forecasting

fits_long_subset <- fits_long[missing, ] # estimate on machine 1 for example

fits_long_2 <- forecast_all(fits_long_subset, var_sets, basefolder = basefolder,
                            n_cores = n_cores)


# all forecasts to 1 file
basefolder <- "~/Documents/bvar_om/estimation/refits/"
all_forecasts <- get_forecasts_from_fit_files(basefolder = basefolder)
write_rds(all_forecasts, paste0(basefolder, "all_forecasts.Rds"))

# go on
all_forecasts <- read_rds(path = paste0(basefolder, "all_forecasts.Rds"))

all_f2 <- attach_y_to_forecasts(all_forecasts, fits_long, rus_macro)
glimpse(all_forecasts)


glimpse(all_f2)
