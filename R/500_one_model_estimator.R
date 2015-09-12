# 500_one_model_estimator.R


source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")

# need to run only once 
source("200_load_after_eviews.R")

#

df <- read_csv("../data/df_2015_final.csv")

# observations used for estimating model
# real number of obs for regression will be T_start - T_end - n_lag + 1
T_start <- 10 
T_end <- 120
variables <- c("ind_prod","cpi","m2")

l_1 <- 1
l_power <- 1 
l_sc <- 1
l_io <- 1
l_const <- 1
l_exo <- 1
keep <- 0 # 0 means only posterior are calculated
n_lag <- 5


num_AR_lags <- 1 # 1/NULL
deltas <- 1
verbose <- FALSE


D <- df[T_start:T_end, variables]

setup <- bvar_conj_setup(D, p=n_lag, 
                         delta = deltas,
                         lambda = c(l_1, l_power, l_sc, l_io, l_const, l_exo), 
                         s2_lag = num_AR_lags)

model <- bvar_conj_estimate(setup=setup, verbose=verbose, keep=keep)
bvar_conj_summary(model)
forecasts <- bvar_conj_forecast(model, out_of_sample = TRUE, h = 3)
