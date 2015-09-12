# 500_one_model_estimator.R


source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")

# need to run only once 
source("200_load_after_eviews.R")

#

# df <- read_csv("../data/df_2015_final.csv")
df <- read_csv("../data/usa_carriero.csv")


# observations used for estimating model
# real number of obs for regression will be T_start - T_end - n_lag + 1
T_start <- 48 
T_end <- 171
variables <- paste0("var",1:3) # colnames(df)

l_1 <- 0.2
l_power <- 1 
l_sc <- 2
l_io <- NA
l_const <- 1000
l_exo <- 1 # not used
carriero_hack <- FALSE

n_lag <- 4
v_prior <- NULL # m+2 if NULL

keep <- 2000 # 0 means only posterior are calculated

num_AR_lags <- 1 # 1/NULL
deltas <- 1
verbose <- FALSE
y_bar_type <- "all" # "all" or "initial"


D <- df[T_start:T_end, variables]

lambda <- c(l_1, l_power, l_sc, l_io, l_const, l_exo)

setup <- bvar_conj_setup(D, p=n_lag, v_prior = v_prior,
                         delta = deltas,
                         lambda = lambda, 
                         s2_lag = num_AR_lags,
                         carriero_hack = carriero_hack)

model <- bvar_conj_estimate(setup=setup, verbose=verbose, keep=keep)
bvar_conj_summary(model)
forecasts <- bvar_conj_forecast(model, out_of_sample = TRUE, h = 1)
forecasts


X_plus <- setup$X_plus # ?
Y_plus <- setup$Y_plus # ?


dummy <- bvar_conj_lambda2dummy(Y_in=D, Z_in=NULL, constant=TRUE, p=n_lag,
                                lambda=lambda, delta=deltas, s2_lag=num_AR_lags,
                                y_bar_type=y_bar_type)

X_plus <- dummy$X_plus # ?
Y_plus <- dummy$Y_plus # ?

X <- setup$X # ok
Y <- setup$Y # ok



