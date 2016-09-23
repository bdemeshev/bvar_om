devtools::install_github("bdemeshev/bvarr")
library("bvarr")
data(Yraw)
priors <- Carriero_priors(Yraw, p = 4, lambdas = c(1,0.2,1,1,1))
model <- bvar_conjugate0(priors = priors)
summary_conjugate(model)

fors <- forecast_conjugate(model, h=2, output="wide")

# compare new code and KK old code results
priors <- KK_code_priors(Yraw, p = 4)
model <- bvar_conjugate0(priors = priors)

model0 <- bvar(Yraw, prior="conjugate")
str(model0)
model0$ALPHA_mean
summary_conjugate(model)

# compare forecasts
model0$Y_pred_mean
forecast_conjugate(model)

model0$SIGMA_mean
# Posterior mean of Sigma (noise covariance) [m = 3 x m =3]
# [,1]         [,2]        [,3]
# [1,]  0.096575439 -0.004421863  0.02823719
# [2,] -0.004421863  0.105771905 -0.09194744
# [3,]  0.028237195 -0.091947442  0.54391353

model0$SIGMA_std
# Posterior sd of Sigma (noise covariance) [m = 3 x m =3]
# [,1]        [,2]       [,3]
# [1,] 0.009440972 0.006966363 0.01623650
# [2,] 0.006966363 0.010347294 0.01772981
# [3,] 0.016236502 0.017729811 0.05356543

# for SIGMA --- ok


# ALPHA mean == Phi mean --- ok
model0$ALPHA_mean
# summary_conjugate(model)
# Posterior mean of Phi (VAR coefficients) [k = 13 x m = 3]:
#   [,1]         [,2]         [,3]
# [1,]  1.491110138  0.008008534  0.537453442
# [2,] -0.263370173  1.285385572 -0.712524637
# [3,] -0.057842066 -0.025195833  0.776555534
# [4,] -0.445884646  0.059046794 -0.748876634
# [5,]  0.188423679 -0.342328595  0.774547281
# [6,]  0.060486424 -0.020736445 -0.030932113
# [7,] -0.080656930 -0.155020939  0.793720496
# [8,] -0.013783310 -0.134118019 -0.339940837
# [9,] -0.005791011  0.095152462  0.104258574
# [10,]  0.035631098  0.124106183 -0.474718459
# [11,]  0.042424828  0.076843821  0.305573177
# [12,]  0.000871978 -0.024545726  0.055467314
# [13,]  0.303753151  0.418185858 -0.003191253

library(MSBVAR)
# 
model <- szbvar(Yraw, p=4, ... )


# test on artificial data --- ok
set.seed(7)
n = 10^5
y1 <- arima.sim(n = n, list(ar=0.5))
y2 <- arima.sim(n = n, list(ar=-0.3))
Yart <- cbind(y1, y2)
priors <- Carriero_priors(Yart, p = 1, lambdas = c(1,0.2,1,1,1))
model <- bvar_conjugate0(Yart, p=1, priors = priors)
summary_conjugate(model)

# test forecast on artificial data 
set.seed(7)
n = 10^5
y1 <- arima.sim(n = n+1, list(ar=0.5))
y2 <- arima.sim(n = n+1, list(ar=-0.3))
Yart <- cbind(y1, y2)
true_y <- tail(Yart,1)

Yart <- head(Yart,-1)
#priors <- KK_code_priors(Yart, p = 4)
priors <- Carriero_priors(Yart, p = 1, lambdas = c(1,0.2,1,1,1))
model <- bvar_conjugate0(Yart, p=1, priors = priors)
# summary_conjugate(model)
forecast_conjugate(model, output="wide", type = "credible")
true_y
tail(Yart,1) * c(0.5,-0.3)

model0 <- bvar(Yart, prior="conjugate", p=4)
model0$Y_pred_mean


model_sz <- szbvar(ts(Yraw), p=2, lambda0=1, lambda1=1,
                   lambda3=1, lambda4=1, lambda5=1,
                   mu5=1, mu6=1)
sz_priors <- SimZha_priors(Yraw, p=2, lambdas=c(1,1,1,1,1,1),
                           mu56=c(1,1))
model_conj <- bvar_conjugate0(priors=sz_priors)
summary_conjugate(model_conj)
str(model_sz)
res <- forecast(model_sz, nsteps = 2)
res[216:217,]

res_c <- forecast_conjugate(model_conj, h=2)
head(res_c,6)

######################

library("MSBVAR")
library("bvarr")
data(Yraw)
priors <- SimZha_priors(Yraw, p=2, lambdas=c(1,1,1,1,1,1),
                           mu56=c(1,1))

model_sz <- szbvar(ts(Yraw), p=2, lambda0=1, lambda1=1,
                   lambda3=1, lambda4=1, lambda5=1,
                   mu5=1, mu6=1)

Y_in=NULL
Z_in=NULL 
constant=TRUE
p=NULL
keep=10000
verbose=FALSE

# run part of function conjugate

Y_my <- Y
X_my <- X

# get inside priors

Y_in <- Yraw
Z_in <- NULL
constant=TRUE
p=2
lambdas=c(1,1,1,1,1,1)
mu56=c(1,1)
VAR_in=c("levels")

# some steps in SimZha priors

head(priors$Y)
head(model_sz$Y)
priors$sigmas_sq
s2

priors$Omega_prior
priors$S_prior
model_sz$S0
