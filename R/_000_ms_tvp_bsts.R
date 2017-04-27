# экспресс-сравнение трёх моделей-пакетов: MSBVAR/bsts/bvarsv
library(MSBVAR)
library(bsts)
library(bvarsv)
library(vars) # VAR
library(forecast) # rw
library(BigVAR)
library(BMR)
library(tidyverse)
library(reshape2)
library(forecastHybrid)
library(caret)

# load data
df <- readr::read_csv("../data/df_2015_final.csv")
df <- dplyr::select(df, -time_y)


lala <- model@OptimalLambda
# BigVAR
Model1 <- constructModel(as.matrix(df), p = 2, struct = "OwnOther", 
                         T1 = 241, T2 = 243, 
                         gran = lala,
                         ownlambdas = TRUE,
                        verbose = TRUE, VARX = list())
model <- cv.BigVAR(Model1)
model


Model1_est <- constructModel(as.matrix(df)[3:243, ], p = 2, struct = "OwnOther", 
                         gran = lala,
                         ownlambdas = TRUE,
                         verbose = TRUE, VARX = list())

model_est <- BigVAR.est(Model1_est)
model_est

a1 <- model@betaPred
a2 <- model_est$B[, , 1]

model@LambdaGrid
model@OptimalLambda
str(model@betaPred)

# granularity = c(grid depth, number of lambda)
# lambda_max is calculated automatically (all beta hat = 0)
# grid depth = lambda_max / lambda_min

plot(model)
SparsityPlot.BigVAR.results(model)

predict(model, n.ahead = 1)

?constructModel

df_matrix <- as.matrix(df)
str(df_matrix)
df_subset <- df_matrix[1:243, ]
model_opt <- constructModel(df_subset, p = 2, struct = "OwnOther",
                         T1 = 1,
                         T2 = 10,
                         ownlambdas = TRUE,
                         gran = model@OptimalLambda,
                         verbose = TRUE, VARX = list())

?BigVAR.est
one_bigvar <- BigVAR.est(model_opt)


all.equal(model@betaPred[, ], one_bigvar$B[, , 1])

# minus one

model_opt_last <- constructModel(df_matrix[1:243, ], p = 2, struct = "OwnOther",
                            ownlambdas = TRUE,
                            gran = model@OptimalLambda,
                            verbose = TRUE, VARX = list())

one_bigvar_last <- BigVAR.est(model_opt_last)
one_bigvar_last$B[1, 1, 1]


for (t_from in 1:4) {
  for (t_to in 238:239) {
    model_opt_last <- constructModel(df_matrix[t_from:t_to, ], p = 2, struct = "OwnOther",
                                     ownlambdas = TRUE,
                                     gran = lala,
                                     verbose = TRUE, VARX = list())
    
    one_bigvar_last <- BigVAR.est(model_opt_last)
    message(t_from, "-", t_to, " b[1, 1] = ", one_bigvar_last$B[1, 1, 1])
  }
}


model@betaPred[1, 1]


all.equal(model@betaPred[, ], one_bigvar_last$B[, , 1])
all.equal(one_bigvar_last$B[, , 1], one_bigvar$B[, , 1])


# MSBVAR


# bsts


# bvarsv

y <- ts(df)
model <- bvar.sv.tvp(y, p = 12, nrep = 50, nburn = 10, k_B = 25)
# big dimension error!

# 3 russian variables
y <- ts(df[, c("ind_prod", "cpi", "ib_rate")])
head(y)
model <- bvar.sv.tvp(y, p = 2, nrep = 1000, nburn = 5000)


# reordering matters!
data("usmacro")

set.seed(42)
y <- usmacro
model <- bvar.sv.tvp(y, p = 2, nrep = 5000, nburn = 1000)

forecast_y1 <- predictive.draws(model, v = 1, h = 1)
forecast_y1

y_b <- usmacro[, c(3, 2, 1)] # reorder

model_b <- bvar.sv.tvp(y_b, p = 2, nrep = 5000, nburn = 1000)

forecast_y1_b <- predictive.draws(model_b, v = 3, h = 1)

y_future <- data_frame(y1 = forecast_y1$y, y1_b = forecast_y1_b$y)
y_future_melt <- melt(y_future)

qplot(data = y_future_melt, geom = "density", x = value, color = variable)

plot(usmacro)
head(usmacro)
nrow(usmacro)
?usmacro

# no movement in A_t

y <- usmacro
model <- bvar.sv.tvp(y, p = 2, nrep = 5000, nburn = 1000, k_S = 0.0000001)
# lower k_S gives error

forecast_y1 <- predictive.draws(model, v = 1, h = 1)
hist(forecast_y1$y)

y <- usmacro[, c(3, 2, 1)] # reorder

model_b <- bvar.sv.tvp(y, p = 2, nrep = 5000, nburn = 1000, k_S = 0.0000001)

forecast_y1_b <- predictive.draws(model, v = 3, h = 1)

hist(forecast_y1$y)
hist(forecast_y1_b$y)


# BMR
y <- usmacro

y <- ts(df)
bvartvptest <- BVARTVP(y, coefprior = NULL, tau = 80,
                      p = 12, keep = 30, burnin = 10, 
                      XiBeta = 4, XiQ = 0.005, gammaQ = NULL, 
                      XiSigma = 1, gammaS = NULL)
predict(bvartvptest, h = 1)
str(bvartvptest)
# VAR


# rw


# forecast hybrid
y <- df$ib_rate
y <- ts(y, frequency = 12, start = c(1990, 1))
model_h <- hybridModel(y, models = "at")
forecast(model_h, h = 12)


