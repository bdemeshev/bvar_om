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

# BigVAR
Model1 <- constructModel(as.matrix(df), p = 12, struct = "OwnOther", gran = c(100, 5),
                        verbose = TRUE, VARX = list())
model <- cv.BigVAR(Model1)
model
plot(model)
SparsityPlot.BigVAR.results(model)

predict(model, n.ahead = 1)

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


