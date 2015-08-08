# read data from Carriero

library("R.matlab")
library("bvarr")
all <- readMat("../data/usa_data.mat")
str(all)

usa_data <- all$usa.data
usa_data
str(usa_data)

# У Карьеро --- независимое Нормальное Уишарта

priors <- Carriero_priors(Y_in = usa_data, p = 4, s2_lag = 1, 
                          lambdas = c(0.2,,))
model <- bvar_conjugate0(priors=priors)
forecasts <- forecast_conjugate(model, h=12)

?Carriero_priors
?forecast_conjugate
?bvar_conjugate0
