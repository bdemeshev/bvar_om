# read data from Carriero

library("R.matlab")
all <- readMat("../data/usa_data.mat")
str(all)

usa_data <- all$usa.data
usa_data
str(usa_data)

priors <- Carriero_priors(Y_in = usa_data, p = 12, s2_lag = 1)
model <- bvar_conjugate0(priors=priors)
forecasts <- forecast_conjugate(model, h=12)

?forecast_conjugate
