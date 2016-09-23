# small scale experiment

# model selection set

library("dplyr")
library("ggplot2")
library("MCS")
library("parallel")


###### Experiment 1: Third alternative irrelevance for MCS???

n <- 100

set.seed(23)
true_y <- rep(0, n) + rnorm(n, sd = 0.0005)
m1 <- rep(0.1, n) + rnorm(n, sd = 0.0005)
m2 <- rep(0.10001, n) + rnorm(n, sd = 0.0005)
m3 <- rep(0.1001, n) + rnorm(n, sd = 0.0005)

m4 <- rep(0.25, n) + rnorm(n, sd = 0.0005)
m5 <- rep(0.25, n) + rnorm(n, sd = 0.0005)
m6 <- rep(0.25, n) + rnorm(n, sd = 0.0005)
m7 <- rep(0.25, n) + rnorm(n, sd = 0.0005)

# MCS with 3 models

forecasts <- cbind(m1, m2, m3)
n_models <- ncol(forecasts)
actual <- matrix(rep(true_y, n_models), ncol = n_models)

loss <- (actual - forecasts) ^ 2
best_models <- MCSprocedure(loss)

# MCS with 7 models

forecasts <- cbind(m1, m2, m3, m4, m5, m6, m7)
n_models <- ncol(forecasts)
actual <- matrix(rep(true_y, n_models), ncol = n_models)

loss <- (actual - forecasts) ^ 2
best_models <- MCSprocedure(loss)

# Moralite
# MCS with 3 models rejects model 3
# MCS with 7 models does not reject model 3 and models with higher loss

###### Experiment 2: Speed as a function of number of models: linear or quadratic

n <- 100



all_n_models <- c(2, 5, 10, 15, 20)
time <- rep(NA, length(all_n_models))

for (i in 1:length(all_n_models)) {
  set.seed(23)
  
  n_models <- all_n_models[i]
  
  true_y <- rep(0, n) + rnorm(n, sd = 0.0005)
  
  forecasts <- matrix(rnorm(n * n_models, sd = 0.0005), ncol = n_models)
  
  actual <- matrix(rep(true_y, n_models), ncol = n_models)
  
  loss <- (actual - forecasts) ^ 2
  
  measure <- system.time(best_models <- MCSprocedure(loss))
  time[i] <- measure[3]
}


qplot( x = all_n_models, y = time)



