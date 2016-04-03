library("dplyr")
library("ggplot2")
library("MCS")
library("parallel")


###### Experiment 3: parallel computing
?MCSprocedure

detectCores()

time <- rep(NA, 8)

for (n_cores in 1:8) {
  message("Testing ", n_cores, " cores :)")
  cluster <- makeCluster(n_cores)

  set.seed(23)

  n <- 100

  n_models <- 15
  true_y <- rep(0, n) + rnorm(n, sd = 0.0005)
  forecasts <- matrix(rnorm(n * n_models, sd = 0.0005), ncol = n_models)
  actual <- matrix(rep(true_y, n_models), ncol = n_models)
  loss <- (actual - forecasts) ^ 2
  best_models <- MCSprocedure(loss, cl = cluster)
  stopCluster(cluster)

  time[n_cores] <- best_models@Info$elapsed.time
}

cat(time)
