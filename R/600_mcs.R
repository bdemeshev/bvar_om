library("MCS")
library("bvarr")
library("dplyr")
library("readr")

help(package = "bvarr")
help(package = "MCS")

data(Yraw)

Ytrain <- Yraw[1:200, ]
Ytest <- Yraw[201:215, ]

setup <- list()
setup[[1]] <- bvar_conj_setup(Ytest, p = 4, 
                      lambda = c(0.2, 1, 1, 1, 100, 100))
setup[[2]] <- bvar_conj_setup(Ytest, p = 2, 
                              lambda = c(0.2, 1, 1, 1, 100, 100))
setup[[3]] <- bvar_conj_setup(Ytest, p = 5, 
                              lambda = c(0.2, 1, 1, 1, 100, 100))

n_models <- length(setup)

model <- list()
pred <- list()

h <- 15

for (i in 1:n_models) {
  model[[i]] <- bvar_conj_estimate(setup[[i]], keep = 100)
  pred[[i]] <- bvar_conj_forecast(model[[i]], h = h, include = "mean")
}



variable_name <- "inflation"

pred_matrix <- matrix(0, nrow = h, ncol = n_models)

for (i in 1:n_models) {
  pred_filtered <- filter(pred[[i]], variable == variable_name)
  pred_matrix[, i] <- pred_filtered$value
}

true_matrix <- matrix(rep(Ytest[, variable_name], n_models), ncol = n_models)

loss_matrix <- (pred_matrix - true_matrix) ^ 2

mcset <- MCSprocedure(loss_matrix, k = 100, 
                      alpha = 0.5, statistic = "Tmax")
# смысл полей:
# _M --- статистика Tmax
# Rank_M -- номер v_M при сортировке по возрастанию
# v_M --- значение статистики 
# MCS_M --- ???
# _R --- статистика TR
# Rank_R -- номер v_R при сортировке по возрастанию
# v_R --- значение статистики
# MCS_R --- ???
# Loss -- это средний лось модели



loss_df <- as.data.frame(loss_matrix)

rownames(loss_df) <- NULL
colnames(loss_df) <- paste0("model_", 1:3)

# write_csv(loss_df, path = "~/Downloads/loss_matrix.csv")
loss_matrix <- read_csv("~/Downloads/loss_matrix.csv")

set.seed(7)
mcset <- MCSprocedure(loss_matrix, statistic = "TR")

data(Loss)
builtin_loss <- Loss
# MCSprocedure(builtin_loss) # too long ~45 min
