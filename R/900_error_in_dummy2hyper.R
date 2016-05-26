library("bvarr")

model <- readRDS("../estimation/bad_model.Rds")

# throws error
mdd <- bvar_conj_mdd(model)

# NaN
prior <- bvar_conj_dummy2hyper(model$Y_plus, model$X_plus)
prior
X_star <- model$X_plus
Y_star <- model$Y_plus



Z_in = NULL 
constant = TRUE
p = 4 
lambda = c(0.2,1,1,1,100,100)
delta = 1
s2_lag = NULL
y_bar_type = c("initial")
v_prior = NULL
carriero_hack = FALSE
