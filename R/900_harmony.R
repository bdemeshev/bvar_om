# нахождение гармонии между формулами априорное-апостериорное и дамми cNIW

library("bvarr")
library("dplyr")

?bvar_conjugate0

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2, lambdas = c(0.2, 1.0, NA, NA, 100.0, 100.0))
X_dummy_cniw <- priors$X_dummy_cniw
Y_dummy_cniw <- priors$Y_dummy_cniw
Omega_prior <- priors$Omega_prior
Phi_prior <- priors$Phi_prior

# Omega_prior^{-1}
sym_inv(Omega_prior) %>% is.diagonal()
crossprod(X_dummy_cniw) %>% is.diagonal()

sym_inv(Omega_prior) %>% diag()
crossprod(X_dummy_cniw) %>% diag()

# Phi_prior 
Phi_prior
solve(crossprod(X_dummy_cniw), crossprod(X_dummy_cniw, Y_dummy_cniw))

model <- bvar_conjugate0(priors = priors, fast_forecast = TRUE)

X_wo_dummy <- attr(model, "data")$X_wo_dummy
Y_wo_dummy <- attr(model, "data")$Y_wo_dummy
X_star <- rbind(X_dummy_cniw, X_wo_dummy)
Y_star <- rbind(Y_dummy_cniw, Y_wo_dummy)

# Phi_post
attr(model,"post")$Phi_post
solve(crossprod(X_star), crossprod(X_star, Y_star))

# Omega_post
attr(model,"post")$Omega_post
sym_inv(crossprod(X_star))

# wo = without dummies 
Phi_wo <- solve(crossprod(X_wo_dummy),crossprod(X_wo_dummy,Y_wo_dummy))
Phi_star <- solve(crossprod(X_star),crossprod(X_star,Y_star))
Phi_cniw <- solve(crossprod(X_dummy_cniw),crossprod(X_dummy_cniw,Y_dummy_cniw))
E_wo <- Y_wo_dummy - X_wo_dummy %*% Phi_wo
E_cniw_wo <- Y_dummy_cniw - X_dummy_cniw %*% Phi_wo
E_star_wo <- Y_star - X_star %*% Phi_wo
E_star <- Y_star - X_star %*% Phi_star
E_cniw <- Y_dummy_cniw - X_dummy_cniw %*% Phi_cniw

# S_prior
priors$S_prior
crossprod(E_cniw)
# no analog?

# S_post
attr(model,"post")$S_post
crossprod(E_star)

crossprod(E_wo)
crossprod(E_cniw_wo)
crossprod(E_star_wo)

