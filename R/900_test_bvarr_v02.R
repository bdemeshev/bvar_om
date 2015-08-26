# test bvarr-v02

devtools::install_github("bdemeshev/bvarr")

library("bvarr")
library("dplyr")

data(Yraw)

old <- Carriero_priors(Yraw, p=2)
X_star <- old$X_dummy_cniw 
Y_star <- old$Y_dummy_cniw 
hyper_prior <- bvar_dummy2hyper(old$Y_dummy_cniw, old$X_dummy_cniw)

d2 <- bvar_conj_lambda2dummy(Yraw,p=2)
hp2 <- bvar_dummy2hyper(d2$Y_cniw, d2$X_cniw)
d2$Y_cniw
d2$X_cniw
# check Omega
hyper_prior$Omega
hp2$Omega
old$Omega_prior


# check S
hyper_prior$S
hp2$S
old$S_prior

# check Phi
hyper_prior$Phi
hp2$Phi
old$Phi_prior


# check hyper -> dummy -> hyper cycle

hp2$Omega
hp2$S
hp2$Phi

dum <- bvar_hyper2dummy(hp2$Omega, hp2$S, hp2$Phi)

dum$X_plus
dum$Y_plus

hp2b <- bvar_dummy2hyper(dum$Y_plus, dum$X_plus)


all.equal(hp2b$Omega, hp2$Omega)
all.equal(hp2b$S, hp2$S)
all.equal(hp2b$Phi, hp2$Phi)

####

data(Yraw)
dummy <- bvar_conj_lambda2dummy(Yraw,p=2)
hyper <- bvar_dummy2hyper(dummy$Y_cniw, dummy$X_cniw)
# these are priors but just for testing let's pretend that they are posteriors
post_sample <- bvar_conj_simulate(v_post=10, hyper$Omega_root, hyper$S, hyper$Phi)
hyper$Omega_root
hyper$S
hyper$Phi

post_sample
