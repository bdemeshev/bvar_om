
source("400_model_funs.R")
source("400_model_lists.R")
source("500_replicate_banbura_function.R")

#### test 1: estimate big model (23 vars, 12 lags)


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- create_var_set_info()
#var_set_info <- var_set_info[-11,]

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_test.csv")

mlist <- read_csv("../estimation/mlist_test.csv")

var_set_info <- create_var_set_info()
model_info <- mlist %>% filter(id == 3)
model_info
estimate_model(model_info, df = df, var_set_info = var_set_info, verbose = TRUE)



#### test 3: RW and WN prediction test in-sample


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- create_var_set_info()

mlist <- create_rwwn_list()
mlist <- estimate_models(mlist, df = df, var_set_info = var_set_info, parallel = "off")
plist <- data.frame(model_id = c(1, 2), h = NA, type = "in-sample")
mlist

pred_info <- plist[1, ]
pred_info
forecast_model(pred_info, mlist, var_set_info = var_set_info, df = df)

pred_info <- plist[2,]
pred_info
forecast_model(pred_info, mlist, var_set_info = var_set_info, df = df)

#### test 4: VAR prediction test in-sample

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- create_var_set_info()

# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list, df = df, var_set_info = var_set_info, parallel = "off")

# forecast VAR
var_forecast_list <- data.frame(model_id = unique(var_list$id), h = NA, type = "in-sample")

pred_info <- var_forecast_list[1, ]
pred_info
forecast_model(pred_info, mlist, df = df, var_set_info = var_set_info)

#### test 5: RW and WN prediction test out-of-sample


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- create_var_set_info()

mlist <- create_rwwn_list()
mlist <- estimate_models(mlist, df = df, var_set_info = var_set_info)
plist <- data.frame(model_id = c(1, 2), h = 12, type = "out-of-sample")
mlist

pred_info <- plist[1, ]
pred_info
res <- forecast_model(pred_info, mlist, df = df, var_set_info = var_set_info)
res

pred_info <- plist[2, ]
pred_info
res <- forecast_model(pred_info, mlist,  df = df, var_set_info = var_set_info)
res

#### test 6: VAR prediction test out-of-sample


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- create_var_set_info()

# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list, parallel = "off", df = df, var_set_info = var_set_info)

# forecast VAR
var_forecast_list <- data.frame(model_id = unique(var_list$id), h = 12, type = "out-of-sample")

pred_info <- var_forecast_list[1, ]
pred_info
res <- forecast_model(pred_info, var_list, df = df, var_set_info = var_set_info)
res

#### test 7: BVAR estimation error


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

bvar_list <- create_bvar_banbura_list()

var_set_info <- create_var_set_info()
model_info <- bvar_list %>% filter(id == 154)
model_info
estimate_model(model_info, df = df, var_set_info = var_set_info, verbose = TRUE)


#### test 8: BVAR forecasting error suspected for n_lag=1, h>1

?forecast_conjugate

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2)
model <- bvar_conjugate0(priors = priors, keep = 1000, fast_forecast = TRUE)
forecast_conjugate(model, h = 5, output = "wide")

#### test 9: BVAR posterior Phi mean: sample vs analytical


?forecast_conjugate

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2)
model <- bvar_conjugate0(priors = priors, keep = 1000, fast_forecast = FALSE)
summary_conjugate(model)
forecast_conjugate(model, h = 5, output = "wide")

#### test 9+: BVAR coefficient/forecast test: fast_forecast=TRUE vs way="cholesky"/"svd"
# everything should be equal in summaries! both mean and sd


data(Yraw)
priors <- Carriero_priors(Yraw, p = 3)
model_chol <- bvar_conjugate0(priors = priors, keep = 10000, fast_forecast = FALSE, 
                              way_omega_post_root = "cholesky")
model_svd <- bvar_conjugate0(priors = priors, keep = 10000, fast_forecast = FALSE,
                             way_omega_post_root = "svd")
model_fast <- bvar_conjugate0(priors = priors, keep = 10000, fast_forecast = TRUE,
                              way_omega_post_root = "svd")

# compare coef mean (sample vs theoretical), coef sd (svd/chol) :
summary_conjugate(model_chol)
summary_conjugate(model_svd)
summary_conjugate(model_fast)

# summary_conjugate(model)
a <- forecast_conjugate(model_chol, h = 12, include = "mean", level = NULL)
b <- forecast_conjugate(model_svd, h = 12, include = "mean", level = NULL)
c <- forecast_conjugate(model_fast, h = 12,  fast_forecast = TRUE)

a <- rename(a, value_chol = value) %>% select(-what)
b <- rename(b, value_svd = value) %>% select(-what)
c <- rename(c, value_fast = value) %>% select(-what)

abc <- left_join(a,b) %>% left_join(c)

abc

## test 10: what is faster chol or sqrtm. 

library("microbenchmark")
library("expm")

H <- matrix(c(2, 1, 1, 2), nrow = 2)

microbenchmark(chol(H), sqrtm(H))

chol(H)
sqrtm(H)
# chol(H) is faster

t(chol(H)) %*% chol(H)
sqrtm(H) %*% sqrtm(H)
