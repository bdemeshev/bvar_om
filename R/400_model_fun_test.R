
source("400_model_funs.R")

#### test 1: estimate big model (23 vars, 12 lags)


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")
#var_set_info <- var_set_info[-11,]

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_test.csv")

mlist <- read_csv("../estimation/mlist_test.csv")
model_info <- mlist %>% filter(id==3)
model_info
estimate_model(model_info)

#### test 2: estimate big model (23 vars, 12 lags) without pairs of variables

df <- read_csv("../data/df_2015_final.csv")

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_test.csv")

mlist <- read_csv("../estimation/mlist_test.csv")
model_info <- mlist %>% filter(id==3)

# тестируем выбрасывание переменной
for (j in 10:31) {
  for (i in (j+1):32) {
    var_set_info <- read_csv("../data/var_set_info.csv")
    var_set_info <- var_set_info[-c(i,j),]
    estimate_model(model_info)
  }
}

# восстановили нормальный список переменных для каждого набора 
var_set_info <- read_csv("../data/var_set_info.csv")

#### test 3: RW and WN prediction test in-sample

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

mlist <- create_rwwn_list()
mlist <- estimate_models(mlist,parallel = parallel)
plist <- data.frame(model_id=c(1,2), h=NA, type="in-sample")
mlist

pred_info <- plist[1,]
pred_info
forecast_model(pred_info, mlist, parallel = "off")

pred_info <- plist[2,]
pred_info
forecast_model(pred_info, mlist, parallel = "off")

#### test 4: VAR prediction test in-sample
df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list,parallel = parallel)

# forecast VAR
var_forecast_list <- data.frame(model_id=unique(var_list$id), h=NA, type="in-sample")

pred_info <- var_forecast_list[1,]
pred_info
forecast_model(pred_info, mlist, parallel = "off")

#### test 5: RW and WN prediction test out-of-sample

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

mlist <- create_rwwn_list()
mlist <- estimate_models(mlist,parallel = parallel)
plist <- data.frame(model_id=c(1,2), h=12, type="out-of-sample")
mlist

pred_info <- plist[1,]
pred_info
res <- forecast_model(pred_info, mlist, parallel = "off")
res

pred_info <- plist[2,]
pred_info
res <- forecast_model(pred_info, mlist, parallel = "off")
res

#### test 6: VAR prediction test out-of-sample
df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list,parallel = parallel)

# forecast VAR
var_forecast_list <- data.frame(model_id=unique(var_list$id), h=12, type="out-of-sample")

pred_info <- var_forecast_list[1,]
pred_info
res <- forecast_model(pred_info, var_list, parallel = "off")
res

#### test 7: BVAR estimation error
df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

bvar_list <- create_bvar_banbura_list()


model_info <- bvar_list %>% filter(id==154)
model_info
estimate_model(model_info, test=TRUE)
estimate_model(model_info)

#### test 8: BVAR forecasting error suspected for n_lag=1, h>1
?forecast_conjugate

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2)
model <- bvar_conjugate0(priors = priors, keep=1000, fast_forecast = TRUE)
forecast_conjugate(model, h=5, output="wide")

#### test 9: BVAR posterior Phi mean: sample vs analytical
?forecast_conjugate

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2)
model <- bvar_conjugate0(priors = priors, keep=1000, fast_forecast = FALSE)
summary_conjugate(model)
forecast_conjugate(model, h=5, output="wide")

#### test 9: BVAR forecast test: fast_forecast TRUE/FALSE

data(Yraw)
priors <- Carriero_priors(Yraw, p = 2)
model <- bvar_conjugate0(priors = priors, keep=1000, fast_forecast = FALSE)
# summary_conjugate(model)
forecast_conjugate(model, h=12, include="mean", level=NULL)
forecast_conjugate(model, h=12,  fast_forecast = TRUE)
