# 500_replicate_banbura.R

# input: "../data/df_2015_final.csv", "../data/var_set_info.csv"
# output: 



library("foreach")

source("400_model_funs.R")
source("500_banbura_funs.R")

parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

##### banbura step 1
# calculate MSFE-0. Estimate RWWN (random walk OR white noise model)

# classify variables into RW and WN
deltas <- delta_i_prior(df, remove_vars = "time_y")
# delta_i_from_ar1(df, remove_vars = "time_y")

deltas <- mutate(deltas, rw_wn = ifelse(delta==1,"rw","wn"))
deltas

# estimate all RW and WN models
rwwn_list <- create_rwwn_list()
rwwn_list <- estimate_models(rwwn_list,parallel = parallel)

# forecast all RW and WN models
rwwn_forecast_list <- data.frame(model_id=c(1,2), h=NA, type="in-sample")
rwwn_forecasts <- forecast_models(rwwn_forecast_list, rwwn_list)

# calculate all MSFE-0
# two ways to calculate MSFE
# a) use all available predictions for each model
# b) use only common available predictions for all models

rwwn_forecasts %>% group_by(model_id) %>% summarise(Tf_start=min(t),Tf_end=max(t))

#### if we prefer way b, just uncomment next two lines
# Tf_common_start <- 2
# rwwn_forecasts <- filter(rwwn_forecasts, t >= Tf_common_start)
####

# build table with corresponding MSFE-0 (RW or WN)

# joining actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

rwwn_forecasts <- rename(rwwn_forecasts, forecast=value)
rwwn_obs <- left_join(rwwn_forecasts, actual_obs, by=c("t","variable"))

rwwn_obs <- mutate(rwwn_obs, sq_error=(forecast-actual)^2)
head(rwwn_obs)

MSFE0_all <- rwwn_obs %>% group_by(variable, model_id) %>% summarise(MSFE=mean(sq_error))
MSFE0_all

# add rw/wn label to MSFE0
rwwn_wlist <- dcast(rwwn_list, id~variable) # wlist = wide list
MSFE0_all <- left_join(MSFE0_all, select(rwwn_wlist, id, type), by=c("model_id"="id") )

MSFE0 <- left_join(deltas, MSFE0_all, by= c("variable"="variable","rw_wn"="type") )

MSFE0

##### banbura step 2
# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list,parallel = parallel)

# calculate MSFE-inf, FIT-inf


##### banbura step 3
# goal: calculate FIT-lam

# create model list to find optimal lambda 
mlist <- create_model_list_banbura()
write_csv(mlist, path = "../estimation/mlist_optimise_banbura.csv")

# estimate models from list
mlist <- read_csv("../estimation/mlist_optimise_banbura.csv")
mlist <- estimate_models(mlist, parallel = parallel, ncpu=ncpu) # status and filename are updated
write_csv(mlist, path = "../estimation/mlist_optimise_banbura.csv")

# calculate MSFE-lam

# calculate FIT-lam

##### banbura step 4
# find optimal lambda


##### banbura step 5
# forecast and evaluate using optimal lambda