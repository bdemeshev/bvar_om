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

deltas

# estimate all RW and WN models
rwwn_list <- create_rwwn_list()
rwwn_list <- estimate_models(rwwn_list,parallel = parallel)

# forecast all RW and WN models

# calculate all MSFE-0

# build table with corresponding MSFE-0 (RW or WN)



##### banbura step 2
# estimate VAR
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