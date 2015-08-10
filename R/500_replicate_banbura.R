# 500_replicate_banbura.R

# input: "../data/df_2015_final.csv", "../data/var_set_info.csv"
# output: 



library("foreach")

source("400_model_funs.R")

parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30

df <- read_csv("../data/df_2015_final.csv")

##### banbura step 1
# calculate MSFE-0. Estimate RWWN (random walk OR white noise model)



# list of RW variables

# this function detects series with unit roots
# for ADF all three types are supported
# for KPSS only trend and constant
# c_0=0 value for stationary series
# c_1=1 value for non-stationary series
# alpha=5%
delta_i_prior <- function(df_ts, varset=colnames(df_ts), test=c("ADF","KPSS"),
                        type=c("trend", "constant", "neither"), c_0=0, c_1=1) {
  df_sel <- df_ts[,varset]
  nvar <- ncol(df_sel)
  stat_res <- stationarity(df_sel, print = FALSE)
  
  type <- match.arg(type) # take first option if not specified
  test <- match.arg(test)

  if (type=="trend") line <- 1
  if (type=="constant") line <- 2
  if (type=="neither") line <- 3
  
  nvar <- ncol(df_sel)
  
  # use time_trend - ADF - 5%
  if(test=="ADF") {
    adf_cr <- stat_res$ADF$`5 Pct`[line]
    adf_obs <- stat_res$ADF[line,1:nvar]  
    result <- ifelse(adf_obs < adf_cr, c_0, c_1)
  }
  
  if(test=="KPSS") {
    kpss_cr <- stat_res$KPSS$`5 Pct`[line]
    kpss_obs <- stat_res$KPSS[line,1:nvar]  
    result <- ifelse(kpss_obs > kpss_cr, c_1, c_0)
  }
  
  return(as.vector(result))
}




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