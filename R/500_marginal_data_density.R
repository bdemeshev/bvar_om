# 500_marginal_data_density.R

# input: "../data/df_2015_final.csv"
# output: 

source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")

# need to run only once 
source("200_load_after_eviews.R")

##################################################################
######################## begin set-up part #######################
##################################################################


df <- read_csv("../data/df_2015_final.csv")

# subset t:
df <- tail(df, -12) # tail(df, -12) # all obs except first twelve

# parallel computing details:
parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30 # number of cores for paralel computing, ignored if parallel=="off"


h_max <- 12 # maximum forecast horizont for VAR and BVAR 

T_available <- nrow(df) # number of observations

# posterior simulation details:
fast_forecast <- TRUE # TRUE = posterior means of coefficients are used for forecast
keep <- 0 # 5000 # number of simulations from posterior (used only if fast_forecast is FALSE)
verbose <- TRUE # turn on/off messages from functions 

num_AR_lags <- 1 # number of lags in AR() model used to estimate sigma^2 
# if num_AR_lags <- NULL then p will be used

# testing mode (less lambdas are estimated, see 400_model_lists.R)
testing_mode <- TRUE

# use wrong formulas from carriero code for dummy cNIW without square root for sigma^2
carriero_hack <- FALSE

num_AR_lags <- NULL # NULL # number of lags in AR() model used to estimate sigma^2 
# if num_AR_lags <- NULL then p will be used

set_delta_by <- "ADF" # "ADF", "KPSS" or "global AR1" or "AR1" or a number

v_prior <- "m+2" # "m+2" / "T_dummy" / константа

################################################

# a lot of models are estimated but only some are reported

desired_n_lag <- 12
desired_h <- c(1,3,6,12)

###############################################
if (parallel=="unix") {
  library("doParallel")
  # library("doMC")
  # registerDoMC(ncpu) # number of CPU cores
}
if (parallel=="windows") {
  library("doSNOW")
}


################################################
# create var_set_info
# describe which variables are included in each set
# dput(colnames(df))

add_3 <- data_frame(var_set="set_03", variable=c("ind_prod", "cpi", "ib_rate"))
add_6 <- data_frame(var_set="set_06", variable=c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))
add_23 <- data_frame(var_set="set_23",variable=c("employment", 
                                                 "ind_prod", 
                                                 "cpi", 
                                                 "ib_rate", "lend_rate", "real_income", 
                                                 "unemp_rate", "oil_price", 
                                                 "ppi", 
                                                 "construction", "real_investment", 
                                                 "wage", 
                                                 "m2", 
                                                 "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                                                 "agriculture", "retail", "gov_balance", 
                                                 "export", 
                                                 "import"
) )

var_set_info <- rbind(add_3,add_6,add_23)

# write_csv(var_set_info,"../data/var_set_info.csv")

##################################################################
######################## end set-up part #########################
##################################################################

####### step 0 (before banbura procedure)
# melting actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value) %>% 
  mutate(variable=as.character(variable))


##### banbura step 1
# calculate msfe-0. Estimate RWWN (random walk OR white noise model)

# set delta (prior hyperparameter)

if (set_delta_by %in% c("ADF","KPSS")) {
  deltas <- delta_i_prior(df, test = set_delta_by,
                          remove_vars = c("time_y","t"), c_0 = 0, c_1 = 1)
} 

if (set_delta_by == "global AR1") {
  deltas <- delta_i_from_ar1(df, remove_vars = c("time_y","t"))
}

if ((set_delta_by == "AR1") | is.finite(set_delta_by) ) { 
  # if "AR1" then will be calculated automatically by bvar_conj_setup
  deltas <- data_frame(
    variable = setdiff(colnames(df),c("t","time_y")),
    delta = set_delta_by)
}


# Important!!!
# value of delta determines prior
# RW/WN determines model for comparison

deltas <- mutate(deltas, 
                 rw_wn = "rw", # always compare with Random Walk
                 variable=as.character(variable))
deltas # 

# estimate all RW and WN models
rwwn_list <- create_rwwn_list()
message("Estimate RW and WN")
rwwn_list <- estimate_models(rwwn_list,parallel = parallel, ncpu=ncpu, verbose=verbose)

# forecast all RW and WN models
rwwn_forecast_list <- data_frame(model_id=c(1,2), h=NA, type="in-sample")
message("Forecast RW and WN")
rwwn_forecasts <- forecast_models(rwwn_forecast_list, rwwn_list,
                                  parallel = parallel, ncpu=ncpu, verbose=verbose)

# calculate all msfe-0
# two ways to calculate msfe (we use b)
# a) use all available predictions for each model
# b) use only common available predictions for all models

rwwn_forecasts %>% group_by(model_id) %>% summarise(Tf_start=min(t),Tf_end=max(t))

plot_forecast(rwwn_forecasts, var_name="m2", mod_id=2)


# build table with corresponding msfe-0 (RW or WN)
msfe0_all <- get_msfe(rwwn_forecasts, actual_obs, 
                      models = dplyr::select(rwwn_list, id, type),
                      msfe_name = "msfe")


# add deltas info to msfe0_all

msfe0 <- left_join(deltas, msfe0_all, by= c("variable"="variable","rw_wn"="type") )

msfe0

##### banbura step 2

# estimate VAR
var_list <- create_var_list()
message("Estimating VAR")
var_list <- estimate_models(var_list,parallel = parallel, ncpu=ncpu, verbose=verbose)

# forecast VAR
var_forecast_list <- data_frame(model_id=var_list$id, h=NA, type="in-sample")
message("Forecasting VAR")
var_forecasts <- forecast_models(var_forecast_list, var_list, 
                                 parallel = parallel, ncpu=ncpu, verbose=verbose)


# calculate msfe-inf
msfe_Inf_info <- get_msfe(var_forecasts, actual_obs, 
                          models = select(var_list, id, type, var_set, n_lag),
                          msfe_name = "msfe_Inf")


msfe_Inf_info %>% head()

##### step 2 (banbura step 6): calculate omsfe for RW/WN/VAR
# look at lists:
var_list
rwwn_list


rwwn_var_unique_list <- bind_rows(var_list, rwwn_list) # %>% mutate_each("as.numeric", n_lag, T_in, T_start)

# every model should be rolled
rwwn_var_out_list <- rolling_model_replicate(rwwn_var_unique_list) %>% 
  mutate(file=paste0(type,"_",id,
                     "_T_",T_start,"_",T_in,"_",
                     var_set,"_lags_",n_lag,".Rds") ) 
message("Estimating rolling rw/wn/var models")
rwwn_var_out_list <- estimate_models(rwwn_var_out_list,parallel = parallel, ncpu=ncpu, verbose=verbose) # takes some minutes

rwwn_var_out_forecast_list <- rwwn_var_out_list %>% rowwise() %>% mutate(model_id=id, 
                                                                         h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")
# without rowwise min will be global and always equal to 1

# forecast all rolling models  
message("Forecasting rolling rw/wn/var models, out-of-sample")
rwwn_var_out_forecasts <- forecast_models(rwwn_var_out_forecast_list, rwwn_var_out_list, progress_bar=TRUE,
                                          parallel = parallel, ncpu=ncpu, verbose=verbose)



omsfe_rwwn_var_table <- get_msfe(rwwn_var_out_forecasts, actual_obs,
                                 models = select(rwwn_var_out_list, id, var_set, n_lag, model_type=type),
                                 plus_group_vars = "model_type",
                                 msfe_name = "omsfe", msfe_type = "out-of-sample")

# select rw or wn omsfe for each variable, for each h
omsfe_selected_rwwn <- left_join(deltas, 
                                 filter(omsfe_rwwn_var_table, model_type %in% c("wn","rw")),
                                 by=c("rw_wn"="model_type","variable"="variable")) %>% 
  select(-var_set, -n_lag, -delta)
omsfe_selected_rwwn %>% glimpse()

omsfe_rwwn_banbura <- omsfe_selected_rwwn %>% 
  # filter(variable %in% filter_variables ) %>%
  rename(model_type=rw_wn) 


#### estimate bvar model and calculate mdd for each t 

# all models for the same moment of time
bvar_list_zero <- create_mdd_list() 

# all models for rolling t
bvar_list <- rolling_model_replicate(bvar_list_zero) %>%
   mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                         var_set,"_lags_",n_lag,
                         "_lams_",round(100*l_1),
                         "_sc_",round(100*l_sc),
                         "_io_",round(100*l_io),
                         "_pow_",round(100*l_power),
                         "_cst_",round(100*l_const),
                         ".Rds") ) 

message("Estimating rolling BVAR")
bvar_list <- estimate_models(bvar_list, parallel = parallel, ncpu=ncpu, verbose=verbose) # status and filename are updated

message("Calculating mdd for rolling BVAR")
bvar_list <- calculate_mdd(bvar_list)

#### for each t find optimal model (maybe we'll add grouping variable)

optimal_by <- c("Tf_start") #, var_set, n_lag

best_bvars <- bvar_list %>% mutate(Tf_start=T_start+n_lag) %>% 
  group_by_(.dots=optimal_by) %>% 
  mutate(mdd_rank=dense_rank(mdd)) %>% filter(mdd_rank == 1) %>% ungroup()

#### check for non-unique optimal mdd
check_uniqueness <- best_bvars %>% group_by_(.dots=optimal_by) %>% 
  summarise(num_of_best_mdd=n())
if (max(check_uniqueness$num_of_best_mdd)>1) {
  warning("**** ACHTUNG ****: non unique optimal mdd")
  check_uniqueness %>% filter(num_of_best_mdd>1)
  message("First in list optimal mdd will be chosen")
  best_bvars <- best_bvars %>% group_by_(.dots=optimal_by) %>%
    mutate(rownum = row_number()) %>% filter(rownum==1) %>% select(-rownum)
}


#### forecast all best BVAR models

best_bvars_forecast_list <- best_bvars %>% rowwise() %>% 
  mutate(model_id=id, h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")

message("Forecasting best BVAR, out-of-sample")
best_bvars_forecasts <- forecast_models(best_bvars_forecast_list, best_bvars, 
                                        parallel = parallel, ncpu=ncpu, verbose=verbose)

# проверить, по каким переменным группировать в get_msfe!!!!!!!!!!!
omsfe_best_mdd <- get_msfe(best_bvars_forecasts, actual_obs,
                             models = select(best_bvars, id, var_set, n_lag),
                             msfe_name = "omsfe", msfe_type = "out-of-sample")
omsfe_best_mdd %>% head()

#### calculate rmsfe
rmsfe_long <- omsfe_best_mdd %>%
  left_join(omsfe_rwwn_banbura %>% select(-model_type, omsfe_rwwn=omsfe), 
            by=c("h","variable")) %>%
  mutate(rmsfe=omsfe/omsfe_rwwn, model_type="rmsfe_mdd") 

# how to automate depending on optimal_by??????????
rmsfe_wide <- rmsfe_long %>% select(-omsfe,-omsfe_rwwn, -var_set) %>% 
  dcast(h+variable+n_lag~model_type, value.var="rmsfe") 


some_rmsfe_wide <- rmsfe_wide %>% filter(h %in% desired_h)

some_rmsfe_wide



