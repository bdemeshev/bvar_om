# 500_replicate_banbura.R

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


# parallel computing details:
parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30 # number of cores for paralel computing, ignored if parallel=="off"

df <- read_csv("../data/df_2015_final.csv")

h_max <- 12 # maximum forecast horizont for VAR and BVAR 

T_available <- nrow(df) # number of observations

# posterior simulation details:
fast_forecast <- TRUE # TRUE = posterior means of coefficients are used for forecast
keep <- 0 # 5000 # number of simulations from posterior (used only if fast_forecast is FALSE)
verbose <- FALSE # turn on/off messages from functions 

# testing mode (less lambdas are estimated, see 400_model_lists.R)
testing_mode <- FALSE

num_AR_lags <- 1 # number of lags in AR() model used to estimate sigma^2 
# if num_AR_lags <- NULL then p will be used

################################################
# create fit_set_info
# describe which msfe ratios are averaged in fit
fit_set_2vars <- data_frame(variable=c("ind_prod","cpi"),fit_set="ind+cpi")
fit_set_3vars <- data_frame(variable=c("ind_prod","cpi","ib_rate"),fit_set="ind+cpi+rate")

fit_set_info <- bind_rows(fit_set_2vars, fit_set_3vars)
fit_set_info

# a lot of models are estimated but only some are reported
desired_fit_set <- "ind+cpi+rate"
desired_n_lag <- 12
desired_h <- c(1,3,6,12)

################################################
# create var_set_info
# describe which variables are included in each set
dput(colnames(df))

add_3 <- data_frame(var_set="set_3", variable=c("ind_prod", "cpi", "ib_rate"))
add_6 <- data_frame(var_set="set_6", variable=c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))
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

var_set_info <- bind_rows(add_3,add_6,add_23)

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

# classify variables into RW and WN
deltas <- delta_i_prior(df, remove_vars = c("time_y","t"), c_0 = 0, c_1 = 1)
# deltas <- delta_i_from_ar1(df, remove_vars = c("time_y","t"))

deltas <- mutate(deltas, rw_wn = ifelse(delta==1,"rw","wn"), variable=as.character(variable))
deltas

# estimate all RW and WN models
rwwn_list <- create_rwwn_list()
message("Estimate RW and WN")
rwwn_list <- estimate_models(rwwn_list,parallel = parallel)

# forecast all RW and WN models
rwwn_forecast_list <- data_frame(model_id=c(1,2), h=NA, type="in-sample")
message("Forecast RW and WN")
rwwn_forecasts <- forecast_models(rwwn_forecast_list, rwwn_list)

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
var_list <- estimate_models(var_list,parallel = parallel)

# forecast VAR
var_forecast_list <- data_frame(model_id=var_list$id, h=NA, type="in-sample")
message("Forecasting VAR")
var_forecasts <- forecast_models(var_forecast_list, var_list)


# calculate msfe-inf
msfe_Inf_info <- get_msfe(var_forecasts, actual_obs, 
                     models = select(var_list, id, type, var_set, n_lag),
                     msfe_name = "msfe_Inf")


msfe_Inf_info %>% head()

# join msfe0
msfe_0_Inf <- left_join(msfe_Inf_info,
                        msfe0 %>% select(msfe,rw_wn,variable) %>% rename(msfe0=msfe),
                        by="variable")

msfe_0_Inf <- msfe_0_Inf %>% select(-model_id,-type) %>% 
  mutate(msfe_ratio=msfe_Inf/msfe0)
msfe_0_Inf


# create empty table
fit_inf_table <- NULL

# cycle all fit_sets:
for (current_fit_set in unique(fit_set_info$fit_set)) {
  fit_variables <- filter(fit_set_info, fit_set==current_fit_set)$variable
  block <- filter(msfe_0_Inf, variable %in% fit_variables) %>% group_by(var_set,n_lag) %>%
    summarise(fit_inf=mean(msfe_ratio)) %>% mutate(fit_set=current_fit_set)
  fit_inf_table <- bind_rows(fit_inf_table,block)
}

fit_inf_table

##### banbura step 3
# goal: calculate fit-lam

# create model list to find optimal lambda 
bvar_list_ <- create_bvar_banbura_list()
# write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# estimate models from list
# bvar_list <- read_csv("../estimation/bvar_list.csv")
message("Estimating BVAR")
bvar_list <- estimate_models(bvar_list_, parallel = parallel, ncpu=ncpu, test=FALSE) # status and filename are updated
# write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# forecast BVAR
bvar_forecast_list <- data_frame(model_id=bvar_list$id, h=NA, type="in-sample")
message("Forecasting BVAR in-sample")
bvar_forecasts <- forecast_models(bvar_forecast_list, bvar_list)
message("Forecasting BVAR in-sample ok")


# calculate msfe-lam
msfe_lam_info <- get_msfe(bvar_forecasts, actual_obs, 
                     models = select(bvar_list, id, type, var_set, n_lag,
                                     l_1, l_const, l_io, l_power, l_sc, n_lag),
                     msfe_name = "msfe_lam")


msfe_lam_info %>% head()

# join msfe0
msfe_0_lam <- left_join(msfe_lam_info,msfe0 
                        %>% select(msfe,rw_wn,variable) %>% rename(msfe0=msfe),
                        by="variable")

msfe_0_lam <- msfe_0_lam %>% select(-model_id,-type) %>% 
  mutate(msfe_ratio=msfe_lam/msfe0)
msfe_0_lam

# calculate fit-lam
# create empty table
fit_lam_table <- NULL

# cycle all fit_sets:
for (current_fit_set in unique(fit_set_info$fit_set)) {
  fit_variables <- filter(fit_set_info, fit_set==current_fit_set)$variable
  block <- filter(msfe_0_lam, variable %in% fit_variables) %>% 
    group_by(var_set,n_lag,l_1,l_const,l_io,l_power,l_sc) %>%
    summarise(fit_lam=mean(msfe_ratio)) %>% mutate(fit_set=current_fit_set)
  fit_lam_table <- bind_rows(fit_lam_table,block)
}

fit_lam_table

##### banbura step 4
# find optimal lambda
# ungroup() is needed! otherwise cannot remove var_set (active group on var_set)
fit_goal <- dplyr::filter(fit_inf_table, var_set=="set_3") %>% ungroup() %>% dplyr::select(-var_set)
fit_goal

fit_lam_table <- left_join(fit_lam_table, fit_goal, by=c("n_lag","fit_set"))
fit_lam_table <- mutate(fit_lam_table, delta_fit = abs(fit_lam-fit_inf))

# for each var_set, n_lag, fit_set the best lambdas are calculated
optimal_by <- c("var_set", "n_lag", "fit_set")
best_lambda <- fit_lam_table %>% group_by_(.dots=optimal_by) %>% 
  mutate(fit_rank=dense_rank(delta_fit)) %>% filter(fit_rank == 1) %>% ungroup()
best_lambda %>% arrange(var_set, n_lag, fit_set)



# check whether best lambda is unique
check_uniqueness <- best_lambda %>% group_by_(.dots=optimal_by) %>% summarise(num_of_best_lambdas=n())
# num of best lambdas should be always one
if (max(check_uniqueness$num_of_best_lambdas)>1) {
  warning("**** ACHTUNG ****: non unique optimal lambdas")
  check_uniqueness %>% filter(num_of_best_lambdas>1)
  message("First in list optimal lambda will be chosen")
  best_lambda <- best_lambda %>% group_by_(.dots=optimal_by) %>%
    mutate(rownum = row_number()) %>% filter(rownum==1) %>% select(-rownum)
  
}


##### banbura step 5: calculate omsfe for bvar
# forecast and evaluate using optimal lambda

# create best models lists with correct time spec
# we need to keep fit_set variable for further comparison
bvar_out_list  <- create_bvar_out_list(best_lambda)
# best_lambda table should have: var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, fit_set

message("Estimate rolling BVAR")
bvar_out_list <- estimate_models(bvar_out_list,parallel = parallel)
#write_csv(bvar_out_list, path = "../estimation/bvar_out_list.csv")

# forecast BVAR
# check structure 
bvar_out_list %>% group_by(var_set, fit_set, n_lag) %>% summarise(n=n())

bvar_out_forecast_list <- bvar_out_list %>% rowwise() %>% mutate(model_id=id, 
                                 h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")
# without rowwise min will be global and always equal to 1

# a lot of forecasts
### WARNING: maybe too many observations!!!
message("Forecasting rolling BVAR, out-of-sample")
bvar_out_forecasts <- forecast_models(bvar_out_forecast_list, bvar_out_list)
message("Forecasting rolling BVAR out-of-sample ok")

omsfe_bvar_table <- get_msfe(bvar_out_forecasts, actual_obs,
                             models = select(bvar_out_list, id, var_set, n_lag, fit_set),
                             plus_group_vars = "fit_set",
                             msfe_name = "omsfe", msfe_type = "out-of-sample")
omsfe_bvar_table %>% head()

##### banbura step 6: calculate omsfe for RW/WN/VAR
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
rwwn_var_out_list <- estimate_models(rwwn_var_out_list,parallel = parallel) # takes some minutes

rwwn_var_out_forecast_list <- rwwn_var_out_list %>% rowwise() %>% mutate(model_id=id, 
  h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")
# without rowwise min will be global and always equal to 1

# forecast all rolling models  
message("Forecasting rolling rw/wn/var models, out-of-sample")
rwwn_var_out_forecasts <- forecast_models(rwwn_var_out_forecast_list, rwwn_var_out_list, progress_bar=TRUE)



omsfe_rwwn_var_table <- get_msfe(rwwn_var_out_forecasts, actual_obs,
                                 models = select(rwwn_var_out_list, id, var_set, n_lag, model_type=type),
                                 plus_group_vars = "model_type",
                                 msfe_name = "omsfe", msfe_type = "out-of-sample")


##### banbura step 7: calculate relative omsfe 

# we need to join:
# omsfe_bvar_table, omsfe_rwwn_var_table, deltas
omsfe_rwwn_var_table
deltas
omsfe_bvar_table

# select rw or wn omsfe for each variable, for each h
omsfe_selected_rwwn <- left_join(deltas, 
                    filter(omsfe_rwwn_var_table, model_type %in% c("wn","rw")),
                    by=c("rw_wn"="model_type","variable"="variable")) %>% 
         select(-var_set, -n_lag, -delta)
omsfe_selected_rwwn %>% glimpse()


omsfe_bvar_table <- ungroup(omsfe_bvar_table) %>% mutate_each("as.numeric", n_lag) 


# replicate banbura table from page 79
filter_variables <- ( fit_set_info %>% filter(fit_set==desired_fit_set) ) $ variable %>% as.character()

omsfe_var_banbura <- ungroup(omsfe_rwwn_var_table) %>% 
  filter(model_type=="var",n_lag==desired_n_lag, variable %in% filter_variables) 

omsfe_bvar_banbura <- omsfe_bvar_table %>% 
  filter(fit_set==desired_fit_set,n_lag==desired_n_lag, variable %in% filter_variables ) %>%
  select(-fit_set) %>% mutate(model_type="bvar")

omsfe_rwwn_banbura <- omsfe_selected_rwwn %>% 
  filter(variable %in% filter_variables ) %>%
  rename(model_type=rw_wn) 

omsfe_var_banbura 
omsfe_bvar_banbura
omsfe_rwwn_banbura 

rmsfe_long <- bind_rows(omsfe_var_banbura,omsfe_bvar_banbura) %>%
  left_join(omsfe_rwwn_banbura %>% select(-model_type, omsfe_rwwn=omsfe), by=c("h","variable")) %>%
  mutate(rmsfe=omsfe/omsfe_rwwn) %>% select(-n_lag)

rmsfe_long

# not automated!!!
rmsfe_wide <- rmsfe_long %>% select(-omsfe,-omsfe_rwwn) %>% 
  dcast(h+variable~var_set+model_type, value.var="rmsfe") %>%
  select(h,variable,set_3_var,set_3_bvar,set_6_var,set_6_bvar,set_23_bvar) 
some_rmsfe_wide <- rmsfe_wide %>% filter(h %in% desired_h)

some_rmsfe_wide
# step 8: optimal VAR selection

# if using SC with lag.max=12 then optimal is 1, no changes are needed
create_best_var_list("SC",12)


