# 500_replicate_banbura.R

# input: "../data/df_2015_final.csv", "../data/var_set_info.csv"
# output: 

source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")

# need to run only once 
source("200_load_after_eviews.R")

##################################################################
######################## begin set-up part #######################
##################################################################


parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30

df <- read_csv("../data/df_2015_final.csv")


T_available <- nrow(df)

fast_forecast <- TRUE # TRUE = posterior means of coefficients are used for forecast
keep <- 5000 # number of simulations from posterior (used only if fast_forecast is FALSE)
verbose <- FALSE # turn on/off messages from functions 
way_omega_post_root <- "svd" # "cholesky" or "svd" # how (Omega_post)^{1/2} is obtained

################################################
# create fit_set_info
# describe which msfe ratios are averaged in fit
fit_set_2vars <- data.frame(variable=c("ind_prod","cpi"),fit_set="ind+cpi")
fit_set_3vars <- data.frame(variable=c("ind_prod","cpi","ib_rate"),fit_set="ind+cpi+rate")

fit_set_info <- rbind(fit_set_2vars, fit_set_3vars)
fit_set_info

# all versions of fit_set are considered, but only desired is reported:
desired_fit_set <- "ind+cpi+rate"

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

var_set_info <- rbind(add_3,add_6,add_23)

# write_csv(var_set_info,"../data/var_set_info.csv")

##################################################################
######################## end set-up part #########################
##################################################################

####### step 0 (before banbura procedure)
# melting actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)


##### banbura step 1
# calculate msfe-0. Estimate RWWN (random walk OR white noise model)

# classify variables into RW and WN
deltas <- delta_i_prior(df, remove_vars = c("time_y","t"), c_0 = 0, c_1 = 1)
# delta_i_from_ar1(df, remove_vars = "time_y")

deltas <- mutate(deltas, rw_wn = ifelse(delta==1,"rw","wn"), variable=as.character(variable))
deltas

# estimate all RW and WN models
rwwn_list <- create_rwwn_list()
message("Estimate RW and WN")
rwwn_list <- estimate_models(rwwn_list,parallel = parallel)

# forecast all RW and WN models
rwwn_forecast_list <- data.frame(model_id=c(1,2), h=NA, type="in-sample")
message("Forecast RW and WN")
rwwn_forecasts <- forecast_models(rwwn_forecast_list, rwwn_list)

# calculate all msfe-0
# two ways to calculate msfe
# a) use all available predictions for each model
# b) use only common available predictions for all models

rwwn_forecasts %>% group_by(model_id) %>% summarise(Tf_start=min(t),Tf_end=max(t))
rwwn_forecasts <- rename(rwwn_forecasts, forecast=value)

plot_forecast(rwwn_forecasts, var_name="m2", mod_id=2)

#### if we prefer way b, just uncomment next two lines
# Tf_common_start <- 2
# rwwn_forecasts <- filter(rwwn_forecasts, t >= Tf_common_start)
####

# build table with corresponding msfe-0 (RW or WN)



# joining actual observations

rwwn_obs <- left_join(rwwn_forecasts, actual_obs, by=c("t","variable"))

rwwn_obs <- mutate(rwwn_obs, sq_error=(forecast-actual)^2)
head(rwwn_obs)

msfe0_all <- rwwn_obs %>% group_by(variable, model_id) %>% summarise(msfe=mean(sq_error))
msfe0_all

# add rw/wn label to msfe0
rwwn_wlist <- dcast(rwwn_list, id~variable) %>% # wlist = wide list
  mutate_each("as.numeric",n_lag,T_start,T_in) 
msfe0_all <- left_join(msfe0_all, dplyr::select(rwwn_wlist, id, type), by=c("model_id"="id") )

msfe0 <- left_join(deltas, msfe0_all, by= c("variable"="variable","rw_wn"="type") )

msfe0

##### banbura step 2

# estimate VAR
var_list <- create_var_list()
message("Estimating VAR")
var_list <- estimate_models(var_list,parallel = parallel)

# forecast VAR
var_forecast_list <- data.frame(model_id=unique(var_list$id), h=NA, type="in-sample")
message("Forecasting VAR")
var_forecasts <- forecast_models(var_forecast_list, var_list)

# joining actual observations
var_forecasts <- rename(var_forecasts, forecast=value)
var_obs <- left_join(var_forecasts, actual_obs, by=c("t","variable"))

var_obs <- mutate(var_obs, sq_error=(forecast-actual)^2)
head(var_obs)

# calculate msfe-inf

msfe_Inf <- var_obs %>% group_by(variable, model_id) %>% summarise(msfe_Inf=mean(sq_error))
msfe_Inf

# join model info
var_wlist <- dcast(var_list, id~variable) %>%
  mutate_each("as.numeric",n_lag,T_start,T_in)
var_wlist %>% glimpse()
msfe_Inf_info <- left_join(msfe_Inf, select(var_wlist, id, type, var_set, n_lag),
                           by=c("model_id"="id"))
msfe_Inf_info
# calculate fit-Inf

# msfe_0_Inf <- rename(msfe_Inf_info, msfe_Inf=msfe)
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
  fit_inf_table <- rbind(fit_inf_table,block)
}

fit_inf_table

##### banbura step 3
# goal: calculate fit-lam

# create model list to find optimal lambda 
bvar_list <- create_bvar_banbura_list()
# write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# estimate models from list
# bvar_list <- read_csv("../estimation/bvar_list.csv")
message("Estimating BVAR")
bvar_list <- estimate_models(bvar_list, parallel = parallel, ncpu=ncpu, test=FALSE) # status and filename are updated
write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# forecast BVAR
bvar_forecast_list <- data.frame(model_id=unique(bvar_list$id), h=NA, type="in-sample")
message("Forecasting BVAR in-sample")
bvar_forecasts <- forecast_models(bvar_forecast_list, bvar_list)

# joining actual observations
bvar_forecasts <- rename(bvar_forecasts, forecast=value)
bvar_obs <- left_join(bvar_forecasts, actual_obs, by=c("t","variable"))

bvar_obs <- mutate(bvar_obs, sq_error=(forecast-actual)^2)
head(bvar_obs)

# calculate msfe-lam

msfe_lam <- bvar_obs %>% group_by(variable, model_id) %>% summarise(msfe_lam=mean(sq_error))
msfe_lam

# join model info
bvar_wlist <- dcast(bvar_list, id~variable) %>%
  mutate_each("as.numeric",n_lag,T_start,T_in)
bvar_wlist %>% select(-file,-type)

msfe_lam_info <- left_join(msfe_lam, 
            select(bvar_wlist, id, type, var_set, n_lag,
                   l_1, l_const, l_io, l_power, l_sc, n_lag),
                           by=c("model_id"="id"))
msfe_lam_info 

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
  fit_lam_table <- rbind(fit_lam_table,block)
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
best_lambda <- fit_lam_table %>% group_by(var_set, n_lag, fit_set) %>% 
  mutate(fit_rank=dense_rank(delta_fit)) %>% filter(fit_rank == 1) %>% ungroup()
best_lambda



# check whether best lambda is unique
check_uniqueness <- best_lambda %>% group_by(var_set, n_lag, fit_set) %>% summarise(num_of_best_lambdas=n())
# num of best lambdas should be always one
if (max(check_uniqueness$num_of_best_lambdas)>1) warning("**** ACHTUNG ****: non unique optimal lambdas")
check_uniqueness %>% filter(num_of_best_lambdas>1)

# choose best_lambda in tie case --- with no sc dummy obs if possible
best_lambda <- best_lambda %>% group_by(var_set, n_lag, fit_set) %>% 
  mutate(sc_na_bonus=is.na(l_sc)) %>% arrange(desc(sc_na_bonus)) %>% 
  mutate(fit_tie_rank=row_number()) %>% 
  filter(fit_tie_rank == 1) %>% ungroup()
best_lambda


##### banbura step 5: calculate omsfe for bvar
# forecast and evaluate using optimal lambda

# create best models lists with correct time spec
# we need to keep fit_set variable for further comparison
bvar_out_list  <- create_bvar_out_list(best_lambda)
# best_lambda table should have: var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, fit_set

message("Estimate rolling BVAR")
bvar_out_list <- estimate_models(bvar_out_list,parallel = parallel)
write_csv(bvar_out_list, path = "../estimation/bvar_out_list.csv")

# forecast VAR
bvar_out_wlist <- dcast(bvar_out_list, id~variable) %>% 
  mutate_each("as.numeric", n_lag, seed, starts_with("l_"), starts_with("T_") )
bvar_out_wlist %>% glimpse()

h_max <- 12
bvar_out_forecast_list <- bvar_out_wlist %>% rowwise() %>% mutate(model_id=id, 
                                 h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")
# without rowwise min will be global and always equal to 1

# a lot of forecasts
### WARNING: maybe to many observations!!!
message("Forecasting rolling BVAR, out-of-sample")
bvar_out_forecasts <- forecast_models(bvar_out_forecast_list, bvar_out_list)

# joining actual observations
bvar_out_forecasts <- rename(bvar_out_forecasts, forecast=value)
bvar_out_obs <- left_join(bvar_out_forecasts, actual_obs, by=c("t","variable"))

bvar_out_obs <- mutate(bvar_out_obs, sq_error=(forecast-actual)^2)
head(bvar_out_obs)

# join models info 
bvar_out_obs <- left_join(bvar_out_obs, select(bvar_out_wlist, id, var_set, n_lag, fit_set),
                          by=c("model_id"="id"))
bvar_out_obs %>% head()

# calculate OMSFE by var_set, h, n_lag, variable
omsfe_bvar_table <- bvar_out_obs %>% group_by(var_set, n_lag, h, variable, fit_set) %>% 
  summarise(omsfe=mean(sq_error))
omsfe_bvar_table %>% head()

##### banbura step 6: calculate omsfe for RW/WN/VAR
# look at lists:
var_wlist
rwwn_wlist
str(rwwn_wlist)

rwwn_var_unique_wlist <- rbind_list(var_wlist, rwwn_wlist) %>% mutate_each("as.numeric", n_lag, T_in, T_start)

# every model should be rolled
rwwn_var_out_wlist <- rolling_model_replicate(rwwn_var_unique_wlist) %>% 
  mutate(file=paste0(type,"_",id,
                     "_T_",T_start,"_",T_in,"_",
                                var_set,"_lags_",n_lag,".Rds") ) 
rwwn_var_out_list <- melt(rwwn_var_out_wlist, id.vars = "id") %>% arrange(id)
message("Estimating rolling rw/wn/var models")
rwwn_var_out_list <- estimate_models(rwwn_var_out_list,parallel = parallel) # takes some minutes

h_max <- 12
rwwn_var_out_forecast_list <- rwwn_var_out_wlist %>% rowwise() %>% mutate(model_id=id, 
  h=min(h_max, T_available - T_start - T_in + 1 ) ) %>%
  select(model_id, h) %>% mutate(type="out-of-sample")
# without rowwise min will be global and always equal to 1

# forecast all rolling models  
message("Forecasting rolling rw/wn/var models, out-of-sample")
rwwn_var_out_forecasts <- forecast_models(rwwn_var_out_forecast_list, rwwn_var_out_list, progress_bar=TRUE)


# joining actual observations
rwwn_var_out_forecasts <- rename(rwwn_var_out_forecasts, forecast=value)
rwwn_var_out_obs <- left_join(rwwn_var_out_forecasts, actual_obs, by=c("t","variable"))

rwwn_var_out_obs <- mutate(rwwn_var_out_obs, sq_error=(forecast-actual)^2)
head(rwwn_var_out_obs)
glimpse(rwwn_var_out_obs)

# join models info 
rwwn_var_out_obs <- left_join(rwwn_var_out_obs, select(rwwn_var_out_wlist, id, var_set, n_lag, model_type=type),
                          by=c("model_id"="id"))
rwwn_var_out_obs %>% head()

# calculate OMSFE by var_set, h, n_lag, variable, model_type
omsfe_rwwn_var_table <- rwwn_var_out_obs %>% group_by(var_set, n_lag, h, variable, model_type) %>% 
  summarise(omsfe=mean(sq_error))
omsfe_rwwn_var_table 

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
# desired_fit_set <- "ind+cpi+rate" # set up in the beginning of file
filter_variables <- ( fit_set_info %>% filter(fit_set==desired_fit_set) ) $ variable %>% as.character()

omsfe_var_banbura <- ungroup(omsfe_rwwn_var_table) %>% 
  filter(model_type=="var",n_lag==12, variable %in% filter_variables) 

omsfe_bvar_banbura <- omsfe_bvar_table %>% 
  filter(fit_set==desired_fit_set,n_lag==12, variable %in% filter_variables ) %>%
  select(-fit_set) %>% mutate(model_type="bvar")

omsfe_rwwn_banbura <- omsfe_selected_rwwn %>% 
  filter(variable %in% filter_variables ) %>%
  rename(model_type=rw_wn) 

omsfe_var_banbura 
omsfe_bvar_banbura
omsfe_rwwn_banbura 

var_bvar_omsfe_banbura_table <- rbind_list(omsfe_var_banbura,omsfe_bvar_banbura) %>%
  left_join(omsfe_rwwn_banbura %>% select(-model_type, omsfe_rwwn=omsfe), by=c("h","variable")) %>%
  mutate(rmsfe=omsfe/omsfe_rwwn) %>% select(-n_lag)

var_bvar_omsfe_banbura_table


all_rmsfe_wide <- var_bvar_omsfe_banbura_table %>% select(-omsfe,-omsfe_rwwn) %>% 
  dcast(h+variable~var_set+model_type, value.var="rmsfe") %>%
  select(h,variable,set_3_var,set_3_bvar,set_6_var,set_6_bvar,set_23_bvar) 
some_rmsfe_wide <- all_rmsfe_wide %>% filter(h %in% c(1,3,6,12))

some_rmsfe_wide
# step 8: optimal VAR selection

# if using SC with lag.max=12 then optimal is 1, no changes are needed
create_best_var_list("SC",12)


