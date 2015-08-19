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

T_available <- nrow(df)

##### banbura step 1
# calculate msfe-0. Estimate RWWN (random walk OR white noise model)

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

# calculate all msfe-0
# two ways to calculate msfe
# a) use all available predictions for each model
# b) use only common available predictions for all models

rwwn_forecasts %>% group_by(model_id) %>% summarise(Tf_start=min(t),Tf_end=max(t))

#### if we prefer way b, just uncomment next two lines
# Tf_common_start <- 2
# rwwn_forecasts <- filter(rwwn_forecasts, t >= Tf_common_start)
####

# build table with corresponding msfe-0 (RW or WN)

# joining actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

rwwn_forecasts <- rename(rwwn_forecasts, forecast=value)
rwwn_obs <- left_join(rwwn_forecasts, actual_obs, by=c("t","variable"))

rwwn_obs <- mutate(rwwn_obs, sq_error=(forecast-actual)^2)
head(rwwn_obs)

msfe0_all <- rwwn_obs %>% group_by(variable, model_id) %>% summarise(msfe=mean(sq_error))
msfe0_all

# add rw/wn label to msfe0
rwwn_wlist <- dcast(rwwn_list, id~variable) %>% # wlist = wide list
  mutate_each("as.numeric",n_lag,T_start,T_in) 
msfe0_all <- left_join(msfe0_all, select(rwwn_wlist, id, type), by=c("model_id"="id") )

msfe0 <- left_join(deltas, msfe0_all, by= c("variable"="variable","rw_wn"="type") )

msfe0

##### banbura step 2

# describe which msfe ratios are averaged in fit

  
fit_set_2vars <- data.frame(variable=c("ind_prod","cpi"),fit_set="ind+cpi")
fit_set_3vars <- data.frame(variable=c("ind_prod","cpi","ib_rate"),fit_set="ind+cpi+rate")

fit_set_info <- rbind(fit_set_2vars, fit_set_3vars)
fit_set_info


# estimate VAR
var_list <- create_var_list()
var_list <- estimate_models(var_list,parallel = parallel)

# forecast VAR
var_forecast_list <- data.frame(model_id=unique(var_list$id), h=NA, type="in-sample")
var_forecasts <- forecast_models(var_forecast_list, var_list)

# joining actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

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
write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# estimate models from list
bvar_list <- read_csv("../estimation/bvar_list.csv")
bvar_list <- estimate_models(bvar_list, parallel = parallel, ncpu=ncpu) # status and filename are updated
write_csv(bvar_list, path = "../estimation/bvar_list.csv")

# forecast BVAR
bvar_forecast_list <- data.frame(model_id=unique(bvar_list$id), h=NA, type="in-sample")
bvar_forecasts <- forecast_models(bvar_forecast_list, bvar_list)

# joining actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

bvar_forecasts <- rename(bvar_forecasts, forecast=value)
bvar_obs <- left_join(bvar_forecasts, actual_obs, by=c("t","variable"))

bvar_obs <- mutate(bvar_obs, sq_error=(forecast-actual)^2)
head(var_obs)

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
# group_by is needed! otherwise cannot remove var_set (active group on var_set)
fit_goal <- dplyr::filter(fit_inf_table, var_set=="set_3") %>% ungroup() %>% dplyr::select(-var_set)
fit_goal

fit_lam_table <- left_join(fit_lam_table, fit_goal, by=c("n_lag","fit_set"))
fit_lam_table <- mutate(fit_lam_table, delta_fit = abs(fit_lam-fit_inf))

# for each var_set, n_lag, fit_set the best lambdas are calculated
best_lambda <- fit_lam_table %>% group_by(var_set, n_lag, fit_set) %>% 
  mutate(fit_rank=dense_rank(delta_fit)) %>% filter(fit_rank==1) %>% ungroup()




# check whether best lambda is unique
check_uniqueness <- best_lambda %>% group_by(var_set, n_lag, fit_set) %>% summarise(num_of_best_lambdas=n())
# num of best lambdas should be always one
if (max(check_uniqueness$num_of_best_lambdas)>1) warning("**** ACHTUNG ****: non unique optimal lambdas")
check_uniqueness


best_lambda %>% arrange(fit_set, n_lag, var_set)

##### banbura step 5: calculate omsfe for bvar
# forecast and evaluate using optimal lambda

# create best models lists with correct time spec
# we need to keep fit_set variable for further comparison
bvar_out_list  <- create_bvar_out_list(best_lambda)
# best_lambda table should have: var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, fit_set

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
message("Forecasting rolling bvar models, out-of-sample")
bvar_out_forecasts <- forecast_models(bvar_out_forecast_list, bvar_out_list)

# joining actual observations
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

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
df <- mutate(df, t=row_number()) 
actual_obs <- melt(df, id.vars="t" ) %>% rename(actual=value)

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

# relative-msfe
rwwn_rel_msfe <- left_join(omsfe_bvar_table %>% rename(omsfe_bvar=omsfe), 
                           omsfe_selected_rwwn %>% rename(omsfe_rwwn=omsfe), 
                           by=c("variable","h")) 
rwwn_rel_msfe %>% filter(fit_set=="ind+cpi+rate",var_set=="set_23",n_lag==6,h==12)

rwwn_rel_msfe %>% glimpse()

# casting VAR models 
var_cast <- ungroup(omsfe_rwwn_var_table) %>%
          filter(model_type=="var") %>% select(-model_type) %>%
          rename(omsfe_var=omsfe) %>% 
          mutate(var_set=plyr::revalue(var_set, c("set_3"="omsfe_var_3","set_6"="omsfe_var_6"))) %>%
  dcast(variable+n_lag+h~var_set, value.var="omsfe_var")
var_cast

# joining omsfe from VAR
rwwn_rel_msfe <- left_join(rwwn_rel_msfe, var_cast, by=c("variable","n_lag","h"))
rwwn_rel_msfe

# calculate rel-msfe
rwwn_rel_msfe <- rwwn_rel_msfe %>% mutate(rmsfe_bvar=omsfe_bvar/omsfe_rwwn,
                                          rmsfe_var_3=omsfe_var_3/omsfe_rwwn,
                                          rmsfe_var_6=omsfe_var_6/omsfe_rwwn)
rwwn_rel_msfe
