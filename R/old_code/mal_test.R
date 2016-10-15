
## ----, echo=FALSE, message=FALSE, include=FALSE, warning=FALSE-----------
message("The user name is '",Sys.info()[7],"'")

if (Sys.info()[7]=="boris") {
  Sys.setenv(X13_PATH ="/usr/local/bin/x13tramo")
  setwd("~/Documents/bvar/mal_art")
}

if (Sys.info()[7]=="Oxana") {
  Sys.setenv(X13_PATH="D:/Oxana/Science/Projects/ProjectXIV/x13as/")
  setwd("D:/Oxana/Science/Projects/ProjectXIV/mal_art/mal_art/")
}


source("mal_art_setup.R")
source("mal_art_fun.R")

df <- fread("data/data_bvar_csv.csv")
# str(df)
# glimpse(df)


## ------------------------------------------------------------------------
mvtsplot(df)


## ------------------------------------------------------------------------
# кроме ib_rate, unempl_rate
df_log <- mutate_each(df, "log", -time, -ib_rate, -unempl_rate)


## ----, "seasonally adjust"-----------------------------------------------
df_omit <- df_log %>% na.omit()
df_notime <- dplyr::select(df_omit,-time)
df_ts <- ts(df_notime, start=c(1995,9), frequency = 12)
to_adjust <- c("empl_manuf", "cpi", "ppi", "industry_prod",
               "real_income", "unempl_rate", "inv_index", "real_wage",
               "construct") # list of series to seasonally adjust
df_sa <- seas_adjust(df_ts, to_adjust )


## ------------------------------------------------------------------------
plot(df_ts[,1:10])
plot(df_sa[,1:10])
plot(df_ts[,11:14])
plot(df_sa[,11:14])


## ----, "unit root table", results='asis'---------------------------------
df_st_test <- df_ts[1:185, ]
non_station <- data_frame(var=colnames(df_ts),
                          kpss_t=bvarm_prior(df_st_test,test="KPSS",type="trend"),
                          kpss_nt=bvarm_prior(df_st_test,test="KPSS",type="constant"),
                          adf_t=bvarm_prior(df_st_test,test="ADF",type="trend"),
                          adf_c=bvarm_prior(df_st_test,test="ADF",type="constant"),
                          adf_0=bvarm_prior(df_st_test,test="ADF",type="neither"))
non_station
pander(non_station)


## ----, "table with priors"-----------------------------------------------
# autodetect using ADF+trend
all_sets <- c("set_3", "set_5", "set_6", "set_14")
set_3 <- data_frame(vars=c("industry_prod", "cpi", "ib_rate"), name="set_3")
set_5 <- data_frame(vars=c(set_3$vars, "ex_rate", "m2"), name="set_5")
set_6 <- data_frame(vars=c(set_5$vars, "oil_price"), name="set_6")
set_14 <- data_frame(vars=colnames(df_ts), name="set_14")
varset_table <- build_varset(list(set_3, set_5, set_6, set_14))
head(varset_table)

# manual correction: all series Random Walks except:
varset_table <- mutate(varset_table, prior = ifelse(var=="ib_rate", 0, 1))


## ----, "best number of lags"---------------------------------------------
lag_table <- NULL
for (i in all_sets) {
  vars <- dplyr::filter(varset_table, set==i)$var
  best_lag <- VARselect(df_ts[,vars])

  lag_table <- rbind(lag_table, best_lag$selection)
}
colnames(lag_table) <- c("AIC","HQ","SC","FPE")
lag_table <- data.frame(lag_table)
lag_table$set <- all_sets
lag_table


## ----, "table of models"-------------------------------------------------
#### describe all models in one table (mod_table)
# p - число лагов
# v, lambda - параметры из работы
# m - число рядов, считается само
# t_for - длина прогнозного периода
# T - длина ряда, считается сама


eps <- 10^(-42) # epsilon # 42 just for fun
# BVARM crashes for HP1 or HP2 exactly zero
# but HP1=HP2=10^(-200) is ok

lambda_all <- c(eps, 0.01, 0.1, 1, 10, 100, Inf)
t_for <- 24
p_all <- 1:5
v_all = c(0.25, 0.5, 0.75, 1)

# names of sets should be real! see build_varset function
bvarm_models <- bvarm_block(lambda=lambda_all, v=v_all, p=p_all, t_for=24,
                            varset=all_sets)
cvar_models <- cvar_block(p=p_all,varset=all_sets)


# mod_table <- readRDS("./estimation/mod_table.Rds")
# mod_table <- NULL # initial (!!!!)

# now one should not manually change mod_table
# only through add_models() function


add_models(bvarm_models)
# rbind(bvarm_models,cvar_models)
# no need to consider cvar_models --- they are just lambda=inf models


### далее идет новый блок, который мы обнаружили после того
## как поняли где лежат оптимальные лямбда

lambda_all <- seq(0.01,2,by=0.01)
bvarm_new <- bvarm_block(lambda=lambda_all, v=1, p=5, t_for=24,
                            varset=c("set_5","set_6","set_14"))
add_models(bvarm_new)



## ----, "estimation phase"------------------------------------------------

estimate_models()

## ----, eval=TRUE, "msfe calculations"------------------------------------
save_msfe <- c("cpi","industry_prod","ib_rate") # msfe to save

if ("msfe_table.Rds" %in% list.files("./estimation/")) { # if msfe_table.Rds exists...
  msfe_table <- readRDS("./estimation/msfe_table.Rds")
} else {
  msfe_table <- NULL # initial (!)
}

mod_table <- readRDS("./estimation/mod_table.Rds")

ids_to_msfe <- 1:nrow(mod_table)
ids_to_msfe <- ids_to_msfe[!ids_to_msfe %in% unique(msfe_table$model_id)]

files <- id2fname(ids_to_msfe)

for (file in files) {
  cat("calculate msfe for: ", save_msfe, "\n")
  cat("file: ", file, "\n")
  model <- readRDS(file=paste0("./estimation/models/",file))

  id <- fname2id(file)
  p <- dplyr::filter(mod_table, model_id==id)$p

  df_for <- forinsample(model)
  df_test <- model$data[-(1:p),]

  msfe_model <- fmeasure(df_for, df_test, measure="mse", vars=save_msfe)
  msfe_add <- data_frame(msfe=msfe_model, var=save_msfe, model_id=id)

  msfe_table <- rbind(msfe_table, msfe_add)
}
saveRDS(msfe_table, "./estimation/msfe_table.Rds")


## ------------------------------------------------------------------------
# mod_zero_inf <- dplyr::filter(mod_table, lambda %in% c(eps,Inf))

msfe_table <- readRDS("./estimation/msfe_table.Rds")

msfe_wide <- dcast(data = msfe_table, model_id~var, value.var="msfe")
head(msfe_wide)

mod_msfe <- dplyr::left_join(x = mod_table, y=msfe_wide)
mod_small <- dplyr::select(mod_msfe, model_id, v, lambda, p, varset, cpi, ib_rate, industry_prod)

mod_0 <- dplyr::filter(mod_small, lambda==eps)
# mod_Inf <- dplyr::filter(mod_small, lambda==Inf)

fit_table <- left_join(x=mod_small, y=mod_0, by=c("p","v","varset"))
fit_table <- mutate(fit_table,
                    fit2 = 1/2*(industry_prod.x/industry_prod.y+cpi.x/cpi.y),
                    fit3 = 1/3*(industry_prod.x/industry_prod.y+cpi.x/cpi.y+ib_rate.x/ib_rate.y))
saveRDS(fit_table, file="./estimation/fit_table.Rds")


## ----, eval=FALSE--------------------------------------------------------
##
## id <- 22 # 22: lambda=0, 28: lambda=Inf
## model <- readRDS(paste0("./estimation/models/",id2fname(id)))
## p <- dplyr::filter(mod_table, model_id==7)$p
##
## df_data <- model$data
## model_var <- VAR(df_data,type = "const")
## # summary(model_var)
## df_for <- forinsample(model)
## df_test <- model$data[-(1:p),]
## fmeasure(df_for, df_test, measure="mse", vars=save_msfe)
## # cbind(df_for,df_test)


## ------------------------------------------------------------------------
fit_table <- readRDS(file="./estimation/fit_table.Rds")

fit_inf <- fit_table %>% dplyr::filter(v==1, lambda.x==Inf, varset=="set_3") %>%
  dplyr::select( p,  fit2, fit3)

# add columns with target fit_inf
fit_joined <- left_join(x=fit_table, y=fit_inf, by="p")

# calculate distance
fit_joined <- mutate(fit_joined, dfit2=abs(fit2.x-fit2.y), dfit3=abs(fit3.x-fit3.y))

# find "best" lambda for each v, p, set
best2 <- fit_joined %>% group_by(v,p,varset) %>% filter(dfit2==min(dfit2))
best3 <- fit_joined %>% group_by(v,p,varset) %>% filter(dfit3==min(dfit3))

# здесь мы обнаружили, что иногда (при больших лямбда, т.е. 10, 100, Inf)
# разница между fit_(lambda,m) и fit_(inf,3) уже ноль
# поэтому там сразу несколько оптимальных фитов

# выберем fit с наибольшим лямбда
best2u <- best2 %>% group_by(v, p, varset) %>% filter(lambda.x==max(lambda.x)) %>%
  dplyr::select(p,v,varset,best_lam2=lambda.x)
best3u <- best3 %>% group_by(v, p, varset) %>% filter(lambda.x==max(lambda.x)) %>%
  dplyr::select(p,v,varset,best_lam3=lambda.x)
best <- left_join(x=best2u, y=best3u, by=c("p","v","varset"))

best <- mutate(best, modsize=mod_size(varset))


## ------------------------------------------------------------------------
qplot(data=best %>% dplyr::filter(modsize>3), x=modsize, y=best_lam2, color=factor(p),
      geom="line") +  facet_grid(~v)


## ------------------------------------------------------------------------
best_v1p5 <- best %>% dplyr::filter(v==1, p==5)


## ----, eval=FALSE, "testing forroling forecast"--------------------------
## forrolling(df_sa, lambda=1, varset="set_3", h=3)


## ------------------------------------------------------------------------
mod_eval_table <- best_v1p5 %>% ungroup() %>% dplyr::select(varset, lambda=best_lam2)
add_0_Inf <- expand.grid(varset=all_sets,lambda=c(eps,Inf)) # lambda should not be zero!
add_0_Inf <- mutate(add_0_Inf, varset=as.character(varset))
mod_eval_table <- rbind(mod_eval_table, add_0_Inf) %>% unique()
h_all <- c(1,3,6)
n_h <- length(h_all)
mod_eval_table <- mod_eval_table[rep(1:nrow(mod_eval_table),each=n_h),]
mod_eval_table$h <- rep(h_all, nrow(mod_eval_table)/n_h)

saveRDS(mod_eval_table, file="./estimation/mod_eval_table.Rds")


## ------------------------------------------------------------------------
save_msfe <- c("cpi","industry_prod","ib_rate") # msfe to safe
eps <- 10^(-42)

mod_eval_table <- readRDS("./estimation/mod_eval_table.Rds")

omsfe_table <- get_omsfe(df_sa, mod_eval_table)
# omsfe_table
saveRDS(omsfe_table, "./estimation/omsfe_table.Rds")




## ------------------------------------------------------------------------
omsfe_table <- readRDS("./estimation/omsfe_table.Rds")
head(omsfe_table)

omsfe_zero <- omsfe_table %>% dplyr::filter(lambda<2*eps) %>% # to play on the safe side
  dplyr::select(-lambda) %>% dplyr::rename(omsfe_zero=omsfe)
omsfe_inf <- omsfe_table %>% dplyr::filter(lambda==Inf) %>%
  dplyr::select(-lambda) %>% dplyr::rename(omsfe_inf=omsfe)
omsfe_join <- left_join(x=omsfe_table, y=omsfe_zero, by=c("var","varset","h") )
omsfe_join <- left_join(x=omsfe_join, y=omsfe_inf, by=c("var","varset","h") )

omsfe_join <- mutate(omsfe_join,
                     rw_romsfe=omsfe/omsfe_zero,
                     var_romsfe=omsfe/omsfe_inf)

romsfe_table <- omsfe_join %>% filter(lambda>2*eps) %>%
  filter( varset=="set_3" | lambda<Inf) %>%
  mutate(modsize=mod_size(varset))
saveRDS(romsfe_table, "./estimation/romsfe_table.Rds")


## ------------------------------------------------------------------------
romsfe_table <- readRDS("./estimation/romsfe_table.Rds")
romsfe_rw <- dplyr::select(romsfe_table, -omsfe_zero,
                              -omsfe_inf, -omsfe, -varset, -var_romsfe, -lambda)
bandura_t_rw <- dcast(data = romsfe_rw, formula = var+h~modsize, value.var = "rw_romsfe")
bandura_t_rw

romsfe_var <- dplyr::select(romsfe_table, -omsfe_zero,
                              -omsfe_inf, -omsfe, -varset, -rw_romsfe, -lambda)
bandura_t_var <- dcast(data = romsfe_var, formula = var+h~modsize, value.var = "var_romsfe")
bandura_t_var
saveRDS(bandura_t_rw, file="./estimation/bandura_t_rw.Rds")
saveRDS(bandura_t_var, file="./estimation/bandura_t_var.Rds")


