# model lists
# this script creates model lists

T_common <- 120 # число наблюдений для которых строится прогнозы внутри выборки
# WN модель будет оцениваться по T_in <- T_common наблюдениям
# RW модель будет оцениваться по T_in <- T_common + 1 наблюдению, чтобы 
# получить ровно T_common прогнозов
# VAR(p), BVAR(p) будут оцениваться по T_in <- T_common + p наблюдениям
p_max <- 12 # для выравнивания первого внутривыборочного прогноза
# T_start <- p_max + 1 - n_lag
# при p_max=n_lag на вход модели будут подаваться наблюдения начиная с первого 


# id variable value
# 1 type "conjugate"/"minnesota"/"independent"/"var"/"theta"/"rw"/"wn" (или не все?)
# 1 n_lag 5 число лагов - подбирается
# 1 T_start 1995.0 номер стартового наблюдения подаваемого на вход функции
# 1 T_in сколько наблюдений подаётся на вход функции
# 1 var_set "3/5/7A/7Б/15" набор эндогенных переменных
# 
# Поля:
# "theta": T_start, T_in, var_set, file, status 
# "var": T_start, T_in, var_set, n_lag, file, status
# "rw": T_start, T_in, var_set, file, status
# "wn": T_start, T_in, var_set, file, status
# "conjugate": T_start, T_in, var_set, n_lag, l0, l1, l3, l4, seed, file, status
# На выходе: сохраняется модель. 
# Там, где не надо (theta), сохраняется что-то типа "ok"
# для rw, wn, theta поле var_set фактические показывает для каких переменных оценивается данная модель
# там имеет смысл указывать самый широкий var_set



# Для построения прогноза:
# model_id, h (количество шагов вперед), type=in-sample/out-of-sample
# функция одного прогноза:
# На выходе: model_id, h, t, переменная (текстом название), значение, type
# функция многих прогнозов:
# 

# this list is for testing purposes only, it is not used in Banbura procedure
create_model_list <- function() {
  # в столбце value получаем тип character
  mlist <- expand.grid(type="conjugate", 
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=12,
                    l_1=c(0.01,0.1,1,2,5,10),
                    l_power=1,
                    l_const=1,
                    l_sc=1,
                    l_io=1,
                    seed=13, # good luck, mcmc
                    status="not estimated")
  mlist <- mlist %>% mutate_each("as.character",type, status, var_set) 
  mlist <- mlist %>% mutate(id=row_number())
  mlist <- mlist %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  mlist <- mlist %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  mlist <- reshape2::melt(mlist, id.vars="id") %>% arrange(id)
  
  return(mlist)
}



create_rwwn_list <- function() {
  # в столбце value получаем тип character
  mlist <- data_frame(type=c("rw","wn"),
                   var_set="set_23",
                   n_lag=c(1,0),
                   T_start=c(p_max,p_max+1), T_in=c(T_common + 1, T_common),
                   status="not estimated")
  mlist <- mlist %>% mutate_each("as.character",type, status, var_set) 
  mlist <- mlist %>% mutate(id=row_number())
  mlist <- mlist %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  mlist <- reshape2::melt(mlist, id.vars="id") %>% arrange(id)
  
  return(mlist)
}



create_var_list <- function() {
  # в столбце value получаем тип character
  mlist <- expand.grid(type="var", 
                    var_set=c("set_3","set_6"), # no set 23 in var
                    n_lag=c(1,6,12),
                    status="not estimated")
  mlist <- mlist %>% mutate_each("as.character",type, status, var_set) 
  mlist <- mlist %>% mutate(id=row_number())
  mlist <- mlist %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  mlist <- mlist %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  mlist <- reshape2::melt(mlist, id.vars="id") %>% arrange(id)
  return(mlist)
}

create_best_var_list <- function(criterion = c("AIC","HQ","SC","FPE"), lag.max = 24) {
  criterion <- match.arg(criterion)
  # в столбце value получаем тип character
  mlist <- expand.grid(type="var", 
                    var_set=c("set_3","set_6"), # no set 23 in var
                    n_lag=NA, # will fill them automatically! :)
                    status="not estimated")
  mlist <- mlist %>% mutate_each("as.character",type, status, var_set) 
  mlist <- mlist %>% mutate(id=row_number())
  
  # here we compute optimal n_lag for each model
  
  for (i in 1:nrow(mlist)) {
    var_set_i <- mlist$var_set[i]
    variables <- filter(var_set_info, var_set == var_set_i)$variable
    D <- df[1:(lag.max+T_common), variables]
    selection_results <- VARselect(D, lag.max = lag.max, type = "const" )
    # let's use AIC # let's use Schwartz Criterion:
    mlist$n_lag[i] <- selection_results$selection[paste0(criterion,"(n)")]
  }
  
  if (max(mlist$n_lag)>p_max) warning("Maximum found lag is greater than global constant p_max!!!")
  
  mlist <- mlist %>% mutate(T_in = T_common + n_lag, T_start = max(mlist$n_lag,p_max) + 1 - n_lag)
  mlist <- mlist %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  mlist_melted <- reshape2::melt(mlist, id.vars="id") %>% arrange(id)
  return(mlist_melted)
}



create_bvar_banbura_list <- function() {
  # в столбце value получаем тип character
  mlist <- expand.grid(type="conjugate", 
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=c(1,6,12),
                    l_1=c(0.01,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.75,1,2,5,Inf),
                    l_power=1,
                    l_const=Inf,
                    # l_sc=1,
                    l_io=c(1,NA), # NA means no initial observations
                    seed=13, # good luck, mcmc
                    status="not estimated")
  
  # add no sc dummy and l_sc=10*l_1
  mlist_sc <- mlist %>% mutate(l_sc=10*l_1)
  mlist_nosc <- mlist %>% mutate(l_sc=NA)
  mlist <- rbind_list(mlist_sc,mlist_nosc)
  
  mlist <- mlist %>% mutate_each("as.character",type, status, var_set) 
  mlist <- mlist %>% mutate(id=row_number())
  mlist <- mlist %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  mlist <- mlist %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  mlist <- reshape2::melt(mlist, id.vars="id") %>% arrange(id)
  
  return(mlist)
}

# add for each model all possible shifts for T_in
rolling_model_replicate <- function(mlist_small) {
  # we need mlist_small with n_lag
  
  # add starting time:
  # T_start should vary from (p_max + 1 - n_lag) to (T_available - T_common - n_lag)
  mlist <- mutate(mlist_small, rownum = row_number(), T_start_min = p_max + 1 - n_lag, 
               T_start_max = T_available - T_common - n_lag,
               T_start_amount = T_start_max - T_start_min + 1)
  replicator <- rep(mlist$rownum, times = mlist$T_start_amount)
  mlist_big <- mlist[replicator,]
  mlist_big <- mlist_big %>% group_by(rownum) %>% 
    mutate(T_start = T_start_min + row_number() - 1)
  
  # remove auxillary variables, create id
  mlist_big <- ungroup(mlist_big) %>% 
    select(-T_start_min, -T_start_max, -T_start_amount, -rownum) %>%
    mutate(status="not estimated", id=row_number())   

  return(mlist_big)
}


create_bvar_out_list <- function(best_lambda) {
  # ungroup and clear junk variables from data.frame:
  # we need to keep fit_set variable for further comparison 
  mlist <- ungroup(best_lambda) %>% select(var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, fit_set)
  mlist <- mutate_each(mlist, "as.numeric", n_lag, l_1, l_const, l_io, l_power, l_sc)
  mlist <- mutate(mlist, status = "not estimated", type = "conjugate", 
               seed=13, T_in = n_lag + T_common )
  
  mlist_big <- rolling_model_replicate(mlist) 
  
  mlist_big <- mlist_big %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # mlist <- mlist %>% mutate_each("as.factor",type, status, var_set) 
  
  # mlist_big <- mlist_big %>% select(-rownum, -T_start_min, -T_start_max, -T_start_amount)
  
  mlist_melted <- reshape2::melt(mlist_big, id.vars="id") %>% arrange(id)
  
  return(mlist_melted)
}
