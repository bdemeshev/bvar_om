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


create_model_list <- function() {
  # в столбце value получаем тип character
  df <- expand.grid(type="conjugate", 
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=12,
                    l_1=c(0.01,0.1,1,2,5,10),
                    l_power=1,
                    l_const=1,
                    l_sc=1,
                    l_io=1,
                    seed=13, # good luck, mcmc
                    status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}



create_rwwn_list <- function() {
  # в столбце value получаем тип character
  df <- data_frame(type=c("rw","wn"),
                   var_set="set_23",
                   T_start=c(p_max,p_max+1), T_in=c(T_common + 1, T_common),
                   status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}



create_var_list <- function() {
  # в столбце value получаем тип character
  df <- expand.grid(type="var", 
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=c(1,6,12),
                    status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  return(df)
}

create_bvar_banbura_list <- function() {
  # в столбце value получаем тип character
  df <- expand.grid(type="conjugate", 
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=c(1,6,12),
                    l_1=c(0.01,0.1,1,2,5,10),
                    l_power=1,
                    l_const=1,
                    l_sc=1,
                    l_io=1,
                    seed=13, # good luck, mcmc
                    status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(T_in = T_common + n_lag, T_start = p_max + 1 - n_lag)
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}



create_bvar_out_list <- function(best_lambda) {
  # ungroup and clear junk variables from data.frame:
  # we need to keep fit_set variable for further comparison 
  df <- ungroup(best_lambda) %>% select(var_set, n_lag, l_1, l_const, l_io, l_power, l_sc, fit_set)
  df <- mutate_each(df, "as.numeric", n_lag, l_1, l_const, l_io, l_power, l_sc)
  df <- mutate(df, status = "not estimated", type = "conjugate", 
               seed=13, T_in = n_lag + T_common, rownum = row_number() )
  
  # add starting time:
  # T_start should vary from (p_max + 1 - n_lag) to (T_available - T_common - n_lag)
  df <- mutate(df, T_start_min = p_max + 1 - n_lag, 
                   T_start_max = T_available - T_common - n_lag,
                   T_start_amount = T_start_max - T_start_min + 1)
  replicator <- rep(df$rownum, times = df$T_start_amount)
  df_big <- df[replicator,]
  df_big <- df_big %>% group_by(rownum) %>% 
    mutate(T_start = T_start_min + row_number() - 1)

  df_big <- ungroup(df_big) %>% mutate(id=row_number())
  df_big <- df_big %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,
                                  "_lams_",round(100*l_1),
                                  "_sc_",round(100*l_sc),
                                  "_io_",round(100*l_io),
                                  "_pow_",round(100*l_power),
                                  "_cst_",round(100*l_const),
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df_big <- df_big %>% select(-rownum, -T_start_min, -T_start_max, -T_start_amount)
  
  df_melted <- reshape2::melt(df_big, id.vars="id") %>% arrange(id)
  
  return(df_melted)
}
