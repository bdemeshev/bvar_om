# model lists
# this script creates model lists

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
                    T_start=1, 
                    T_in=100,
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
                   T_start=1, T_in=120,
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


