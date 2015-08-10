library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
library("dplyr")


source("400_one_model.R")


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
# model_id, h (количество шагов вперед), type=in-sample/out-of-sample, status, file
# На выходе: model_id, h, t, переменная, значение, type, status, file 







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
                    l_kron=1,
                    seed=13, # good luck, mcmc
                    status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,"_lams_",
                                  round(100*l0),"_",round(100*l1),"_",
                                  round(100*l3),"_",round(100*l4),
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}



create_model_list_banbura <- function() {
  # в столбце value получаем тип character
  df <- expand.grid(type="conjugate", 
                    T_start=1, 
                    T_in=100,
                    var_set=c("set_3","set_6","set_23"),
                    n_lag=12,
                    l0=c(0.01,0.1,1,2,5,10),
                    l1=1,
                    l3=1,
                    l4=1,
                    seed=13, # good luck, mcmc
                    status="not estimated")
  df <- df %>% mutate_each("as.character",type, status, var_set) 
  df <- df %>% mutate(id=row_number())
  df <- df %>% mutate(file=paste0(type,"_",id,"_T_",T_start,"_",T_in,"_",
                                  var_set,"_lags_",n_lag,"_lams_",
                                  round(100*l0),"_",round(100*l1),"_",
                                  round(100*l3),"_",round(100*l4),
                                  ".Rds") ) 
  # df <- df %>% mutate_each("as.factor",type, status, var_set) 
  
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}



estimate_models <- function(mlist, parallel = parallel, 
                            no_reestimation=TRUE, ncpu=4, test=FALSE, do_log=FALSE) {
  mlist_todo <- mlist
  model_ids <- unique(mlist$id)

  
  if (no_reestimation) {
    temp <- dplyr::filter(mlist, variable=="status") %>% filter(!value=="estimated")
    model_ids <- temp$id
    mlist_todo <- filter(mlist, id %in% model_ids)
  }  
  
  requested_packages <- c("mvtnorm","MHadaptive","MCMCpack","bvarr")
  
  if (parallel=="windows") {
    library("doSNOW")
    cl <- makeCluster(ncpu, outfile="") # number of CPU cores
    registerDoSNOW(cl)
    
    foreach(i=model_ids, .packages=requested_packages) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, do_log=do_log, test=test)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
    stopCluster(cl)
  }
  
  if (parallel=="unix") {
    library("doMC")
    registerDoMC(ncpu) # number of CPU cores

    foreach(i=model_ids, .packages=requested_packages) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, do_log=do_log, test=test)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
  }
  
  if (parallel=="off") {
    for (i in model_ids)  {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, do_log=do_log, test=test)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
  }
  
  
  return(mlist) # statuses are updated
}




# function to make forecast of one model for one dataset
# pred_info - one line describing desired forecast
# mlist - list of all models
forecast_model <- function(pred_info, mlist, parallel = parallel, 
                            ncpu=4, test=FALSE, do_log=FALSE) {
  
  # get info about model:
  model_id <- pred_info$model_id 
  model_info <- mlist %>% filter(id==model_id)
  model <- readRDS(model_info$file)
  
  if (model_info$type=="wn") {
    if (pred_info$type=="in-sample") {
      
    }
    
  }

  if (model_info$type=="rw") {
    
    
  }
  
  if (model_info$type=="conjugate") {
    
    
  }
  
  

  
  
    
  
  
  
  
  return(plist) # ???
}





# function to make forecasts of many model for many datasets
forecast_models <- function(plist, mlist, parallel = parallel, 
                            ncpu=4, test=FALSE, do_log=FALSE) {
  
  
  
  
  return(plist) # ???
}





