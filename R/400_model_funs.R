library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
library("dplyr")


source("400_one_model.R")


# id variable value
# 1 type "conjugate"/"minnesota"/"independent"/"var"/"theta"/"rw" (или не все?)
# 1 n_lag 5 число лагов - подбирается
# 1 T_in сколько наблюдений подаётся на вход функции
# 1 T_start 1995.0 номер стартового наблюдения подаваемого на вход функции
# 1 var_set "3/5/7A/7Б/15" набор эндогенных переменных
# 
# Перебираемые параметры:
# "theta" нет
# "var" нет 
# "rw" нет



create_model_list <- function() {
  df <- expand.grid(type="mcmc", chain_no=1:5, status="not estimated", 
                    sample=c("west","east","all"), W=c("border","road"))
  df <- df %>% mutate_each("as.factor",type, status, sample, W) %>% mutate(seed=chain_no+42)
  df <- df %>% mutate(file=paste0(type,"_",W,"_",sample,"_",chain_no,".Rds"), id=row_number())
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}


estimate_models <- function(mlist, parallel = parallel, 
                            no_reestimation=TRUE, ncpu=4) {
  mlist_todo <- mlist
  model_ids <- unique(mlist$id)

  
  if (no_reestimation) {
    temp <- dplyr::filter(mlist, variable=="status") %>% filter(!value=="estimated")
    model_ids <- temp$id
    mlist_todo <- filter(mlist, id %in% model_ids)
  }  
  
  
  if (parallel=="windows") {
    library("doSNOW")
    cl <- makeCluster(ncpu, outfile="") # number of CPU cores
    registerDoSNOW(cl)
    
    foreach(i=model_ids, .packages=c("mvtnorm","MHadaptive","MCMCpack")) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, log=TRUE)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
    stopCluster(cl)
  }
  
  if (parallel=="unix") {
    library("doMC")
    registerDoMC(ncpu) # number of CPU cores

    foreach(i=model_ids, .packages=c("mvtnorm","MHadaptive","MCMCpack")) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, log=TRUE)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
  }
  
  if (parallel=="off") {
    for (i in model_ids)  {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, log=TRUE)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
  }
  
  
  return(mlist) # statuses are updated
}



