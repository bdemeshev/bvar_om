library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
library("dplyr")


source("400_one_model.R")

create_model_list <- function() {
  df <- expand.grid(type="mcmc", chain_no=1:5, status="not estimated", 
                    sample=c("west","east","all"), W=c("border","road"))
  df <- df %>% mutate_each("as.factor",type, status, sample, W) %>% mutate(seed=chain_no+42)
  df <- df %>% mutate(file=paste0(type,"_",W,"_",sample,"_",chain_no,".Rds"), id=row_number())
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}


estimate_models <- function(mlist, parallel = parallel, 
                            log_file="../estimation/estimation_log.txt",
                            no_reestimation=TRUE) {
  mlist_todo <- mlist
  model_ids <- unique(mlist$id)

  
  if (no_reestimation) {
    temp <- dplyr::filter(mlist, variable=="status") %>% filter(!value=="estimated")
    model_ids <- temp$id
    mlist_todo <- filter(mlist, id %in% model_ids)
  }  
  
  log_con <- file(log_file, open="a")
  
  if (parallel=="windows") {
    library("doSNOW")
    cl <- makeCluster(10, outfile="") # number of CPU cores
    registerDoSNOW(cl)
    
    foreach(i=model_ids, .packages=c("mvtnorm","MHadaptive","MCMCpack")) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, log_con)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
    stopCluster(cl)
  }
  
  if (parallel=="off") {
    for (i in model_ids)  {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, log_con)
      mlist$value[(mlist$id==i)&(mlist$variable=="status")] <- status
    }
  }
  
  close(log_con)
  
  return(mlist) # statuses are updated
}

do_log <- function(log_con,...) {
  msg <- paste0(...,"\n")
  cat(msg,file=log_con)
}



