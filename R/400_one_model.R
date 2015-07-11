

estimate_model <- function(model_info, 
                           log=FALSE ) {
  minfo <- reshape2::dcast(model_info, id~variable)
  model_full_path <- paste0("../estimation/models/",minfo$file)
  set.seed(minfo$seed)  # wish your good luck, MCMC
  
  if (log) {
    log_con <- file(paste0("../estimation/log_model_",minfo$id,".csv"), open="a")
    cat(file=log_con, Sys.time, ", model ", minfo$id,", started")
  }
  
  # ...
  
  
  
  saveRDS(model, model_full_path)
  if (log) cat(file=log_con, Sys.time, ", model ", minfo$id,", finished")
  
  status <- "estimated"
  
  return(status)
}

