
# возможно функции надо передавать df
# параллельные вычисления?
estimate_model <- function(model_info, 
                           log=FALSE, test=FALSE ) {
  minfo <- reshape2::dcast(model_info, id~variable)
  model_full_path <- paste0("../estimation/models/",minfo$file)
  
  seed <- as.numeric(minfo$seed)
  set.seed(seed)  # wish your good luck, MCMC
  
  if (log) {
    log_con <- file(paste0("../estimation/log_model_",minfo$id,".csv"), open="a")
    cat(file=log_con, Sys.time, ", model ", minfo$id,", started")
  }
  
  
  if (minfo$type=="conjugate") {
    # select variables
    v_set <- minfo$var_set
    variables <- filter(var_set_info, var_set==v_set) $ variable
    
    # select observations
    T_start <- as.numeric(minfo$T_start)
    T_in <- as.numeric(minfo$T_in)
    T_end <- T_start+T_in-1
    
    D <- df[T_start:T_end, variables]
    
    # create priors
    l0 <- as.numeric(minfo$l0)
    l1 <- as.numeric(minfo$l1)
    l3 <- as.numeric(minfo$l3)
    l4 <- as.numeric(minfo$l4)
    n_lag <- as.numeric(minfo$n_lag)
    
    priors <- Carriero_priors( D, p=n_lag,    # p=n_lag, for 23 only 5 
                lambdas=c(l0,l1,l3,l4) )
    
    # priors$X_dummy <- NULL
    # priors$Y_dummy <- NULL
    # estimate model
    if (test) model <- bvar_conjugate0(priors = priors, verbose =TRUE, keep=100) 
    if (!test) model <- bvar_conjugate0(priors = priors)
  }
  
  
  
  saveRDS(model, model_full_path)
  if (log) cat(file=log_con, Sys.time, ", model ", minfo$id,", finished")
  
  status <- "estimated"
  
  return(status)
}

