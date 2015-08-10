
# возможно функции надо передавать df
# параллельные вычисления?
estimate_model <- function(model_info, 
                           do_log=FALSE, test=FALSE ) {
  minfo <- reshape2::dcast(model_info, id~variable)
  model_full_path <- paste0("../estimation/models/",minfo$file)
  
  
  
  if (do_log) {
    log_con <- file(paste0("../estimation/log_model_",minfo$id,".csv"), open="a")
    cat(file=log_con, Sys.time, ", model ", minfo$id,", started")
  }
  
  # select observations
  T_start <- as.numeric(minfo$T_start)
  T_in <- as.numeric(minfo$T_in)
  T_end <- T_start+T_in-1
  
  # select variables
  v_set <- minfo$var_set
  variables <- filter(var_set_info, var_set==v_set) $ variable
  
  D <- df[T_start:T_end, variables]
  
  if (minfo$type=="conjugate") {
    
    seed <- as.numeric(minfo$seed)
    set.seed(seed)  # wish your good luck, MCMC
    
    
    
    # create priors
    l_1 <- as.numeric(minfo$l_1)
    l_power <- as.numeric(minfo$l_power)
    l_sc <- as.numeric(minfo$l_sc)
    l_io <- as.numeric(minfo$l_io)
    l_exo <- 1 # does not matter as we don't have exo variables
    n_lag <- as.numeric(minfo$n_lag)
    
    priors <- Carriero_priors( D, p=n_lag,    # p=n_lag, for 23 only 5 
                lambdas=c(l_1, l_power, l_sc, l_io, l_const, l_exo) )
    
    # priors$X_dummy <- NULL
    # priors$Y_dummy <- NULL
    # estimate model
    if (test) model <- bvar_conjugate0(priors = priors, verbose =TRUE, keep=100) 
    if (!test) model <- bvar_conjugate0(priors = priors)
  }
  
  
  if (minfo$type=="wn") {
    model_vector <- c(t(apply(D, MARGIN=2, mean))) # just sample means of variables
    model <- data.frame(variables=variables)
    model$coef <- model_vector
  }
  
  if (minfo$type=="rw") {
    # just mean growth rate
    model_vector <- c(t((tail(D,1)-head(D,1))/ (nrow(D)-1)))
    model <- data.frame(variables=variables)
    model$coef <- model_vector
  }
  
  
  saveRDS(model, model_full_path)
  if (do_log) {
    cat(file=log_con, Sys.time, ", model ", minfo$id,", finished")
    close(log_con)
  }
  
  status <- "estimated"
  
  return(status)
}

