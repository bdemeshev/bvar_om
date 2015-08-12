library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
library("dplyr")
library("bvarr")

source("400_model_lists.R")

########################################################################
################ estimating functions
########################################################################

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
  # as.numeric is needed because in mlist all values have the same type (most general, character)
  T_start <- as.numeric(minfo$T_start)
  T_in <- as.numeric(minfo$T_in) # number of supplied observations
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
    l_const <- as.numeric(minfo$l_const)
    # l_kron <- as.numeric(minfo$l_kron)
    l_io <- as.numeric(minfo$l_io)
    l_exo <- 1 # does not matter as we don't have exo variables
    n_lag <- as.numeric(minfo$n_lag)
    
    fast_forecast <- TRUE # later we may add as option
    
    priors <- Carriero_priors( D, p=n_lag,    # p=n_lag, for 23 only 5 
                               lambdas=c(l_1, l_power, l_sc, l_io, l_const, l_exo) )
    
    # priors$X_dummy <- NULL
    # priors$Y_dummy <- NULL
    # estimate model
    
    if (test) { # testing mode
      verbose <- TRUE
      keep <- 100
    } else { # normal mode
      verbose <- FALSE
      keep <- 10000
    }
    
    model <- bvar_conjugate0(priors = priors, 
                             verbose=verbose, keep=keep, 
                             fast_forecast = fast_forecast) 
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

########################################################################
################ forecasting functions
########################################################################

# function to make forecast of one model for one dataset
# pred_info - one line describing desired forecast
# mlist - list of all models
forecast_model <- function(pred_info, mlist, parallel = parallel, 
                            ncpu=4, test=FALSE, do_log=FALSE) {
  
  
  # get info about model:
  model_id <- pred_info$model_id 
  model_info <- mlist %>% filter(id==model_id)
  minfo <- reshape2::dcast(model_info, id~variable)
  model_full_path <- paste0("../estimation/models/",minfo$file)
  model <- readRDS(model_full_path)
  
  T_start <- as.numeric(minfo$T_start) # starting observation for model estimation 
  T_in <- as.numeric(minfo$T_in) # number of observation supplied for model estimation
  
  # select variables
  v_set <- minfo$var_set
  variables <- filter(var_set_info, var_set==v_set) $ variable
  n_vars <- length(variables)
  D <- df[,variables]
  
  answer <- NULL
  
  if (minfo$type=="wn") {
    if (pred_info$type=="in-sample") {
      Tf_start <- T_start
      Tf_length <- T_in
      h <- Tf_length
      Tf_end <- Tf_start + Tf_length - 1
      
      value <- model$coef[rep(1:n_vars,Tf_length)]
      t <- rep(Tf_start:Tf_end, each=n_vars)
      answer <- data.frame(value=value, t=t, variable = rep(model$variables, Tf_length), h=NA)
    }
  }

  if (minfo$type=="rw") {
    Tf_start <- T_start + 1 # is not possible to forecast first observation
    Tf_length <- T_in - 1
    h <- Tf_length
    Tf_end <- Tf_start + Tf_length - 1
    
    ## multivariate analog of simple idea: y_first + (1:h)*Delta_y
    y_first <- c(t(D[1,])) # c(t()) is a transformation of data.frame row into a vector
    value <- y_first + rep(1:Tf_length, each=n_vars) * model$coef[rep(1:n_vars,Tf_length)] 
    t <- rep(Tf_start:Tf_end, each=n_vars)
    answer <- data.frame(value=value, t=t, variable = rep(model$variables, Tf_length), h=NA)
  }
  
  if (minfo$type=="conjugate") {
    
    
  }
  
  
  answer$model_id <- model_id
  # answer$h <- pred_info$h # WRONG! max_h is requested, but all h will be computed!!!!
  answer$type <- pred_info$type
  
  return(answer) 
}





# function to make forecasts of many model for many datasets
forecast_models <- function(plist, mlist, parallel = parallel, 
                            ncpu=4, test=FALSE, do_log=FALSE) {
  
  
  
  
  return(plist) # ???
}





