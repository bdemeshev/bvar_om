library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
#library("data.table")
library("reshape2")
library("vars")
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
  
  
  if (test) message("Estimating model for file ", minfo$file, "...")
  
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
    
    priors <- Carriero_priors( D, p=n_lag,    # p=n_lag
                               lambdas=c(l_1, l_power, l_sc, l_io, l_const, l_exo) )
    
    # priors$X_dummy <- NULL
    # priors$Y_dummy <- NULL
    # estimate model
    
    
    verbose <- FALSE
    keep <- 10000
    
    if (test) { # testing mode
      verbose <- TRUE
      keep <- 100
    } 
    
    model <- bvar_conjugate0(priors = priors, 
                             verbose=verbose, keep=keep, 
                             fast_forecast = TRUE) 
    status <- "estimated"
  }
  
  
  if (minfo$type=="var") {
    n_lag <- as.numeric(minfo$n_lag)
    model <- VAR(D, p=n_lag, type="const") 
    status <- "estimated"
    
    if (sum(is.na(unlist(coef(model))))>0) {
      message("Can't estimate VAR model id = ",minfo$id,", file =",minfo$file)
      status <- "failed"
    }
  }
  
  if (minfo$type=="wn") {
    model_vector <- c(t(apply(D, MARGIN=2, mean))) # just sample means of variables
    model <- data.frame(variables=variables)
    model$coef <- model_vector
    status <- "estimated"
  }
  
  if (minfo$type=="rw") {
    # just mean growth rate
    model_vector <- c(t((tail(D,1)-head(D,1))/ (nrow(D)-1)))
    model <- data.frame(variables=variables)
    model$coef <- model_vector
    status <- "estimated"
  }
  
  
  saveRDS(model, model_full_path)
  if (do_log) {
    cat(file=log_con, Sys.time, ", model ", minfo$id,", finished")
    close(log_con)
  }
  
  
  
  return(status)
}




estimate_models <- function(mlist, parallel = c("off","windows","unix"), 
                            no_reestimation=TRUE, ncpu=4, test=FALSE, do_log=FALSE,
                            progress_bar=TRUE) {
  
  parallel <- match.arg(parallel)
  
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
    
    if (length(model_ids)>0) {
      if (progress_bar) pb <- txtProgressBar(min = 1, max = length(model_ids), style = 3)
      for (i in 1:length(model_ids))  {
        model_info <- mlist_todo %>% dplyr::filter(id==model_ids[i])
        status <- estimate_model(model_info, do_log=do_log, test=test)
        mlist$value[(mlist$id==model_ids[i])&(mlist$variable=="status")] <- status
        
        if (progress_bar) setTxtProgressBar(pb, i)
      }
      if (progress_bar) close(pb)
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
  
  ###### in-sample forecasts
  
  if ((minfo$type=="wn") & (pred_info$type=="in-sample")) {
      Tf_start <- T_start
      Tf_length <- T_in
      Tf_end <- Tf_start + Tf_length - 1
      
      value <- model$coef[rep(1:n_vars,Tf_length)]
      t <- rep(Tf_start:Tf_end, each=n_vars)
      answer <- data.frame(value=value, t=t, 
                           variable = rep(variables, Tf_length), h=NA)
  }

  if ((minfo$type=="rw") & (pred_info$type=="in-sample")) {
    Tf_start <- T_start + 1 # is not possible to forecast first observation
    Tf_length <- T_in - 1
    Tf_end <- Tf_start + Tf_length - 1
    
    ## multivariate analog of simple idea: y_first + (1:h)*Delta_y
    y_first <- c(t(D[1,])) # c(t()) is a transformation of data.frame row into a vector
    value <- y_first + rep(1:Tf_length, each=n_vars) * model$coef[rep(1:n_vars,Tf_length)] 
    t <- rep(Tf_start:Tf_end, each=n_vars)
    answer <- data.frame(value=value, t=t, 
                         variable = rep(variables, Tf_length), h=NA)
  }
  
  if ((minfo$type=="conjugate") & (pred_info$type=="in-sample")) {
    n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + n_lag # is not possible to forecast first n_lag observation
    Tf_length <- T_in - n_lag
    Tf_end <- Tf_start + Tf_length - 1
    
    predictions <- forecast_conjugate(model, 
                  fast_forecast = TRUE, out_of_sample = FALSE)
    answer <- mutate(predictions, t=h+Tf_start-1, h=NA) %>% select(-what) 
  }
  
  if ((minfo$type=="var") & (pred_info$type=="in-sample")) {
    n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + n_lag # is not possible to forecast first n_lag observation
    Tf_length <- T_in - n_lag
    Tf_end <- Tf_start + Tf_length - 1
    
    predictions <- data.frame(fitted(model))
    predictions$t <- Tf_start:Tf_end
    answer <- melt(predictions, id.vars = "t")
    answer$h <- NA
  }
  
  ###### out-of-sample forecasts
  
  if ((minfo$type=="conjugate") & (pred_info$type=="out-of-sample")) {
    # n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + T_in # where forecast starts
    Tf_length <- pred_info$h
    Tf_end <- Tf_start + Tf_length - 1
    
    predictions <- forecast_conjugate(model, h = pred_info$h, include="mean", level=NULL, 
                                      fast_forecast = TRUE)
    answer <- mutate(predictions, t=h+Tf_start-1) %>% select(-what) 
  }
  
  if ((minfo$type=="wn") & (pred_info$type=="out-of-sample")) {
    Tf_start <- T_start + T_in # where forecast starts
    Tf_length <- pred_info$h
    Tf_end <- Tf_start + Tf_length - 1
    
    value <- model$coef[rep(1:n_vars,Tf_length)]
    t <- rep(Tf_start:Tf_end, each=n_vars)
    answer <- data_frame(value=value, t=t, 
                         variable = rep(variables, Tf_length), h=t-Tf_start+1)
  }
  
  if ((minfo$type=="rw") & (pred_info$type=="out-of-sample")) {
    Tf_start <- T_start + T_in # where forecast starts
    Tf_length <- pred_info$h
    Tf_end <- Tf_start + Tf_length - 1
    
    ## multivariate analog of simple idea: y_last + (1:h)*Delta_y
    y_last <- c(t(D[Tf_start-1,])) # c(t()) is a transformation of data.frame row into a vector
    value <- y_last + rep(1:Tf_length, each=n_vars) * model$coef[rep(1:n_vars,Tf_length)] 
    t <- rep(Tf_start:Tf_end, each=n_vars)
    answer <- data_frame(value=value, t=t, 
                         variable = rep(variables, Tf_length), h=t-Tf_start+1)
  }
  
  if ((minfo$type=="var") & (pred_info$type=="out-of-sample")) {
    n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + T_in # where forecast starts
    Tf_length <- pred_info$h
    Tf_end <- Tf_start + Tf_length - 1
    
    answer <- predict(model, n.ahead = pred_info$h)$fcst %>% melt() %>%
      filter(Var2=="fcst")
    
    # in cases (h==1) and (h>1) predict+melt returns different Var1
    if (pred_info$h==1) {
      answer <- answer %>% select(variable=L1, value) %>%
         mutate(h=1, t=h+Tf_start-1)
    } else { # case (h>1)
      answer <- answer %>% select(h=Var1, variable=L1, value) %>%
         mutate(t=h+Tf_start-1)
    }
    answer <- mutate(answer, value=as.numeric(value))
  }
  
  
  answer$model_id <- model_id
  # easy answer$h <- pred_info$h is WRONG! max_h is requested, but all h will be computed!!!!
  answer$type <- pred_info$type
  
  return(answer) 
}





# function to make forecasts of many model for many datasets
forecast_models <- function(plist, mlist, parallel = c("off","windows","unix"), 
                            ncpu=4, test=FALSE, do_log=FALSE,
                            progress_bar=TRUE) {
  
  parallel <- match.arg(parallel)
  
  answer <- NULL
  if (parallel=="off") {
    # fast version using rbindlist from data.table
    all_data <- list()
    if (progress_bar) pb <- txtProgressBar(min = 1, max = nrow(plist), style = 3)
    for (i in 1:nrow(plist)) {
      all_data[[i]] <- forecast_model(plist[i,],mlist=mlist, parallel = parallel,
                                      ncpu=ncpu, test=test, do_log=do_log)
      if (progress_bar) setTxtProgressBar(pb, i)
    }
    if (progress_bar) close(pb)
    answer <- rbind_all(all_data)
  }
  return(answer) 
}





