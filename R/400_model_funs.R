# 400_model_funs.R

# library("foreach") # parallel processing
library("readr") # reading csv files
# library("MHadaptive")
library("MCMCpack") # IW
# library("mvtnorm") # multivariate normal 
# library("data.table")
library("reshape2")
library("dplyr")
library("ggplot2")

library("vars")
library("bvarr")


########################################################################
################ estimating functions
########################################################################

# maybe we need to pass actual df for parallel computing
estimate_model <- function(minfo, 
                           do_log=FALSE, verbose=FALSE ) {
  model_full_path <- paste0("../estimation/models/",minfo$file)
  
  
  if (verbose) message("Estimating: ", minfo$file, "...")
  
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
    
    deltas_table_part <- filter(deltas, variable %in% variables)
    
    
    

    setup <- bvar_conj_setup(D, p=n_lag, v_prior=v_prior,
                delta = deltas_table_part$delta,
                lambda = c(l_1, l_power, l_sc, l_io, l_const, l_exo), 
                s2_lag = num_AR_lags, carriero_hack=carriero_hack)

    
    # priors$X_dummy <- NULL
    # priors$Y_dummy <- NULL
    # estimate model
    
    model <- bvar_conj_estimate(setup=setup, verbose=verbose, keep=keep)
  
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
    model <- data_frame(variables=variables)
    model$coef <- model_vector
    status <- "estimated"
  }
  
  if (minfo$type=="rw") {
    # just mean growth rate
    model_vector <- c(t((tail(D,1)-head(D,1))/ (nrow(D)-1)))
    model <- data_frame(variables=variables)
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
                            progress_bar=TRUE, verbose=FALSE) {
  start_time <- Sys.time()
  
  parallel <- match.arg(parallel)
  
  mlist_todo <- mlist
  if (no_reestimation) mlist_todo <- dplyr::filter(mlist, !status=="estimated")
  
  model_ids <- mlist$id
  

  requested_packages <- c("mvtnorm","MHadaptive","MCMCpack","bvarr")
  
  if (parallel=="windows") {
    library("doSNOW")
    cl <- makeCluster(ncpu, outfile="") # number of CPU cores
    registerDoSNOW(cl)
    
    foreach(i=model_ids, .packages=requested_packages) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, do_log=do_log, test=test, verbose=verbose)
      mlist$status[mlist$id==i] <- status
    }
    stopCluster(cl)
  }
  
  if (parallel=="unix") {
    library("doMC")
    registerDoMC(ncpu) # number of CPU cores

    foreach(i=model_ids, .packages=requested_packages) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info, do_log=do_log, test=test, verbose=verbose)
      mlist$status[mlist$id==i] <- status
      
    }
  }
  
  if (parallel=="off") {
    
    if (length(model_ids)>0) {
      if (progress_bar) pb <- txtProgressBar(min = 1, max = length(model_ids), style = 3)
      for (i in 1:length(model_ids))  {
        model_info <- mlist_todo %>% dplyr::filter(id==model_ids[i])
        status <- estimate_model(model_info, do_log=do_log, verbose=verbose)
        mlist$status[mlist$id==i] <- status
        
        if (progress_bar) setTxtProgressBar(pb, i)
      }
      if (progress_bar) close(pb)
    }
  }
  
  end_time <- Sys.time()
  message("Time elapsed: ", round(end_time - start_time,1)," ", attr(end_time - start_time,"units"))
  
  return(mlist) # statuses are updated
}

########################################################################
################ forecasting functions
########################################################################

# function to make forecast of one model for one dataset
# pred_info - one line describing desired forecast
# mlist - list of all models
forecast_model <- function(pred_info, mlist, parallel = parallel, 
                            ncpu=4, test=FALSE, do_log=FALSE, verbose=FALSE) {
  
  
  # get info about model:
  model_id <- pred_info$model_id 
  minfo <- mlist %>% filter(id==model_id) 
  model_full_path <- paste0("../estimation/models/",minfo$file)
  if (verbose) message("Forecasting: ",minfo$file)
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
      answer <- data_frame(value=value, t=t, 
                           variable = rep(variables, Tf_length), h=NA)
  }

  if ((minfo$type=="rw") & (pred_info$type=="in-sample")) {
    Tf_start <- T_start + 1 # is not possible to forecast first observation
    Tf_length <- T_in - 1
    Tf_end <- Tf_start + Tf_length - 1
    
    ## multivariate analog of simple idea: y_lagged + Delta_y
    y_lagged <- c(t(D[(Tf_start-1):(Tf_end-1),])) # c(t()) is a transformation of data.frame row into a vector
    value <- y_lagged + model$coef[rep(1:n_vars,Tf_length)] 
    t <- rep(Tf_start:Tf_end, each=n_vars)
    answer <- data_frame(value=value, t=t, 
                         variable = rep(variables, Tf_length), h=NA)
  }
  
  if ((minfo$type=="conjugate") & (pred_info$type=="in-sample")) {
    n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + n_lag # is not possible to forecast first n_lag observation
    Tf_length <- T_in - n_lag
    Tf_end <- Tf_start + Tf_length - 1
    
    predictions <- bvar_conj_forecast(model, 
                  fast_forecast = TRUE, out_of_sample = FALSE)
    # in-sample forecasts are one-step predictions, so fast_forecast is TRUE
    answer <- mutate(predictions, t=h+Tf_start-1, h=NA) %>% select(-what) 
  }
  
  if ((minfo$type=="var") & (pred_info$type=="in-sample")) {
    n_lag <- as.numeric(minfo$n_lag)
    Tf_start <- T_start + n_lag # is not possible to forecast first n_lag observation
    Tf_length <- T_in - n_lag
    Tf_end <- Tf_start + Tf_length - 1
    
    predictions <- as.data.frame(fitted(model))
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
    
    predictions <- bvar_conj_forecast(model, h = pred_info$h, include="mean", level=NULL, 
                                      fast_forecast = fast_forecast)
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
  
  answer <- mutate(answer, variable=as.character(variable)) %>% 
    rename(forecast=value)
  
  return(answer) 
}





# function to make forecasts of many model for many datasets
forecast_models <- function(plist, mlist, parallel = c("off","windows","unix"), 
                            ncpu=4, do_log=FALSE,
                            progress_bar=TRUE, verbose=verbose) {
  start_time <- Sys.time()
  
  parallel <- match.arg(parallel)
  
  answer <- NULL
  if (parallel=="off") {
    all_data <- list()
    if (progress_bar) pb <- txtProgressBar(min = 1, max = nrow(plist), style = 3)
    for (i in 1:nrow(plist)) {
      all_data[[i]] <- forecast_model(plist[i,],mlist=mlist, parallel = parallel,
                                      ncpu=ncpu, verbose=verbose, do_log=do_log)
      if (progress_bar) setTxtProgressBar(pb, i)
    }
    if (progress_bar) close(pb)
    answer <- data.table::rbindlist(all_data, use.names = TRUE) %>% as.data.frame()
  }
  
  end_time <- Sys.time()
  message("Time elapsed: ", round(end_time - start_time,1)," ", attr(end_time - start_time,"units"))
  return(answer) 
}

# forecasted_models -- data frame with "t", "forecast", "variable" and "model_id"
# var_name -- name of variable to forecast
# mod_id --- id of models 
plot_forecast <- function(forecasted_models, var_name, mod_id ) {
  forecasts_filtered <- filter(forecasted_models, variable==var_name, model_id==mod_id)
  actual_filtered <- filter(actual_obs, variable==var_name)
  data <- left_join(actual_filtered, forecasts_filtered, by="t")
  ggplot(data=data, aes(x=t))+geom_line(aes(y=actual)) + geom_line(aes(y=forecast), col="red")
}


