seas_adjust <- function(df_ts, to_adjust=colnames(df_ts))  {
  df_sa <- df_ts # make a copy
  for (i in 1:ncol(df_ts)) {
    if (names(df_notime)[i] %in% to_adjust) {
      # cat(names(df_notime)[i], i,"\n") 
      df_sa[,i] <- final(seas(df_ts[,i]))
    }
  }
  return(df_sa)
}


fname2id <- function(fname) {
  id <- str_match(pattern="[0-9]+",fname) %>% as.numeric()
  return(id)
}

id2fname0 <- function(id) {
  fname <- paste0(c("mod_id_", rep(0,4-nchar(id)), id, ".Rds"), collapse = "")
  return(fname)
} 

id2fname <- Vectorize(id2fname0)


fmeasure <- function(df_for, df_test, measure="mape", vars=colnames(df_test)) {
  answer <- NULL
  for (i in vars) {
    answer[i] <- ftsa::error(forecast = df_for[,i], 
                             true = df_test[,i], method=measure)
  }  
  return(answer)
}

# fmeasure(df_for, df_test, "mape")

#### describe all models in one table (mod_table)
# p - число лагов 
# v, lambda - параметры из работы
# m - число рядов, считается само
# t_for - длина прогнозного периода
# T - длина ряда, считается сама

bvarm_block <- function(lambda=0.5,p=4,t_for=24,v=1,H3=Inf,varset) {
  new_block <- expand.grid(lambda=lambda, v=v, p=p, t_for=t_for, varset=varset, H3=H3 )
  new_block <- mutate(new_block, varset=as.character(varset),
                      H1 = lambda^2, H2 = v* H1,
                      type="BVARM", estimated=FALSE, T=NA, m=NA)
  
  #model_id <- nrow(mod_table)+(1:nrow(new_block))
  new_block$model_id <- NA
  return(new_block)
}

cvar_block <- function(p=1:4,t_for=24,varset) {
  new_block <- expand.grid(p=p, t_for=t_for, varset=varset)
  new_block <- mutate(new_block, varset=as.character(varset), H1=NA, H2=NA, H3=NA, 
                      v = NA, lambda = NA, type="CVAR", estimated=FALSE, T=NA, m=NA)
  
  #model_id <- nrow(mod_table)+(1:nrow(new_block))
  new_block$model_id <- NA
  return(new_block)
}

# this function detects series with unit roots
# for ADF all three types are supported
# for KPSS only trend and constant
# c_0=0 value for stationary series
# c_1=1 value for non-stationary series
# alpha=5%
bvarm_prior <- function(df_ts, varset=colnames(df_ts), test=c("ADF","KPSS"),
                        type=c("trend", "constant", "neither"), c_0=0, c_1=1) {
  df_sel <- df_ts[,varset]
  nvar <- ncol(df_sel)
  stat_res <- stationarity(df_sel, print = FALSE)
  
  type <- type[1] # take first option
  if (type=="trend") line <- 1
  if (type=="constant") line <- 2
  if (type=="neither") line <- 3
  
  nvar <- ncol(df_sel)
  
  # use time_trend - ADF - 5%
  if(test[1]=="ADF") {
    adf_cr <- stat_res$ADF$`5 Pct`[line]
    adf_obs <- stat_res$ADF[line,1:nvar]  
    result <- ifelse(adf_obs < adf_cr, c_0, c_1)
  }
  
  if(test[1]=="KPSS") {
    kpss_cr <- stat_res$KPSS$`5 Pct`[line]
    kpss_obs <- stat_res$KPSS[line,1:nvar]  
    result <- ifelse(kpss_obs > kpss_cr, c_1, c_0)
  }
  
  return(as.vector(result))
}


build_varset <- function(set_list) {
  varset_table <- NULL
  for (set in set_list) {
    set_vars <- set$vars
    set_name <- set$name
    set_priors <- bvarm_prior(df_ts, set$vars)
    varset_block <- data_frame(var=set_vars, prior = set_priors, set=set_name)
    varset_table <- rbind(varset_table, varset_block)
  }  
  return(varset_table)
}



# работает с BVARM и CVAR 
forinsample <- function(model, t_for=NULL ) {
  beta <- model$Beta 
  
  data <- model$data # original data
  for_data <- data # here we will store forecasts
  
  intercept <- model$constant # if intercept is present
  m <- ncol(data) # number of variables
  t_tot <- nrow(data)
  
  p <- (nrow(beta) - intercept)/m # number of lags  
  if (is.null(t_for)) t_for <- nrow(data) - p # number of forecasted periods 
  
  t_start <- t_tot - t_for + 1 # first forecasted period
  
  for (i in t_start:t_tot) {
    if(intercept) y_f <- beta[1, ] else y_f <- rep(0,m) # add constants if needed    
    
    # !!!!! Error in names in BVARM
    # needs t() !
    for (l in 1:p) y_f <- y_f + t(beta[1 + intercept + (((l-1)*m):(l*m - 1)),]) %*% data[i-l,]    
    for_data[i, ] <- y_f
  }
  return(for_data[t_start:t_tot, ])
}


varset2vars <- function(varset) {
  var_table <- dplyr::filter(varset_table,set==varset)
  return(var_table$var)
}

varset2priors <- function(varset) {
  var_table <- dplyr::filter(varset_table,set==varset)
  return(var_table$prior)
}





forrolling <- function(df_sa, lambda, varset, h, max_h=6, t_win=120, t_for=24, p=5, v=1, H3=Inf) {
  
  var_table <- dplyr::filter(varset_table,set==varset)
  vars <- var_table$var
  coefprior <- var_table$prior
  
  cat("Variables: ",paste0(vars,collapse = ", "),"\n")
  
  df_sel <- df_sa[,vars]
  
  t_tot <- nrow(df_sel)
  t_start_f <- t_tot - t_for + 1 # here all forecasts begins
  
  df_test <- df_sel[t_start_f:t_tot,]
  df_for <- df_test # later we'll fill forecasted values point by point
  
  for (j in t_start_f:t_tot) {
    t_end_est <- j - h
    t_start_est <- t_end_est - t_win -1
    message("forecasting ",j,", estimation from=",t_start_est,", to=",t_end_est)
    
    model <- BVARM(df_sel[t_start_est:t_end_est,], p = p, coefprior = coefprior,
                   HP1 = lambda^2, HP2 = v*lambda^2, HP3=H3,
                   irf.periods = 1 )
    obj_for <- BMR::forecast(model, periods=h, shocks=TRUE, plot=FALSE)
    
    for_one_obs <- tail(matrix(obj_for$PointForecast, ncol=ncol(df_sel)),1)
    # если не поставить в строчке выше matrix, то при прогнозе  h=1
    # результатом будет вектор и tail отрежет не прогнозы для последнего периода для m переменных
    # а только одно число, что бяка
    
    colnames(for_one_obs) <- vars
    df_for[j-t_start_f+1,] <- for_one_obs    
  }
  return(df_for)
}



get_omsfe <- function(df_sa, mod_eval_table, t_win=120, t_for=24) {
  
  omsfe_table <- NULL
  for (j in 1:nrow(mod_eval_table)) {
    
    lambda <- mod_eval_table$lambda[j]
    varset <- mod_eval_table$varset[j]
    h <- mod_eval_table$h[j]
    
    message("-----------")
    message("lambda=",lambda,", varset=",varset,", h=",h, ", total:",j,"/",nrow(mod_eval_table))
    message("-----------")
    
    df_for <- forrolling(df_sa, varset=varset, t_win = t_win, t_for = t_for, 
                         lambda=lambda, h=h) 
    
    vars <- varset2vars(varset)
    df_test <- tail(df_sa[,vars], t_for)
    
    omsfe_model <- fmeasure(df_for, df_test, measure="mse", vars=save_msfe)
    omsfe_add <- data_frame(omsfe=omsfe_model, var=save_msfe, h=h, lambda=lambda, varset=varset)
    
    omsfe_table <- rbind(omsfe_table, omsfe_add)
  }
  return(omsfe_table)
}


mod_size <- function(varset) {
  ans <- as.numeric(str_match(varset,pattern = "[0-9]+"))
  return(ans)
}


add_models <- function(add_block) {
  if ("mod_table.Rds" %in% list.files("./estimation/")) { # if mod_table.Rds exists...
    mod_table <- readRDS("./estimation/mod_table.Rds")
  } else {
    mod_table <- NULL # initial (!)
  }
  
  # I am doing fake rbind to (possibly) change the structure of mod_table and add_block
  mod_table <- rbind_list(mod_table, head(add_block,0))
  add_block <- rbind_list(add_block, head(mod_table,0))
  
  # remove already estimated, or non-estimated but specified in mod_table models
  # from add_block
  j_by <- colnames(mod_table)
  j_by <- j_by[!j_by %in% c("estimated","model_id","T","m")]
  add_block_filtered <- anti_join(add_block, mod_table, by = j_by )
  
  mod_table <- rbind_list(mod_table, add_block_filtered)
  
  mod_table$model_id <- 1:nrow(mod_table)
  saveRDS(mod_table, file="./estimation/mod_table.Rds") 
}

estimate_models <- function() {
  mod_table <- readRDS("./estimation/mod_table.Rds")
  
  non_estimated <- which(mod_table$estimated==FALSE)
  n_2estimate <- length(non_estimated)
  message("I will estimate ", n_2estimate, " models.")
  message("List of ids: ",non_estimated)
  
  
  for (i in non_estimated) {
    
    message("Estimating model ",i," out of ",n_2estimate)
    id <- mod_table$model_id[i]
    mod_line <- mod_table[i,]
    message("type=", mod_line$type, ", H1=",mod_line$H1,", H2=",mod_line$H2,
        ", H3=",mod_line$H3, ", lambda=",mod_line$lambda)
    
    
    var_table <- filter(varset_table,set==mod_line$varset)
    vars <- var_table$var
    message("Variables: ",paste0(vars,collapse = ", "))
    
    
    model <- NULL
    df_for <- NULL
    obj_for <- NULL
    mape <- NA
    
    df_sel <- df_sa[,vars]
    
    m <- ncol(df_sel) # the number of time series
    T <- nrow(df_sel) # total length 
    t_for <- mod_line$t_for
    p <- mod_line$p
    
    mod_table$m[i] <- m
    mod_table$T[i] <- T
    
    
    df_train <- head(df_sel, T-t_for)
    df_test <- tail(df_sel, t_for)
    
    
    set.seed(42)
    if (mod_line$type=="BVARM") {
      HP1 <- mod_line$H1
      HP2 <- mod_line$H2
      HP3 <- mod_line$H3
      coefprior <- var_table$prior
      
      model <- BVARM(df_train, p=p, coefprior=coefprior, HP1=HP1, HP2=HP2, HP3=HP3, irf.periods = 1)
      mod_table$estimated[i] <- TRUE
      # obj_for <- BMR::forecast(model, periods=t_for, shocks=TRUE, plot=FALSE)
      # df_for <- obj_for$PointForecast
      # colnames(df_for) <- vars
      # mape <- mean(fmeasure(df_for, df_test))
    }
    
    if (mod_line$type=="CVAR") {
      model <- CVAR(df_train, p=p)
      mod_table$estimated[i] <- TRUE
      # obj_for <- BMR::forecast(model,periods=t_for ,plot=FALSE,confint=0.95)
      # df_for <- obj_for$PointForecast
      # colnames(df_for) <- vars
      # mape <- mean(fmeasure(df_for, df_test))
    }
    # ouput:
    # model --- estimated model
    # obj_for --- forecast with some junk
    # df_for --- dataset with point forecasts (t_for x m)
    # df_test --- save it for simplicity (t_for x m)
    
    # mod_table$mape[i] <- mape # save mape
    
    
    saveRDS(model, file=paste0("./estimation/models/",id2fname(id)))  
    saveRDS(mod_table, file="./estimation/mod_table.Rds") 
  } 
}





