# translation to R by Boris Demeshev
# error with constant = FALSE in the original script is probably corrected
# based on the code by Gary Koop
# http://personal.strath.ac.uk/gary.koop/bayes_matlab_code_by_koop_and_korobilis.html

Yraw <- read.table("~/Documents/bvarr/bvar_matlab/BVAR_Analytical/Yraw.dat")

# used for mlag2 function
makeblock <- function(X,i,p) {
  Xblock <- X[(p+1-i):(nrow(X)-i),] # get useful lines of X
  Xblock <- as.matrix(Xblock) # assure X is a matrix, not a vector
  Xblock <- rbind(matrix(0,nrow=p,ncol=ncol(X)),Xblock)  # append p zero lines at top
  return(Xblock)
}

mlag2 <- function(X,p) {
  X <- as.matrix(X)
  
  # we need to bind horizontally p blocks
  Xlag <- matrix(nrow=nrow(X),ncol=0) # create empty matrix with correct number of raws
  for (i in 1:p) Xlag <- cbind(Xlag,makeblock(X,i,p)) # bind blocks horizontally
  return(Xlag)  
}


# constant = TRUE/FALSE;        TRUE: if you desire intercepts, FALSE: otherwise 
# p = 2;  Number of lags on dependent variables
# forecasting = TRUE;      TRUE: Compute h-step ahead predictions, FALSE: no prediction
# forecast_method = "direct" Direct forecasts, "iterated"  Iterated forecasts
# h = 4;    Number of forecast periods
# prior "diffuse" "Minnesota" "conjugate" Natural conjugate Prior
# a_bar = c(0.5,0.5,10^2)  # Hyperparameters on the Minnesota variance of alpha

# for debug:
forecasting <- TRUE
constant <- FALSE
forecast_method<-"direct"
h<-4
prior<-"Minnesota"
p<-2
a_bar <- c(0.5, 0.5, 10^2)

bvar_analyt <- function(Yraw,constant=TRUE,p=2,forecasting=TRUE,
                        forecast_method="direct",h=4,
                        prior="Minnesota", a_bar=c(0.5, 0.5, 10^2)) {
  
  Yraw <- as.matrix(Yraw)

  ######### sanity check
  if (forecasting & h<=0) 
    stop("You have set forecasting, but the forecast horizon h<=0")
  if (forecasting & !(forecast_method %in% c("direct","iterated")))
    stop("forecast_method should be 'direct' or 'iterated'")
  
  ######### Data handling
  
  
  #  Get initial dimensions of dependent variable
  Traw <- nrow(Yraw)
  M <- ncol(Yraw)
  

  if (forecasting & forecast_method=="direct") {     
    Y1 <- Yraw[(h+1):Traw,]
    Y2 <- Yraw[2:(Traw-h),]
    Traw <- Traw - h - 1
  }
  if (forecasting & forecast_method=="iterated") {      
    Y1 <- Yraw
    Y2 <- Yraw
  }
  if (!forecasting) {
    Y1 <- Yraw
    Y2 <- Yraw    
  }
  
  # enerate lagged Y matrix. This will be part of the X matrix
  Ylag <-  mlag2(Y2,p); # Y is [T x M]. ylag is [T x (Mp)]
  
  # Now define matrix X which has all the R.H.S. variables (constant, lags of
  # the dependent variable and exogenous regressors/dummies)
  
  X1 <- Ylag[(p+1):Traw,]
  if (constant)
    X1 <- cbind(1, X1) # add column of ones to the left of Ylag 
  
     
  
  # Get size of final matrix X
  Traw3 <- nrow(X1)
  K <- ncol(X1)
  
  # Create the block diagonal matrix Z
  Z1 <- kronecker(diag(M),X1)
  
  # Form Y matrix accordingly
  # Delete first "LAGS" rows to match the dimensions of X matrix
  Y1 <- Y1[(p+1):Traw,] # This is the final Y matrix used for the VAR
  
  # Traw was the dimension of the initial data. T is the number of actual 
  # time series observations of Y and X
  T <- Traw - p
    
  ########## FORECASTING SET-UP:
  # Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
  if (forecasting&forecast_method=="direct")  {
    Y <- head(Y1,-1)                # Y is equal to Y1 without 1 last observation
    X <- head(X1,-1)
    Z <- kronecker(diag(M),X)
    T <- T - 1
  }
  if (forecasting&forecast_method=="iterated")  {              
    Y <- head(Y1,-h)                # Y is equal to Y1 without h last observation
    X <- head(X1,-h)
    Z <- kronecker(diag(M),X)
    T <- T - h
  }
  if (!forecasting) {
    Y <- Y1
    X <- X1
    Z <- Z1
  }
  
  
  #--------------------------------PRIORS------------------------------------
  # First get Ordinary Least Squares (OLS) estimators
  A_OLS <- solve(t(X)%*%X)%*% t(X) %*%Y # This is the matrix of regression coefficients
  a_OLS <- as.vector(A_OLS)     # This is the vector of coefficients, i.e. it holds
  # that a_OLS = vec(A_OLS)
  SSE <- t(Y - X%*%A_OLS)%*%(Y - X%*%A_OLS)
  SIGMA_OLS <- SSE/(T-K)

  #-----------------Prior hyperparameters for bvar model
  
  
  if ( prior == "diffuse" ) {} 
  # I guess there is nothing to specify in this case!
  # Posteriors depend on OLS quantities
  
  if( prior == "conjugate")  { # Normal-Wishart (natural conjugate)
  # Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
    A_prior <- matrix(0,K,M) ### original matlab: = 0*ones(K,M)  
    a_prior <- as.vector(A_prior)
    V_prior <- 10*diag(K)
    # Hyperparameters on solve(SIGMA) ~ W(v_prior,solve(S_prior))
    v_prior <- M
    S_prior <- diag(M)
    inv_S_prior <- solve(S_prior)
  }
  
  if ( prior == "Minnesota" ) {
    
    A_prior <- matrix(0,K,M) # matlab: 0*ones(K,M);   
    a_prior <- as.vector(A_prior) # matlab: A_prior(:);
    
   
    
    # Now get residual variances of univariate p-lag autoregressions. Here
    # we just run the AR(p) model on each equation, ignoring the constant
    # and exogenous variables 
    # (if they have been specified for the original VAR model)
    sigma_sq <- rep(0,M) # matlab: zeros(M,1); # vector to store residual variances
    for (i in 1:M) {
      # Create lags of dependent variable in i-th equation
      Ylag_i <- mlag2(Yraw[,i],p)
      Ylag_i <- Ylag_i[(p+1):Traw,]
      # Dependent variable in i-th equation
      Y_i <- Yraw[(p+1):Traw,i]
      # OLS estimates of i-th equation
      alpha_i <- solve(t(Ylag_i) %*% Ylag_i) %*% t(Ylag_i) %*% Y_i
      sigma_sq[i] <- (1/(T-p+1))*sum((Y_i - Ylag_i %*% alpha_i)^2)
    }
    
    
    # DEBUG ONLY
    # sigma_sq <- 1:M
    
    # Now define prior hyperparameters.
    # Create an array of dimensions K x M, which will contain the K diagonal
    # elements of the covariance matrix, in each of the M equations.
    V_i <- matrix(0,K,M)
    
    # index in each equation which are the own lags
    ind <- matrix(0,M,p)
    for (i in 1:M)
       ind[i,] <- constant+seq(from=i,by=M,len=p) # matlab, i:M:K 
       # but i:M:K sometime produces sequence longer than p
    
    for (i in 1:M) { # for each i-th equation
    for (j in 1:K) {   # for each j-th RHS variable
      if (constant) {
        if (j==1) V_i[j,i] <- a_bar[3]*sigma_sq[i] # variance on constant                
        if (j %in% ind[i,]) 
          V_i[j,i] <- a_bar[1]/(ceiling((j-1)/M)^2) # variance on own lags  
        if ((j>1)& (!(j %in% ind[i,]))) {         
          #for (kj in 1:M) {
          #  if (j %in% ind[kj,]) ll <- kj                               
          #}
          ll <- (j-1) %% M # remainder of division by M
          if (ll==0) ll <- M # replace remainder==0 by M
          
          V_i[j,i] <- (a_bar[2]*sigma_sq[i])/((ceiling((j-1)/M)^2)*sigma_sq[ll]) 
        }
      }
      if (!constant) {
          if (j %in% ind[i,])
             V_i[j,i] <- a_bar[1]/(ceiling(j/M)^2) # variance on own lags
             # !!!! was (j-1)/M in the original script
          if (!(j %in% ind[i,])) {
                 
            ll <- j %% M # remainder of division by M
            if (ll==0) ll <- M # replace remainder==0 by M
                        
            V_i[j,i] <- (a_bar[2]*sigma_sq[i])/((ceiling(j/M)^2)*sigma_sq[ll])            
            # !!!! was (j-1)/M in the original script
            
          }
      }
    }} # i,j cycle
    
    
    # Now V is a diagonal matrix with diagonal elements the V_i
    V_prior <- diag(as.vector(V_i))  # this is the prior variance of the vector a  
    
    # SIGMA is equal to the OLS quantity
    SIGMA <- SIGMA_OLS
    
  } 
  
  
  #============================ POSTERIORS ==================================
  #--------- Posterior hyperparameters of ALPHA and SIGMA with Diffuse Prior
  if (prior == "diffuse") {
    # Posterior of alpha|Data ~ Multi-T(kronecker(SSE,solve(X'X)),alpha_OLS,T-K)
    V_post <- solve(t(X) %*% X)
    a_post <- a_OLS
    A_post <- matrix(as.vector(a_post),nrow=K) # matlab: reshape(a_post,K,M); # ?????????????
                                           
    # posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
    S_post <- SSE
    v_post <- T-K
                                           
    # Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean <- a_post
    alpha_var <- (1/(v_post - M - 1))*kronecker(S_post,V_post)
                                           
  }
  
  #--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
  if (prior == "Minnesota") {
    # ******Get all the required quantities for the posteriors 

    V_post <- solve( solve(V_prior) + kronecker(solve(SIGMA),t(X) %*% X) )
    a_post <- V_post %*% ( solve(V_prior) %*% a_prior + kronecker(solve(SIGMA),t(X) %*% X) %*% a_OLS )
    A_post <- matrix(as.vector(a_post),nrow=K) # matlab: reshape(a_post,K,M)
  
    # In this case, the mean is a_post and the variance is V_post
    alpha_mean <- a_post
    alpha_var <- V_post # added
  }
  #--------- Posterior hyperparameters of ALPHA and SIGMA with Normal-Wishart Prior
  if (prior == "conjugate") {
    # ******Get all the required quantities for the posteriors       
    # For alpha
    V_post <- solve( solve(V_prior) + t(X) %*% X )
    A_post <- V_post %*% ( solve(V_prior) %*% A_prior + t(X) %*% X %*% A_OLS )
    a_post <- as.vector(A_post) # matlab: A_post(:)
  
    # For SIGMA
    S_post <- SSE + S_prior + t(A_OLS) %*% t(X) %*% X %*% A_OLS + 
      t(A_prior) %*% solve(V_prior) %*% A_prior - 
      t(A_post) %*% ( solve(V_prior) + t(X) %*% X ) %*% A_post
    v_post <- T + v_prior
    # Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean <- a_post
    alpha_var <- (1/(v_post - M - 1))*kronecker(S_post,V_post) 
  }
  
  #==========================================================================
  
  X_tplus1 <- c(Y[T,], X[T,2:(M*(p-1)+1)])
  if (constant) X_tplus1 <- c(1, X_tplus1)
  # was always with 1 in the original script
  
  # As a point forecast use predictive mean 
  
  #print(length(X_tplus1))
  #print(dim(A_post))
  
  Pred_mean <- X_tplus1 %*% A_post
  
  
  
  return(list(alpha_mean=alpha_mean,alpha_var=alpha_var,Pred_mean=Pred_mean))
}


bvar_analyt(Yraw,constant=TRUE)

# 1. check with constant = FALSE !!!!!!!!!!!! ??? works but not sure whether it is correct
# 2. where is alpha_var for Minnesota??? --- V_post! (ok)
# 3. separate forecasting from estimation
# 4. use msbvar syntax


