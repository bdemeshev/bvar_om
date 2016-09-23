% BVAR_FULL.m 
% This code replicates the results from the 1st empirical illustration 
% in Koop and Korobilis (2009).
%
% You can chose 6 different priors. For some priors, analytical results are
% available, so Monte Carlo Integration is used. For other priors, you need
% to use the Gibbs sampler. For Gibbs sampler models I take a number of
% 'burn-in draws', so that I keep only the draws which have converged.
%
% The specification of the prior hyperparmeters are in the file
% prior_hyper.m. See there for details.
%
% The convention used here is that ALPHA is the K x M matrix of VAR coefficients,
% alpha is the KM x 1 column vector of vectorized VAR coefficients, i.e.
% alpha = vec(ALPHA), and SIGMA is the M x M VAR covariance matrix.
%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models using posterior simulation. Dependent on your choice of forecasting,
% the VAR model is:
%
% In this code we provide direct (as opposed to iterated) forecasts
% Direct h-step ahead foreacsts:
%     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
%
% so that in this case there are also p lags of Y (from 0 to p-1).
%
% In any of the two cases, the model is written as:
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A).
%--------------------------------------------------------------------------
% NOTES: The code sacrifices efficiency for clarity. It follows the
%        theoretical equations in the monograph and the manual.
%
% AUTHORS: Gary Koop and Dimitris Korobilis
% CONTACT: dikorombilis@yahoo.gr
%--------------------------------------------------------------------------

clear all;
clc;
randn('seed',2); %#ok<RAND>
rand('seed',2); %#ok<RAND>

%------------------------------LOAD DATA-----------------------------------
% Load Quarterly US data on inflation, unemployment and interest rate, 
% 1953:Q1 - 2006:Q3
load Yraw.dat;
%Yraw = Yraw(29:189,:);
% or Simulate data from a simple VAR Data Generating process
%[Yraw] = bvardgp();

% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.

%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model


constant = 1;        % 1: if you desire intercepts, 0: otherwise 
p = 4;               % Number of lags on dependent variables
forecasting = 1;     % 1: Compute h-step ahead predictions, 0: no prediction

repfor = 50;         % Number of times to obtain a draw from the predictive 
                     % density, for each generated draw of the parameters                     
h = 1;               % Number of forecast periods
impulses = 1;        % 1: compute impulse responses, 0: no impulse responses
ihor = 24;           % Horizon to compute impulse responses

% Set prior for BVAR model:
prior = 6;  % prior = 1 --> Diffuse ('Jeffreys') (M-C Integration)
            % prior = 2 --> Minnesota            (M-C Integration)
            % prior = 3 --> Normal-Wishart       (M-C Integration)           
            % prior = 4 --> Independent Normal-Wishart      (Gibbs sampler)
            % prior = 5 --> SSVS in mean-Wishart            (Gibbs sampler)
            % prior = 6 --> SSVS in mean-SSVS in covariance (Gibbs sampler)
            

% Gibbs-related preliminaries
nsave = 10000;         % Final number of draws to save
nburn = 2000;         % Draws to discard (burn-in)
% For models using analytical results, there are no convergence issues (You
% are not adviced to change the next 3 lines)
if prior ==1 || prior == 2 || prior == 3
    nburn = 0*nburn;
end
ntot = nsave + nburn;  % Total number of draws
it_print = 2000;       % Print on the screen every "it_print"-th iteration
%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw M] = size(Yraw);

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0 % Check for wrong (incompatible) user input
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end
    % Now create VAR specification according to forecast method
   
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
   
else
   Y1 = Yraw;
   Y2 = Yraw;
end


% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y2,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies).
% Note that in this example I do not include exogenous variables (other macro
% variables, dummies, or trends). You can load a file with exogenous
% variables, call them, say W, and then extend variable X1 in line 133, as:
%            X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:) W(p+1:Traw,:)];
% and line 135 as:
%            X1 = [Ylag(p+1:Traw,:)  W(p+1:Traw,:)];
if constant
    X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
else
    X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
end

% Get size of final matrix X
[Traw3 K] = size(X1);

% Create the block diagonal matrix Z
Z1 = kron(eye(M),X1);

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if forecasting==1
    Y_pred = zeros(nsave*repfor,M); % Matrix to save prediction draws
    PL = zeros(nsave,1);            % Matrix to save Predictive Likelihood
    
    % Direct forecasts, we only need to keep the last observation for evaluation
        Y = Y1(1:end-1,:);          
        X = X1(1:end-1,:);
        Z = kron(eye(M),X);
        T = T - 1;
   
else % if no prediction is present, keep all observations
    Y = Y1;
    X = X1;
    Z = Z1;
end

%========= IMPULSE RESPONSES SET-UP:
% Create matrices to store forecasts

if impulses == 1;
    imp_infl = zeros(nsave,M,ihor); % impulse responses to a shock in inflation
    imp_une = zeros(nsave,M,ihor);  % impulse responses to a shock in unemployment
    imp_int = zeros(nsave,M,ihor);  % impulse responses to a shock in interest rate
    bigj = zeros(M,M*p);
    bigj(1:M,1:M) = eye(M);
end

%-----------------------------PRELIMINARIES--------------------------------
% First get ML estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of parameters, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);   % Sum of squared errors
SIGMA_OLS = SSE./(T-K+1);

% Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;     % This is the single draw from the posterior of alpha
ALPHA = A_OLS;     % This is the single draw from the posterior of ALPHA
SSE_Gibbs = SSE;   % This is the SSE based on each draw of ALPHA
SIGMA = SIGMA_OLS; % This is the single draw from the posterior of SIGMA
IXY =  kron(eye(M),(X'*Y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% translation line %%%%%%%%%%%%%%%%%%%%%
% Storage space for posterior draws
alpha_draws = zeros(nsave,K*M);   % save draws of alpha
ALPHA_draws = zeros(nsave,K,M);   % save draws of alpha
SIGMA_draws = zeros(nsave,M,M);   % save draws of ALPHA

%-----------------Prior hyperparameters for bvar model
% load file which sets hyperparameters for chosen prior
prior_hyper;
%-------------------- Prior specification ends here
    
%========================== Start Sampling ================================