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
nsave = 10;         % Final number of draws to save
nburn = 5;         % Draws to discard (burn-in)
% For models using analytical results, there are no convergence issues (You
% are not adviced to change the next 3 lines)
if prior ==1 || prior == 2 || prior == 3
    nburn = 0*nburn;
end
ntot = nsave + nburn;  % Total number of draws
it_print = 1;       % Print on the screen every "it_print"-th iteration
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


% Storage space for posterior draws
alpha_draws = zeros(nsave,K*M);   % save draws of alpha
ALPHA_draws = zeros(nsave,K,M);   % save draws of alpha
SIGMA_draws = zeros(nsave,M,M);   % save draws of ALPHA

%-----------------Prior hyperparameters for bvar model
% load file which sets hyperparameters for chosen prior
prior_hyper;
%-------------------- Prior specification ends here
    
%========================== Start Sampling ================================
%==========================================================================
tic;
disp('Number of iterations');
for irep = 1:ntot  %Start the Gibbs "loop"
    if mod(irep,it_print) == 0 % print iterations
        disp(irep);
        toc;
    end
    
    %--------- Draw ALPHA and SIGMA with Diffuse Prior
    if prior == 1
        % Posterior of alpha|SIGMA,Data ~ Normal
        V_post = kron(SIGMA,inv(X'*X));
        alpha = a_OLS + chol(V_post)'*randn(K*M,1);% Draw alpha
        ALPHA = reshape(alpha,K,M); % Create draw of ALPHA       
        
        % Posterior of SIGMA|Data ~ iW(SSE_Gibbs,T-K) 
        SIGMA = inv(wish(inv(SSE_Gibbs),T-K));% Draw SIGMA
        
    %--------- Draw ALPHA and SIGMA with Minnesota Prior
    elseif prior == 2
        %Draw ALPHA
        for i = 1:M
            V_post = inv( inv(V_prior((i-1)*K+1:i*K,(i-1)*K+1:i*K)) + inv(SIGMA(i,i))*X'*X );
            a_post = V_post*(inv(V_prior((i-1)*K+1:i*K,(i-1)*K+1:i*K))*a_prior((i-1)*K+1:i*K,1) + inv(SIGMA(i,i))*X'*Y(:,i));
            alpha((i-1)*K+1:i*K,1) = a_post + chol(V_post)'*randn(K,1); % Draw alpha
        end
        ALPHA = reshape(alpha,K,M); % Create draw in terms of ALPHA
        
        % SIGMA in this case is a known matrix, whose form is decided in
        % the prior (see prior_hyper.m)
        
    %--------- Draw ALPHA and SIGMA with Normal-Wishart Prior
    elseif prior == 3
        % ******Get all the required quantities for the posteriors       
        V_post = inv( inv(V_prior) + X'*X );
        A_post = V_post*(inv(V_prior)*A_prior + X'*X*A_OLS);
        a_post = A_post(:);
    
        S_post = SSE + S_prior + A_OLS'*X'*X*A_OLS + A_prior'*inv(V_prior)*A_prior - A_post'*(inv(V_prior) + X'*X)*A_post;
        v_post = T + v_prior;
    
        % This is the covariance for the posterior density of alpha
        COV = kron(SIGMA,V_post);
    
        % Posterior of alpha|SIGMA,Data ~ Normal
        alpha = a_post + chol(COV)'*randn(K*M,1);  % Draw alpha
        ALPHA = reshape(alpha,K,M); % Draw of ALPHA
        
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
        SIGMA = inv(wish(inv(S_post),v_post));% Draw SIGMA
        
    %--------- Draw ALPHA and SIGMA with Independent Normal-Wishart Prior
    elseif prior == 4
        VARIANCE = kron(inv(SIGMA),eye(T));
        V_post = inv(V_prior + Z'*VARIANCE*Z);
        a_post = V_post*(V_prior*a_prior + Z'*VARIANCE*Y(:));
        alpha = a_post + chol(V_post)'*randn(n,1); % Draw of alpha
        
        ALPHA = reshape(alpha,K,M); % Draw of ALPHA
        
        % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
        v_post = T + v_prior;
        S_post = S_prior + (Y - X*ALPHA)'*(Y - X*ALPHA);
        SIGMA = inv(wish(inv(S_post),v_post));% Draw SIGMA
        
    %--------- Draw ALPHA and SIGMA using SSVS prior 
    elseif prior == 5 || prior == 6
        % Draw SIGMA
        if prior == 5 % Wishart
            % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
            v_post = T + v_prior;
            S_post = inv_S_prior + (Y - X*ALPHA)'*(Y - X*ALPHA);
            SIGMA = inv(wish(inv(S_post),v_post));% Draw SIGMA
        elseif prior == 6 % SSVS
            % Draw psi|alpha,gamma,omega,DATA from the GAMMA dist.
            % Get S_[j] - upper-left [j x j] submatrices of SSE
            % The following loop creates a cell array with elements S_1,
            % S_2,...,S_j with respective dimensions 1x1, 2x2,...,jxj
            S=cell(1,M);
            for kk_2 = 1:M                                       
                S{kk_2} = SSE_Gibbs(1:kk_2,1:kk_2);   
            end
            % Set also SSE =(s_[i,j]) & get vectors s_[j]=(s_[1,j] , ... , s_[j-1,j])
            s=cell(1,M-1);
            for kk_3 = 2:M
                s{kk_3 - 1} = SSE_Gibbs(1:(kk_3 - 1),kk_3);
            end
            % Parameters for Heta|omega ~ N_[j-1](0,D_[j]*R_[j]*D_[j]), see eq. (15)
            % Create and update h_[j] matrix
            % If omega_[ij] = 0 => h_[ij] = kappa0, else...
            hh=cell(1,M-1);
            for kk_4 = 1:M-1                  
                omeg = cell2mat(omega(kk_4));
                het = cell2mat(hh(kk_4));
                for kkk = 1:size(omeg,1)           
                    if omeg(kkk,1) == 0
                        het(kkk,1) = kappa_0;
                    else                        
                        het(kkk,1) = kappa_1;               
                    end
                end
                hh{kk_4} = het;
            end            
            % D_j = diag(hh_[1j],...,hh_[j-1,j])
            D_j=cell(1,M-1);
            for kk_5 = 1:M-1           
                D_j{kk_5} = diag(cell2mat(hh(kk_5)));
            end
            % Now create covariance matrix D_[j]*R_[j]*D_[j], see eq. (15)
            DD_j=cell(1,M-1);
            for kk_6 = 1:M-1
                DD = cell2mat(D_j(kk_6));
                DD_j{kk_6} = (DD*DD);
            end
            % Create B_[i] matrix
            B=cell(1,M);
            for rr = 1:M           
                if rr == 1
                    B{rr} = b_i + 0.5*(SSE(rr,rr));
                elseif rr > 1
                    s_i = cell2mat(s(rr-1));
                    S_i = cell2mat(S(rr-1));
                    DiDi = cell2mat(DD_j(rr-1));
                    rr;
                    disp(size(s_i));
                    disp(size(S_i));
                    disp(size(DiDi));
                    disp('-----');
                    B{rr} = b_i + 0.5*(SSE_Gibbs(rr,rr) - s_i'*inv(S_i + inv(DiDi))*s_i);
                end
            end
            % Now get B_i from cell array B, and generate (psi_[ii])^2
            B_i = cell2mat(B);
            psi_ii_sq = zeros(M,1);
            for kk_7 = 1:M	              
                psi_ii_sq(kk_7,1) = gamm_rnd(1,1,(a_i + 0.5*T),B_i(1,kk_7));
            end
            
            % Draw eta|psi,phi,gamma,omega,DATA from the [j-1]-variate
            % NORMAL dist.
            eta = cell(1,M-1);
            for kk_8 = 1:M-1       
                s_i = cell2mat(s(kk_8));
                S_i = cell2mat(S(kk_8));
                DiDi = cell2mat(DD_j(kk_8));
                miu_j = - sqrt(psi_ii_sq(kk_8+1))*(inv(S_i + inv(DiDi))*s_i);
                Delta_j = inv(S_i + inv(DiDi));
                
                eta{kk_8} = miu_j + chol(Delta_j)'*randn(kk_8,1);
            end
           
            % Draw omega|eta,psi,phi,gamma,omega,DATA from BERNOULLI dist.
            omega_vec = []; %temporary vector to store draws of omega
            for kk_9 = 1:M-1       
                omeg_g = cell2mat(omega(kk_9));
                eta_g = cell2mat(eta(kk_9));
                for nn = 1:size(omeg_g)  % u_[ij1], u_[ij2], see eqs. (32 - 33)                          
                    u_ij1 = (1./kappa_0)*exp(-0.5*((eta_g(nn))^2)./((kappa_0)^2))*q_ij;
                    u_ij2 = (1./kappa_1)*exp(-0.5*((eta_g(nn))^2)./((kappa_1)^2))*(1-q_ij);
                    ost = u_ij1./(u_ij1 + u_ij2);
                    omeg_g(nn,1) = bernoullirnd(ost);
                    omega_vec = [omega_vec ; omeg_g(nn,1)]; %#ok<AGROW>
                end
                omega{kk_9} = omeg_g; %#ok<AGROW>
            end
            
            % Create PSI matrix from individual elements of "psi_ii_sq" and "eta"
            PSI_ALL = zeros(M,M);
            for nn_1 = 1:M  % first diagonal elements
                PSI_ALL(nn_1,nn_1) = sqrt(psi_ii_sq(nn_1,1));   
            end
            for nn_2 = 1:M-1 % Now non-diagonal elements
                eta_gg = cell2mat(eta(nn_2));
                for nnn = 1:size(eta_gg,1)
                    PSI_ALL(nnn,nn_2+1) = eta_gg(nnn);
                end
            end
            % Create SIGMA
            SIGMA = inv(PSI_ALL*PSI_ALL');        
        end % END DRAWING SIGMA 
            
        % Draw alpha              
        % Hyperparameters for alpha|gamma ~ N_[m](0,D*D)
        h_i = zeros(n,1);   % h_i is tau_0 if gamma=0 and tau_1 if gamma=1
        for nn_3 = 1:n
            if gammas(nn_3,1) == 0               
                h_i(nn_3,1) = tau_0(nn_3);
            elseif gammas(nn_3,1) == 1        
                h_i(nn_3,1) = tau_1(nn_3);       
            end
        end
        D = diag(h_i'*eye(n)); % Create D. Here D=diag(h_i) will also do
        DD = D*D;   % Prior covariance matrix for Phi_m
        isig=inv(SIGMA);
        psi_xx = kron(inv(SIGMA),(X'*X));
        V_post = inv(psi_xx + inv(DD));
%        a_post = V_post*((psi_xx)*a_OLS + (inv(DD))*a_prior);
        
         visig=isig(:);
       a_post = V_post*(IXY*visig + inv(DD)*a_prior);
        alpha = a_post + chol(V_post)'*randn(n,1); % Draw alpha
        
        alpha = a_post + chol(V_post)'*randn(n,1); % Draw alpha
   
        ALPHA = reshape(alpha,K,M); % Draw of ALPHA
        
        % Draw gamma|phi,psi,eta,omega,DATA from BERNOULLI dist.    
        for nn_6 = 1:n
            u_i1 = (1./tau_0(nn_6))*exp(-0.5*(alpha(nn_6)./(tau_0(nn_6)))^2)*p_i;           
            u_i2 = (1./tau_1(nn_6))*exp(-0.5*(alpha(nn_6)./(tau_1(nn_6)))^2)*(1-p_i);
            gst = u_i1./(u_i1 + u_i2);
            gammas(nn_6,1) = bernoullirnd(gst); %#ok<AGROW>
        end
        
        % Save new Sum of Squared Errors (SSE) based on draw of ALPHA  
        SSE_Gibbs = (Y - X*ALPHA)'*(Y - X*ALPHA);
    
    end
    % =============Estimation ends here
    
    
    % ****************************|Predictions, Responses, etc|***************************
    if irep > nburn               
        %=========FORECASTING:
        if forecasting==1
          
                Y_temp = zeros(repfor,M);
                % compute 'repfor' predictions for each draw of ALPHA and SIGMA
                for ii = 1:repfor
                    X_fore = [Y(T,:) X(T,2:M*(p-1)+1)];
                    if constant == 1
                      X_fore = [1 X_fore];
                    end
                    % Forecast of T+1 conditional on data at time T
                    Y_temp(ii,:) = X_fore*ALPHA + randn(1,M)*chol(SIGMA);
                end
                % Matrix of predictions
                Y_pred(((irep-nburn)-1)*repfor+1:(irep-nburn)*repfor,:) = Y_temp;
                % Predictive likelihood
                PL(irep-nburn,:) = mvnpdf(Y1(T+1,:),X(T,:)*ALPHA,SIGMA);
                if PL(irep-nburn,:) == 0
                    PL(irep-nburn,:) = 1;
                end
          
        end % end forecasting
        %=========Forecasting ends here
        
        %=========IMPULSE RESPONSES:
        if impulses==1
            % ------------Identification code I:
            Bv = zeros(M,M,p);
            for i_1=1:p
                alpha_index = (i_1-1)*M + 1:i_1*M;
                if constant==1
                    alpha_index = alpha_index + 1
                end
                Bv(:,:,i_1) = ALPHA(alpha_index,:);               
            end
                        
            % st dev matrix for structural VAR
            shock = chol(SIGMA)';
            d = diag(diag(shock));
            shock = inv(d)*shock;
            
            [responses]=impulse(Bv,shock,ihor);
            
            % Restrict to policy shocks
            responses1 = squeeze(responses(:,1,:));
            responses2 = squeeze(responses(:,2,:));
            responses3 = squeeze(responses(:,3,:));
            
            imp_infl(irep-nburn,:,:) = responses1;
            imp_une(irep-nburn,:,:) = responses2;
            imp_int(irep-nburn,:,:) = responses3;
                    


        end
               
        %----- Save draws of the parameters
        alpha_draws(irep-nburn,:) = alpha;
        ALPHA_draws(irep-nburn,:,:) = ALPHA;
        SIGMA_draws(irep-nburn,:,:) = SIGMA;
        if prior == 5 || prior == 6
            gamma_draws(irep-nburn,:) = gammas; %#ok<AGROW>
            if prior == 6
                omega_draws(irep-nburn,:) = omega_vec; %#ok<AGROW>
            end
        end

    end % end saving results
       
end %end the main Gibbs for loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% translation line %%%%%%%%%%%%%%%%%%%%%
%====================== End Sampling Posteriors ===========================

%Posterior mean of parameters:
ALPHA_mean = squeeze(mean(ALPHA_draws,1)); %posterior mean of ALPHA
SIGMA_mean = squeeze(mean(SIGMA_draws,1)); %posterior mean of SIGMA

%Posterior standard deviations of parameters:
ALPHA_std = squeeze(std(ALPHA_draws,1)); %posterior std of ALPHA
SIGMA_std = squeeze(std(SIGMA_draws,1)); %posterior std of SIGMA

%or you can use 'ALPHA_COV = cov(alpha_draws,1);' to get the full
%covariance matrix of the posterior of alpha (of dimensions [KM x KM] )

if prior == 5 || prior == 6   
    % Find average of restriction indices Gamma
    gammas = mean(gamma_draws,1);
    gammas_mat = reshape(gammas,K,M);
    if prior == 6
        % Find average of restriction indices Omega
        omega = mean(omega_draws,1)';
        omega_mat = zeros(M,M);
        for nn_5 = 1:M-1
            ggg = omega(((nn_5-1)*(nn_5)/2 + 1):(nn_5*(nn_5+1)/2),:);
            omega_mat(1:size(ggg,1),nn_5+1) = ggg;
        end
    end
end


% mean prediction and log predictive likelihood
if forecasting == 1
    Y_pred_mean = mean(Y_pred,1); % mean prediction
    Y_pred_std = std(Y_pred,1);   % std prediction
    log_PL = mean((log(PL)),1);

    %This are the true values of Y at T+h:
 
        true_value = Y1(T+1,:);
    
    %(subsequently you can easily also get MSFE and MAFE)
 
    %======PLOT posterior predictive
    figure
    bars = 3000;
    subplot(3,1,1)
    hist(Y_pred(:,1),bars);
    title('Inflation')
    text(Y_pred_mean(:,1),max(hist(Y_pred(:,1),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,1)) ', std = ' num2str(Y_pred_std(:,1))],...       
        'HorizontalAlignment','left')
    subplot(3,1,2)
    hist(Y_pred(:,2),bars);
    title('Unemployment')
    text(Y_pred_mean(:,2),max(hist(Y_pred(:,2),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,2)) ', std = ' num2str(Y_pred_std(:,2))],...
        'HorizontalAlignment','left')
    subplot(3,1,3)
    hist(Y_pred(:,3),bars);
    title('Interest Rate')
    text(Y_pred_mean(:,3),max(hist(Y_pred(:,3),bars)),['\leftarrow mean = '...
        num2str(Y_pred_mean(:,3)) ', std = ' num2str(Y_pred_std(:,3))],...
        'HorizontalAlignment','left')    
end



% You can also get other quantities, like impulse responses
if impulses==1;
    % Set quantiles from the posterior density of the impulse responses
    qus = [.1, .5, .90];
    imp_resp_infl = squeeze(quantile(imp_infl,qus));       
    imp_resp_une = squeeze(quantile(imp_une,qus));
    imp_resp_int = squeeze(quantile(imp_int,qus));
    
    %======= Plot impulse responses
    figure
    set(0,'DefaultAxesColorOrder',[0 0 0],...
        'DefaultAxesLineStyleOrder','--|-|--')
    subplot(3,3,1)
    plot(squeeze(imp_resp_infl(:,1,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Inflation')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,2)
    plot(squeeze(imp_resp_une(:,1,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Inflation')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,3)
    plot(squeeze(imp_resp_int(:,1,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Interest Rate, Shock to Inflation')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,4)
    plot(squeeze(imp_resp_infl(:,2,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Unemployment')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,5)
    plot(squeeze(imp_resp_une(:,2,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Unemployment')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,6)
    plot(squeeze(imp_resp_int(:,2,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Interest Rate, Shock to Unemployment')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,7)
    plot(squeeze(imp_resp_infl(:,3,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Inflation, Shock to Interest Rate')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,8)
    plot(squeeze(imp_resp_une(:,3,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Unemployment, Shock to Interest Rate')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
    subplot(3,3,9)
    plot(squeeze(imp_resp_int(:,3,:))')
    hold;
    plot(zeros(1,ihor),'-')
    title('Response of Interest Rate, Shock to Interest Rate')
    xlim([1 ihor])
    set(gca,'XTick',0:4:ihor)
end



% Clear screen and print elapsed time
clc;
toc;

% Print some directions to the user
disp('Please find the means and variances of the VAR parameters in the vectors')
disp('ALPHA_mean and ALPHA_std for the VAR regression coefficients, and ')
disp('SIGMA_mean and SIGMA_std for the VAR covariance matrix. The predictive')
disp('mean and standard deviation are in Y_pred_mean and Y_pred_std, respectively.')
disp('The log Predictive Likelihood is given by variable log_PL. The true value')
disp('of y(t+h) is given in the variable true_value. For example the mean squared')
disp('forecast error can be obtained using the command')
disp('                MSFE = (Y_pred_mean - true_value).^2')
disp('If you are using the SSVS prior, you can get the averages of the restriction')
disp('indices $\gamma$ and $\omega$. These are in the variables gammas_mat and omegas_mat') 


if prior == 1; name = 'diffuse';
elseif prior == 2; name = 'minnesota';
elseif prior == 3; name = 'normal_wishart' ;
elseif prior == 4; name = 'indep_normal_wishart';
elseif prior == 5; name = 'SSVS_wishart';
elseif prior == 6; name = 'SSVS_full';
end    
save(sprintf('%s_%s.mat','results',name),'-mat');