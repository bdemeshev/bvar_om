% Define hyperparameters under different priors
% I indicate the exact line in which the prior hyperparameters can be found


if prior == 1 % Diffuse
    % I guess there is nothing to specify in this case!
    % Posteriors depend on OLS quantities
elseif prior == 2 % Minnesota-Whishart
    % Prior mean on VAR regression coefficients
    A_prior = [zeros(1,M); 0.9*eye(M); zeros((p-1)*M,M)];  %<---- prior mean of ALPHA (parameter matrix) 
    a_prior = A_prior(:);               %<---- prior mean of alpha (parameter vector)
    
    % Minnesota Variance on VAR regression coefficients
    % First define the hyperparameters 'a_bar_i'
    a_bar_1 = 0.5;
    a_bar_2 = 0.5;
    a_bar_3 = 10^2;
    
    % Now get residual variances of univariate p_MIN-lag autoregressions. Here
    % we just run the AR(p) model on each equation, ignoring the constant
    % and exogenous variables (if they have been specified for the original
    % VAR model)
    p_MIN = 6;
    sigma_sq = zeros(M,1); % vector to store residual variances
    for i = 1:M
        % Create lags of dependent variables   
        Ylag_i = mlag2(Yraw(:,i),p);
        Ylag_i = Ylag_i(p_MIN+1:Traw,:);
        X_i = [ones(Traw-p_MIN,1) Ylag_i];
        Y_i = Yraw(p_MIN+1:Traw,i);
        % OLS estimates of i-th equation
        alpha_i = inv(X_i'*X_i)*(X_i'*Y_i);
        sigma_sq(i,1) = (1./(Traw-p_MIN))*(Y_i - X_i*alpha_i)'*(Y_i - X_i*alpha_i);
    end
    % Now define prior hyperparameters.
    % Create an array of dimensions K x M, which will contain the K diagonal
    % elements of the covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    
    % index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if constant==1
                if j==1 % if there is constant, use this code
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil((j-1)/M)^2); % variance on own lags           
                    % Note: the "ceil((j-1)/M)" command finds the associated lag 
                    % number for each parameter
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil((j-1)/M)^2)*sigma_sq(ll,1));           
                end
            else   % if no constant is defined, then use this code
                if find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil(j/M)^2); % variance on own lags
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil(j/M)^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    V_prior = diag(V_i(:));  % this is the prior variance of the vector alpha
    
    %NOTE: No prior for SIGMA. SIGMA is simply a diagonal matrix with each
    %diagonal element equal to sigma_sq(i). See Kadiyala and Karlsson (1997)
    SIGMA = diag(sigma_sq);
    
elseif prior == 3 % Normal-Whishart
    % Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
    A_prior = 0*ones(K,M);   %<---- prior mean of ALPHA (parameter matrix)
    a_prior = A_prior(:);    %<---- prior mean of alpha (parameter vector)
    V_prior = 10*eye(K);     %<---- prior variance of alpha
    
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior = M+1;                 %<---- prior Degrees of Freedom (DoF) of SIGMA
    S_prior = eye(M);      %<---- prior scale of SIGMA
    inv_S_prior = inv(S_prior);    
    
elseif prior == 4  % Independent Normal-Wishart
    n = K*M; % Total number of parameters (size of vector alpha)
    a_prior = 0*ones(n,1);   %<---- prior mean of alpha (parameter vector)
    V_prior = 10*eye(n);     %<---- prior variance of alpha
    
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior = M+1;             %<---- prior Degrees of Freedom (DoF) of SIGMA
    S_prior = eye(M);            %<---- prior scale of SIGMA
    inv_S_prior = inv(S_prior);
    
elseif prior == 5 || prior == 6 % SSVS on alpha, Wishart or SSVS on SIGMA    
    n = K*M; % Total number of parameters (size of vector alpha)
    % mean of alpha
    a_prior = zeros(n,1);
    
    % This is the std of the OLS estimate ofalpha. You can use this to 
    % scale tau_0 and tau_1 (see below) if you want.
    sigma_alpha = sqrt(diag(kron(SIGMA,inv(X'*X))));
    % otherwise, set ' sigma_alpha = ones(n,1); '
    
    % SSVS variances for alpha
    tau_0 = 0.1*sigma_alpha;   % Set tau_[0i], tau_[1i]
    tau_1 = 10*sigma_alpha;
    
    % Priors on SIGMA
    if prior == 6 %SSVS on SIGMA
        % SSVS variances for non-diagonal elements of SIGMA     
        kappa_0 = 0.1; % Set kappa_[0ij], kappa_[1ij]
        kappa_1 = 6;
   
        % Hyperparameters for diagonal elements of SIGMA (Gamma density)
        a_i = .01;
        b_i = .01;
    elseif prior == 5  %Wishart on SIGMA
        % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
        v_prior = M+1;             %<---- prior Degrees of Freedom (DoF) of SIGMA
        S_prior = eye(M);            %<---- prior scale of SIGMA  
        inv_S_prior = inv(S_prior);
    end
    
    % Hyperparameters for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
    p_i = .5;
    
    % Hyperparameters for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
    q_ij = .5;
    
    % Initialize Gamma and Omega vectors
    gammas = ones(n,1);       % vector of Gamma
    omega=cell(1,M-1);
    for kk_1 = 1:(M-1)
        omega{kk_1} = ones(kk_1,1);	% Omega_j
    end
    
    % Set space in memory for some vectors that we are using in SSVS
    gamma_draws = zeros(nsave,n); % vector of gamma draws
    omega_draws = zeros(nsave,.5*M*(M-1)); %vector of omega draws    
else
    error('Wrong choice of prior, prior = 1-6 only')
end