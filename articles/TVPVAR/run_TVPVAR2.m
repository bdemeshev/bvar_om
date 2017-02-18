
% Bayesian TVP-VAR à la Primiceri (2005), corrected as in Del Negro and
% Primiceri (2015)
%
% Note: This code replicates the results of Del Negro and Primiceri (2015)
% as outlined in my notes available under www.bkolb.eu/codes/TVPVAR.zip.
% The differences to Del Negro and Primiceri are mainly a different
% drawing of the covariance matrices, a differentiation across elements of
% A(t), and the use of a burn-in. For further information, see my note and
% the mentioned papers.
%
% This code is based on code by Haroon Mumtaz, Dimitris Korobilis, and
% Marco Del Negro and Giorgio Primiceri, as referenced in my note. All
% errors are my own.
%
% By B. Kolb, October 2015
%
% Some comments:
% - there are references to Primiceri (2005), indicated e.g. as "p845".
% - code is hard-wired for one constant, i.e. number of states
%   nB=1+nlag*ny, and for block diagonalisation à la Primiceri (2005). Any
%   other changes should be done in Part A, the rest of the code runs by
%   itself.


% Housekeeping
clear; clc; close all force;


%% Part A: Settings
% Choose some parameters
nlag = 2;       % # lags in the VAR
T0   = 40;      % size of training sample
nrep = 500;     % # iterations for Gibbs sampler (Primiceri: 110000)
burn = 100;     % burn-in period                 (Primiceri: 109000)
thin = 10;      % thinning parameter; set=1 to use all draws
scal = 3.5e-04; % scale down cov matrix of small training sample
Aind = [4 7 8]; % Index for A(t): in 3x3 case, a_21(t) is 4th element etc.
kQ   = 0.01;
kS   = 0.1;
kW   = 0.01;
TQ   = T0;
TS   = [2 3];
TW   = 4;
cbar = 0.001; % offset constant used in drawing Sigma(t), Primiceri p. 846

% Parameters of the seven normals to be mixed to approximate log
% Chi-squared distribution (see Primiceri p. 846f)
qj  = [.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mj  = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] ...
    - 1.2704;  % means already adjusted by - 1.2704
v2j = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];

% For plotting
hor  = 40;      % impulse response horizon
time = 1963.5 : 0.25 : 2001.5; % vector of dates


% -------------------  No changes below! -----------------------------


%% Part B: DATA
load ('data'); % load data from Del Negro and Primiceri (2005)
Y = data;

ny = size(Y,2);      % # variables
X  = [ ones(size(Y,1),1),  makelag(Y,1), makelag(Y,2) ];
nB = size(X,2)*ny;   % # states
nA = ny*(ny-1)/2;    % # lower-triangular elements
Y  = Y(nlag+1:end,:);
X  = X(nlag+1:end,:);


%% Part C: Priors (from OLS estimation on training sample)
% 1: Priors for B(t) and cov(B(t)) by OLS
y0   = Y(1:T0,:);
x0   = X(1:T0,:);
Bpr  = (x0'*x0)\x0'*y0;        % prior for B(t)

r0   = y0-x0*Bpr;              % residuals
Om0  = (r0'*r0)/T0;            % covariance matrix of residuals, Omega(0)
VBpr = kron(Om0,eye(nB/ny)/(x0'*x0));  % prior covariance matrix of B(t)

% 2: Priors for A(t) by OLS on residuals
% use Primiceri formula: A(t)*Omega(t)*A(t)'=Sigma(t)*Sigma(t)'
yA   = reshape(r0,T0*ny,1); % LHS variable in A(t) regression
aux  = kron(eye(ny),r0);
xA   = aux(:,Aind);         % RHS variable in A(t) regression
Apr  = (xA'*xA)\xA'*yA;     % prior for A(t)

eA   = reshape(yA - xA*Apr,T0,ny);
OmA  = (eA'*eA)/T0;
VAx  = (xA'*xA)\xA'*(kron(OmA,eye(T0)))*xA/(xA'*xA);
VApr = zeros(nA);           % prior covariance matrix of A(t)
m    = 1;
for j=1:ny-1
    VApr(m:m+j-1, m:m+j-1) = VAx(m:m+j-1, m:m+j-1);
    m = m + j;
end

% 3. Priors for H(t)
Rpr  = diag(eA'*eA/length(eA));
Hpr  = 0.5*log(Rpr);            % prior for H(t)
VHpr = eye(ny);                 % prior covariance matrix of H(t)

% remove initial sample
Y = Y(T0+1:end,:);
X = X(T0+1:end,:);
T = size(X,1);

B_s  = zeros(T,nB,nrep);
A_s  = zeros(T,nA,nrep+1);
H_s  = zeros(T,ny,nrep+1);
Q_s  = zeros(nB,nB,nrep+1);
S1_s = zeros(1,1,nrep+1);
S2_s = zeros(2,2,nrep+1);
W_s  = zeros(ny,ny,nrep+1);

%% Part D: Initialisation
% initialise Q
covQ = TQ*VBpr*kQ^2;
Q    = iwishrnd(covQ,TQ);

% initialise S
covS2  = TS(1)*VApr(1,1)*kS^2;
S2     = iwishrnd(covS2,TS(1));
S(1,1) = S2;

covS3      = TS(2)*VApr(2:3,2:3)*kS^2;
S3         = iwishrnd(covS3,TS(2));
S(2:3,2:3) = S3;

% initialise W
covW = TW*VHpr*kW^2;
W    = iwishrnd(covW,TW);

% initialise H and Sigma(t)
H    = NaN(T,ny);
Sig  = NaN(2,2,T);
H(1,:)      = mvnrnd(Hpr,VHpr,1);
Sig(:,:,1)  = diag(exp(2*H(1,2:3)));
for t=2:T
    H(t,:)      = mvnrnd(H(t-1,:),W,1);
    Sig(:,:,t)  = diag(exp(2*H(t,2:3)));
end

% initialise A
A      = NaN(T,nA);
A(1,:) = mvnrnd(Apr,4*VApr,1);
for t=2:T
    A(t,:) = mvnrnd(A(t-1,:),S,1);
end

VBpr = 4*VBpr;


%% Part E: Initialise (note difference to p.77)
A_s(:,:,1)   = A;
H_s(:,:,1)   = H;
Q_s(:,:,1)   = Q;
S2_s(:,:,1)  = S2;
S3_s(:,:,1)  = S3;
S_s(:,:,1)   = S;
W_s(:,:,1)   = W;


%% Part F: Gibbs sampling algorithm
tic;
wbGS = waitbar(0,'Running the Gibbs sampler');

for m = 1:nrep
    waitbar(m/nrep)
    
    
    %% Stage 1: Draw beta(t) and Q
    % Step 1.1: Get draws for beta and VAR residuals by backward recursions
    % a) get Omega from   Omega(t) = A(t)^-1 * H(t)*H(t)' * (A(t)')^-1  :
    Omega = NaN(ny,ny,T);
    Alt   = kron(ones(T,1),eye(ny)); % lower-triangular version of A
    for t = 1:T
        Alt((t-1)*ny+2,1) = -A(t,1); % fill alpha_21 elements
        Alt((t-1)*ny+3,1) = -A(t,2); % fill alpha_31 elements
        Alt((t-1)*ny+3,2) = -A(t,3); % fill alpha_32 elements
        Omega(:,:,t) = Alt((t-1)*ny+1:t*ny,:) \ diag(exp(2*H(t,:))) ...
            / Alt((t-1)*ny+1:t*ny,:)';
    end
    
    % b) run Kalman filter and Carter-Kohn recursion algorithm together
    [B, resB] = KalmanCarterKohn(Y, X, reshape(Bpr,nB,1), VBpr, nlag, Omega, Q, 0);
    % for debugging/testing:
    % p00=VBpr; bt00=reshape(Bpr,nB,1); Rtot=Omega; t=1; chck=0;
    
    % Step 1.2: Sample Q from the IW distribution (given beta draws)
    errorB = diff(B);
    scaleQ = errorB'*errorB + covQ;
    Q      = iwishrnd(scaleQ,T-1+TQ); % QB: why -1?
    
    
    %% Stage 2: Draw A(t) and S
    % Step 2.1: Get draws for A and VAR residuals by backward recursions
    XA  = zeros(T,nA,2);       % regressors for A(t)
    Sig = zeros(nA-1,nA-1,T);  % measurement error for A(t)
    for t = 1:T
        XA(t,1,1) = resB(t,1);
        XA(t,2,2) = resB(t,1);
        XA(t,3,2) = resB(t,2);
        Sig(:,:,t) = diag(exp(2*H(t,2:3)));
    end
    
    % Step 2.2: Get draw for alpha_21, alpha_31 and alpha_32 from backward
    % recursion
    [A, resA] = KalmanCarterKohn(resB(:,2:3), XA, Apr, 4*VApr, 1, Sig, S,0);
    % for debugging/testing:
    % Y = resB(:,2:3); X = XA; beta00 = Apr; p00 = 4*VApr; nlag = 1; Rall = Sig; Q = S; t = 1;
    
    % Step 2.3: Draw S, the covariance of A(t), from inverse Wishart
    errorA = diff(A);
    sqerrA = errorA'*errorA;
    
    scaleS2 = sqerrA(1,1) + covS2;
    S2      = iwishrnd(scaleS2, T-1+TS(1));
    
    scaleS3 = sqerrA(2:3,2:3) + covS3;
    S3      = iwishrnd(scaleS3, T-1+TS(2));
    
    S(1,1)     = S2;
    S(2:3,2:3) = S3;
    
    
    %% Stage 3: Draw H(t)
    % Step 3.1: Get y**(t)
    % First, create A(t) as in (A.3), p. 845
    At = repmat(eye(ny),T,1);
    for  ii = 1:T
        At((ii-1)*ny+2,1) = - A(ii,1); % fill alpha_21 elements
        At((ii-1)*ny+3,1) = - A(ii,2); % fill alpha_31 elements
        At((ii-1)*ny+3,2) = - A(ii,3); % fill alpha_32 elements
    end
    % Second, create y*(t) as in (A.3), p845 and y**(t) as in (A.4), p846
    ys = NaN(T,ny);
    for t = 1:T
        ys(t,:) = At((t-1)*ny+1:(t-1)*ny+ny,:)*resB(t,:)';
    end
    yss = log(ys.^2 + cbar);
    
    % Step 3.2: Use Kim et al. (1998) procedure for "mixture of normals
    % approximation of the log chi-squared distribution", Primiceri p846
    s = zeros(T,ny);
    for i=1:ny
        % sample S from a 7-point discrete density (see Kim et al., 1998)
        q = repmat(qj,T,1).*normpdf( repmat(yss(:,i),1,7), 2*repmat(H(:,i),1,7) + repmat(mj,T,1), repmat(sqrt(v2j),T,1)  );
        q = q./repmat(sum(q,2),1,7);
        f = (q + makelag(q',1)' + makelag(q',2)' + makelag(q',3)' ...
            + makelag(q',4)' + makelag(q',5)' + makelag(q',6)');
        z1 = unidrnd(100000,T,1)/100000;
        z2 = sum(z1*ones(1,7) >= f, 2) + 1;
        z3 = z2 > 7;
        s(:,i) = z2 - z3;
    end
    
    % Step 3.3: Get draw for alpha_21, alpha_31 and alpha_32 from backward
    % recursion
    covH = zeros(ny,ny,T); % set up covariance matrix
    for t = 1:T
        covH(:,:,t) = diag(v2j(s(t,:)));
    end
    
    [H, resH] = KalmanCarterKohn(yss-mj(s), ones(T,1)*2, Hpr, VHpr, 1, covH, W, 0);
    % for debugging/testing:
    % Y=yss-mj(s); X=ones(T,1)*2; beta00=Hpr; p00=VHpr; nlag=1; Rall=covH; Q=W; t=1; checkstable=0;
    
    % Step 3.4: Draw S, the covariance of A(t), from inverse Wishart
    errorH = diff(H);
    scaleW = errorH'*errorH + TW*VHpr*kW^2;
    W      = iwishrnd(scaleW,T-1+TW);
    
    
    %% Save output from Gibbs sampler
    B_s(:,:,m)     = B;   % draws of B(t)
    Q_s(:,:,m+1)   = Q;   % draws of Q from IW
    A_s(:,:,m+1)   = A;   % draws of A(t)
    S2_s(:,m+1)    = S2;  % draws of S2 from IW
    S3_s(:,:,m+1)  = S3;  % draws of S3 from IW
    H_s(:,:,m+1)   = H;   % draws of H(t)
    W_s(:,:,m+1)   = W;   % draws of W from IW
    
    
end % end Gibbs sampler
close(wbGS)
toc;


%% Plot figures as in DP
% see Del Negro and Primiceri (2013), p ix

rt.B =     B_s(:, :, thin:thin:end);
rt.A =     A_s(:, :, thin:thin:end);
rt.H = exp(H_s(:, :, thin:thin:end));

rt.M = size(rt.B, 3);
rt.N = round(burn/thin);

rt.Bm = mean(rt.B(:, :, rt.N:end),3);
rt.Am = mean(rt.A(:, :, rt.N:end),3);
rt.Hm = mean(rt.H(:, :, rt.N:end),3);


% Figure 1: Stoch vola - redoing Figure 1 in Del Negro and Primiceri (2015)
figure('name','Posterior mean and (16,84)th percentiles of the standard deviation of residuals')
subplot(3,1,1); plot(time,rt.Hm(:,1)); hold on; title('(a)')
plot(time,prctile(squeeze(rt.H(:,1,rt.N:end))',16),'--r')
plot(time,prctile(squeeze(rt.H(:,1,rt.N:end))',84),'--r')
axis([1962 2002 0 .8])

subplot(3,1,2); plot(time,rt.Hm(:,2)); hold on; title('(b)')
plot(time,prctile(squeeze(rt.H(:,2,rt.N:end))',16),'--r')
plot(time,prctile(squeeze(rt.H(:,2,rt.N:end))',84),'--r')
axis([1962 2002 0 1.2])

subplot(3,1,3); plot(time,rt.Hm(:,3)); hold on; title('(c)')
plot(time,prctile(squeeze(rt.H(:,3,rt.N:end))',16),'--r')
plot(time,prctile(squeeze(rt.H(:,3,rt.N:end))',84),'--r')
axis([1962 2002 0 5])
