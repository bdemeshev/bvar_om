
function [bt, res] = KalmanCarterKohn(Y,X,bt00,p00,nlag,Rtot,Q,chck)

% Run first the Kalman filter and then the Carter-Kohn backward recursions
% to obtain smoothed estimates of states from the following state-space
% model:
%
%    Y(t) = X(t)*bt(t)   + e(t),  e(t)~N(0,R(t))    (observation equation)
%   bt(t) =      bt(t-1) + v(t),  v(t)~N(0,Q)        (transition equation)
%
% Usage: [bt, res] = KalmanCarterKohn(Y,X,bt00,p00,nlag,Rtot,Q,chck)
%
% Inputs: Y    - data (T x ny)
%         X    - regressors (T x ns), where ns=(1+nlag*ny)*ny
%         bt00 - initial values for bt (1 x ns)
%         p00  - initial value for covariance matrix of bt (ns x ns)
%         Rtot - covariance matrix of measurement error, R(t) (ns x T)
%         Q    - covariance matrix of transition equation error
%         chck - (opt.) set to 1 to allow only stable solutions
%                         (default = 0)
%
% Output: bt - estimated states (T x ns)
%         res  - smoothed residuals (T x ny)
%
% Based on code by Haroon Mumtaz (example3.m), see
% http://www.bankofengland.co.uk/education/Documents/ccbs/...
% technical_handbooks/Coding/code.zip
% By B. Kolb, Oct. 2015


%% Step 1: Set up matrices for the Kalman Filter
[T, ny] = size(Y);
ns      = size(bt00,1);
bt_tt   = NaN(T,ns);
ptt     = NaN(T,ns,ns);

bt11    = bt00;
p11     = p00;

if nargin == 7
    chck = 0;
end


%% Step 2: Run Kalman Filter
for t=1:T
    if ndims(Rtot) == 3   % i.e. for B(t) and A(2:3,2:3,t)
        R = Rtot(:,:,t);
    elseif ismatrix(Rtot) % i.e. for A(1,1,t)
        R = Rtot(:,t);
    end
    
    if ismatrix(X)       % i.e. for B(t)
        x = kron(eye(ny),X(t,:));
        y = Y(t,:);
    elseif ndims(X) == 3 % i.e. for A(t)
        x = squeeze(X(t,:,:))';
        y = Y(t,:);
    end
    
    bt10 = bt11;             % prediction
    p10 = p11 + Q;           % prediction
    eta = y - (x*bt10)';     % prediction error
    K = p10*x'/(x*p10*x'+R); % Kalman gain
    bt11 = (bt10+K*eta');    % updating
    p11 = p10 - K*x*p10;     % updating
    ptt(t,:,:) = p11;        % saving
    bt_tt(t,:) = bt11;       % saving
end


%% Step 3: Backward recursion (Carter-Kohn algorithm)
% (get mean and variance of states, and use these to take draws)
stable = -1;
while stable < 0 % stability not yet satisfied
    
    bt = zeros(T,ns); % saves draw of the state variable
    res = zeros(T,ny);
    instable = zeros(T,1);
    
    % Start at period T
    p00     = squeeze(ptt(T,:,:));
    bt(T,:) = mvnrnd(bt_tt(T,:)',p00/2+p00'/2,1);  % ensures symmetric p00
    
    if ismatrix(X)       % residuals for B(t) regression
        res(T,:)  = Y(T,:)-X(T,:)*reshape(bt(T,:),size(X,2),ny);
    elseif ndims(X) == 3 % residuals for A(t) regression
        res(T,:)  = Y(T,:)-(squeeze(X(T,:,:))'*bt(T,:)')';
    end
    
    if chck == 1
        instable(T) = stabilityB(bt(T,:)',ny,nlag);  % instable=0: stable
    else
        instable = 0;
    end
    
    % Periods T-1 to 1
    for t = T-1:-1:1
        pt = squeeze(ptt(t,:,:));
        % update bt for information in bt[t+1], eq 8.16 p193 in Kim/Nelson
        bm = bt_tt(t:t,:) + (pt/(pt+Q)*(bt(t+1:t+1,:)-bt_tt(t,:))')';
        pm = pt - pt/(pt+Q)*pt;              % update covariance of bt
        bt(t:t,:) = mvnrnd(bm,pm/2+pm'/2,1); % ensures pm is symmetric
        
        if ismatrix(X)       % residuals for B(t) regression
            res(t,:)=Y(t,:)-X(t,:)*reshape(bt(t:t,:),size(X,2),ny);
        elseif ndims(X) == 3 % residuals for A(t) regression
            res(T,:)  = Y(t,:)-(squeeze(X(t,:,:))'*bt(t,:)')';
        end
        
        if chck == 1
            instable(t) = stabilityB(bt(t,:)',ny,nlag); % stable if 0
        end
    end
    
    if sum(instable)==0 % not instable
        stable=1;
    end
end


end