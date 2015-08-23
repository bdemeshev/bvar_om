function yhat = get_pathsVAR(N, HORZ, T, L, Y, beta2, sigma)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


% Draw N(0,1) innovations for variance and mean equation:
csigma=chol(sigma);
uu = randn(HORZ+L,N);
% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.

%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y(T-L+1:T,:);

for fi=L+1:HORZ+L
    
    xhat=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
    end
    xhat=[xhat 1];
    
   
    yhat(fi,:) = xhat*reshape(beta2,N*L+1,N) + uu(fi,:)*csigma;
    
end
