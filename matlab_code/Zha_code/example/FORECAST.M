function yhat = forecast(Bh,phi,nn)
% yhat = forecast(Bh,phi,nn)
%
% Forecast: unconditional forecating: yhat = forecast(Bh,phi,nn)
%       y_hat(t+h) = c + x_hat(t+h-1)*Bh, X: 1*k; Bh: k*nvar; y_hat: 1*nvar
%  where Bh: the (posterior) estimate of B;
%        phi: the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                (last period plus lags before the beginning of forecast);
%        nn: [nvar,lags,forep], forep: forecast periods;
%        yhat: forep*nvar.

% ** setup
nvar = nn(1);
lags = nn(2);
forep = nn(3);
tcwc = nvar*lags;     % total coefficients without constant

% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
% **       where phi = x(t+h-1) with last column being constant
yhat = zeros(forep,nvar);
for k=1:forep
   yhat(k,:) = phi*Bh;
   phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
   phi(1,1:nvar) = yhat(k,:);
end
