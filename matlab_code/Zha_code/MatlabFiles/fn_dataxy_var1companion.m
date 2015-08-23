function [Yt, Ytm1] = fn_dataxy_var1companion(nvar,lags,z)
% function [Yt, Ytm1] = fn_dataxy_companion(nvar,lags,z)   
%
% nvar:  number of endogenous variables.
% lags: the maximum length of lag
% z: T*nvar matrix of raw or original data (no manipulation involved)
%       with sample size including lags and with exogenous variables other than a constant.
%       Order of columns: (1) nvar endogenous variables; (2) (nexo-1) exogenous variables;
%                         (3) constants will be automatically put in the last column.
% -------------------
% phi:  X; fss-by-k; column: [nvar for 1st lag, ..., nvar for last lag, const term]
% y:    Y: fss-by-nvar where T=fss
% ncoef: number of coefficients in *each* equation. RHS coefficients only, nvar*lags+nexo
%
% Tao Zha, March 2014.
% Copyright (C) 2014- Tao Zha
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


%*** original sample dimension without dummy prior
nSample = size(z,1); % the sample size (including lags, of course)
sb = lags+1;   % original beginning without dummies
ncoef = nvar*lags+1;  % number of coefficients in *each* equation, RHS coefficients only.

if indxDummy   % prior dummy prior
   %*** expanded sample dimension by dummy prior
   ndobs=nvar+1;         % number of dummy observations
   fss = nSample+ndobs-lags;

   %
   % **** nvar prior dummy observations with the sum of coefficients
   % ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
   % **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, exo var, const]
   % **    Now, T=T+ndobs -- added with "ndobs" dummy observations
   %
   phi = zeros(fss,ncoef);
   %* constant term
   const = ones(fss,1);
   const(1:nvar) = zeros(nvar,1);
   phi(:,ncoef) = const;      % the first nvar periods: no or zero constant!
   %* other exogenous (than) constant term
   phi(ndobs+1:end,ncoef-nexo+1:ncoef-1) = z(lags+1:end,nvar+1:nvar+nexo-1);
   exox = zeros(ndobs,nexo);
   phi(1:ndobs,ncoef-nexo+1:ncoef-1) = exox(:,1:nexo-1);
            % this = [] when nexo=1 (no other exogenous than constant)

   xdgel = z(:,1:nvar);  % endogenous variable matrix
   xdgelint = mean(xdgel(1:lags,:),1); % mean of the first lags initial conditions
   %* Dummies
   for k=1:nvar
      for m=1:lags
         phi(ndobs,nvar*(m-1)+k) = xdgelint(k);
         phi(k,nvar*(m-1)+k) = xdgelint(k);
         % <<>> multiply hyperparameter later
      end
   end
   %* True data
   for k=1:lags
      phi(ndobs+1:fss,nvar*(k-1)+1:nvar*k) = xdgel(sb-k:nSample-k,:);
      % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, exo var, const]
      %                     Thus, # of columns is nvar*lags+nexo = ncoef.
   end
   %
   % ** Y with "ndobs" dummies added
   y = zeros(fss,nvar);
   %* Dummies
   for k=1:nvar
      y(ndobs,k) = xdgelint(k);
      y(k,k) = xdgelint(k);
      % multiply hyperparameter later
   end
   %* True data
   y(ndobs+1:fss,:) = xdgel(sb:nSample,:);

   phi(1:nvar,:) = 1*mu(5)*phi(1:nvar,:);    % standard Sims and Zha prior
   y(1:nvar,:) = mu(5)*y(1:nvar,:);      % standard Sims and Zha prior
   phi(nvar+1,:) = mu(6)*phi(nvar+1,:);
   y(nvar+1,:) = mu(6)*y(nvar+1,:);

   [xq,xr]=qr(phi,0);
   xtx=xr'*xr;
   xty=phi'*y;
   [yq,yr]=qr(y,0);
   yty=yr'*yr;
   Bh = xr\(xr'\xty);   % xtx\xty where inv(X'X)*(X'Y)
   e=y-phi*Bh;
else
   fss = nSample-lags;
   %
   % ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
   % **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, exo var, const]
   %
   phi = zeros(fss,ncoef);
   %* constant term
   const = ones(fss,1);
   phi(:,ncoef) = const;      % the first nvar periods: no or zero constant!
   %* other exogenous (than) constant term
   phi(:,ncoef-nexo+1:ncoef-1) = z(lags+1:end,nvar+1:nvar+nexo-1);
            % this = [] when nexo=1 (no other exogenous than constant)

   xdgel = z(:,1:nvar);  % endogenous variable matrix
   %* True data
   for k=1:lags
      phi(:,nvar*(k-1)+1:nvar*k) = xdgel(sb-k:nSample-k,:);
      % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, exo var, const]
      %                     Thus, # of columns is nvar*lags+nexo = ncoef.
   end
   %
   y = xdgel(sb:nSample,:);

   [xq,xr]=qr(phi,0);
   xtx=xr'*xr;
   xty=phi'*y;
   [yq,yr]=qr(y,0);
   yty=yr'*yr;
   Bh = xr\(xr'\xty);   % xtx\xty where inv(X'X)*(X'Y)
   e=y-phi*Bh;
end
