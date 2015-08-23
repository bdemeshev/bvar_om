function [xtx,xty,yty,fss,phi,y,ncoef,xr,Bh,e] = dataxy(nvar,lags,z,mu,indxDummy,nexo)
% [xtx,xty,yty,fss,phi,y,ncoef,xr,Bh,e] = dataxy(nvar,lags,z,mu,indxDummy,nexo)
%
%    Export arranged data matrices for future estimation for DM models.
%    See Wagonner and Zha's Gibbs sampling paper.
%
% nvar:  number of endogenous variables.
% lags: the maximum length of lag
% nhp:  number of haperparameters
% z: T*(nvar+(nexo-1)) matrix of raw or original data (no manipulation involved)
%                 with sample size including lags.
% mu: 6-by-1 vector of hyperparameters (the following numbers for Atlanta Fed's forecast)
%       mu(1): overall tightness and also for A0;  (0.57)
%       mu(2): relative tightness for A+;  (0.13)
%       mu(3): relative tightness for the constant term;  (0.1)
%       mu(4): tightness on lag decay;  (1)
%       mu(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
%         Here, only mu(5) and mu(6) are used for dummy observations
% indxDummy = 1;  % 1: add dummy observations to the data; 0: no dummy added.
% nexo:  number of exogenous variables.  Besides constant terms, we have nexo-1 exogenous variables.
% -------------------
% xtx:  X'X: k-by-k where k=ncoef
% xty:  X'Y: k-by-nvar
% yty:  Y'Y: nvar-by-nvar
% fss:  T: sample size excluding lags.  With dummyies, fss=nSample-lags+ndobs.
% phi:  X; T-by-k; column: [nvar for 1st lag, ..., nvar for last lag, const]
% y:    Y: T-by-nvar where T=fss
% ncoef: number of coefficients in *each* equation. RHS coefficients only, nvar*lags+ncoef
% xr:  the economy size (ncoef-by-ncoef) in qr(phi) so that xr=chol(X'*X) or xr'*xr=X'*X
% Bh: ncoef-by-nvar estimated reduced-form parameter; column: nvar;
%           row: ncoef=[nvar for 1st lag, ..., nvar for last lag, exogenous variables, const]
% e:  estimated residual e = y -x*Bh,  T-by-nvar
%
% Tao Zha, February 2000
% Copyright (C) 1997-2012 Tao Zha
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


if nargin == 5
   nexo=1;    % default for constant term
elseif nexo<1
   error('We need at least one exogenous term so nexo must >= 1')
end

%*** original sample dimension without dummy prior
nSample = size(z,1); % the sample size (including lags, of course)
sb = lags+1;   % original beginning without dummies
ncoef = nvar*lags+nexo;  % number of coefficients in *each* equation, RHS coefficients only.

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
