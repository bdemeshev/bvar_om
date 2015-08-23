function [phi,y,fss,xtx,xty,yty,Bols] = xydata(xdgel,lags,mu)
% [phi,y,xtx,xty,yty,Bols,fss] = xydata(xdgel,lags,mu)
%  Construct matrices X, Y, X'X, etc. so that y = X*B + U, where y is T-by-1, X T-by-ncoef,
%     and B is ncoef-by-1.
%
% xdgel: the general matrix of the original data (no manipulation involved)
%                     with sample size including lags
% lags:  the lag length in the AR(p) process
% mu: 2-by-1 vector of hyperparameters for dummy observations (FRA's parameter settings)
%       mu(1): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(2): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
%---------------
% phi:  X: T-by-ncoef where ncoef=nvar*lags + constant and T=fss
% y:  Y: T-by-nvar where T=fss
% fss: effective sample size (including dummies if mu is provided) exclusing all lags
% xtx:  X'X (Data)
% xty:  X'Y (Data)
% yty:  Y'Y (Data)
% Bols:  ncoef-by-1: OLS estimate of B
% Copyright (c) September 1999 Tao Zha
%
% Copyright Kilian and Zha
%function [Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
%             idmat0,idmatpp] = szasbvar(idfile,q_m,lags,xdgel,mu)
%  Now version: [Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
%             idmat0,idmatpp] = szasbvar(idfile,q_m,lags,xdgel,mu)
%
%    Estimating the Bayesian VAR of Sims and Zha with asymmetric prior (as)
%
% idfile:  Identification filename with rows corresponding to equations.
%              But all the output will be transposed to columns
%              corresponding to equations.
% q_m:  quarter or month
% lags: the maximum length of lag
% nhp:  number of haperparameters
% xdgel: the general matrix of the original data (no manipulation involved)
%                    with sample size including lags
% mu: 6-by-1 vector of hyperparameters (the following numbers for Atlanta Fed's forecast)
%       mu(1): overall tightness and also for A0;  (0.57)
%       mu(2): relative tightness for A+;  (0.13)
%       mu(3): relative tightness for the constant term;  (0.1)
%       mu(4): tightness on lag decay;  (1)
%       mu(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
% --------------------
% --------------------
% Gb: cell(nvar,1). Each cell, postmultiplied by A0(:,i), is used to compute A+(:,i)
%               where A+ is k-by-m, where k--ncnoef, m--nvar, if asymmetric prior
% Sbd: cell(nvar,1). Sbd=diag(S(1), ..., S(m)).  Already divided by 'fss' for the
%        posterior of a0 or A0(:) in Waggoner and Zha when one has asymmetric prior.
%        Note,"bd" stands for block diagonal.
% Bh: reduced form B, k(mp+1)-by-m, column--equation. Bh=NaN if asymmetric prior.
%     In the form of y=X*B+U where B=[B1|B2| ...|Bp|C]
%       where Bi: m-by-m (i=1,..,p--lag length), C: 1-by-ml, constant.
% SpH:  divided by T, the final S in the Waggoner-Zha exponential part of p(A0|y)
%               SpH=NaN if asymmetric prior
% fss:  in-sample-size (for forecasting). Including dummy observations
% ndobs: number of dummy observations, ndobs=nvar+1
% phi:  X in the form of y = X*B+U. Row: nSmaple-lags+ndobs. Column: ncoef
% y:    y in the form y = X*B+U. Including the dummy observations too,
%         T-lags+ndobs-by-nvar.
% nvar: number of variables
% ncoef: number of coefficients in *each* equation. RHS coefficients only, nvar*lags+1
% xxhpc: chol(X'X+inv(H_p_tilde)): upper triangular but its transpose
%                                      is lower triagular Choleski
% a0indx: the index number for free parameters in A0.  Column meaning equation
% na0p:   number of free parameters in A0
% idmat0:  identification matrix for A0 with asymmetric prior;  column -- equation.
% idmatpp: asymmetric prior variance for A+ without constant; column -- equation.
%
% Revisions by CAS 9/27/96:  If we transmit idmat, rather than a0indx, we can scale
%   the elements of idmat0 or idmatp so that it carries information about relative
%   prior variances.
% Revsions by TZ, 10/13/96:  Efficiency improvement by streamlining the previous
%   code according to the general setup in Sims and Zha's IER and according
%   to Zha's notes Forecast (0) (p.3) and Forecast (1) (p.9).
%
% Copyright (c) September 1999 by Tao Zha
%
%  NOTE: "nSample" is something I deleted as an input recently, so it may not
%    be compatible with old programs, 10/15/98 by TAZ.
%
%  NOTE2: I put "mu" as an input arg and delete "nhp" as an input arg, so it may
%    not be compatible with old programs, 03/06/99 by TAZ
% Copyright (C) 1997-2012 Christopher A. Sims and Tao Zha 
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


[nSample,nvar]=size(xdgel);   % sample size including lags and number of variables
ncoef=nvar*lags+1;  % number of coefficients in *each* equation, RHS coefficients only,
                    % with only constant
ndobs=nvar+1;    % number of dummy observations
sb = lags+1;   % the beginning of original sample without dummies

if nargin==3     % activate the option of dummy observations
   fss = nSample+ndobs-lags;
                    % Effective sample size including dummies but exclusing all lags

   %
   % **** nvar prior dummy observations with the sum of coefficients
   % ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
   % **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, const]
   % **    Now, T=T+ndobs -- added with "ndobs" dummy observations
   %
   phi = zeros(fss,ncoef);
   const = ones(fss,1);
   const(1:nvar) = zeros(nvar,1);
            % the first nvar periods: no or zero constant (designed for dummies)!
   phi(:,ncoef) = const;
   %
   xdgelint = mean(xdgel(1:lags,:),1); % mean of the first lags initial conditions
   %**** Dummies for phi
   for k=1:nvar
      for m=1:lags
         phi(ndobs,nvar*(m-1)+k) = xdgelint(k);
         phi(k,nvar*(m-1)+k) = xdgelint(k);
         % <<>> multiply hyperparameter later
      end
   end
   %* True data for phi
   for k=1:lags
      phi(ndobs+1:fss,nvar*(k-1)+1:nvar*k) = xdgel(sb-k:nSample-k,:);
      % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
      %                     Thus, # of columns is nvar*lags+1 = ncoef.
   end
   %
   % ** Y with "ndobs" dummies added
   y = zeros(fss,nvar);
   %* Dummies for y
   for k=1:nvar
      y(ndobs,k) = xdgelint(k);
      y(k,k) = xdgelint(k);
      % multiply hyperparameter later
   end
   %* True data for y
   y(ndobs+1:fss,:) = xdgel(sb:nSample,:);
   %
   %*** Dummies once again (finfal version)
   phi(1:nvar,:) = 1*mu(1)*phi(1:nvar,:);    % standard Sims and Zha prior
   y(1:nvar,:) = mu(1)*y(1:nvar,:);      % standard Sims and Zha prior
   %
   phi(nvar+1,:) = mu(2)*phi(nvar+1,:);
   y(nvar+1,:) = mu(2)*y(nvar+1,:);
else      % dummy observations
   fss = nSample-lags;
            % Effective sample size with no dummies but exclusing all lags

   %
   % **** nvar prior dummy observations with the sum of coefficients
   % ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
   % **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, const]
   % **    Now, T=T+ndobs -- added with "ndobs" dummy observations
   %
   phi = zeros(fss,ncoef);
   const = ones(fss,1);
   phi(:,ncoef) = const;
   %* True data for phi
   for k=1:lags
      phi(1:fss,nvar*(k-1)+1:nvar*k) = xdgel(sb-k:nSample-k,:);
      % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
      %                     Thus, # of columns is nvar*lags+1 = ncoef.
   end
   %
   %* True data for y
   y(1:fss,:) = xdgel(sb:nSample,:);
end
%*** data input for the posterior density
[xq,xr]=qr(phi,0);
xtx=xr'*xr;   % X'X
xty=phi'*y;  % X'Y
yty = y'*y;   % Y'Y

Bols = xr\(xr'\xty);  % inv(X'X)*X'Y

