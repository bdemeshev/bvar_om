function yhat = fn_forecastfixe(Bh,A0h,phi,nn,Estr,nexo,Xfexo)
% yhat = fn_forecastfixe(Bh,A0h,phi,nn,Estr,nexo,Xfexo)
%   Forecasts conditional on fixed structural shocks
%       y_hat(t+h) = c + x_hat(t+h-1)*Bh + Estr(t+h,:), X: 1*k; Bh: k*nvar; y_hat: 1*nvar
%
% Bh: k-by-nvar, the (posterior) estimate of B;
% A0h: nvar-by-nvar, columns correponding to equations
% phi: the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%          (last period plus lags before the beginning of forecast).
% nn: [nvar,lags,nfqm], nfqm: forecast periods (months or quarters).
% Estr:  nfqm-by-nvar, each column corresponds to strcutural shocks from a particular
%            source such as MS;
% nexo:  number of exogenous variables.  The constant term is the default setting.
%              Besides this term, we have nexo-1 exogenous variables.
% Xfexo:  nfqm-by-nexo-1 vector of exoenous variables in the forecast horizon where
%         nfqm:  number of forecast periods.
%-----------------
% yhat: nfqm*nvar matrix of forecasts.
%
% See also fn_forecast.m and fn_forecastsim.m.
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

% ** setup
nvar = nn(1);
lags = nn(2);
nfqm = nn(3);
tcwx = nvar*lags;  % total coefficeint without exogenous variables

if nexo>1
   if (nfqm > size(Xfexo,1))
      disp(' ')
      warning('Make sure the forecast horizon in the exogenous variable matrix Xfexo > forecast periods')
      disp('Press ctrl-c to abort')
      pause
   elseif ((nexo-1) ~= size(Xfexo,2))
      disp(' ')
      warning('Make sure that nexo matchs the exogenous variable matrix Xfexo')
      disp('Press ctrl-c to abort')
      pause
   end
end

% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
% **       where phi = x(t+h-1) with last column being constant
Ures = Estr/A0h;
yhat = zeros(nfqm,nvar);
for k=1:nfqm
   yhat(k,:) = phi*Bh + Ures(k,:);
   %*** lagged endogenous variables
   phi(1,nvar+1:tcwx) = phi(1,1:tcwx-nvar);
   phi(1,1:nvar) = yhat(k,:);
   %*** exogenous variables excluding constant terms
   if (nexo>1)
      phi(1,tcwx+1:tcwx+nexo-1) = Xfexo(k,:);
   end
end
