function yhat = fn_forecast_consts(Bh,nn)
%yhat = fn_impulse2shock(Bh,shocks,nn)  WARNING: Not sure what this function is trying to do.  T. Zha, 2013, Nov 13.
%  yhat: nsteps X nvar matrix of responses
%-----------------
%  Bh is the estimated reduced form coefficient in the form
%       Y(T*nvar) = XB + U, X: T*k (may include all exogenous terms), B: k*nvar.
%       The matrix form and dimension are the same as "Bh" from the function "sye.m";
%       Rows: 1st lag (with nvar variables) to lags (with nvar variables) + const = k.
%       Note: columns correspond to equations.
%  permshocks: Permanent shocks stored as a 1Xnvar vector of reduced-form shocks.
%  nn is the numbers of inputs [nvar, lags, # of forecast steps].
%
%  Written by Tao Zha.
% Copyright (c) 2013 by Tao Zha
% Copyright (C) 2013 Tao Zha
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

nvar = nn(1);
lags = nn(2);
nfqm = nn(3);
ntc = nvar*lags + 1;  % total number of coefficeints (with constant only).
tcwx = nvar*lags;  % total number of coefficients without constant terms.

%*** reconstruct x(t) for y(t+h) = x(t+h-1)*B
%***       where phi = x(t+h-1) with last column being constant
phi = ones(1,ntc);
for (ki=1:lags)
   phi((1+(ki-1)*nvar):(ki*nvar)) = Bh(end,:);
end   
yhat = zeros(nfqm,nvar);
for k=1:nfqm
   yhat(k,:) = phi*Bh;
   %*** lagged endogenous variables
   phi(1,nvar+1:tcwx) = phi(1,1:tcwx-nvar);
   phi(1,1:nvar) = yhat(k,:);
end
