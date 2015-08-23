function yhat = fn_impulseconstants(Bh,const_matrix,nn)
%yhat = fn_impulseconstants(Bh,const_matrix,nn)
%  yhat: nsteps X nvar matrix of responses
%-----------------
%  Bh is the estimated reduced form coefficient in the form
%       Y(T*nvar) = XB + U, X: T*k (may include all exogenous terms), B: k*nvar.
%       The matrix form and dimension are the same as "Bh" from the function "sye.m";
%       Column: 1st lag (with nvar variables) to lags (with nvar variables) + const = k.
%       Note: columns correspond to equations.
%  permshocks: Permanent shocks stored as a 1Xnvar vector of reduced-form shocks.
%  nn is the numbers of inputs [nvar,lags,# of steps of impulse responses].
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
tcwx = nvar*lags;  % total coefficeint without exogenous variables

%*** reconstruct x(t) for y(t+h) = x(t+h-1)*B
%***       where phi = x(t+h-1) with last column being constant
phi = zeros(1,tcwx);
Bhwx = Bh(1:tcwx,:);
yhat = zeros(nfqm,nvar);
for k=1:nfqm
   yhat(k,:) = phi*Bhwx + const_matrix(k,:);
   phi(1,nvar+1:tcwx) = phi(1,1:tcwx-nvar);
   phi(1,1:nvar) = yhat(k,:);
end
