function yhat = fshock(Bh,phi,shock,nn)
% fshock: unconditionally forecasts supplied with reduced-form shocks
%           yhat = fshock(Bh,phi,shock,nn)
%
%  y_hat(t+h) = c + x_hat(t+h-1)*Bh, X: 1*k; Bh: k*nvar; y_hat: 1*nvar
%  where Bh: the (posterior) estimate of B;
%        phi: the 1-by-(nvar*lags+1) initial data matrix where k=nvar*lags+1
%                (last period plus lags before the beginning of forecast);
%        shock: forep-by-nvar reduced-form residuals, where forep: forecast periods;
%        nn: [nvar,lags,forep];
%        yhat: forep-by-nvar.
%
% For structural shocks, see forefixe.m
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

% ** setup
nvar = nn(1);
lags = nn(2);
forep = nn(3);
tcwc = nvar*lags;     % total coefficients without constant

% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
% **       where phi = x(t+h-1) with last column being constant
yhat = zeros(forep,nvar);
for k=1:forep
   yhat(k,:) = phi*Bh + shock(k,:);
   phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
   phi(1,1:nvar) = yhat(k,:);
end
