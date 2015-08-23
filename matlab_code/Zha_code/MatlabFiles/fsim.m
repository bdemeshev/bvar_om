function yhat = fsim(Bh,A0h,phi,nn)
% yhat = fsim(Bh,A0h,phi,nn)
%    Unconditional forecasts with simulated shocks (which do not depend on reduced-form or structural shocks).
%       y_hat(t+h) = c + x_hat(t+h-1)*Bh, X: 1*k; Bh: k*nvar; y_hat: 1*nvar
%  where Bh: the (posterior) estimate of B;
%        A0h: has columns correpond to equations
%        phi: the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                (last period plus lags before the beginning of forecast);
%        nn: [nvar,lags,forep], forep: forecast periods;
%        yhat: forep*nvar.
%
% 3/22/99.  Revised so that A0h has columns corr. to equations.  Previous
%     may use A0h' (papers before Sims and Zha IER paper).  So please double
%     check.
%
% See forecasterr.m when A0h_in (instead of A0h) is used.
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

Ures = randn(forep,nvar)/A0h;    % Unconditional forecast
         % Now, forep-by-nvar -- ready for forecasts
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks

yhat = zeros(forep,nvar);
for k=1:forep
   yhat(k,:) = phi*Bh + Ures(k,:);
   phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
   phi(1,1:nvar) = yhat(k,:);
end
