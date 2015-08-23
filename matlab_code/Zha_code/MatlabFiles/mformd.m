function [phi,YtY] = mformd(z,nn)
% mformd: arrange matrix form data: [Y,YtY] = mformd(z,nn) as structural form
%                   YA = E, Y: T*(nvar*lags+nvar+1)
%
%    where z is the (T+lags)-by-(nvar+1) raw data matrix (nvar of variables + constant);
%          nn is the numbers of inputs [nvar,lags,sample period (total)];
%          phi: Y as in the structural form YA = E, Y: T*(nvar*lags+nvar+1)
%          YtY:  Y'Y: (nvar*lags+nvar+1)*(nvar*lags+nvar+1)
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

% ** setup of orders and lengths **
nvar = nn(1);
lags = nn(2);
sp = nn(3);    % sample period

ess = sp-lags;  % effective sample size
sb = lags+1;   % sample beginning
sl = sp;       % sample last period
ncoe = nvar*lags + nvar + 1;     % with constant and contemporaneous data

% ** construct Y as in YA = E where phi = Y **
x = z(:,1:nvar);
C = z(:,nvar+1);
phi = zeros(ess,ncoe);
phi(:,1:nvar) = x(sb:sl,:);
phi(:,ncoe) = C(1:ess);
for k=1:lags, phi(:,nvar*k+1:nvar*(k+1)) = x(sb-k:sl-k,:); end
% column: T; row: [nvar for 0th lag, nvar for 1st lag,
%                                ..., nvar for last lag, const]
%                     Thus, # of columns is nvar*lags+nvar+1 = ncoef.
% ** YtY, residuals **
[u d v]=svd(phi,0); %trial
vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
YtY=vd*vd';        % YtY = phi'*phi;      % X'X, k*k (ncoe*ncoe)