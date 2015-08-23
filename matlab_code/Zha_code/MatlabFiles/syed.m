function [Bh,e,xtx,phi,y] = syed(z,nn)
% syed: estimate a system of equations: [Bh,e,xtx,phi,y] = syed(z,nn)
%               Y((T-lags)*nvar) = XB + u, X: (T-lags)*k, B: k*nvar.
%    where z is the T*(nvar+ndt) raw data matrix (nvar of variables +
%                               ndt -- number of deterministic terms);
%          nn is 5 inputs [auindx, ndt, nvar,lags,sample period (total)];
%               auindx = 0 (no autoregressive) and 1 (autoregressive);
%               total -- including lags, etc.
%          Bh: the estimated B;  column: nvar;  row: [nvar for 1st lag, ...,
%                         nvar for last lag, deterministic terms (ndt)]
%          e:  estimated residual e = y -xBh,  (T-lags)*nvar
%          xtx:  X'X
%          phi:  X;  column: [nvar for 1st lag, ...,
%                        nvar for last lag, deterministic terms (ndt)]
%          y:    Y
%
%          See also "sye.m".
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
auindx = nn(1);
ndt = nn(2);   % # of deterministic terms including constant
nvar = nn(3);  % # of endogenous variables
lags = nn(4);
sp = nn(5);    % sample period

ess = sp-lags;  % effective sample size
sb = lags+1;   % sample beginning
sl = sp;       % sample last period

% ** construct X for Y = X*B + U where phi = X **
x = z(:,1:nvar);
C = z(:,nvar+1:nvar+ndt);  % C = [] when ndt=0
%
if auindx == 0
   ncoe = ndt;     % with deterministic terms
   phi = zeros(sp,ncoe);     % preallocating
   y = x;
else
   y = x(sb:sl,:);
   ncoe = nvar*lags + ndt;     % with deterministic terms
   phi = zeros(ess,ncoe);     % preallocating
   for k=1:lags, phi(:,nvar*(k-1)+1:nvar*k) = x(sb-k:sl-k,:); end
end
%
if length(C) == 0
   phi(:,ncoe-ndt+1:ncoe) = C;
else
   phi(:,ncoe-ndt+1:ncoe) = C(1:sp,:);   % perhaps, it should have been be C(sb:sp,:).  2/24/00
end
%
% row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag,
%                          deterministic terms (ndt)]
%   Thus, # of columns is nvar*lags+ndt = ncoe.
% ** estimate: B, XTX, residuals **
[u d v]=svd(phi,0); %trial
%xtx = phi'*phi;      % X'X, k*k (ncoe*ncoe)
vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
dinv = 1./diag(d);    % inv(diag(d))
vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
xtx=vd*vd';
xtxinv = vdinv*vdinv';
%xty = phi'*y;        % X'Y
uy = u'*y; %trial
xty = vd*uy; %trial
%Bh = xtx\xty;        %inv(X'X)*(X'Y), k*m (ncoe*nvar).
Bh = xtxinv*xty;
%e = y - phi*Bh;      % from Y = XB + U, e: (T-lags)*nvar
e = y - u*uy;
