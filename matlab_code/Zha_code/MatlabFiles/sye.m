function [Bh,e,xtx,xty,phi,y,ncoe,xr] = sye(z,lags)
% Now [Bh,e,xtx,xty,phi,y,ncoe,xr] = sye(z,lags)
% Old: [Bh,e,xtx,xty,phi,y,ncoe,Sigu,xtxinv] = sye(z,lags)
%
%    Estimate a system of equations in the form of  y(T*nvar) = XB + u,
%          X--phi: T*k, B: k*nvar; where T=sp-lags, k=ncoe,
%
% z: (T+lags)-by-(nvar+1) raw data matrix (nvar of variables + constant).
% lags: number of lags
%--------------------
% Bh: k-by-nvar estimated reduced-form parameter; column: nvar;
%           row: k=ncoe=[nvar for 1st lag, ..., nvar for last lag, const]
% e:  estimated residual e = y -xBh,  T-by-nvar
% xtx:  X'X: k-by-k
% xty:  X'Y: k-by-nvar
% phi:  X; T-by-k; column: [nvar for 1st lag, ..., nvar for last lag, const]
% y:    Y: T-by-nvar
% ncoe: number of coeffcients per equation: nvar*lags + 1
% xr:  the economy size (k-by-k) in qr(phi) so that xr=chol(X'*X)
% Sigu: e'*e: nvar-by-nvar. Note, not divided (undivided) by "nobs"
% xtxinv: inv(X'X): k-by-k
%
% See also syed.m (allowing for more predetermined terms) which has not
%        been yet updated as "sye.m".
%
%  Note, "lags" is something I changed recently, so it may not be compatible
%       with old programs, 10/15/98 by TAZ.
%
% Revised, 5/2/99.  Replaced outputs Sigu and xtxinv with xr so that previous
%                programs may be incompatible.

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
[sp,nvar] = size(z);   % sp: sample period T include lags
nvar = nvar-1;     % -1: takes out the counting of constant

ess = sp-lags;  % effective sample size
sb = lags+1;   % sample beginning
sl = sp;       % sample last period
ncoe = nvar*lags + 1;     % with constant

% ** construct X for Y = X*B + U where phi = X **
x = z(:,1:nvar);
C = z(:,nvar+1);
phi = zeros(ess,ncoe);
phi(:,ncoe) = C(1:ess);
for k=1:lags, phi(:,nvar*(k-1)+1:nvar*k) = x(sb-k:sl-k,:); end
% row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
%                     Thus, # of columns is nvar*lags+1 = ncoef.
% ** estimate: B, XTX, residuals **
y = x(sb:sl,:);
%
%**** The following, though stable, is too slow *****
%  [u d v]=svd(phi,0); %trial
%  %xtx = phi'*phi;      % X'X, k*k (ncoe*ncoe)
%  vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
%  dinv = 1./diag(d);    % inv(diag(d))
%  vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
%  xtx=vd*vd';
%  xtxinv = vdinv*vdinv';
%  %xty = phi'*y;        % X'Y
%  uy = u'*y; %trial
%  xty = vd*uy; %trial
%  %Bh = xtx\xty;        %inv(X'X)*(X'Y), k*m (ncoe*nvar).
%  Bh = xtxinv*xty;
%  %e = y - phi*Bh;      % from Y = XB + U, e: (T-lags)*nvar
%  e = y - u*uy;
%**** The following, though stable, is too slow *****

%===== (Fast but perhaps less accurate) alternative to the above =========
[xq,xr]=qr(phi,0);
xtx=xr'*xr;
xty=phi'*y;
Bh = xr\(xr'\xty);
e=y-phi*Bh;
%===== (Fast but perhaps less accurate) alternative to the above =========


%* not numerically stable way of computing e'*e
%Sigu = y'*y-xty'*Bh;
%Sigu = y'*(eye(ess)-phi*xtxinv*phi')*y;    % probablly better way, commented out
                                            % by TZ, 2/28/98.  See following
%Sigu = y'*(eye(ess)-u*u')*y;    % I think this is the best, TZ, 2/28/98
                % Note, u*u'=x*inv(x'x)*x.
