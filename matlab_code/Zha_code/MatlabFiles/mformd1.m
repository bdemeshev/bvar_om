function [phi,y,ncoe] = mformd1(z,lags)
%
% [phi,y,ncoe] = mformd1(z,lags) where size(z,1)=T+lags.
%    This file is a subset of sye.m, but only designed to arrange matrix form data as
%             y(T*nvar) = XB + u, X--phi: T*k, B: k*nvar; where T=sp-lags, k=ncoe.
%    Note T>=0.  If T=0, y is empty, but phi (X) is still well defined.
%
% z: (T+lags)-by-(nvar+1) raw data matrix (nvar of variables + constant).
%             e.g., C = ones(nSample,1); z=[xdget C];
% lags: number of lags
% phi:  (T-lags)-by-k X; where k=ncoe=[nvar for 1st lag, ..., nvar for last lag, const]
% y:    Y: (T-lags)-by-nvar
% ncoe: number of coeffcients per equation: k=nvar*lags + 1
%
% See also sye.m
%
% October 1998 Tao Zha.  Revised, 03/13/99
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
if ess<0
   warning('The length of the data must be >= lags')
   disp('Check yrEnd and qmEnd in parac.m to make sure')
   disp('Press ctrl-c to abort')
   disp(' ')
   pause
elseif ess==0
   x = z(:,1:nvar);
   C = z(:,nvar+1);
   phi = zeros(1,ncoe);
   phi(1,ncoe) = C(1);
   for k=1:lags, phi(1,nvar*(k-1)+1:nvar*k) = x(sb-k:sl-k+1,:); end
   % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
   %                     Thus, # of columns is nvar*lags+1 = ncoef.
   y = x(sb:sl,:);
else
   x = z(:,1:nvar);
   C = z(:,nvar+1);
   phi = zeros(ess,ncoe);
   phi(:,ncoe) = C(1:ess);
   for k=1:lags, phi(:,nvar*(k-1)+1:nvar*k) = x(sb-k:sl-k,:); end
   % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
   %                     Thus, # of columns is nvar*lags+1 = ncoef.
   y = x(sb:sl,:);
end