function g = a0lhgh(x0,s,nobs,nvar,a0indx);
%function g = a0lhgh(x0,s,nobs,nvar,a0indx);
%      computes the hessian with analytical gradient provided
%  x0 (parameter vector),
%  s (covariance matrix of innovations): note divided by "nobs"
%  nobs (no of obs),
%  nvar (no of variables),
%  a0indx (matrix indicating the free parameters in A0, and each column in A0 corresponds
%                    to an equation)
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

%%a = zeros(nvar*nvar,1);
a = zeros(nvar);
g = zeros(nvar*nvar,1);
badg = 0;
a(a0indx) = x0;
b = -nobs*inv(a') + nobs*s*a;
b = reshape(b,nvar*nvar,1);
g = b(a0indx);