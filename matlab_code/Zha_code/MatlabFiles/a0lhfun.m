function of = a0lhfun(x,s,nobs,nvar,a0indx)
% of = a0lhfun(x,s,nobs,nvar,a0indx) -- negative logLH
%          Note: columns correpond to equations
% general program to setup A0 matrix and compute the likelihood
% requires
%  x (parameter vector),
%  s (covariance matrix of innovations): note already divided by "nobs"
%  nobs (no of obs),
%  nvar (no of variables),
%  a0indx (matrix indicating the free parameters in A0, and each column in A0 corresponds
%                    to an equation)
% written by Eric Leeper
% Copyright (C) 1997-2012 Eric Leeper
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


a0 = zeros(nvar);
a0(a0indx) = x;
% Note: each column in a0 corresponds to an equation!!
%
%%ada = chol(a0'*a0);
%%ada = log(abs(diag(ada)));
%%ada = sum(ada);
% **  TZ, 10/15/96, the above two lines can be improved by the following three lines
[a0l,a0u] = lu(a0);
%ada=diag(abs(a0u));
%ada=sum(log(ada));
ada = sum(log(abs(diag(a0u))));

%
%tra = trace(a0'*s*a0);
tra = reshape(s,nvar*nvar,1)'*reshape(a0*a0',nvar*nvar,1);
%if ada == 0
%   of = 1.0;
%   else
%   of = -nobs*ada + nobs*.5*tra;
%end
% commented out by T.Z. for there appears no sense of using "if-end"


of = -nobs*ada + nobs*.5*tra;