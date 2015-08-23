function of = a0freefun(b,Ui,nvar,n0,fss,H0inv)
% of = a0freefun(b,Ui,nvar,n0,fss,H0inv)
%
%    Negative logPosterior function for squeesed A0 free parameters, which are b's in the WZ notation
%        Note: columns correspond to equations
%
% b: sum(n0)-by-1 vector of A0 free parameters
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
% nvar:  number of endogeous variables
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
% fss:  nSample-lags (plus ndobs if dummies are included)
% H0inv: cell(nvar,1).  In each cell, posterior inverse of covariance inv(H0) for the ith equation,
%           resembling old SpH in the exponent term in posterior of A0, but not divided by T yet.
%----------------
% of:  objective function (negative logPosterior)
%
% Tao Zha, February 2000
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

b=b(:);

A0 = zeros(nvar);
n0cum = cumsum(n0);
tra = 0.0;
for kj = 1:nvar
   if kj==1
      bj = b(1:n0(kj));
      A0(:,kj) = Ui{kj}*bj;
      tra = tra + 0.5*bj'*H0inv{kj}*bj;   % negative exponential term
   else
      bj = b(n0cum(kj-1)+1:n0cum(kj));
      A0(:,kj) = Ui{kj}*bj;
      tra = tra + 0.5*bj'*H0inv{kj}*bj;   % negative exponential term
   end
end

[A0l,A0u] = lu(A0);

ada = -fss*sum(log(abs(diag(A0u))));    % negative log determinant of A0 raised to power T

of = ada + tra;
