function [g,badg] = fn_a0freegrad(b,Ui,nvar,n0,fss,H0inv)
% [g,badg] = a0freegrad(b,Ui,nvar,n0,fss,H0inv)
%    Analytical gradient for a0freefun.m in use of csminwel.m.  See Dhrymes's book.
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
%---------------
% g: sum(n0)-by-1 analytical gradient for a0freefun.m
% badg: 0, the value that is used in csminwel.m
%
% Tao Zha, February 2000.  Revised, August 2000
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


b=b(:);  n0 = n0(:);

A0 = zeros(nvar);
n0cum = [0;cumsum(n0)];
g = zeros(n0cum(end),1);
badg = 0;

%*** The derivative of the exponential term w.r.t. each free paramater
for kj = 1:nvar
   bj = b(n0cum(kj)+1:n0cum(kj+1));
   g(n0cum(kj)+1:n0cum(kj+1)) = H0inv{kj}*bj;
   A0(:,kj) = Ui{kj}*bj;
end
B=inv(A0');

%*** Add the derivative of -Tlog|A0| w.r.t. each free paramater
for ki = 1:sum(n0)
   n = max(find( (ki-n0cum)>0 ));  % note, 1<=n<=nvar
   g(ki) = g(ki) - fss*B(:,n)'*Ui{n}(:,ki-n0cum(n));
end
