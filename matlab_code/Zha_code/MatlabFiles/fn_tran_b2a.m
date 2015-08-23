function A0 = fn_tran_b2a(b,Ui,nvar,n0)
% A0 = fn_tran_b2a(b,Ui,nvar,n0)
%    Transform free parameters b's to A0.  Note: columns correspond to equations
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
%----------------
% A0: nvar-by-nvar, contempareous matrix (columns correspond to equations)
%
% Tao Zha, February 2000.  Revised, August 2000.
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

b=b(:);  n0=n0(:);
A0 = zeros(nvar);
n0cum = [0; cumsum(n0)];
for kj = 1:nvar
   A0(:,kj) = Ui{kj}*b(n0cum(kj)+1:n0cum(kj+1));
end
