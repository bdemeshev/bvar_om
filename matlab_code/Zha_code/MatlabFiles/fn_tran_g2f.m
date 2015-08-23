function F = fn_tran_g2f(g,Vi,nvar,ncoef,np)
%    Transform free parameters g's to F (A+).  Note: columns correspond to equations
%    See Waggoner and Zha's ``A Gibbs sampler for structural VARs''
%
% g: sum(np)-by-1 stacked vector of all free lagged parameters A+.
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
% nvar:  number of endogeous variables
% ncoef:  number of original lagged variables per equation
% np: nvar-element vector, ith element represents the number of free A+ parameters in ith equation
%---------------
% F: ncoef-by-nvar matrix of original lagged parameters A+.  Column corresponding to equation.
%
% August 2000, Tao Zha.
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

g=g(:); np=np(:);
npcum = [0;cumsum(np)];
F = zeros(ncoef,nvar); % ncoef: maximum original lagged parameters per equation
for kj=1:nvar
   F(:,kj) = Vi{kj}*g(npcum(kj)+1:npcum(kj+1));
end
