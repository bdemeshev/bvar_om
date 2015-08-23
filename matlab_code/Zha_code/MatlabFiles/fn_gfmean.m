function [Fmat,gvec] = fn_gfmean(b,P,Vi,nvar,ncoef,n0,np)
% [Fmat,gvec] = fn_gfmean(b,P,Vi,nvar,ncoef,n0,np)
%
%    Mean of free lagged parameters g and original lagged parameters F, conditional on comtemporaneous b's
%    See Waggoner and Zha's Gibbs sampling
%
% b: sum(n0)-element vector of mean estimate of A0 free parameters
% P: cell(nvar,1).  In each cell, the transformation matrix that affects the posterior mean of A+ conditional on A0.
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
% nvar:  number of endogeous variables
% ncoef:  number of original lagged variables per equation
% n0: nvar-element vector, ith element represents the number of free A0 parameters in ith equation
% np: nvar-element vector, ith element represents the number of free A+ parameters in ith equation
%---------------
% Fmat: ncoef-by-nvar matrix of original lagged parameters A+.  Column corresponding to equation.
% gvec: sum(np)-by-1 stacked vector of all free lagged parameters A+.
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


b=b(:); n0=n0(:); np=np(:);

n0cum = [0;cumsum(n0)];
npcum = [0;cumsum(np)];
gvec = zeros(npcum(end),1);
Fmat = zeros(ncoef,nvar); % ncoef: maximum original lagged parameters per equation

if ~(length(b)==n0cum(end))
   error('Make inputs n0 and length(b) match exactly')
end

for kj=1:nvar
   bj = b(n0cum(kj)+1:n0cum(kj+1));
   gj = P{kj}*bj;
   gvec(npcum(kj)+1:npcum(kj+1)) = gj;
   Fmat(:,kj) = Vi{kj}*gj;
end
