function [Fmat,gvec] = gfmean(b,P,Vi,nvar,ncoef,n0,np)
% [Fmat,gvec] = gfmean(b,P,Vi,nvar,ncoef,n0,np)
%
%    Mean of free lagged parameters g and original lagged parameters F, conditional on comtemporaneous b's
%    See Waggoner and Zha's Gibbs sampling
%
% b: sum(n0)-by-1 vector of A0 free parameters
% P: cell(nvar,1).  In each cell, posterior linear transformation for random walk prior for the ith equation % tld: tilda
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
% nvar:  number of endogeous variables
% ncoef:  number of original lagged variables per equation
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
% np: nvar-by-1, ith element represents the number of free A+ parameters in ith equation
%---------------
% Fmat: 0, the value that is used in csminwel.m
% gvec: sum(n0)-by-1 analytical gradient for a0freefun.m
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


n0cum = cumsum(n0);
npcum = cumsum(np);
gvec = zeros(npcum(end),1);
Fmat = zeros(ncoef,nvar); % ncoef: maximum original lagged parameters per equation

if ~(length(b)==n0cum(end))
   error('Make inputs n0 and length(b) match exactly')
end

for kj=1:nvar
   if kj==1
      bj = b(1:n0(kj));
      gj = P{kj}*bj;
      gvec(1:np(kj)) = gj;
      Fmat(:,kj) = Vi{kj}*gj;
   else
      bj = b(n0cum(kj-1)+1:n0cum(kj));
      gj = P{kj}*bj;
      gvec(npcum(kj-1)+1:npcum(kj)) = gj;
      Fmat(:,kj) = Vi{kj}*gj;
   end
end
