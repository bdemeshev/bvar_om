function [P,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Ui,Vi)
% [P,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Ui,Vi)
%
%    Exporting deterministic (no random prior) Bayesian posterior matrices with linear restrictions
%    See Waggoner and Zha's Gibbs sampling paper
%
% xtx:  X'X: k-by-k where k=ncoef
% xty:  X'Y: k-by-nvar
% yty:  Y'Y: nvar-by-nvar
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.  Imported from dnrprior.m.
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k (ncoef) is a total number of RHS variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.  Imported from dnrprior.m.
%-----------------
% P: cell(nvar,1).  In each cell, the transformation matrix that affects the posterior mean of A+ conditional on A0.
% H0inv: cell(nvar,1).  In each cell, posterior inverse of covariance matrix (i.e., the exponential
%           term is a_i'*H0*a_i) for the ith equation (not divided by T yet).  It resembles
%           old SpH in the exponent term in posterior of A0, but not divided by T yet.
% Hpinv: cell(nvar,1).  In each cell, posterior inverse of covariance matrix Hp (A+) for the
%           ith equation.
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


nvar = size(yty,1);

P = cell(nvar,1); % tld: tilda
H0inv = cell(nvar,1);  % posterior inv(H0), resemble old SpH, but not divided by T yet.
Hpinv = cell(nvar,1);  % posterior inv(Hp).

for n=1:nvar       % one for each equation
   Hpinv{n} = Vi{n}'*xtx*Vi{n};
   P1 = Vi{n}'*xty*Ui{n};
   P{n} = Hpinv{n}\P1;
   H0inv{n} =  Ui{n}'*yty*Ui{n} - P1'*(Hpinv{n}\P1);
end

