function [P,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0invtld,Hpinvtld,Ui,Vi)
% [P,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0tld,Hptld,Ui,Vi)
%
%    Exporting random (i.e., random prior) Bayesian posterior matrices with linear restrictions
%    See Waggoner and Zha's Gibbs sampling paper
%
% xtx:  X'X: k-by-k where k=ncoef
% xty:  X'Y: k-by-nvar
% yty:  Y'Y: nvar-by-nvar
% Ptld: cell(nvar,1), transformation matrix that affects the (random walk) prior mean of A+ conditional on A0.
% H0invtld: cell(nvar,1), transformed inv covaraince for free parameters in A0(:,i).
% Hpinvtld: cell(nvar,1), transformed inv covaraince for free parameters in A+(:,i);
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
%      In other words, the posterior mean (of g_i) = P{i}*b_i where g_i is a column vector of free parameters
%        of A+(:,i)) given b_i (b_i is a column vector of free parameters of A0(:,i)).
% H0inv: cell(nvar,1).  Not divided by T yet.  In each cell, inverse of posterior covariance matrix H0.
%          The exponential term is b_i'*inv(H0)*b_i for the ith equation where b_i = U_i*a0_i.
%          It resembles old SpH or Sbd in the exponent term in posterior of A0, but not divided by T yet.
% Hpinv: cell(nvar,1).  In each cell, posterior inverse of covariance matrix Hp (A+) for the free parameters
%          g_i = V_i*A+(:,i) in the ith equation.
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
   Hpinv{n} = Vi{n}'*xtx*Vi{n} + Hpinvtld{n};
   P1 = Vi{n}'*xty*Ui{n} + Hpinvtld{n}*Ptld{n};
   P{n} = Hpinv{n}\P1;
   H0inv{n} =  Ui{n}'*yty*Ui{n} + H0invtld{n} + Ptld{n}'*Hpinvtld{n}*Ptld{n} ...
                   - P1'*P{n};     %P{n} = (Hpinv{n}\P1);
end
