function [Tinv,UT,VHphalf,PU,VPU] = fn_gibbsrvar_setup(H0inv, Ui, Hpinv, Pmat, Vi, nvar, fss)
% [Tinv,UT,VHphalf,PU,VPU] = fn_gibbsrvar_setup.m(H0inv, Ui, Hpinv, Pmat, Vi, fss, nvar)
%    Global setup outside the Gibbs loop to be used by fn_gibbsvar().
%    Reference: "A Gibbs sampler for structural VARs" by D.F. Waggoner and T. Zha,  ``
%               Journal of Economic Dynamics & Control (JEDC) 28 (2003) 349-366.
%    See Note Forecast (2) pp. 44-51, 70-71, and Theorem 1 and Section 3.1 in the WZ JEDC paper.
%
% H0inv: cell(nvar,1).  Not divided by T yet.  In each cell, inverse of posterior covariance matrix H0.
%          The exponential term is b_i'*inv(H0)*b_i for the ith equation where b_i = U_i*a0_i.
%          It resembles old SpH or Sbd in the exponent term in posterior of A0, but not divided by T yet.
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.  Imported from dnrprior.m.
% Hpinv: cell(nvar,1).  In each cell, posterior inverse of covariance matrix Hp (A+) for the free parameters
%          g_i = V_i*A+(:,i) in the ith equation.
% Pmat: cell(nvar,1).  In each cell, the transformation matrix that affects the posterior mean of A+ conditional on A0.
%      In other words, the posterior mean (of g_i) = Pmat{i}*b_i where g_i is a column vector of free parameters
%        of A+(:,i)) given b_i (b_i is a column vector of free parameters of A0(:,i)).
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k (ncoef) is a total number of RHS variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.  Imported from dnrprior.m.
% nvar:  number of endogenous variables or rank of A0.
% fss: effective sample size (in the exponential term) = nSample - lags + ndobs (ndobs = # of dummy observations
%        is set to 0 when fn_rnrprior_covres_dobs() is used where dummy observations are included as part of the explicit prior.
%-------------
% Tinv: cell(nvar,1).  In each cell, inv(T_i) for T_iT_i'=S_i where S_i is defined on p.355 of the WZ JEDC paper.
% UT: cell(nvar,1).   In each cell, U_i*T_i.
% VHphalf: cell(nvar,1).  In each cell, V_i*sqrt(Hp_i).
% PU: cell(nvar,1).  In each cell, Pmat{i}*U_i where Pmat{i} = P_i defined in (13) on p.353 of the WZ JEDC paper.
% VPU: cell(nvar,1).  In each cell, V_i*P_i*U_i
%
% Written by Tao Zha, September 2004.
% Copyright (c) 2004 by Waggoner and Zha
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



%--- For A0.
Tinv = cell(nvar,1);   % in each cell, inv(T_i) for T_iT_i'=S_i where S_i is defined on p.355 of the WZ JEDC paper.
UT = cell(nvar,1);  % in each cell, U_i*T_i.
%--- For A+.
VHphalf = cell(nvar,1);  % in each cell, V_i*sqrt(Hp_i).
PU = cell(nvar,1);  % in each cell, Pmat{i}*U_i where Pmat{i} = P_i defined in (13) on p.353 of the WZ JEDC paper.
VPU = cell(nvar,1);  % in each cell, V_i*P_i*U_i
%
for ki=1:nvar
   %--- For A0.
   Tinv{ki} = chol(H0inv{ki}/fss);   % Tinv_i'*Tinv_i = inv(S_i) ==> T_i*T_i' = S_i where S_i = H0inv{i}/fss is defined on p.355 of the WZ JEDC paper.
   UT{ki} = Ui{ki}/Tinv{ki};    % n-by-qi: U_i*T_i in (14) on p. 255 of the WZ JEDC paper.
   %--- For A+.
   VHphalf{ki} = Vi{ki}/chol(Hpinv{ki}); % where chol(Hpinv_i)*chol(Hpinv_i)'=Hpinv_i.
   PU{ki} = Pmat{ki}*Ui{ki}';
   VPU{ki} = Vi{ki}*PU{ki};
end
