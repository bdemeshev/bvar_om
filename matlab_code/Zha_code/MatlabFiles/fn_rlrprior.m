function [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Ui,Vi,Pi,H0multi,Hpmulti,nvar)
% [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Ui,Vi,Pi,H0multi,Hpmulti,nvar)
%
%    Exporting random Bayesian prior with linear restrictions
%    See Waggoner and Zha's Gibbs sampling paper
%
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
% Pi: ncoef-by-nvar matrix for the ith equation under random walk.  Same for all equations
% H0multi: nvar-by-nvar-by-nvar; H0 for different equations under asymmetric prior
% Hpmulti: ncoef-by-ncoef-by-nvar; H+ for different equations under asymmetric prior
% nvar:  number of endogenous variables
% --------------------
% Ptld: cell(nvar,1).  The prior mean of g_i is Ptld{i}*b_i;
% H0invtld: cell(nvar,1). Transformed inv covaraince for b_i, the free parameters in A0(:,i);
% Hpinvtld: cell(nvar,1). Transformed inv covaraince for g_i, the free parameters in A+(:,i);
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

Ptld = cell(nvar,1); % tld: tilda
H0invtld = cell(nvar,1);  % H0 for different equations under linear restrictions
Hpinvtld = cell(nvar,1);  % H+ for different equations under linear restrictions

for n=1:nvar       % one for each equation
   Hpinvtld{n} = Vi{n}'*(Hpmulti(:,:,n)\Vi{n});
   Ptld{n} = (Hpinvtld{n}\Vi{n}')*(Hpmulti(:,:,n)\Pi)*Ui{n};
   H0invtld{n} = Ui{n}'*(H0multi(:,:,n)\Ui{n}) + Ui{n}'*Pi'*(Hpmulti(:,:,n)\Pi)*Ui{n} ...
                      - Ptld{n}'*Hpinvtld{n}*Ptld{n};
end
