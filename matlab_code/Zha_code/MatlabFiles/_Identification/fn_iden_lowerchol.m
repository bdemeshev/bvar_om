function [Ui,Vi,n0,np] = fn_iden_lowerchol(lags,nvar,nexo)
%  [Ui,Vi,n0,np] = fn_iden_lowerchol(lags,nvar,nexo)
%  Lowercholeski decompositions with columns corresponding to equations.  Automatically applies to all variables.
%
%    Exporting orthonormal matrices for the deterministic linear restrictions (equation by equation)
%    See Waggoner and Zha's Gibbs sampling paper.
%
% q_m:  quarter or month
% lags: the maximum length of lag
% nvar:  number of endogeous variables
% nexo:  number of exogenous variables.  If nexo is not supplied, nexo=1 as default for a constant
% indxC0Pres: index for cross-A0-A+ restrictions.  if 1: cross-A0-and-A+ restrictions; 0: idfile is all we have
%                Example for indxOres==1: restrictions of the form P(t) = P(t-1).
%                These restrictions have to be manually and carefully keyed in.
%-----------------
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
% np: nvar-by-1, ith element represents the number of free A+ parameters in ith equation
%
% Copyright (C) 2000-2012 Tao Zha
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


Ui = cell(nvar,1);  % initializing for contemporaneous endogenous variables
Vi = cell(nvar,1);  % initializing for lagged and exogenous variables
n0 = zeros(nvar,1); % ith element represents the number of free A0 parameters in ith equation
np = zeros(nvar,1); % ith element represents the number of free A+ parameters in ith equation

if (nargin==2)
   nexo = 1;  % 1: constant as default where nexo must be a nonnegative integer
elseif (nargin==3)
   indxC0Pres = 0;  % default is no cross-A0-and-A+ restrictions.
end

k = lags*nvar+nexo;  % maximum number of lagged and exogenous variables in each equation

Qi = zeros(nvar,nvar,nvar);   % for nvar contemporaneous equations
Ri = zeros(k,k,nvar);    % for nvar lagged and exogenous equations
  % Row corresponds to equation. 0 means no restriction.
  %                              1 means exclusion restriction such that the corresponding parameter is restricted to 0.

%nfvar = 6;   % number of foreign (Granger causing) variables
%nhvar = nvar-nfvar;  % number of home (affected) variables.

identity=eye(nvar);
for i=2:nvar
    Qi(1:i-1,:,i)=identity(1:i-1,:);
end

for n=1:nvar   %  initializing loop for each equation
   Ui{n} = null(Qi(:,:,n));
   Vi{n} = null(Ri(:,:,n));
   n0(n) = size(Ui{n},2);
   np(n) = size(Vi{n},2);
end


