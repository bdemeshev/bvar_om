function agrates_percent = fn_annualizedgrowthrates(BhatCom, Bhat, nvar,lags, q_m, n_eigclose1)
%agrates_percent = fn_annualizedgrowthrates(BhatCom, Bhat, nvar,lags, q_m, n_eigclose1)
%
% BhatCom: nvar*lags X nvar*lags companion matrix to strack the VAR(p) to the VAR(1) with z_t = BhatC z_{t-1}.
% Bhat:  ncoef-by-nvar where ncoef=nvar*lags+nexo and nvar is the number of endogenous variables.
%    Columns corresponds to equations with
%    ncoef=[nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
%                       ..., nvar coef in the last lag, and nexo coefficients.
%    Note that entries in the rows of Bhat that > nvar*lags are irrelevant.
% nvar: number of endogenous variables.
% lags: number of lags.
% q_m: 12 (monthly) or 4 (quarterly).
% n_eigclose1: number of eigenvalues close to 1. 
%---
% agrates:  annualized growth rates.
%
% See fn_annualizedgrowthrates2trends.m for better calculation of trend growth rates.
%
% Copyright (C) 1997-2013 Daniel Waggoner and Tao Zha
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

if (size(BhatCom,1) ~= (nvar*lags))
   error('fn_annualizedgrowthrates(): Check the dimension of BhatCom against nvar and lags');
end   
[eigU, eigD] = eig(BhatCom);
[eigDsort, indxascend] = sort(diag(eigD),1,'ascend');
eigUsort = eigU(:,indxascend);

companionC = zeros(nvar*lags,1);
companionC(1:nvar) = (Bhat(end,:))';
gbetas = eigUsort\companionC;

pickbetas = gbetas;
pickbetas(1:(nvar*lags-n_eigclose1)) = 0;
growthrates0 = eigUsort*pickbetas;
agrates_percent = growthrates0(1:nvar)*q_m*100;