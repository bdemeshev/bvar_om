function [trends_ave, trends_subsample, relcon_stationary_ave, relcon_stationary_subsample] = fn_annualizedgrowthrates2trends(BhatCom, Bhat, phi, nvar,lags, q_m, n_eigclose1)
%[trends_ave, trends_subsample] = fn_annualizedgrowthrates2trends(BhatCom, Bhat, phi, nvar,lags, q_m, n_eigclose1)
%
% BhatCom: nvar*lags X nvar*lags companion matrix to strack the VAR(p) to the VAR(1) with z_t = BhatC z_{t-1}.
% Bhat:  ncoef-by-nvar where ncoef=nvar*lags+nexo and nvar is the number of endogenous variables.
%    Columns corresponds to equations with
%    ncoef=[nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
%                       ..., nvar coef in the last lag, and nexo coefficients.
%    Note that entries in the rows of Bhat that > nvar*lags are irrelevant.
% phi:  the Tsub X (nvar*lags+1) subsample of right-hand-side variables in the VAR.
% nvar: number of endogenous variables.
% lags: number of lags.
% q_m: 12 (monthly) or 4 (quarterly).
% n_eigclose1: number of eigenvalues close to 1. 
% flag_display: 1: displaying the relative contribution of stationary part to growth (should be << 1); 0: no such a display.
%---
% trends_subsample: Tsub X nvar matrix of trends.
% trends_ave:    nvar X 1 vector of average trends.
% relcon_stationary_subsample: Tsub X nvar matrix of relative contributions of stationary part to growth (should be << 1).
% relcon_stationary_ave:  nvar X 1 vector of average relative contributions of stationary part to growth (should be << 1).
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
nsubsample = size(phi,1);

[eigU, eigD] = eig(BhatCom);
[eigdsort, indxascend] = sort(diag(eigD),1,'ascend');
eigUsort = eigU(:,indxascend);

companionC = zeros(nvar*lags,1);
companionC(1:nvar) = (Bhat(end,:))';
gbetas = eigUsort\companionC;

trends_subsample = zeros(nsubsample,nvar);
relcon_stationary_subsample = zeros(nsubsample,nvar);
for (ti=1:nsubsample)
   phi_ini = phi(ti,:);
   galphas = eigUsort\(phi_ini(1:(end-1)))';
   factors = (eigdsort-1.0) .* galphas + gbetas;
   %factors = gbetas;
   pickfactors = factors;
   pickfactors(1:(nvar*lags-n_eigclose1)) = 0;
   growthrates0 = eigUsort*pickfactors;
   agrates_percent = growthrates0(1:nvar)*q_m*100;
   %trends_subsample(ti,:) = real(agrates_percent');
   trends_subsample(ti,:) = agrates_percent';
   
   %--- Report the relative contribution of stationary part to growth (should be << 1).
   zerofactors = factors;
   zerofactors(((nvar*lags-n_eigclose1)+1):end) = 0;   
   zerogrowth = eigUsort*zerofactors;
   relcon_stationary = real(zerogrowth(1:nvar)) ./ real(growthrates0(1:nvar));
   relcon_stationary_subsample(ti,:) = relcon_stationary';
   if (0)   
      %disp('------- factors: ----------')
      %abs(factors)

      %disp('------- eigUsort*zerofactors (should be zero relative to eigUsort*pickfactors): ----------')
      %real(zerogrowth(1:nvar))
      
      %disp('------- eigUsort*pickfactors: ----------')
      %real(growthrates0(1:nvar))
      
      disp('------- (eigUsort*zerofactors)/(eigUsort*pickfactors) (should be << 1): ----------')
      real(zerogrowth(1:nvar)) ./ real(growthrates0(1:nvar))      
   end   
end
%--- Relative contribution of stationary part to growth (should be << 1):
relcon_stationary_ave = (mean(relcon_stationary_subsample))';
%--- Accurate (best) growth rates:
trends_ave = (mean(trends_subsample))';
