function [BhatCom, C_Com] = fn_companionform(Bhat,nvar,lags)
%BhatC = fn_companionform(Bhat,nvar,lags)
%
%    Using eigenvalues to find the inverse of all roots associated with the VAR proceess:
%          y_t' = C + y_{t-1}'*B_1 + ... + Y_{t-p}'*B_p + u_t'.
%    where columns correspond to equations.  See also Judge (1), pp.753-755 where rows correspond to equations.
% Bhat:  ncoef-by-nvar where ncoef=nvar*lags+nexo and nvar is the number of endogenous variables.
%    Columns corresponds to equations with
%    ncoef=[nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
%                       ..., nvar coef in the last lag, and nexo coefficients.
%    Note that entries in the rows of Bhat that > nvar*lags are irrelevant.
% nvar: number of endogenous variables.
% lags: number of lags.
%-------
% BhatCom: nvar*lags X nvar*lags companion matrix to strack the VAR(p) to the VAR(1) with z_t = BhatC z_{t-1}.
% C_Com:  nvar*lags X 1 companion constant vector.
%
% Tao Zha, September 2013, Revised, March 2014
% Copyright (C) 1997-2014 Tao Zha
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


if size(Bhat,1)<nvar*lags
   disp(' ')
   warning('Make sure that Bhat has at least nvar*lags rows')
   return
end

%--------- Strack the VAR(p) to the VAR(1) with z_t = Az_{t-1}.
%
A1 = diag(ones(nvar*(lags-1),1));
A2 = [A1 zeros(nvar*(lags-1),nvar)];
A = [Bhat(1:nvar*lags,:)'; A2];
BhatCom = A;

C_Com = cat(1, (Bhat(end,:))', zeros(nvar*(lags-1),1));
%rootsinv=eig(A);
