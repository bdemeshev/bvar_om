function [YCom_t, YCom_tm1] = fn_data_companionform(nvar,lags,ydata)
%[YCom_t, YCom_tm1] = fn_data_companionform(nvar,lags,ydata)
%
% nvar: number of endogenous variables.
% lags: number of lags.
% ydata: nvar X (fss+lags) matrix of raw or original data (no manipulation involved)
%       with sample size including lags.
%-------
% YCom_t: nvar*lags X fss companion vector at time t
% YCom_tm1: nvar*lags X fss companion vector at time t-1.
%
% Tao Zha, March 2014
% Copyright (C) 2014- Tao Zha
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


%--------- Strack the VAR(p) to the VAR(1) with YCom_t = C_Com + BhatCom YCom_{t-1} + UCom_{t}.
nSample = size(ydata,2);
sb = lags+1;   % original beginning without dummies
fss = nSample-lags;
ncoefwoc = nvar*lags;
%
% ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
% **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, exo var, const]
%
YCom_t = zeros(ncoefwoc,fss);
YCom_tm1 = zeros(ncoefwoc,fss);
for k=1:lags
   YCom_t(nvar*(k-1)+1:nvar*k,:) = ydata(:,(sb-k+1):(nSample-k+1));
   YCom_tm1(nvar*(k-1)+1:nvar*k,:) = ydata(:,sb-k:nSample-k);
   % row: [nvar for 1st lag, ..., nvar for last lag]; column: fss; 
   %        Thus, # of rows is nvar*lags.                 
end



