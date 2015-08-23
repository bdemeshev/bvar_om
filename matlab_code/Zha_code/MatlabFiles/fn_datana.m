function [yactyrge,yactyre,yactqmyge,yactqmge,yactqme] = fn_datana(xdatae,q_m,vlistlog,vlistper,Byrqm,Eyrqm)
% [yactyrge,yactyre,yactqmyge,yactqmge,yactqme] = fn_datana(Byrqm,Eyrqm,xdatae,q_m,vlistlog,vlistper)
%
%      Generate prior period, period-to-period-in-last-year, and annual growth rates
%      For annual rates, works for both calendar and any annual years, depending on Byrqm and Eyrqm
%      If monthly data, we haven't got time to get quarterly growth rates yet.
%
% xdatae:  all data (logged levels or interest rates/100, some of which may be NaN) with the first
%          2 columns indicating years and periods.
% q_m: quarter or month period
% vlistlog: sublist for logged variables
% vlistper: sublists for percent variables
% Byrqm: [year quarter(month)]: beginning year and period.  If Byqm(2)~=1, we don't get
%     calendar annual rates.  In other words, the first column of yactyge (which
%     indicates years) does not mean calendar years.  Byqm(2) must be specified; in other
%     words, it must be not set to 0 as in, say, fn_dataext.
% Eyrqm: [year period]:  end year and period.  Eyqm(2) must be specified; in other words, it
%     must be not set to 0 as in, say, fn_dataext.
%    NOTE: if no inputs Byrqm and Eyrqm are specified, all growth rates begin at xdatae(1,1:2).
%----------
% yactyrge: annual growth rates with dates in the first 2 columns.
% yactyre:  annual average logged level with dates in the 1st 2 columns.
% yactqmyge: period-to-period-in-last-year annual growth rates with dates in the first 2 columns.
% yactqmge:  prior-period annualized growth rates with dates in the first 2 columns.
% yactqme:  data (logged levels or interest rates/100) with dates in the first 2 columns.
%           Same as xdatae but with Brow:Erow.
%
% Tao Zha, April 2000.
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


if size(xdatae,1)<2*q_m
   error('We need at least two years of xdatae to get annual rates.  Check xdatae!!')
end

if nargin==4
   Brow=1; Erow=length(xdatae(:,1));
   nyr = floor((Erow-Brow+1)/q_m);
   yrsg = [xdatae(q_m+1,1):xdatae(q_m+1,1)+nyr-2]';    % for annual growth later on
else
   if Byrqm(2)<1 | Eyrqm(2)<1
      error('This function requires specifying both years and months (quarters) in Byrqm and Eyrqm')
   end

   Brow = min(find(xdatae(:,1)==Byrqm(1)));
   if isempty(Brow)
      error('Byrqm is outside the date range of xdatae(:,1:2)!')
   end
   nadt=Byrqm(2)-xdatae(Brow,2);
   if nadt<0
      error('Byrqm is outside the date range indicated by xdatae(:,1:2)!')
   end
   Brow=Brow+nadt;
   %
   Erow = min(find(xdatae(:,1)==Eyrqm(1)));
   if isempty(Brow)
      error('Eyrqm is outside the date range of xdatae(:,1:2)!')
   end
   nadt2=Eyrqm(2)-xdatae(Erow,2);
   if nadt<0
      error('Eyrqm is outside the date range indicated by xdatae(:,1:2)!')
   end
   Erow=Erow+nadt2;

   nyr = floor((Erow-Brow+1)/q_m);
   yrsg = [Byrqm(1)+1:Byrqm(1)+nyr-1]';    % for annual growth later on, which will
                %   start at Byrqm(1) instead of Byrqm(1)+1
end
%
yactqme = xdatae(Brow:Erow,:);   % with dates
yactqm = yactqme(:,3:end);   % only data

%======== prior period change (annaluized rate)
yactqmg = yactqm(2:end,:);    % start at second period to get growth rate
yactqmg(:,vlistlog) = (yactqm(2:end,vlistlog) - yactqm(1:end-1,vlistlog)) .* q_m;
                             % monthly, 12*log(1+growth rate), annualized growth rate
yactqmg(:,vlistlog) = 100*(exp(yactqmg(:,vlistlog))-1);
yactqmg(:,vlistper) = 100*yactqmg(:,vlistper);
yactqmge = [yactqme(2:end,1:2) yactqmg];

%======== change from the last year
yactqmyg = yactqm(q_m+1:end,:);    % start at the last-year period to get growth rate
yactqmyg(:,vlistlog) = (yactqm(q_m+1:end,vlistlog) - yactqm(1:end-q_m,vlistlog));
yactqmyg(:,vlistlog) = 100*(exp(yactqmyg(:,vlistlog))-1);
yactqmyg(:,vlistper) = 100*yactqmyg(:,vlistper);
yactqmyge = [yactqme(q_m+1:end,1:2) yactqmyg];

%======== annual growth rates
nvar = length(xdatae(1,3:end));
ygmts = yactqm(1:nyr*q_m,:);   % converted to the multiplication of q_m
ygmts1 = reshape(ygmts,q_m,nyr,nvar);
ygmts2 = sum(ygmts1,1) ./ q_m;
ygmts3 = reshape(ygmts2,nyr,nvar);  % converted to annual average series
%
yactyrg = ygmts3(2:end,:);    % start at the last-year period to get growth rate
yactyrg(:,vlistlog) = ygmts3(2:end,vlistlog) - ygmts3(1:end-1,vlistlog);
                          % annual rate: log(1+growth rate)
yactyrg(:,vlistlog) = 100*(exp(yactyrg(:,vlistlog))-1);
yactyrg(:,vlistper) = 100*yactyrg(:,vlistper);
yactyrge = [yrsg zeros(nyr-1,1) yactyrg];
yrsg1=[yrsg(1)-1:yrsg(end)]';
yactyre = [yrsg1 zeros(nyr,1) ygmts3];
