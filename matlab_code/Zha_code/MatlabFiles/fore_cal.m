function [yforelml,yforemgml,yforeqgml,yforeCalygml] = fore_cal(yfore,xdata,nvar,...
              nSample,nSampleCal,forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,...
				  vlistper,lmqyIndx)
% Converting oringinal forecast "yfore" to series by calendar years and series
%                  with growth rates (annualized for monthly and quarterly series).
%
%function [yforelml,yforemgml,yforeqgml,yforeCalygml] = fore_cal(yfore,xdata,nvar,...
%              nSample,nSampleCal,forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,...
%				  vlistper,mqyIndx)
%
% yfore:   oringal forecast series, all logged except R, U, etc.
% xdata:   oringal data set beyond the sample into forecast horizon
%                           until yrFin:qmFin, all logged except R, U, etc.
% nvar:  number of variables
% nSample:  sample size (including lags or initial periods)
% nSampleCal:  sample size in terms of calendar years
% forep:   forecast periods (monthly)
% forepq:  forecast periods (quarterly)
% forepy:  forecast periods (yearly)
% q_m:   quarterly or monthly for the underlying model
% qmEnd:  last q_m before out-of-sample forecasting
% vlist:  a list of variables
% vlistlog: sub list of variables that are in log
% vlistper: sub list of variables that are in percent
% lmyqIndx:  4-by-1 1 or 0's. Index for level, mg, qg, and yg; 1: yes; 0: no
%----------------
% yforelml:  monthly level forecast (in percent for R, U, etc.) with the same size as "yfore" in log
% yforemgml:   monthly growth (at annual rates), in percent
% yforeqgml:  ML forecast: quarterly growth (at annual rates), in percent
% yforeCalygml:  ML forecast: annual growth (by calendar years), in percent
%              forepy-by-nvar
%
% Copyright (c) March 1998 by Tao Zha
% Revision, October 1998.  Added lmyqIndx so previous programs may not be compatible.
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


%=================================================
%  Making everything presentable at FOMC
%=================================================
%
%%%%
%$$$  Monthly, prior quarter, and year-over-year change
%$$$           Out-of-sample forecasts
%%%%


if length(lmqyIndx)~=4
	warning('lmqyIndx must be a 4-by-1 vector containing 1 or 0')
	return
end



%---------------------
%  Actual data
%---------------------
yact = xdata(nSample-q_m+1:nSample,:);   % the latest (not calender) year
yactCal = xdata(nSampleCal-q_m+1:nSampleCal,:);   % the lastest calendar year



%-----------------------------------------
%   Converted to monthly level
%-----------------------------------------
if lmqyIndx(1)
	yforelml=yfore;
	yforelml(:,vlistlog)=exp(yfore(:,vlistlog));       % mode, all levels
	yforelml(:,vlistper)=100*yfore(:,vlistper);       % mode, all levels
else
	yforelml = NaN;
end



%-----------------------------------------
%   Converted to monthly growth
%-----------------------------------------
if lmqyIndx(2)
	yactm = zeros(1,length(vlist));  % the latest month prior to forecasting
	yforem = zeros(1+forep,length(vlist));   % including the 1 month prior to forecasting
	%
	yactm = yact(length(yact(:,1)),:);      % last month prior to forecasting
	yforem(1,:) = yactm;                   % last month prior to forecasting
	yforem(2:forep+1,:) = yfore(1:forep,:);  % monthly forecasts, all logged.
	%
	%@@ monthly growth rate (annualized)
	yforemg = yforem(2:forep+1,:);
	yforemg(:,vlistlog) = ...
	               ( yforem(2:forep+1,vlistlog) - yforem(1:forep,vlistlog) ) .* q_m;
	                      % monthly change (annualized), 12*log(1+growth rate)
	yforemgml=yforemg;
	yforemgml(:,vlistlog) = 100*(exp(yforemg(:,vlistlog))-1);
	                      % monthly growth rate (annualized)
	yforemgml(:,vlistper) = 100*yforemg(:,vlistper);   % monthly growth rate (annualized)
else
	yforemgml=NaN;
end



%-----------------------------------------
%   Converted to quarterly
%-----------------------------------------
if lmqyIndx(3)
	yactQ = xdata(nSample-mod(qmEnd,3)-q_m+1:nSample-mod(qmEnd,3),:);
	                              % the latest actual quarter
	yactq = zeros(4,length(vlist));  % the latest 4 quarters prior to forecasting
	yforeq = zeros(4+forepq,length(vlist));
	                              % including the 4 quarters prior to forecasting
	%
	qT1 = length(yactQ(:,1))/3;
	if qT1 ~= 4
		error('Must hae 4 actual quarters to compute quarterly change for forecasting!')
	end
	for i = 1:qT1
	   i1 = 1+3*(i-1);
	   i2 = 3*i;
	   yactq(i,:) = sum(yact(i1:i2,:)) ./ 3;
	end
	yforeq(1:4,:) = yactq;
	%
	yforeq_1 = sum(xdata( nSample-mod(qmEnd,3)+1:nSample,:),1);
	yforeq(5,:) = ( yforeq_1 + sum(yfore(1:3-mod(qmEnd,3),:),1) ) ./ 3;
				% 1st quarterly forecast which may contain actual data within quarter
	         % note, dimension "1" in sum(ytem,1) is necessary because when
	         %   qmEnd=1, sum(ytem,1) will give us a nomber, not a vector.  But this
				%   is not what we want.
	for i = 2:forepq
	   i1 = 3-mod(qmEnd,3) + 1+3*(i-2);
	   i2 = 3-mod(qmEnd,3) + 3*(i-1);
	   yforeq(4+i,:) = sum(yfore(i1:i2,:)) ./ 3;
	end
	%
	%@@ prior quarter growth rate (annualized)
	yforeqg = yforeq(5:forepq+4,:);
	%yforeqg(:,vlistlog) = (100*(yforeq(5:qT2+4,vlistlog)-yforeq(1:qT2,vlistlog))) ...
	%                       ./ yforeq(1:qT2,vlistlog);   % year-over-year
	%yforeqg(:,vlistlog) = 100*((yforeq(5:qT2+4,vlistlog)./yforeq(4:qT2+3,vlistlog)).^4 - 1 );
	                             % prior quarter
	yforeqg(:,vlistlog) = ...
	              ( yforeq(5:forepq+4,vlistlog) - yforeq(4:forepq+3,vlistlog) ) .* 4;
	                             % prior quarter, 4*log(1+growth rate)
	yforeqgml=yforeqg;
	yforeqgml(:,vlistlog) = 100*(exp(yforeqg(:,vlistlog))-1);
	                             % quarterly growth rate (annualized)
	yforeqgml(:,vlistper) = 100*yforeqg(:,vlistper);
	                             % quarterly growth rate (annualized)
else
	yforeqgml = NaN;
end



%-----------------------------------------
%   Converted to calendar years
%-----------------------------------------
if lmqyIndx(4)
	yactCaly = zeros(1,length(vlist));  % the latest calendar year
	yforeCaly = zeros(1+forepy,length(vlist));
	                         % including the calendar year prior to forecasting
	%
	yT1 = length(yactCal(:,1))/q_m;
	if yT1 ~= 1
		error('yT1 Beginings or ends of monthly and calendar series are not the same!')
	end
	for i = 1:yT1
	   i1 = 1+q_m*(i-1);
	   i2 = q_m*i;
	   yactCaly(i,:) = sum(yactCal(i1:i2,:)) ./ q_m;
	end
	yforeCaly(1,:) = yactCaly;
	%
	%@@ initial monthly actual data for calendar years
	if qmEnd == q_m
	   yforeCaly_1 = 0;
	else
	   ytem = xdata(nSampleCal+1:nSampleCal+qmEnd,:);
	   %ytem(:,vlistlog) = exp(ytem(:,vlistlog));
	   yforeCaly_1 = sum(ytem,1);
	         % note, dimension "1" in sum(ytem,1) is necessary because when
	         %   qmEnd=1, sum(ytem,1) will give us a nomber, not a vector.  But this
				%   is not what we want.
	end
	%
	if qmEnd == q_m
	   for i = 1:forepy
	      i1 = 1+q_m*(i-1);
	      i2 = q_m*i;
	      yforeCaly(1+i,:) = sum(yfore(i1:i2,:)) ./ q_m;
	   end
	else
	   yforeCaly(2,:) = (yforeCaly_1+sum(yfore(1:q_m-qmEnd,:),1)) ./ q_m;
				% note, dimension "1" in sum(yfore,1) is necessary because when
	         %   q_m-qmEnd=1, sum(yfore) will give us a number, not a vector.  But this
				%   is not what we want.
	   for i = 2:forepy
	      i1 = q_m-qmEnd+1+q_m*(i-2);
	      i2 = q_m-qmEnd+q_m*(i-1);
	      yforeCaly(1+i,:) = sum(yfore(i1:i2,:)) ./ q_m;
	   end
	end
	%
	%@@ year-over-year growth rate
	yforeCalyg = yforeCaly(2:forepy+1,:);
	yforeCalyg(:,vlistlog) = ...
	                  yforeCaly(2:forepy+1,vlistlog) - yforeCaly(1:forepy,vlistlog);
	                             % year-over-year, log(1+growth rate)
	yforeCalygml=yforeCalyg;
	yforeCalygml(:,vlistlog) = 100*(exp(yforeCalyg(:,vlistlog))-1);
	                                                         % annaul growth rate
	yforeCalygml(:,vlistper) = 100*yforeCalyg(:,vlistper);   % annaul growth rate
else
	yforeCalygml = NaN;
end



%ndraws
%
% posterior mean of out-of-sample forecasts
%yfore1 = (1.0/ndraws)*yfore1;       % mean
%yforeqg1 = (1.0/ndraws)*yforeqg1;       % mean
%yforeCalyg1 = (1.0/ndraws)*yforeCalyg1;       % mean
%
%@@@ correlation of inflation and U
%IUcov = zeros(forepq,1);   % forepq quarters
%IUbeta = zeros(forepq,2);   % forepq quarters
%for i = 1:forepq
%   indI = 4*forepq+i;
%   indU = 5*forepq+i;
%   %junk = corrcoef(yforeqgw(:,indI),yforeqgw(:,indU));
%   %IUcov(i,1) = junk(1,2);
%   xjnk = [100*yforeqgw(:,indU) ones(ndraws,1)];
%   yjnk = yforeqgw(:,indI);
%   junk = (xjnk'*xjnk)\(xjnk'*yjnk);
%   IUbeta(i,:) = junk';
%end


%%
% *** .68 and .90 probability bands of out-of-sample forecasts
%yforel = zeros(forep*nvar,1);     % preallocating
%yforeh = zeros(forep*nvar,1);     % preallocating
%yforeqgl = zeros(forepq*nvar,1);     % preallocating
%yforeqgh = zeros(forepq*nvar,1);     % preallocating

%clear yfores yforepgs yforeCalygs
%*** write out final results
%yforeml = reshape(yforeml,forep,nvar);
%yforeqgml = reshape(yforeqgml,forepq,nvar);
%yforeCalygml = reshape(yforeCalygml,forepy,nvar);