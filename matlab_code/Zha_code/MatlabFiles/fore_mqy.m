function [yforelml,yforemgml,yforeqgml,yforeCalygml,yactl,yactmg,yactqg,yactCalyg,yactlog] = ...
                      fore_mqy(yfore,xdata,nvar,nSample,nSampleCal,...
                      forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,vlistper,...
                      lmqyIndx,indxCal)
% [yforelml,yforemgml,yforeqgml,yforeCalygml,yactl,yactmg,yactqg,yactCalyg,yactlog] = ...
%                     fore_mqy(yfore,xdata,nvar,nSample,nSampleCal,...
%                     forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,vlistper,...
%                     lmqyIndx,indxCal)
%
% Converting oringinal forecast "yfore" to series by calendar years and series
%            with growth rates (annualized for monthly and quarterly series) and
%            with corresponding actual data (ONLY valid when we have actual data.
%            Later should add an index to allow an option out 3/23/99).
% 3/23/99 I guess that this program is more general than fore_cal because it allows
%    for non-Calendar annual growth rate calculation if nSample before "forep" and
%    and q_m before "vlist" are used.
%
% yfore:   oringal forecast series, all logged except R, U, etc.
% xdata:   oringal data set, up to the period before (peudo-) out-of-sample forecasting,
%                 all logged except R, U, etc.
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
% indxCal:  1: calendar years; 0: (NOT calendar) years
%----------------
% yforelml:  forep-by-nvar, ymonthly ML forecast, all at level, in percent
% yforemgml:  ML forecasts, monthly growth (at annual rates), in percent
% yforeqgml:  ML forecast: quarterly growth (at annual rates), in percent
% yforeCalygml:  forepq-by-nvar ML forecast: annual growth, in percent
% yactl:   forep-by-nvar, actual, all at level, in percent
% yactmg:
% yactqg:
% yactCalyg: forepq-by-nvar, actual, annual growth, in percent
% yactlog:  forep-by-nvar, actual, all log except R and U.
%
% Copyright (c) March 1998 by Tao Zha
% Revised, 3/23/99.  Added "lmqyIndx" so that previous programs may not be compatible.
% Revised, 4/27/99.  Added "indxCal" so that previous programs may not be compatible.

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
%

if length(lmqyIndx)~=4
	warning('lmqyIndx must be a 4-by-1 vector containing 1 or 0')
	return
end



%---------------------
%  Actual data
%---------------------
yact = xdata(nSample-q_m+1:nSample,:);   % the latest (not calender) year
yactCal = xdata(nSampleCal-q_m+1:nSampleCal,:);
        % the lastest calendar year (always earlier than the lastest (not calendar) year



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
   yforemg(:,vlistlog) = ( yforem(2:forep+1,vlistlog) - yforem(1:forep,vlistlog) ) .* q_m;
                        % monthly change (annualized), 12*log(1+growth rate)
   yforemgml=yforemg;
   yforemgml(:,vlistlog) = 100*(exp(yforemg(:,vlistlog))-1);  % monthly growth rate (annualized)
   yforemgml(:,vlistper) = 100*yforemg(:,vlistper);   % monthly growth rate (annualized)
else
	yforemgml=NaN;
end



%-----------------------------------------
%   Converted to quarterly
%-----------------------------------------
if lmqyIndx(3)
   yactQ = xdata(nSample-mod(qmEnd,3)-q_m+1:nSample-mod(qmEnd,3),:);   % the latest actual quarter
   yactq = zeros(4,length(vlist));  % the latest 4 quarters prior to forecasting
   yforeq = zeros(4+forepq,length(vlist));   % including the 4 quarters prior to forecasting
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
   yforeqg(:,vlistlog) = ( yforeq(5:forepq+4,vlistlog) - yforeq(4:forepq+3,vlistlog) ) .* 4;
                              % prior quarter, 4*log(1+growth rate)
   yforeqgml=yforeqg;
   yforeqgml(:,vlistlog) = 100*(exp(yforeqg(:,vlistlog))-1);  % quarterly growth rate (annualized)
   yforeqgml(:,vlistper) = 100*yforeqg(:,vlistper);    % quarterly growth rate (annualized)
else
	yforeqgml = NaN;
end


%-----------------------------------------
%   Converted to annual years (may not be calendar)
%-----------------------------------------
if lmqyIndx(4)
   yactCaly = zeros(1,length(vlist));  % the latest calendar year
   yforeCaly = zeros(1+forepy,length(vlist)); % including the calendar year prior to forecasting
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
   yforeCalyg(:,vlistlog) = yforeCaly(2:forepy+1,vlistlog) - yforeCaly(1:forepy,vlistlog);
                              % year-over-year, log(1+growth rate)
   yforeCalygml=yforeCalyg;
   yforeCalygml(:,vlistlog) = 100*(exp(yforeCalyg(:,vlistlog))-1);   % annaul growth rate
   yforeCalygml(:,vlistper) = 100*yforeCalyg(:,vlistper);   % annaul growth rate
else
	yforeCalygml = NaN;
end



%----------------------------------------------------------------------
% --------------------- Actual -----------------------------------------
%----------------------------------------------------------------------
if indxCal
   Rem_qm = q_m-qmEnd;   % remaining months within the first calendar year
   Rem_y = floor((forep-Rem_qm)/q_m);   % remaining calendar years after the first one
   yact = xdata(nSample-2*q_m+1:nSample+Rem_qm+Rem_y*q_m,:);
            % 2*q_m (not calendar) months + months ended at last calendar year
else
   yact = xdata(nSample-q_m+1:nSample+forep,:);   % q_m (not calendar) months + forep
end
mT1 = length(yact(:,1));
yactlog = yact(q_m+1:mT1,:);


%-----------------------------------------
%   Converted to monthly level
%-----------------------------------------
if lmqyIndx(1)
   yactl=yactlog;
   yactl(:,vlistlog)=exp(yactlog(:,vlistlog));       % mode, all levels
   yactl(:,vlistper)=100*yactlog(:,vlistper);       % mode, all levels
else
   yactl=NaN;
end

%-----------------------------------------
%   Converted to monthly growth
%-----------------------------------------
if lmqyIndx(2)
   yactmg = yactl;
   yactmg(:,vlistlog) = (yact(q_m+1:mT1,vlistlog) - yact(q_m:mT1-1,vlistlog)) .* 12;
                              % monthly, 12*log(1+growth rate), annualized growth rate
   yactmg(:,vlistlog) = 100*(exp(yactmg(:,vlistlog))-1);
   yactmg(:,vlistper) = 100*yactmg(:,vlistper);
else
   yactmg=NaN;
end

%-----------------------------------------
%   Converted to quarterly
%-----------------------------------------
if lmqyIndx(3)
   warning(' ')
   disp('Not worked out yet.  4/27/99')
   disp('Press ctrl-c to abort!')
   pause

   yactQ = yact(10:mT1,:);      % October: beginning of the 4th quarter
   yactq = zeros(1+forepq,length(vlist));
   qT1 = 1+forepq;
   %
   for i = 1:qT1
      i1 = 1+3*(i-1);
      i2 = 3*i;
      yactq(i,:) = sum(yactQ(i1:i2,:)) ./ 3;
   end
   %
   yactqg = yactq(2:qT1,:);    % 4+1=5 where 4 means 4 quarters
   yactqg(:,vlistlog) = (yactq(2:qT1,vlistlog) - yactq(1:qT1-1,vlistlog)) .* 4;
                                 % quarterly, 4*log(1+growth rate), annualized growth rate
   yactqg(:,vlistlog) = 100*(exp(yactqg(:,vlistlog))-1);
   yactqg(:,vlistper) = 100*yactqg(:,vlistper);
else
   yactqg=NaN;
end


%-----------------------------------------
%   Converted to annual years (may not be calendar)
%-----------------------------------------
if lmqyIndx(4)
   if indxCal
      yT1 = 1+forepy;
      yactCaly = zeros(yT1,length(vlist));  % the latest 1+forepy calendar years
      for i = 1:yT1
         i1 = ( Rem_qm+q_m*(qmEnd==q_m) ) + 1+q_m*(i-1);
         i2 = ( Rem_qm+q_m*(qmEnd==q_m) ) + q_m*i;
         yactCaly(i,:) = sum(yact(i1:i2,:)) ./ q_m;
      end
      %
      yactCalyg = yactCaly(2:yT1,:);    % 1+1=2 where 1 means 1 year
      yactCalyg(:,vlistlog) = (yactCaly(2:yT1,vlistlog) - yactCaly(1:yT1-1,vlistlog));
                              % annual rate: log(1+growth rate)
      yactCalyg(:,vlistlog) = 100*(exp(yactCalyg(:,vlistlog))-1);
      yactCalyg(:,vlistper) = 100*yactCalyg(:,vlistper);
   else
      yT1 = 1+forepy;
      for i = 1:yT1
         i1 = 1+q_m*(i-1);
         i2 = q_m*i;
         yactCaly(i,:) = sum(yact(i1:i2,:)) ./ q_m;
      end
      %
      yactCalyg = yactCaly(2:yT1,:);    % 1+1=2 where 1 means 1 year
      yactCalyg(:,vlistlog) = (yactCaly(2:yT1,vlistlog) - yactCaly(1:yT1-1,vlistlog));
                              % annual rate: log(1+growth rate)
      yactCalyg(:,vlistlog) = 100*(exp(yactCalyg(:,vlistlog))-1);
      yactCalyg(:,vlistper) = 100*yactCalyg(:,vlistper);
   end
else
   yactCalyg=NaN;
end