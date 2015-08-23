function [yact,yactmg,yactqg,yactCalyg,yactCal] = datactcon(xinput)
% [yact,yactmg,yactqg,yactCalyg,yactCal] = datactcon(xinput)
%   Data-actual-conditions.  Take the raw data xdata to mlog, mg, qg, and yg.
%   Note yact(mlog), yactqg, and yactCalyg corr. to myears, qyears, and yearsCal
%
% xdata=xinput{1}: data, all logged except R, U, etc.
% nSample=xinput{2}: sample size (including lags or initial periods)
% nSampleCal=xinput{3}: converting nSample to sample marked by calendar years;
% actup=xinput{4}: periods for actual data (monthly)
% actupq=xinput{5}: periods for actual data (quarterly)
% actupy=xinput{6}: periods for actual data (calendar yearly)
% vlist=xinput{7}: list of variables
% vlistlog=xinput{8}: sub list of variables that are in log
% vlistper=xinput{9}: sub list of variables that are in percent
% q_m=xinput{10}: month or quarter in the model
% forep=xinput{11}: forecast periods (monthly)
% forepq=xinput{12}: quarters in the forecast period
% forepy=xinput{13}: calendar years in the forecast period
% Psuedo=xinput{14}: 1: Psuedo out-of-sample; 0: real time out-of-sample
% qmEnd=xinput{15}:  last month (or quarter) of the sample
%-------------
% yact:  (actup+forep*Psuedo)-by-nvar; log(y) except R, ect. (monthly).
%        Begin: nSample-actup+1; End: nSample+forep*Psuedo.  Match dates myears
% yactmg: (size(yact,1)-1)-by-nvar; month-to-month annualized growth for actual data.
%        Begin: nSample-actup+2; End: nSample+forep*Psuedo.  Match dates "myears(2:end)"
% yactqg: prior-quarter annualized growth for actual data
%        Match dates "qyears."
% yactCalyg: annual growth for actual data
%        Match dates "yearsCal."
% yactCal:  weirdo -- seldom used.  Same as yact but ends at the end of the
%           last calendar year.  If Psuedo, it includes some actual data in
%           the 1st forecast calendar year.  I haven't found the use of this weirdo.
%
% Copyright (c) March 1998 by Tao Zha
% Revised, October 1998
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


xdata=xinput{1}; nSample=xinput{2}; nSampleCal=xinput{3};
actup=xinput{4}; actupq=xinput{5}; actupy=xinput{6}; vlist=xinput{7};
vlistlog=xinput{8}; vlistper=xinput{9}; q_m=xinput{10}; forep=xinput{11};
forepq=xinput{12}; forepy=xinput{13}; Psuedo=xinput{14}; qmEnd=xinput{15};

%--------------------------------------
%  Monthly log, acutal data
%  yact the way it is and yactCal ends at the end of calendar year
%----------------------------------------
yact = xdata(nSample-actup+1:nSample+forep*Psuedo,:);
	                         % actup (not calendar) months + forep
yactCal = xdata(nSampleCal-actupy*q_m+1:nSampleCal+forepy*q_m*Psuedo,:);
	                         % calendar years




%---------------------------------------------------------------
% Actual monthly growth rate, Prior quarter, and calendar yearly change
%---------------------------------------------------------------
%
%@@@ Monthly change (annaluized rate)
%
yactm = yact;
mT1 = length(yact(:,1));
%
yactmg = yactm(2:mT1,:);    % start at second month to get growth rate
yactmg(:,vlistlog) = (yactm(2:mT1,vlistlog) - yactm(1:mT1-1,vlistlog)) .* q_m;
                             % monthly, 12*log(1+growth rate), annualized growth rate
yactmg(:,vlistlog) = 100*(exp(yactmg(:,vlistlog))-1);
yactmg(:,vlistper) = 100*yactmg(:,vlistper);

%
%@@@ Prior Quarter change (annaluized rate)
%
yactQ = xdata(nSample-mod(qmEnd,3)-actupq*3+1:nSample-mod(qmEnd,3)+forepq*3*Psuedo,:);
yactq = zeros(actupq+forepq*Psuedo,length(vlist));
%
qT1 = length(yactQ(:,1))/3;
qT1a = length(yactq(:,1));
if qT1 ~= qT1a
	warning('line #')
	error('qT1: Beginings or ends of monthly and quarterly series do not match!')
end
%
for i = 1:qT1
   i1 = 1+3*(i-1);
   i2 = 3*i;
   yactq(i,:) = sum(yactQ(i1:i2,:)) ./ 3;
end
%
yactqg = yactq(5:qT1,:);    % 4+1=5 where 4 means 4 quarters
yactqg(:,vlistlog) = (yactq(5:qT1,vlistlog) - yactq(4:qT1-1,vlistlog)) .* 4;
                              % quarterly, 4*log(1+growth rate), annualized growth rate
yactqg(:,vlistlog) = 100*(exp(yactqg(:,vlistlog))-1);
yactqg(:,vlistper) = 100*yactqg(:,vlistper);

%
%@@@ Calendar year-over-year change
%
yactCaly = zeros(actupy+forepy*Psuedo,length(vlist));
                       % past "actupy" calendar years + forepy*Psuedo
yT1 = length(yactCal(:,1))/q_m;
yT1a = length(yactCaly(:,1));
if yT1 ~= yT1a
	warning('Line #')
	error('yT1: Beginings or ends of monthly and quarterly series are not the same!')
end
%
for i = 1:yT1
   i1 = 1+q_m*(i-1);
   i2 = q_m*i;
   yactCaly(i,:) = sum(yactCal(i1:i2,:)) ./ q_m;
end
%
yactCalyg = yactCaly(2:yT1,:);    % 1+1=2 where 1 means 1 year
yactCalyg(:,vlistlog) = (yactCaly(2:yT1,vlistlog) - yactCaly(1:yT1-1,vlistlog));
                          % annual rate: log(1+growth rate)
yactCalyg(:,vlistlog) = 100*(exp(yactCalyg(:,vlistlog))-1);
yactCalyg(:,vlistper) = 100*yactCalyg(:,vlistper);
