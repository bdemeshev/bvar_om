function [yactCalyg,yforeCalygml,yAg,yFg,yactqg,yforeqgml,qAg,qFg,...
                      yactmg,yforemgml,mAg,mFg,yact,yactCal] = fore_gh(xinput)
% Plot graphics
%    [yactCalyg,yforeCalygml,yAg,yFg,yactqg,yforeqgml,qAg,qFg,...
%                      yactmg,yforemgml,mAg,mFg,yact,yactCal] = fore_gh(xinput)
%
% xdata=xinput{1}: data, all logged except R, U, etc.
% nvar=xinput{2}: number of variables
% nSample=xinput{3}: sample size (including lags or initial periods)
% nSampleCal=xinput{4}: converting nSample of calendar years;
% yforelml=xinput{5}: monthly ML forecast, all in level, in percent
% yforemgml=xinput{6}: monthly growth -- annualized rate, in percent
% yforeqgml=xinput{7}: prior-quarter growth -- annualized rate, in percent
% yforeCalygml=xinput{8}: annual rate, expressed in percentage point, in percent
% actup=xinput{9}: periods for actual data (monthly)
% actupq=xinput{10}: periods for actual data (quarterly)
% actupy=xinput{11}: periods for actual data (calendar yearly)
% vlist=xinput{12}: list of variables
% vlistlog=xinput{13}: sub list of variables that are in log
% vlistper=xinput{14}: sub list of variables that are in percent
% q_m=xinput{15}: month or quarter in the model
% forep=xinput{16}: forecast periods (monthly)
% ylab=xinput{17}: labels for the y-axis
% forepq=xinput{18}: quarters in the forecast period
% forepy=xinput{19}: calendar years in the forecast period
% Psuedo=xinput{20}: 1: Psuedo out-of-sample; 0: real time out-of-sample
% Graphfore=xinput{21}: 1: graphics in this file; 0: no graphics in the execution.
% yearsCal=xinput{22}:  calendar years matching output "yactCalyg"
% qmEnd=xinput{23}:  last month (or quarter) of the sample
%-------------
% yactCalyg: annual growth for actual data
% yforeCalygml: annaul growth for ML forecasts, expressed in percent
% yAg: length of yactCalyg (actual data)
% yFg: lenght of yforeCalyg (forecasts)
%
% yactqg: prior-quarter annualized growth for actual data
% yforeqgml: prior-quarter annaulized growth for ML forecasts, expressed in percent
% qAg: length of yactqg (actual data)
% qFg: lenght of yforeqg (forecasts)
%
% yactmg: month-to-month annualized growth for actual data
% yforemgml: month-to-month annaulized growth for ML forecasts, expressed in percent
% mAg: length of yactqg (actual data)
% mFg: lenght of yforeqg (forecasts)
%
% yact:  log(y) except R, ect. (monthly).  If Psuedo, yact includes "forep" months.
%             i.e., (actup+forep)-by-nvar
% yactCal:  same as yact but ends at the end of the last calendar year.  If Psuedo,
%     (actup+forep)-by-nvar, which include some actual data in the 1st calendar year.
%
% Copyright (c) March 1998 by Tao Zha
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


xdata=xinput{1}; nvar=xinput{2}; nSample=xinput{3}; nSampleCal=xinput{4};
yforelml=xinput{5}; yforemgml=xinput{6}; yforeqgml=xinput{7}; yforeCalygml=xinput{8};
actup=xinput{9}; actupq=xinput{10}; actupy=xinput{11}; vlist=xinput{12};
vlistlog=xinput{13}; vlistper=xinput{14}; q_m=xinput{15}; forep=xinput{16};
ylab=xinput{17}; forepq=xinput{18}; forepy=xinput{19}; Psuedo=xinput{20};
Graphfore=xinput{21}; yearsCal=xinput{22}; qmEnd=xinput{23};

% ========= plot the graphics with actural vs forecast vs Full Bayesian =========
%
%*** actual data
if Psuedo
   yact = xdata(nSample-actup+1:nSample+forep,:);   % actup (not calendar) months + forep
   yactCal = xdata(nSampleCal-actupy*q_m+1:nSampleCal+forepy*q_m,:);  % calendar years
else
   yact = xdata(nSample-actup+1:nSample,:);    % actup (not calendar) months
   yactCal = xdata(nSampleCal-actupy*q_m+1:nSampleCal,:);   % caledar years
end
%



%%%%---------------------------------------------------------------
%$$$ Actual monthly growth rate, Prior quarter, and year-over-year change
%%%%---------------------------------------------------------------

%
%@@@ Monthly change (annaluized rate)
%
yactm = yact;
mT1 = length(yact(:,1));
%
yactmg = yactm(2:mT1,:);    % start at second month to get growth rate
yactmg(:,vlistlog) = (yactm(2:mT1,vlistlog) - yactm(1:mT1-1,vlistlog)) .* 12;
                             % monthly, 12*log(1+growth rate), annualized growth rate
yactmg(:,vlistlog) = 100*(exp(yactmg(:,vlistlog))-1);
yactmg(:,vlistper) = 100*yactmg(:,vlistper);
mAg = length(yactmg(:,1));


%
%@@@ Prior Quarter change (annaluized rate)
%
if Psuedo
	yactQ = xdata(nSample-mod(qmEnd,3)-actupq*3+1:nSample-mod(qmEnd,3)+forepq*3,:);
   yactq = zeros(actupq+forepq,length(vlist));
else
	yactQ = xdata(nSample-mod(qmEnd,3)-actupq*3+1:nSample-mod(qmEnd,3),:);
   yactq = zeros(actupq,length(vlist));
end
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
qAg = length(yactqg(:,1));


%
%@@@ Calendar year-over-year change
%
if Psuedo
   yactCaly = zeros(actupy+forepy,length(vlist));  % past "actupy" calendar years + forepy
else
   yactCaly = zeros(actupy,length(vlist));  % past "actupy" calendar years
end
yT1 = length(yactCal(:,1))/q_m;
yT1a = length(yactCaly(:,1));
if yT1 ~= yT1a
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
yAg = length(yactCalyg(:,1));


%----------------------------------------------------------------------
%=================  Graphics of actual vs forecast ===================
%----------------------------------------------------------------------
yactl(:,vlistlog)=exp(yact(:,vlistlog));
yactl(:,vlistper)=100*yact(:,vlistper);  % convert to levels, in percent
%
mFg = length(yforemgml(:,1));
qFg = length(yforeqgml(:,1));
yFg = length(yforeCalygml(:,1));

if Psuedo
   t1=1:actup+forep;
else
   t1=1:actup;
end
t2=actup+1:actup+forep;
%
if Graphfore
	figure           % The graphics below are expressed in level (monthly)
	for i = 1:nvar
	   subplot(nvar/2,2,i)
   	plot(t1,yactl(:,i),t2,yforelml(:,i),'--')
   	%title('solid-actual, dotted-forecast');
   	%title(eval(['forelabel']));
   	%ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
	end
end

t1 = 1:mAg;
if Psuedo
   t2 = mAg-mFg+1:mAg;
else
   t2 = mAg+1:mAg+mFg;
end
%
if Graphfore
	figure   % The graphics below are monthly growth rates except R and U, at annualized rates
	for i = 1:nvar
	   subplot(nvar/2,2,i)
	   plot(t1,yactmg(:,i),t2,yforemgml(:,i),'--')
	   %title('solid-actual, dotted-forecast');
	   %xlabel(eval(['forelabel']));
	   %ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
	end
end

t1 = 1:qAg;
if Psuedo
   t2 = qAg-qFg+1:qAg;
else
   t2 = qAg+1:qAg+qFg;
end
%
if Graphfore
	figure  % The graphics below are prior quarter growth rate except R and U, at annualized rates
	for i = 1:nvar
	   subplot(nvar/2,2,i)
	   plot(t1,yactqg(:,i),t2,yforeqgml(:,i),'--')
	   %title('quarterly growth, unconditional');
	   %xlabel(eval(['forelabel']));
   	%ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
	end
end


%t1 = 1:yAg;
t1 = yearsCal;
yAgCal=yearsCal(length(yearsCal));
if Psuedo
   t2 = yAgCal-yFg+1:yAgCal;
else
   t2 = yAgCal+1:yAgCal+yFg;
end
%
if Graphfore
	figure
	% The graphics below are over-the-year growth rate except R and U, annually
	for i = 1:nvar
	   subplot(nvar/2,2,i)
   	plot(t1,yactCalyg(:,i),t2,yforeCalygml(:,i),'--')
	   %title('solid-actual, dotted-forecast');
	   %xlabel(eval(['forelabel']));
	   %ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
	end
end