function [HDratio,HDvalue,yhatn,yforen] = history2(HDindx,A0,Bhml,phil,nn,...
            actup,Estrexa,xdata,nvar,nSample,nSampleCal,forep,forepq,forepy,q_m,...
					 qmEnd,vlist,vlistlog,vlistper,lmqyIndx)
%  [HDratio,HDvalue,yhatn,yforen] = history2(HDindx,A0,Bhml,phil,nn,actup,...
%                xdata,nvar,nSample,nSampleCal,forep,forepq,forepy,q_m,...
%					 qmEnd,vlist,vlistlog,vlistper,lmqyIndx)
%
%  (Out-of-sample) historical decompostions: alternative approach to history.m in 3 aspects
%         (1) compute the percentage (not the value itself)
%         (2) NOT cumulative HD (THIS is NOT correct, TAZ, 03/17/99)
%         (3) conditional on time t (not time 0).
%         (4) condtional and uncondtional forecasts
%
% HDindx:  k-by-1 cell where k groups of shocks; index number for particular shocks
% A0h:  nvar-by-nvar for y(t)A0 = X*A+ + E, where column means equation
% Bh: the (posterior) estimate of B for y(t) = X*Bh + E*inv(A0)
% phil: the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% nn: [nvar,lags,forep], forep: forecast periods (monthly)
% actup:  actual data periods (monthly) before the beginning of forecasts
% xdata:   oringal data set, up to the period before (peudo-) out-of-sample forecasting,
%                 all logged except R, U, etc.
% Estrexa:  forep-by-nvar -- backed out structural shocks for out-of-sample forecast
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
% lmyqIndx:  4-by-1 1 or 0's. Index for m log, m level, qg, and yg; 1: yes; 0: no;
%             if lmyqIndx(1)==1, both monthly log and monthly level are returned
%--------------
% HDratio:  5-by-k cells, where k groups (see HDindx) of shocks associated with decomposition;
%           5: monthly log, monthly level, mg, qg, and calendar yg in this order;
%           each cell is forep(q)(y)-by-nvar
% HDvalue:  same dimension as HDratio but with the values (differences between
%                 conditional and unconditional forecasts
% yhatn:  same dimension as HDratio but with conditional forecasts
% yforen: 5-by-1 cells and each cell is forep(q)(y)-by-nvar unconditional forecasts
%
% October 1998 by Tao Zha.
% Last change 3/19/99 on the dimension of "Estrexa" so that previous programs may not be
%        compatible.
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

nyout = 1+nargout('fore_cal');
              % +1 because we want the original logged variables as well.
[nHD,jnk] = size(HDindx);
yhat=cell(nHD,1);
yhatn = cell(nyout,nHD);
         % row: output such as l, mg, qg, or yg; column: how many decomps



%--------------------------------------------------------------------------
%  Unconditional point forecasts for log, monthly levle, mg, qg, and yg.
%--------------------------------------------------------------------------
yforeml = forecast(Bhml,phil,nn);
yforen = cell(nyout,1);
yforen{1} = yforeml;
oputstr = '[';   % output string
for n=2:nyout
	if n<nyout
		oputstr = [oputstr 'yforen{' num2str(n) '},'];
	else
		oputstr = [oputstr 'yforen{' num2str(n) '}'];
	end
end
oputstr =[oputstr ']'];
%
eval([oputstr ' = fore_cal(yforeml,xdata,nvar,nSample,nSampleCal,forep,'...
                'forepq,forepy,q_m,qmEnd,vlist,vlistlog,vlistper,lmqyIndx);']);



%--------------------------------------------------------------------------
%  Conditional forecasts on specified paths of shocks
%      for log, monthly levle, mg, qg, and yg.
%--------------------------------------------------------------------------
for k=1:nHD
	Estr = zeros(forep,nvar);     % out of sample
   Estr(:,HDindx{k}) = Estrexa(:,HDindx{k});    % out of sample, MS
	yhat{k} = forefixe(A0,Bhml,phil,nn,Estr);
	yhatn{1,k} = yhat{k};         % original logged variables
	%
	oputstr = '[';   % output string
	for n=2:nyout
		if n<nyout
			oputstr = [oputstr 'yhatn{' num2str(n) ',k},'];
		else
			oputstr = [oputstr 'yhatn{' num2str(n) ',k}'];
		end
	end
	oputstr =[oputstr ']'];
	%
	eval([oputstr ' = fore_cal(yhat{k},xdata,nvar,nSample,nSampleCal,forep,'...
	                'forepq,forepy,q_m,qmEnd,vlist,vlistlog,vlistper,lmqyIndx);']);
end



%--------------------------------------------------------------------------
%  Historical decompositions: both values and raitos (%)
%--------------------------------------------------------------------------
HDratio=cell(nyout,nHD);
HDvalue=HDratio;
for n=2:nyout
	if lmqyIndx(n-1)
		if (n-1==1)
	  		yhatotal = zeros(size(yhatn{1,k}));
	  		for k=1:nHD
	  			HDvalue{1,k} = yhatn{1,k}-yforen{1};
	  			yhatotal = yhatotal + abs(HDvalue{1,k});
	           		     % a sum of absolute values.  One can also use square
	  		end
	  		%
	  		for k=1:nHD
	  			HDratio{1,k} = abs(HDvalue{1,k})*100 ./ yhatotal;
	  		end
		end
		%
		yhatotal = zeros(size(yhatn{n,k}));
		for k=1:nHD
			HDvalue{n,k} = yhatn{n,k}-yforen{n};
			yhatotal = yhatotal + abs(HDvalue{n,k});
         		     % a sum of absolute values.  One can also use square
		end
		%
		for k=1:nHD
			HDratio{n,k} = abs(HDvalue{n,k})*100 ./ yhatotal;
		end
	end
end


%** check.  This first may not be zeros, but the second (monthly log) must be zeros.
%HDvalue{5,1}+HDvalue{5,2}-...
%      (yactCalyg(size(yactCalyg,1)-3:size(yactCalyg,1),:)-yforeCalygml)
%HDvalue{1,1}+HDvalue{1,2}-...
%      (yact(size(yact,1)-forep+1:size(yact,1),:)-yforen{1,1})