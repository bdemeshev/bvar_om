function gforerr1(yrEnd,qmEnd,q_m,forepy,yhatCalygml,yhatCalygl1,...
              yhatCalygh1,xlab,ylab,forelabel,conlab,keyindx,rnum,cnum)

%
%   One (.68) error band: (calendar year) forecasts only (no actual)
%
% yrEnd:  the end year for the sample
% qmEnd:  last q_m before out-of-sample forecasting
% q_m:   quarterly or monthly for the underlying model
% yearsCal: calendar years including actual and forecast data
% yactCalyg:  actual calendar annual growth data
% yhatCalygml:  conditional (calendar annual) forecasts on particular shocks
% yhatCalygl1:  0.68 lower band
% yhatCalygh1:  0.68 upper band
% forepy:  forecast periods (yearly)
% xlab:  label for structural equations
% ylab:  label for the variables
% forelabel:  title label for as of time of forecast
% conlab:  label for what conditions imposed; e.g., conlab = 'All bar MS shocks inspl'
% keyindx:  index for the variables to be graphed
% rnum:  number of rows in subplot
% cnum:  number of columns in subplot
%-------------
% No output argument for this graph file
%
% 03/19/99

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


%*** Begining period
fbegm = qmEnd+1;   % +1 because the first month out of sample
if fbegm > q_m    % the first forecast month into the next year
	fbegm = 1;
	fbegy = yrEnd+1;   % the first forecast year for the graph
else
	fbegy = yrEnd;    % the first forecast year for the graph
end
%
%  abegy = min(yearsCal);   % the actual beginning year for the graph
%  abegm = 1;
%  ifbegy = find(yearsCal==fbegy-1);   % index for the year before the first
%                 % forecast first year in actual data "yactCalyg"

%*** Final period
ffiny = fbegy+forepy-1;   % the forecast final year
%afiny = max(yearsCal);    % the actual final year



figure

nvar = size(yhatCalygml,2);      %  Pcm,M2,FFR,GdP, CPI, U
%keyindx = 1:nvar;      %  Pcm,M2,FFR,GdP, CPI, U
%keyindx = [4:nvar 3 2];      %  GdP, CPI, U, FFR, M2

hornum = cell(length(fbegy:ffiny),1);    % horizontal number
count=0;
for k=fbegy:ffiny
	count=count+1;
	jnk=num2str(k);
	hornum{count}=jnk(3:4);   % e.g., with '1990', we have '90'
end

count=0;
%
for i = keyindx
   count = count+1;
   subplot(rnum,cnum,count)
	%
   plot(fbegy:ffiny,yhatCalygml(:,i),'-',...
           fbegy:ffiny,yhatCalygl1(:,i),'--',fbegy:ffiny,yhatCalygh1(:,i),'--' )
   set(gca,'XTick',fbegy:ffiny)
	set(gca,'XTickLabel',char(hornum))
   if i<=cnum
   	title(forelabel)
	elseif i>=nvar-1
		xlabel(conlab)
	end
	%
	grid
	ylabel(char(ylab(i)))
end