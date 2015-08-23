function fn_mseriesgraph(ydate,keyindx,rnum,cnum,q_m,xlab,ylab,tlab)
% fn_mseriesgraph(ydate,keyindx,rnum,cnum,q_m,xlab,ylab,tlab)
%   Graph actual series or point forecasts or both (annual or monthly or quarterly)
%     with multiple series in one subplot.
%
% ydate:  T-by-n+2-by-h series data with dates in the first 2 columns in the order of year and month (or quarter).
%           T: the length of an individual series.
%           n+2:  n series plus the first 2 columns indicating dates.
%           h:    the number of sereies on the 3rd dimension, put in the same subplot as in each of the n series.
%                  Thus, the first 2 columns in the 3rd dimension can be NaN while
%                  some elements are also allowed to be NaN if no data are available.
% keyindx:  index for the variables to be graphed
% rnum:  number of rows in subplot
% cnum:  number of columns in subplot
% q_m:  if 4 or 12, quarterly or monthly data
% xlab:  1-by-1 x-axis label for what conditions imposed; e.g., conlab = 'All bar MS shocks inspl'
% ylab:  string array for the length(keyindx)-by-1 variables
% tlab:  1-by-1 title label for (e.g., as of time of forecast)
%-------------
% No output argument for this graph file
%  See fn_foregraph.m, fn_forerrgraph.m.
%
% Tao Zha, September 2000
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


vyrs = ydate(:,1);     % vector column
hornum = cell(length(vyrs),1);    % horizontal year (number)
count=0;
for k=vyrs'
   count=count+1;
   jnk=num2str(k);
   if length(jnk)==4   % years
      hornum{count}=jnk(3:4);   % e.g., with '1990', we have '90'
   else
      hornum{count}=jnk;   % e.g., with '1', '2', ...
   end
end

if rnum*cnum<length(keyindx)
   disp(' ')
   warning('The total number of subplots must be greater that the number of selected series for graphing!')
   disp('Press ctrl-c to abort')
   pause
end

count=0;
for i = keyindx
   count = count+1;
   subplot(rnum,cnum,count)
   plot(ydate(:,1)+ydate(:,2)/q_m,squeeze(ydate(:,2+i,:)),...
        ydate(:,1)+ydate(:,2)/q_m,zeros(length(vyrs),1),'-')

   if (ydate(1,2)==0) & (length(num2str(ydate(1,1)))==4)      % only for annual growth rates (not for, say, monthly annualized rates)
      set(gca,'XLim',[vyrs(1) vyrs(end)])
      set(gca,'XTick',vyrs)
      set(gca,'XTickLabel',char(hornum))
   end

   if i==keyindx(1)
      title(tlab)
   elseif i>=length(keyindx)   %i>=length(keyindx)-1
      xlabel(xlab)
   end
   %
   grid
   ylabel(char(ylab(i)))
end
