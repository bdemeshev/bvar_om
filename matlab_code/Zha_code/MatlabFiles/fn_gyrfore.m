function fn_gyrfore(yfore,yacte,keyindx,rnum,cnum,ylab,forelabel,conlab)
%
%   Graph annual (or calendar year) forecasts vs actual (all data from "msstart.m")
%
% yfore:  actual and forecast annual growth data with dates.
% yacte:  actual annual growth data with dates.
% keyindx:  index for the variables to be graphed
% rnum:  number of rows in subplot
% cnum:  number of columns in subplot
% ylab:  label for the variables
% forelabel:  title label for as of time of forecast
% conlab:  label for what conditions imposed; e.g., conlab = 'All bar MS shocks inspl'
%-------------
% No output argument for this graph file
%
% Tao Zha, March 2000
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


vyrs = yfore(:,1);
hornum = cell(length(vyrs),1);    % horizontal year (number)
count=0;
for k=vyrs'
   count=count+1;
   jnk=num2str(k);
   hornum{count}=jnk(3:4);   % e.g., with '1990', we have '90'
end

count=0;
for i = keyindx
   count = count+1;
   subplot(rnum,cnum,count)
   %
   plot(yacte(:,1),yacte(:,2+i),vyrs,yfore(:,2+i),'--')
   set(gca,'XLim',[vyrs(1) vyrs(end)])
   set(gca,'XTick',vyrs)
   set(gca,'XTickLabel',char(hornum))
   if i==keyindx(1)
      title(forelabel)
   elseif i>=length(keyindx)  %i>=length(keyindx)-1
      xlabel(conlab)
   end
   %
   grid
   ylabel(char(ylab(i)))
end
