function scaleout = fn_mimfgraph(imfe,nvar,q_m,imstp,xlab,ylab,tlab,xTick)
%     mimfgraph:  multiple impulse functions plotted in subplots put in one chart.
%
% imfe:  imstp-by-n^2+2-by-h series data with dates in the first 2 columns in the order of year and month (or quarter).
%           imstp: the number of impulse response steps.
%           nvar^2+2:  n series plus the first 2 columns indicating dates.  For the last nvar^2 columns,
%               the order is: nvar responses to the 1st shock, ..., nvar responses to the last shock.
%           h: The number of sereies on the 3rd dimension, put in the same subplot as in each of the n series.
%              The first 2 columns in the 3rd dimension can be NaN while some other columns
%                  are also allowed to be NaN if no data are available.
%              The last 2 columns in the 3rd dimension must be used for error bands if "area"
%                  is used to plot these bands.
%  nvar: number of variables
%  q_m:  monthly (12) or quarterly (4)
%  imstp:  number of steps of impulse responses
%  xlab,ylab,tlab:   x-axis, y-axis, and title labels
%  xTick:  optional.  Eg: [12 24 36].
%---------------
%  scaleout: column 1 represents maximums; column 2 minimums.  Rows: nvar variables.
%
%  See ftd_mseriesgraph.m and fn_mseriesgraph.m
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


if nargin < 7, xTick = []; end

if nvar^2~=size(imfe,2)-2
   disp(' ')
   warning('The number of columns in imfe must be nvar^2+2 (columns of dates)')
   disp('Press ctrl-c to abort')
   pause
end



t = 1:imstp;
temp1=zeros(nvar,1);
temp2=zeros(nvar,1);
minval=zeros(nvar,1);   % for each variable to nvar shocks
maxval=zeros(nvar,1);   % for each variable to nvar shocks
for i = 1:nvar   % variable
   for j = 1:nvar   % shock
      tmpimf = squeeze(imfe(:,2+(j-1)*nvar+i,:));
      temp1(j) = min(min(tmpimf));
      temp2(j) = max(max(tmpimf));
	end
   minval(i)=min(temp1);
   maxval(i)=max(temp2);
end

scaleout = [minval(:) maxval(:)];

%--------------
%  Column j: Shocks 1 to N; Row i: Variable responses to
%-------------

rowlabel = 1;
for i = 1:nvar       % variable
   columnlabel = 1;

   if minval(i)<0
      if maxval(i)<=0
         yt=[minval(i) 0];
      else
         yt=[minval(i) 0 maxval(i)];
      end
   else % (minval(i) >=0)
      if maxval(i) > 0
         yt=[0 maxval(i)];
      else % (identically zero responses)
         yt=[-1 0 1];
      end
   end

   scale=[1 imstp minval(i) maxval(i)];
   for j = 1:nvar      % shock
      k1=(i-1)*nvar+j;
      k2=(j-1)*nvar+i;
      subplot(nvar,nvar,k1)
      if 0   % Plot area for error bands at the last two 3-rd dimensions.
         area(imfe(:,1)+imfe(:,2)/q_m,squeeze(imfe(:,2+k2,end)),-100,'EdgeColor','none','FaceColor','y')   % yellow
         hold on
         area(imfe(:,1)+imfe(:,2)/q_m,squeeze(imfe(:,2+k2,end-1)),-100,'EdgeColor','none','FaceColor',[1 1 1])   % white
         set(gca,'ColorOrder',[0 0 0]); % turn the color off and set it to black
         set(gca,'LineStyleOrder', '-|-.|--|:');  % cycle through the newly defined LineSytleOrder
         plot(t,squeeze(imfe(:,2+k2,1:end-2)),t,zeros(length(imfe(:,k2)),1),'-');
         set(gca,'Layer','top')
         hold off
      else   % No error bands plotted
         %=== set color to black and cycle through the newly defined LineSytleOrder
         set(0,'DefaultAxesColorOrder',[0 0 0], ...
                  'DefaultAxesLineStyleOrder','-|-.|--|:')
         %set(gca,'ColorOrder',[0 0 0]); % turn the color off and set it to black
         %set(gca,'LineStyleOrder', '-|-.|--|:');  % cycle through the newly defined LineSytleOrder
         plot(t,squeeze(imfe(:,2+k2,1:end)),t,zeros(length(imfe(:,k2)),1),'-');
      end
      if 1 % Get legends
         if (i==2) & (j==3)  % | (i==6)
            legend('Const','S1','S2','S3',0)
            %  legend('Actual Data','Absent policy shocks',0)
         end
      end

      set(gca,'XTick',xTick)
      set(gca,'YTick',yt)
      grid

      axis(scale);   % put limits on both axes.
      %set(gca,'YLim',[minval(i) maxval(i)])   % put the limit only on the y-axis
      if isempty(xTick)  %1     % no numbers on axes
         set(gca,'XTickLabel',' ');
         set(gca,'YTickLabel',' ');
      else   % put numbers on both axes
         if i<nvar
            set(gca,'XTickLabelMode','manual','XTickLabel',[])
         end
         if j>1
            set(gca,'YTickLabel',' ');
         end
      end

      if rowlabel == 1
         %title(['x' num2str(j)])
         %title(eval(['x' num2str(j)]))
         title(char(xlab(j)))
      end
      if columnlabel == 1
         %ylabel(['x' num2str(i)])
         %ylabel(eval(['x' num2str(i)]))
         ylabel(char(ylab(i)))
      end
      columnlabel = 0;
   end
   rowlabel = 0;
end

subtitle(tlab)
