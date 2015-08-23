function scaleout = fn_multigraphn_3g_shadedbands(imfn,...
                         xlab,ylab,tstring,XTick,YTickIndx,scaleIndx,nrowg,ncolg)
%Searching <<>> for ad hoc and specific changes.
%
%Stacking n sets of impulse responses in one graph.  See fn_multigraph2_ver2.m.
%   imfn: row: "nstp" time horizon (in the graphics),
%         column: "nrow "variables such as responses (row in the graphics),
%         3rd D: across "ncol" different situations such as shocks (column in the graphics),
%         4th D: low band, estimate, high band (in each graph).
%   xlab:   x-axis labels on the top
%   ylab:   y-axis labels on the left
%   tstring:  string for time (e.g., month or quarter) -- x-axis labels at the bottom
%   YTickIndx:  1: enable YTick; 0: disable;
%     To get a better picture, it is sometimes better to set YtickIndx to zero.
%   scaleIndx:  1: enable scale along Y-axis; 0: disable
%   nrowg: number of rows in the graph
%   ncolg: number of columns in the graph
%
%  See imrgraph, imcerrgraph, imrerrgraph
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

nstp = size(imfn,1);
nrow = size(imfn,2);
ncol = size(imfn,3);
nmodels = size(imfn,4);
t = 1:nstp;
treverse = fliplr(t);

if (nargin < 9)
  nrowg = nrow;
  ncolg = ncol;
end
if (nrowg*ncolg < nrow*ncol)
   nrowg
   ncolg
   nrow
   ncol

   error('fn_multigraphn_3g_shadedbands.m: nrowg*ncolg must be greater than nrow*ncol')
end


tempmax=zeros(ncol,1);
tempmin=zeros(ncol,1);
maxval=zeros(nrow,1);
minval=zeros(nrow,1);
for i = 1:nrow
   for j = 1:ncol
      tempmax(j) = -realmax;
      tempmin(j) = realmax;
      for k=1:nmodels
         jnk = max(imfn(:,i,j,k));
         tempmax(j) = max([jnk tempmax(j)]);
         %
         jnk = min(imfn(:,i,j,k));
         tempmin(j) = min([jnk tempmin(j)]);
      end
	end
   maxval(i)=max(tempmax);
   minval(i)=min(tempmin);
end

scaleout = [maxval(:) minval(:)];

%--------------
%  Column j: Shock 1 to N; Row i: Responses to
%-------------
%figure


rowlabel = 1;
for i = 1:nrow      % column: from top to bottom
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


   scale=[1 nstp minval(i) maxval(i)];
   
   for j = 1:ncol        % row: from left to right
      k1=(i-1)*ncol+j;
      subplot(nrowg,ncolg,k1)

      %set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
      %   'DefaultAxesLineStyleOrder','-|--|:')
      %set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
      %---<<>>
      %set(0,'DefaultAxesColorOrder',[0 0 0],...
      %   'DefaultAxesLineStyleOrder','-.|-.|-|--|-.|-*|--o|:d')

      nseries = zeros(nstp, nmodels);
      for k=1:nmodels
         nseries(:,k) = imfn(:,i,j,k);
      end     
      
      %---<<>> 
      fill([t treverse],[nseries(:,[3])' fliplr(nseries(:,[1])')],[0.8,0.8,0.8],'EdgeColor','none');
      %plot(t,nseries(:,[1 3]), '-.k','LineWidth',1.0); %,'Color','k');
      hold on
      plot(t,nseries(:,[2]), '-k','LineWidth',1.5);
      %plot(t,nseries(:,[4]), '--k','LineWidth',1.6);
      hold off
      
      %set(gca,'LineStyleOrder','-|--|:|-.')
      %set(gca,'LineStyleOrder',{'-*',':','o'})
      grid;
      if scaleIndx
         axis(scale);
      end
      %set(gca,'YLim',[minval(i) maxval(i)])
      %
		set(gca,'XTick',XTick)
      if YTickIndx
         set(gca,'YTick',yt)
      end

      if i<nrow
        set(gca,'XTickLabelMode','manual','XTickLabel',[])
     end
      %set(gca,'XTickLabel',' ');
      if (scaleIndx) && (j>1)
         set(gca,'YTickLabel',' ');
      end
      if rowlabel == 1
         %title(['x' num2str(j)])
         %title(eval(['x' num2str(j)]))
         if (~isempty(xlab)), title(char(xlab(j))), end
      end
      if columnlabel == 1
         %ylabel(['x' num2str(i)])
         %ylabel(eval(['x' num2str(i)]))
			if (~isempty(ylab)), ylabel(char(ylab(i))), end
      end
      if (i==nrow)  && (~isempty(tstring))
         xlabel(tstring)
      end
      columnlabel = 0;
   end
   rowlabel = 0;
end


%Order of line styles and markers used in a plot.
%This property specifies which line styles and markers to use and in what order
%when creating multiple-line plots. For example,set(gca,'LineStyleOrder', '-*|:|o')sets LineStyleOrder to solid line with asterisk
%marker, dotted line, and hollow circle marker. The default is (-), which specifies
%a solid line for all data plotted. Alternatively, you can create a cell array
%of character strings to define the line styles:set(gca,'LineStyleOrder',{'-*',':','o'})MATLAB supports four line styles, which you can specify any number of
%times in any order. MATLAB cycles through the line styles only after using
%all colors defined by the ColorOrder property. For example,
%the first eight lines plotted use the different colors defined by ColorOrder with
%the first line style. MATLAB then cycles through the colors again, using the
%second line style specified, and so on.You can also specify line style and color directly with the plot and plot3 functions
%or by altering the properties of theline or
%lineseries objects after creating the graph. High-Level Functions and LineStyleOrderNote that, if the axes NextPlot property is set
%to replace (the default), high-level functions like plot reset
%the LineStyleOrder property before determining the line
%style to use. If you want MATLAB to use a LineStyleOrder that
%is different from the default, set NextPlot to replacechildren. Specifying a Default LineStyleOrderYou can also specify your own default LineStyleOrder.
%For example, this statementset(0,'DefaultAxesLineStyleOrder',{'-*',':','o'})
%creates a default value for
