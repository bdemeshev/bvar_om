function scaleout = fn_multigraph1(imf1,nrow,ncol,xlab,ylab,XTick,YTickIndx,scaleIndx)
%Plot only one set of impulse responses in one graph.  See fn_multigraph2.m.
%imf1 -- each has 3 dimensions.  Row: horizon; column: 'nrow' variables; 3rd dim: 'ncol' situations such as regimes 1 and 2.
%   imf: row: "nstp" time horizon (in the graphics), column: "nrow "variables (row in
%             the graphics), 3rd D: across "ncol" different situations (column in
%             the graphics)
%   nrow:  # of rows for the graphics
%   ncol:  # of columns for the graphics
%   YTickIndx:  1: enable YTick; 0: disable
%   scaleIndx:  1: enable scale along Y-axis; 0: disable
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

nstp = size(imf1,1);
t = 1:nstp;
nrow;
ncol;

temp1=zeros(ncol,1);
temp2=zeros(ncol,1);
maxval=zeros(nrow,1);
minval=zeros(nrow,1);
for i = 1:nrow
   for j = 1:ncol
      %jnk1=max(firsth(:,i,j));
      %jnk2=max(firstl(:,i,j));
      %jnk3=max(firsth1(:,i,j));
      %jnk4=max(imf2(:,i,j));
      jnk5=max(imf1(:,i,j));

      temp1(j)=max([jnk5]);
      %
      %jnk1=min(firstl(:,i,j));
      %jnk2=min(firsth(:,i,j));
      %jnk3=min(firstl1(:,i,j));
      %jnk4=min(imf2(:,i,j));
      jnk5=min(imf1(:,i,j));

      temp2(j)=min([jnk5]);
	end
   maxval(i)=max(temp1);
   minval(i)=min(temp2);
end

scaleout = [maxval(:) minval(:)];

%--------------
%  Column j: situation 1 to N; Row i: variables (ie responses to)
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
      subplot(nrow,ncol,k1)
      %plot(t,imf1(:,i,j)) %,t,imf2(:,i,j),'--')
                %t,[firstl(:,i,j) firsth(:,i,j)],':');
      plot(t,imf1(:,i,j),':*b');
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
      if j>1
         set(gca,'YTickLabel',' ');
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
      if i==nrow
         xlabel('Quarters')
      end
      columnlabel = 0;
   end
   rowlabel = 0;
end
