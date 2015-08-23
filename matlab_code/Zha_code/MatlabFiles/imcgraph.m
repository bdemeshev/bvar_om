function scaleout = imcgraph(imf,nvar,imstp,xlab,ylab,indxGimfml)
%   scaleout = imcgraph(imf,nvar,imstp,xlab,ylab,indxGimfml)
%     imcgraph: impulse, c (column: shock 1 to N), graph
%     Graph the ML point impulse response
%
%     imf:  impulse responses, column (responses to 1st shock, responses to 2nd shock
%               etc), row (impusle steps),
%     nvar: number of variables
%     imstp:  step of impulse responses
%     xlab,ylab:   labels
%     indxGimfml:  1, graph; 0, no graph
%
%  NOTE: I added "indxGimfml" so this function may not be compatible with programs
%            older than 03/06/99, TZ
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

t = 1:imstp;
temp1=zeros(nvar,1);
temp2=zeros(nvar,1);
maxval=zeros(nvar,1);
minval=zeros(nvar,1);
for i = 1:nvar
	for j = 1:nvar
		temp1(j)=max(imf(:,(j-1)*nvar+i));
		temp2(j)=min(imf(:,(j-1)*nvar+i));
	end
   maxval(i)=max(temp1);
   minval(i)=min(temp2);
end

scaleout = [maxval(:) minval(:)];

%--------------
%  Column j: Shock 1 to N; Row i: Responses to
%-------------
if indxGimfml
   figure(1)
   rowlabel = 1;
   for i = 1:nvar
      columnlabel = 1;
      for j = 1:nvar
         k1=(i-1)*nvar+j;
         k2=(j-1)*nvar+i;
         subplot(nvar,nvar,k1)
         plot(t,imf(:,k2),t,zeros(length(imf(:,k2)),1),':');
         set(gca,'YLim',[minval(i) maxval(i)])
         set(gca,'XTickLabel',' ');
         set(gca,'YTickLabel',' ');
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
end