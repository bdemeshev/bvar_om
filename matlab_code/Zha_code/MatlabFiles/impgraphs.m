function impgraphs(response,collabel,rowlabel)
%function impgraphs(response,collabel,rowlabel)
% collabel and rowlabel are cell arrays with string contents
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

XTick=[12 24 36 48]; %modify this line for other time units.
%XTickLabel=[] % Use this to suppress labeling.
[nrow,ncol,nt]=size(response);
for irow=1:nrow
	smin=min(min(response(irow,:,:)));
	smax=max(max(response(irow,:,:)));
	if smin<0
		if smax<=0
			yt=[smin 0];
		else
			yt=[smin 0 smax];
		end
	else % (smin >=0)
		if smax > 0
			yt=[0 smax];
		else % (identically zero responses)
			yt=[-1 0 1];
		end
	end
	for i=1:length(yt)
		if yt(i)==0.0
			if length(yt)==3 & -yt(1)<.2*yt(3)
				ytc(i)={''};
			else
				ytc(i)={'      0'};
			end
		else
			ytc(i)={sprintf('%2.4f',yt(i))};
		end
	end
	scale=[1 nt smin smax];
	for icol=1:ncol
		subplot(nrow,ncol,(irow-1)*ncol+icol);
		plot(squeeze(response(irow,icol,:)));grid
		axis(scale);
		if irow==1
			title(collabel{icol});
		end
		if icol==1
			ylabel(rowlabel{irow});
		end
		set(gca,'XTick',XTick)
		set(gca,'YTick',yt)
		if irow<nrow
			set(gca,'XTickLabelMode','manual','XTickLabel',[])
		end
		if icol>1
			set(gca,'YTickLabelMode','manual','YTickLabel',[])
		else
			set(gca,'YTickLabelMode','manual','YTickLabel',char(ytc))
		end
	end
end