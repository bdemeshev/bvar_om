function noargout = ghistd(qmEnd,yrEnd,forepy,q_m,HDv,lendlab,titlab,ylab,timelab,rnum,cnum)
%
%    Plot historical decompositions (annually at this time), no argument out
%
% qmEnd: the end month before the forecast
% yrEnd: the end year before the forecast
% forepy: the number of forecast calendar years
% q_m: 12 or 4 (monthly or quarterly)
% HDv: nHD-by-1 cell where nHD is number of decompositions;
%      each cell is forepy-by-nvar
% lendlab: label for the lengend
% titlab: label for the title
% ylab: label for the variables
% timelab:  label for the time of forecast and indication of insample or out-of-sample
% rnum:  row number for the subplot (e.g., 3)
% cnum:  column number for the subplot (e.g., 2)
%---------------
% no argument out
%
% October 1998 Tao Zha
%% Copyright (C) 1997-2012 Tao Zha
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

%  To be done:
%    lengend does not work well
%    colormap should be replaced by some non-color shading
%    title: need to have only one.

fbegm = qmEnd+1;   % +1 because the first month out of sample
if fbegm > q_m    % the first forecast month into the next year
	fbegm = 1;     % reset
	fbegy = yrEnd+1;   % the first forecast year for the graph
else
	fbegy = yrEnd;    % the first forecast year for the graph
end
ffiny = fbegy+forepy-1;   % the forecast final year

figure
%cnum = 2;
nvar = size(HDv{1},2);
rnum = nvar/cnum;
if ~(rnum==fix(rnum))
	warning('Make that subplot is divided properly')
	disp('Check rnum as well as cnum')
	return
end

hornum = cell(length(fbegy:ffiny),1);    % horizontal number
count=0;
for k=fbegy:ffiny
	count=count+1;
	jnk=num2str(k);
	hornum{count}=jnk(3:4);   % e.g., with '1990', we have '90'
end


for k=1:nvar
	bardata = [];
	for n=1:size(HDv,1)
		bardata = [bardata HDv{n}(:,k)];
	end

   subplot(rnum,cnum,k)
	bar(fbegy:ffiny,bardata,0.75,'group')
	set(gca,'XTickLabel',char(hornum))
	if k==1
		legend(char(lendlab),1)
	end
	colormatrix = [0.2 0.2 0.2;0.8 0.8 0.8];
	colormap(colormatrix)
	if k<=cnum
		title(char(titlab))
	end
	if k>((rnum-1)*cnum)
		xlabel(timelab)
	end
	grid
	ylabel(char(ylab{k}))
end