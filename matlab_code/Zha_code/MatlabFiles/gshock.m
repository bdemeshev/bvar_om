function ExactSmyear = gshock(begy,begm,finy,finm,xlab,idfile,tlinput)
% ExactSmyear = gshock(begy,begm,finy,finm,xlab)
%  Plot the structural shocks for the range given by the inputs
%
% begy:  the beginning year for the graph
% begm:  the begiining month for the graph
% finy:  the end year for the graph
% finm:  the end month for the graph
% xlab:  label for each structural equation
%-----------
% ExactSmyear:  time (myears) and all structural shocks
%  if no nargout, plot graphics.
%
% October 1998 by Tao A. Zha
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



eval(['load ' idfile '.mat']);

%
ibegy = find(myears==begy);
ifiny = find(myears==finy);

ibegy1 = find(myears==begy+1);

%
if isempty(ibegy) & isempty(ibegy1)
	warning('Either ibegy or begm is out of the range of myears')
	disp('Print myears to see the range')
	disp('Or change actuap to make the range longer and save new results to Xshock.mat')
	pause
elseif isempty(ibegy) & (begm<myears(1))
	warning('Either ibegy or begm is out of the range of myears')
	disp('Print myears to see the range')
	disp('Or change actuap to make the range longer and save new results to Xshock.mat')
	pause
elseif isempty(ibegy)
   ibegy = -myears(1)+2;  % +2 is needed for the following -1 in
	                        % bar(Estrexa(ibegy+begm-1:ifiny+finm-1,i))
									%so that it 2-1=1 so that +1 is what we want
end

%
if isempty(ifiny)
	warning('Either ifiny or finm is out of the range of myears')
	disp('Print myears to see the range')
	disp('Or change actuap to make the range longer and save new results to Xshock.mat')
	pause
elseif (ifiny+finm-1>size(myears,1))
	warning('Either ifiny or finm is out of the range of myears')
	disp('Print myears to see the range')
	disp('Or change actuap to make the range longer and save new results to Xshock.mat')
	pause
end


if nargout==0
	t=length(Estrexa(ibegy+begm-1:ifiny+finm-1,1));
	for i=1:size(Estrexa,2)
	   figure
	   %plot(1:t, Estrexa(ibegy+begm-1:ifiny+finm-1,i));
		bar(Estrexa(ibegy+begm-1:ifiny+finm-1,i))
		title([ tlinput ' Monthly Structural Shocks' ])
		ylabel(char(xlab(i)))
		xlabel([ 'From ' num2str(begy) ':' num2str(begm) ' to ' ...
							  num2str(finy) ':' num2str(finm) ])
		grid
	end
else
   ExactSmyear = [myears(ibegy+begm-1:ifiny+finm-1) Estrexa(ibegy+begm-1:ifiny+finm-1,:)];
end
