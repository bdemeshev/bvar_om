function actual = gactual(begy,begm,begq,finy,finm,finq,ylab,mqyIndx,idfile)
%
% actual = gactual(begy,begm,begq,finy,finm,finq,ylab,mqyIndx,idfile)
%    Actual data according to mqyIndx; graph if no output argument is specified.
%
% begy:  the beginning year for the graph
% begm:  the begiining month for the graph
% begq:  the begiining quarter for the graph
% finy:  the end year for the graph
% finm:  the end month for the graph
% finq:  the end quarter for the graph
% ylab:  label for each variable
% mqyIndx:  1-by-3 index of monthly log, qg, and calendar yg. 1: yes; 0: no.
%             Only one at a time
% idfile:  import xinample.mat in the LZ paper under \condz
%--------------
% actual:  T-by-nvar+1; 1st column is the dates; 2nd-7th columns: actual data of
%              of nvar variables with mlog, qg, or calendar yg, one at a time,
%              depending on mqyIndx
% graph if no output argument is specified.
%
%
%   When insampleg.m is run, xinsample.mat is loaded for a long
%                                           history (e.g.,1960-1998)
%   When mspoint.m is run,  outxshock.mat from msstart is loaded, which varies
%                               with parac.m.
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

if length(find(mqyIndx))>1
	warning('To get the number (not graph) out, only one at a time')
	disp('Make sure that mqyIndx has only 1 in a vector of 1 and 0s')
	disp(' ')
	disp('Press Enter to continue or Ctrl to abort')
	pause
end

if mqyIndx(1)==1     % monthly log
	ibegy = find(myears==begy);
	ifiny = find(myears==finy);

	ibegy1 = find(myears==begy+1);
	%
	if isempty(ibegy) & isempty(ibegy1)
      warning('Either ibegy or begm is out of the range of myears in xinsample.mat')
		disp('Print myears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
      disp('Press ctrl-c to abort now')
      pause
	elseif isempty(ibegy) & (begm<myears(1))
      warning('Either ibegy or begm is out of the range of myears in xinample.mat')
		disp('Print myears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
      disp('Press ctrl-c to abort now')
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
		return
	elseif (ifiny+finm-1>size(myears,1))
		warning('Either ifiny or finm is out of the range of myears')
		disp('Print myears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	end


	if nargout==0
		t=length(yact(ibegy+begm-1:ifiny+finm-1,1));
		for i=1:size(yact,2)
		   figure
		   plot(1:t, yact(ibegy+begm-1:ifiny+finm-1,i),'*:');
			%bar(yact(ibegy+begm-1:ifiny+finm-1,i))
			title('Monthly Log')
			ylabel(char(ylab(i)))
			xlabel([ 'From ' num2str(begy) ':' num2str(begm) ' to ' ...
								  num2str(finy) ':' num2str(finm) ])
			grid
   	end
	end

	actual = [myears(ibegy+begm-1:ifiny+finm-1) yact(ibegy+begm-1:ifiny+finm-1,:)];
end


if mqyIndx(2)==1     % quarterly
	ibegy = find(qyears==begy);
	ifiny = find(qyears==finy);

	ibegy1 = find(qyears==begy+1);
	%
	if isempty(ibegy) & isempty(ibegy1)
		warning('Either ibegy or begq is out of the range of qyears')
		disp('Print qyears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	elseif isempty(ibegy) & (begq<qyears(1))
		warning('Either ibegy or begq is out of the range of qyears')
		disp('Print qyears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	elseif isempty(ibegy)
      ibegy = -qyears(1)+2;  % +2 is needed for the following -1 in
		                        % bar(Estrexa(ibegy+begq-1:ifiny+finm-1,i))
										%so that it 2-1=1 so that +1 is what we want
	end
	%
	if isempty(ifiny)
		warning('Either ifiny or finq is out of the range of qyears')
		disp('Print qyears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	elseif (ifiny+finq-1>size(qyears,1))
		warning('Either ifiny or finq is out of the range of qyears')
		disp('Print qyears to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	end


	if nargout==0
		t=length(yactqg(ibegy+begq-1:ifiny+finq-1,1));
		for i=1:size(yactqg,2)
	   	figure
	   	plot(1:t, yactqg(ibegy+begq-1:ifiny+finq-1,i),'*:');
			title('Quarter-to-quarter Growth Rate')
			ylabel(char(ylab(i)))
			xlabel([ 'From ' num2str(begy) ':' num2str(begq) ' to ' ...
								  num2str(finy) ':' num2str(finq) ])
			grid
   	end
	end

	actual = [qyears(ibegy+begq-1:ifiny+finq-1) yactqg(ibegy+begq-1:ifiny+finq-1,:)];
end


if mqyIndx(3)==1     % calendar year
	ibegy = find(yearsCal==begy);
	ifiny = find(yearsCal==finy);
	%
	if isempty(ibegy)
		warning('Either ibegy or begq is out of the range of yearsCal')
		disp('Print yearsCal to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	end
	%
	if isempty(ifiny)
		warning('Either ifiny or finq is out of the range of yearsCal')
		disp('Print yearsCal to see the range')
		disp('Or change actuap to make the range longer and save new results to Xshock.mat')
		return
	end

	if nargout==0
		t=length(yactCalyg(ibegy:ifiny,1));
		for i=1:size(yactCalyg,2)
	   	figure
	   	plot(1:t, yactCalyg(ibegy:ifiny,i),'*:');
			title('Calendar Annual Average Growth Rate')
			ylabel(char(ylab(i)))
			xlabel([ 'From 19' num2str(begy) ' to 19' ...
								  num2str(finy) ])
			grid
   	end
	end

	actual = [yearsCal(ibegy:ifiny) yactCalyg(ibegy:ifiny,:)];
end
