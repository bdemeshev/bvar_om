function startup()
% Remembers a previously set startup directory to set the directory, invoke startd.  E.g.,
%  sd=cd % Set sd to the current directory
%  startd(sd)
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

%----------------------------------

%path(path,'d:\Program Files\MATLAB\R2007a\work\cstz')
%path(path,'c:\softwdisk\matlabr12\toolbox\cstz\cmexfiles\csminwelfinal\csminwelmex')
%path(path,'C:\ZhaData\TZCcode')
%path(path,'C:\Program Files\Intel\MKL\ia32\lib')
%path(path,'C:\Program Files\Intel\MKL\include')
%path(path,'d:\matlabr12\toolbox\cstz\lzpaper2')
%path(path,'d:\matlabr12\toolbox\cstz\rvarcode')
%path(path,'C:\Program Files\MATLAB\R2006b\work\cstz')
%path(path,'/Users/tzha/ZhaData/Git/TZcode/MatlabFiles') 
%path(path,'/Users/tzha/ZhaData/Git/TZcode/MatlabFiles/_MSV')
if exist('startdir0.mat')==2
	load /Users/tzha/ZhaGit/TZcode/MatlabFiles/startdir0
	cd(sd)
end
format compact

fn_reset_ini_seed(0); %Reset the random seed to clockcycle.  Can be reoverwritten by fn_reset_ini_seed(number);



