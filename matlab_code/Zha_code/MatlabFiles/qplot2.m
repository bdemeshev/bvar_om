function dummyout = qplot2(datavector,inityear,initquart,quarter_label)
%QPLOT Quarterly data plot
%QPLOT(DATA,YEAR,QUARTER,LABEL) plots quarterly data against years
%and quarters begining with YEAR:QUARTER. Quarters are not labeled if
%LABEL=0.
%
%e.g. if GDP is a vector containing quarterly observations
%on gdp beginning in the third quarter of 1985, the command
%QPLOT(GDP,85,3) plots the observations against years and quarters.
%You may specify the year as 85 or 1985 but using 85 looks better.
%
%Setting LABEL to the values 5 or 10 will cause only every fifth
%or every tenth year to be labeled beginning with the first decade
%after YEAR.

% QPLOT was written by Clark A. Burdick of the research
% department of the Federal Reserve Bank of Atlanta.
% Original: June 18, 1997
% Last Modified: July 23, 1997

% Copyright (C) 1997-2012 Clark A. Burdick
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


% TO BE DONE:
% 1. Add better bullet proofing
% 2. Cleaning and streamlining
% 3. Better string manipulation
% 4. Options for XTick


if nargin == 1
   warning('Proper usage is QPLOT(data,year,quarter,label')
   inityear = input('What is the initial year?  ');
   initquart = input('What is the initial quarter?  ');
   quarter_label = input('Enter 1 for labeled quarters, 0 for unlableled quarters: ');
end

if nargin == 3
   quarter_label = 1;
end

nquart = length(datavector);
nyear = ceil(nquart/4);

years = inityear + [0:nyear];
if quarter_label == 5
   years = years.*[round(years/5) == years/5];
end
if quarter_label == 10
   years = years.*[round(years/10) == years/10];
end

xtickvec = (kron(years,[1 0 0 0])+kron(ones(1,nyear+1),[0 2 3 4]))';

if quarter_label ~= 1
   xtickvec = num2str(xtickvec);
   xtickvec = strrep(xtickvec,{' 0'},{' '});
   xtickvec = strrep(xtickvec,{' 2'},{' '});
   xtickvec = strrep(xtickvec,{' 3'},{' '});
   xtickvec = strrep(xtickvec,{' 4'},{' '});
end


plot(datavector);
set(gca,'XTick',[1:nquart]);
set(gca,'XTickLabel',xtickvec(initquart:initquart+(nquart-1)));