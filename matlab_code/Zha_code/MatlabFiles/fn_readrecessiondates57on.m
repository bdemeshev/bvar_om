function rec_dates = fn_readrecessiondates57on()
% Recession dates.  For some reasons, Matlab code will only take pairs of rec_dates when using rectangle.
% Note that -1 in the following is for the graph lie in the correct position.
%   For example, 1961+(1-1)/12 means Jan 1961 plotted right at 1961.
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

rec_dates=[
    1957+(8-1)/12 1958+(4-1)/12
    1960+(4-1)/12 1961+(2-1)/12
    1969+(12-1)/12 1970+(11-1)/12
    1973+(11-1)/12 1975+(3-1)/12
    1980+(1-1)/12 1980+(7-1)/12
    1981+(7-1)/12 1982+(11-1)/12
    1990+(7-1)/12 1991+(3-1)/12
    2001+(3-1)/12 2001+(11-1)/12
    2007+(12-1)/12 2009+(6-1)/12];

