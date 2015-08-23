function fn_olotRecessionshades(rec_dates, AxisDat)
% Recession shades.  For some reasons, Matlab code will only take pairs of rec_dates when using rectangle.
% Inputs: rec_dates -- must be in pairs for rectangle to work.
%         AxisDat -- example:  AxisDat=axis([1976 2010 -.2 .15]);
%                      1976 will automatically cut off the dates before 1976 in the following lines.
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

for RecessionNumber=1:length(rec_dates);
    rectangle('Position', [rec_dates(RecessionNumber,1) AxisDat(3) (rec_dates(RecessionNumber,2)-rec_dates(RecessionNumber,1)) (AxisDat(4)-AxisDat(3))], 'FaceColor','y');
end
