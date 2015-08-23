function reset_ini_seed(seednumber)
%After Matlab starts, run this function to reset the seed number; otherwise, Matlab automatically resets to the same initial seed.
% seednumber:  0 -- random state reset to the clock time;
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

if seednumber
   randn('state',seednumber);
   rand('state',seednumber);
else
   randn('state',fix(100*sum(clock)));
   rand('state',fix(100*sum(clock)));
end
