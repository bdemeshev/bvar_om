function fn_reset_ini_seed(seednumber)
%fn_reset_ini_seed(seednumber) or preferablly use dynare_xxx/matlab/set_dynare_seed.m or RandStream.getDefaultStream for a stand-alone product
%
% After Matlab starts, run this function to reset the seed number; otherwise, Matlab automatically resets to the same initial seed.
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

if (seednumber==0)
   seednumber = sum(100*clock);     % if 0, random state at each clock time
end   
rng_stream = RandStream('mt19937ar','seed',seednumber);
RandStream.setGlobalStream(rng_stream);
%RandStream.setDefaultStream(rng_stream);  %This is compatible with older
%versions of Matlab, say before the 2012 release.  

%--- Old random number generator -- please do not use.
%if seednumber
%   randn('state',seednumber);
%   rand('state',seednumber);
%else
%   randn('state',fix(100*sum(clock)));
%   rand('state',fix(100*sum(clock)));
%end
