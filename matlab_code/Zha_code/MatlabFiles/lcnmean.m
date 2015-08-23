function cutmean = lcnmean(x,pnIndx)
%
% cutmean = lcnmean(x,pnIndx)
%
%   The mean of the truncated normal with the lower bound x if pnIndx>0
%            or the upper bound x if pnIndx<0.  Zha Forecast (2) p.27
%
% x:   the lower bound x if pnIndx > 0; the upper bound x if pnIndx < 0
% pnIndx:  1: truncated for posivie mean (tight); 0 truncated for negative mean (loose)
%------------
% cutmean: the backed-out mean for the truncated normal
%
% October 1998 Tao Zha
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

if pnIndx
	cutmean = exp(-x^2/2)/(sqrt(2*pi)*(1-cdf('norm',x,0,1)));
else
	cutmean = -exp(-x^2/2)/(sqrt(2*pi)*cdf('norm',x,0,1));
end
