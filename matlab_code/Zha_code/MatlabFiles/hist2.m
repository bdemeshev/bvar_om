function [no,xo,binwidth] = hist2(y,x)
%HIST2  Histogram.
%   N = HIST2(Y) bins the elements of Y into 10 equally spaced containers
%   and returns the number of elements in each container.  If Y is a
%   matrix, HIST2 works down the columns.
%
%   N = HIST2(Y,M), where M is a scalar, uses M bins.
%
%   N = HIST2(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X.
%
%   [N,X] = HIST2(...) also returns the position of the bin centers in X.
%
%   [N,X,BW] = HIST2(...) also returns the position of the bin centers in X and the width of the bin.
%
%   HIST2(...) without output arguments produces a histogram bar plot of
%   the results.

%   J.N. Little 2-06-86
%   Revised 10-29-87, 12-29-88 LS
%   Revised 8-13-91 by cmt, 2-3-92 by ls.
%   Copyright (c) 1984-97 by The MathWorks, Inc.
%   $Revision: 5.10 $  $Date: 1997/04/08 05:22:50 $

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

if nargin == 0
    error('Requires one or two input arguments.')
end
if nargin == 1
    x = 10;
end
if min(size(y))==1, y = y(:); end
if isstr(x) | isstr(y)
    error('Input arguments must be numeric.')
end
[m,n] = size(y);
if length(x) == 1
    miny = min(min(y));
    maxy = max(max(y));
    binwidth = (maxy - miny) ./ x;
    xx = miny + binwidth*(0:x);
    xx(length(xx)) = maxy;
    x = xx(1:length(xx)-1) + binwidth/2;
else
    xx = x(:)';
    miny = min(min(y));
    maxy = max(max(y));
    binwidth = [diff(xx) 0];
    xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
    xx(1) = miny;
    xx(length(xx)) = maxy;
end
nbin = length(xx);
nn = zeros(nbin,n);
for i=2:nbin
    nn(i,:) = sum(y <= xx(i));
end
nn = nn(2:nbin,:) - nn(1:nbin-1,:);
if nargout == 0
    bar(x,nn,'hist');
else
  if min(size(y))==1, % Return row vectors if possible.
    no = nn';
    xo = x;
  else
    no = nn;
    xo = x';
  end
end