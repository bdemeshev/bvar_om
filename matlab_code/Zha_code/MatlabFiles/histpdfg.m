function [n,z] = histpdfg(seqdraws,binnum,titstr,xlabstr,ylabstr)
% [n,z] = histpdfg(seqdraws,binnum,titstr,xlabstr,ylabstr)
%     Plot (if no nargout) and export (if nargout) pdf by scaling histogram
%                  of an unsorted sequence of draws
%
% seqdraws:   an unordered sequence of draws
% binnum:  the total number of bins in histogram, with 10 being a starting value
% titstr:  title string; if [], no title
% xlabstr:  xlabel string; if [], no xlab
% ylabstr:  ylabel string; if [], no ylab
%--------
% n: a vector of the number of elements in each bin or container
% z: a vector of the position of the center of the bin
% plot p.d.f. graph if no nargout
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

if (nargin==2), titstr=[]; xlabstr=[]; ylabstr=[]; end

[n,z,bw]=hist2(seqdraws,binnum);
n=(n/length(seqdraws))/bw;    % make it p.d.f.
%bar(z,n)
if nargout
   z=z';
   n=n';
else
   %figure
   plot(z,n)
   title(titstr)
   xlabel(xlabstr)
   ylabel(ylabstr)
end
