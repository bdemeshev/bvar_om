function [xpd,xc,xp,w,bw] = fn_histwpdfg(s,nbin,gIdx,w,xlabstr,ylabstr,titstr)
% [xpd,xc,xp,w,bw] = fn_histwpdfg(s,nbin,gIdx,w,xlabstr,ylabstr,titstr)
%   Generate probabilities and plot scaled densitys, given the unsorted draws and weights
%
% s:  ndraws-by-n draws where ndraws is # of draws and n is # of series
% nbin:  the number of bins (the maximum is ndraws)
% gIdx:  1 if plotting the pdf; 0 if no graph
% w (optional): ndraws-by-n (unscaled) weights where ndraws is # of draws and n is # of series
% xlabstr (optional): xlabel string
% ylabstr (optional):  ylabel string
% titstr (optional):  title string
%-------------
% xpd: nbin-by-n: density or pdf (not probability) in the centered bin on x-axis.
% xc: nbin-by-n: the position of centered bin on x-axis, from top to bottom.
%                All columns are identical
% xp: nbin-by-n: probability (not density) in the centered bin on x-axis.
%                 NOTE: sum(xp) must be 1 for each column
% w: ndraws-by-n scaled weights so that sum(w)=1 for each column
% bw: bandwidth
%
% August 1999 by Tao Zha
% July 18 2000.  Place w after gIdx so that previous programs need modifications accordingly.
% Oct 1 2000.  Change the order of xpd, xc, and xp so that previous programs need modifications accordingly.
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

if (nargin<=4), titstr=[]; xlabstr=[]; ylabstr=[]; end
[m,n] = size(s);
if (nargin<=3)
   w=ones(m,n)/m;
else
   %** normalized to 1 for the probability
   wsum = repmat(sum(w), [m 1]);
   w = w ./ wsum;
end

%if min(size(s))==1, s=s(:); w=w(:); end % making sure they are column vectors

%*** the position of the center of the bin on the x-axis
mins = min(min(s));
maxs = max(max(s));
bw = (maxs-mins)/nbin;   % binwidth for x-axis
x = mins + bw*(0:nbin);
x(end) = maxs;      % in case nbin is not an integer
xc = x(1:end-1) + bw/2;  % the position of the center of the x-bin
xc = xc';                % nbin-by-1 row vector: from top to bottom
xc = repmat(xc,[1 n]);   % nbin-by-n, same for each column


%*** the probability at each bin on the x-axis
nbin = nbin+1;   % + 1 to get the difference for getting probability xp
nn = zeros(nbin,n);
for i=2:nbin
   for k=1:n
      xidx = find(s(:,k) <= x(i));   % index for the positions
      nn(i,k) = sum(w(xidx,k));
   end
end
xp = nn(2:nbin,:) - nn(1:nbin-1,:);
if (bw<eps)
   xpd=Inf*ones(size(xp));
else
   xpd = xp/bw;    % the density, NOT the probability as xp
end

if gIdx
   plot(xc,xpd)
   %set(gca,'XLim',[-5 5])   % put the limit only on the y-axis
   title(titstr), xlabel(xlabstr), ylabel(ylabstr);
end
