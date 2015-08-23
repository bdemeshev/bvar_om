function [xc,yc,z,zn] = fn_histwpdfg2D(s,n,w,ditp,xlabstr,ylabstr)
% [xc,yc,z,zn] = histwpdfg2D(s,w,n,ditp,xlabstr,ylabstr)
%     Given uneven weights and 2D draws, generate 2D pdf and probabilites by scaling 2D histogram
%
% s: ndraws-by-2 matrix where ndraws is # of draws and 2 variables
% n: 1-by-2 vector -- # of bins (integers in general) for both x-axis and y-axis
% w: ndraws-by-1 vector of (uneven) weights. Optional.  If [] or no entry, even weights are given.
% ditp: Postive number -- degrees of bicubic interpolation. Bigger ditp is, more points for z with
%          finer the (x,y) plance are interpolated.  Optional.  If [] or no entry, no interpolation.
% xlabstr:  xlabel string.  Optional.  If [] or no entry, no x label.
% ylabstr:  ylabel string.  Optional.  If [] or no entry, no y label.
%---------------------------------
% xc: the position of centered bin on x-axis, from left to right.  All rows are identical
% yc: the position of centered bin on y-axis, from top to bottom.  All columns are identical
% z: size(xc,2)-by-size(yc,1) -- the pdf value in each rectangular cell
% zn: size(xc,2)-by-size(yc,1) -- the probability value in each rectangular cell
% if nargout==0, plot 2D p.d.f. graphics
%
% January 1999 by Tao Zha.  Revised, March 2002.
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

if (size(s,2)~=2)
   error('The size of 1st argument must be *-by-2 (i.e., 2 columns)')
end

if (nargin<3), w=[]; ditp=[]; xlabstr=[]; ylabstr=[]; end
if (nargin<4), ditp=[]; xlabstr=[]; ylabstr=[]; end
if (nargin<5), xlabstr=[]; ylabstr=[]; end
if (nargin<6), ylabstr=[]; end

if (length(n(:))~=2)
   error('2nd argument must have exactly 2 elements')
end

if (min(n)<3)
   error('2nd argument -- bin size -- must have at least 3 for each axis')
end

ndraws=size(s,1);

%** normalized to 1 for the probability
if isempty(w)
   w = ones(ndraws,1)/ndraws;
else
   w=w(:);  % make sure it's a column vector
   w = w/sum(w);
end


%*** x-axis
minx = min(min(s(:,1)));
maxx = max(max(s(:,1)));
h1 = (maxx-minx)/n(1);   % binwidth for x-axis
x = minx + h1*(0:n(1));
xlen = length(x);    % n(1)+1
x(xlen) = maxx;      % in case n(1) is not an integer
xc = x(1:xlen-1) + h1/2;  % the position of the center of the x-bin
                          % 1-by-n(1) row vector: from left to right
xc = repmat(xc,[n(2) 1]);

%*** y-axis
miny = min(min(s(:,2)));
maxy = max(max(s(:,2)));
h2 = (maxy-miny)/n(2);   % binwidth for y-axis
y = miny + h2*(0:n(2));
ylen = length(y);
y(ylen) = maxy;      % in case n(2) is not an integer
yc = y(1:ylen-1) + h2/2;  % the position of the center of the y-bin
yc = yc(:);               % n(2)-by-1 column vector: from top to bottom
yc = repmat(yc,[1 n(1)]);


zn = zeros(n(2),n(1));  % the reverse of x and y is designed for mesh.
                        % see meshgrid to understand this.
tic
for draws=1:ndraws
   k1 = floor((s(draws,1)-minx)/h1)+1;
   k2 = floor((s(draws,2)-miny)/h2)+1;
   %
   %  if k1==0
   %     k1=1;
   %  end
   %  if k2==0;
   %     k2=1;
   %  end
   %
   if k1>n(1)
      k1=n(1);
   end
   if k2>n(2)
      k2=n(2)
   end
   zn(k2,k1) = zn(k2,k1)+w(draws);   % probability in each rectangular cell
end
timeloop=toc;
disp(['Loop time -- minutes ' num2str(timeloop/60)])

z=zn/(h1*h2);   % converted to the height of the p.d.f.


%*** scaled 2D histogram or p.d.f.
%figure(1)
%*** interpolation
if isempty(ditp)
   colormap(zeros(1,3));    % Set color to black
   mesh(xc,yc,z)
   title('Scaled histogram or p.d.f.')
   xlabel(xlabstr)
   ylabel(ylabstr)
else
   [xi,yi]=meshgrid(minx:h1/ditp:maxx,miny:h2/ditp:maxy);
   zi = interp2(xc,yc,z,xi,yi,'bicubic');
   figure(30)
   colormap(zeros(1,3));    % Set color to black
   mesh(xi,yi,zi)   % or mesh
   title('Interpolated p.d.f.')
   xlabel(xlabstr)
   ylabel(ylabstr)
end
