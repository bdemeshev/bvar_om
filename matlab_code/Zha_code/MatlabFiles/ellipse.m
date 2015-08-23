function [x,y,z] = ellipse(sigma,mu,k)
%ELLIPSE Creates an ellipsoid
%[X,Y,Z] = ELLIPSE(SIGMA,MU,K) Creates an arbitrary
%ellipsoid in 2 or 3 dimensions. The ellipsoid represents
%the solutions to the quadratic equation:
%
%             (x-MU)'*SIGMA*(x-MU) = K
%
%SIGMA is a positive definate matrix of size 2x2 or 3x3.
%SIGMA determines the shape of the ellipsoid.
%MU is a vector conformable with SIGMA; either 2x1 or 3x1.
%MU determines the location of the ellipsoid.
%K is a real constant.
%K determines the size of the ellipsoid.
%The vector x=[X;Y] in 2D or x=[X,Y,Z] in 3D.
%
%If no output aruments are specified, the resulting
%ellipsoid is plotted. One can create these plots manually
%using PLOT(X,Y) in the 2 dimensional case and
%using MESH(X,Y,Z) in the 3 dimensional case.

% QPLOT was written by Clark A. Burdick of the research
% department of the Federal Reserve Bank of Atlanta.
% Original: August 19, 1997
% Last Modified: August 19, 1997
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
%      I'm open to suggestions.

n=50;  % <<>>
theta=pi*(-n:2:n)/n;
phi=(pi/2)*(-n:2:n)'/n;

R=chol(sigma);
rho=sqrt(k);

if length(mu) == 2
   X = rho*cos(theta);
   Y = rho*sin(theta);

   for i = 1:length(theta)
      circlevec = [X(i);Y(i)];
      ellipsevec = R\circlevec;     % T.A. Zha, 8/16/98
      x(i) = mu(1) + ellipsevec(1);
      y(i) = mu(2) + ellipsevec(2);
   end
   if nargout == 0
      plot(x,y)
   end
   z='This was only a 2 dimensional ellipse';
end

if length(mu) == 3
   X = rho*cos(phi)*cos(theta);
   Y = rho*cos(phi)*sin(theta);
   Z = rho*sin(phi)*ones(size(theta));

   for i = 1:length(phi)
      for j = 1:length(theta)
         spherevec = [X(i,j); Y(i,j); Z(i,j)];
         ellipsevec = R\spherevec;    % T.A. Zha, 8/16/98
         x(i,j) = mu(1) + ellipsevec(1);
         y(i,j) = mu(2) + ellipsevec(2);
         z(i,j) = mu(3) + ellipsevec(3);
      end
   end
   if nargout == 0
      mesh(x,y,z)
   end
end