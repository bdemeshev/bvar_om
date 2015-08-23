function H = bfgsi(H0,dg,dx)
% H = bfgsi(H0,dg,dx)
% dg is previous change in gradient; dx is previous change in x;
% 6/8/93 version that updates inverse hessian instead of hessian
% itself.
% Copyright by Christopher Sims 1996.  This material may be freely
% reproduced and modified.

% Copyright (C) 1996 Christopher Sims
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
dispIndx = 1;   % 1: turn on all the diplays on the screen; 0: turn off (Added by T. Zha)

if size(dg,2)>1
   dg=dg';
end
if size(dx,2)>1
   dx=dx';
end
Hdg = H0*dg;
dgdx = dg'*dx;
if (abs(dgdx) >1e-12)
   H = H0 + (1+(dg'*Hdg)/dgdx)*(dx*dx')/dgdx - (dx*Hdg'+Hdg*dx')/dgdx;
else
   if dispIndx
      disp('bfgs update failed.')
      disp(['|dg| = ' num2str(sqrt(dg'*dg)) '|dx| = ' num2str(sqrt(dx'*dx))]);
      disp(['dg''*dx = ' num2str(dgdx)])
      disp(['|H*dg| = ' num2str(Hdg'*Hdg)])
   end
   H=H0;
end
save H.dat H
