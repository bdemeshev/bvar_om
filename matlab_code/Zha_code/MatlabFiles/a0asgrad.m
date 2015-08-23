function [g,badg] = a0asgrad(x,s,nobs,nvar,a0indx);
% Computes analytical gradient for use in csminwel.m
%      function [g,badg] = a0asgrad(x,s,nobs,nvar,a0indx);
%
%  x (parameter vector),
%  s (diag(S1,...,Sm)): note, as in "a0lhfun", already divided by "nobs"
%  nobs (no of obs),
%  nvar (no of variables),
%  a0indx (matrix indicating the free parameters in A0, and each column in A0 corresponds
%                    to an equation)
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

a0 = zeros(nvar);
%%g = zeros(nvar*nvar,1);     % 4/27/97: not necessary.
badg = 0;

a0(a0indx) = x;

b1=-nobs*inv(a0');
b2 = zeros(nvar,nvar);
for i=1:nvar
   b2(:,i) = nobs*s{i}*a0(:,i);
end

%b = -nobs*inv(a') + nobs*s*a;
b = b1(:)+b2(:);
g = b(a0indx);