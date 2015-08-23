function a = clmonq(q)
% function a = monq(q)
%  Find the monthly AR coefficient for interpolated data given
%    an estimate of the quarterly AR coefficient
%    Written by E.M. Leeper
%
% Copyright (C) 1997-2012 E. M. Leeper and Tao Zha
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
% First take the quarterly AR coefficient, q, and then
%     seek the root, a, that solves
%        q = (a^5 + 2 a^4 + 3 a^3 + 2 a^2 + a) / (2 a^2 + 4 a + 3)  (1)

% if n = numerator poly and d = denominator poly in (1)
n = [1 2 3 2 1 0];
d = [2 4 3];
ln = length(n);
ld = length(d);
pad = ln - ld;
dd = [zeros(1,pad) d];
qdd = q.*dd;
newpol = n - qdd;
r = roots(newpol);
for i = 1:length(r)
 if imag(r(i)) == 0
   a = r(i);
 end
end
