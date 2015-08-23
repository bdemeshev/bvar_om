function zstar = zroot(B)
% zroot(B);  find roots of a matrix polynomial
% B:  the 3 dimensional array B, e.g. B(:,:,1) is B_1, B(:,:,2) is B_2 and
%     B(:,:,p) is B_p.    Note:  rows correpond to equations.  See Judge (1), pp. 763-764.
%   For the system of equations:
%          y_t = B_1*y_{t-1} + ... + B_p*Y_{t-p} + u_t.
%   The corresponding matrix polynomial is:
%         det(eye(m)-B_1*z-...-B_p*z^p) = 0
% zroot(B) solves for the roots of this above polynomial.  The matrices B_1, B_2,....B_p
%        are passed through B.  The degree of the matrix polynomial is determined
%        by the size of the third dimension of B.
% Each of the B(:,:,i) must be an m by m matrix.

% ZROOT was written by Clark A. Burdick of the research
% department of the Federal Reserve Bank of Atlanta.
% Original: July 31, 1997
% Last Modified: July 31, 1997

% TO BE DONE:
% Better bullet proofing
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


[a b c] = size(B);
syms z;
thepoly = eye(size(B(:,:,1)))-B(:,:,1)*z;
for i=2:c
    thepoly = thepoly-B(:,:,i)*z^i;
end
zstar=roots(sym2poly(det(thepoly)));

