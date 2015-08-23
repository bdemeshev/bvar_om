function R = chol2(A)
% R = chol2(A)
%
%     Returns an upper triangular R such that % R * R' is a factorization of
%           a symmetric, positive definite A
%
% Written by Tao Zha, July 1996
% Copyright (C) 1996-2012 Tao Zha
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

%* The following two lines give another expression of the same result
%%function [R,p] = chol2(A)
%%[R,p] = chol(fliplr(flipud(A)));

R = chol(fliplr(flipud(A)));
R = fliplr(flipud(R))';
