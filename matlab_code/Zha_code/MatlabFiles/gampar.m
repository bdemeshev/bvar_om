function f = gampar(ab, XLO, XUP, PLO, PUP);

% The function takes as inputs the parameters ab=[a, b] of the gamma
% distribution, the bounds of the support [XLO, XUP], the the corresponding
% probability of the bounds [PLO, PUP] and returns the residual value f.
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

a = abs(ab(1)); b = abs(ab(2));  % abs() is used to allow continuous search for the fsolve function find_gampar.m.
f1 = PLO - gamcdf(XLO, a, 1/b);  %Note that in Matlab, it is 1/b, NOT b.
f2 = PUP - gamcdf(XUP, a, 1/b);  %Note that in Matlab, it is 1/b, NOT b.
% Equivalently, one can use the gaminv function to define the zero
% conditions:
% f1 = XLO - gaminv(PLO, a, b);
% f2 = XUP - gaminv(PUP, a, b);

f = [f1, f2];
