function [a, b, XLO, XUP] = find_gampar(XLO, XUP, PLO, PUP, a0, b0);

% This function takes as inputs the bounds [XLO, XUP] in the support of the
% gamma distribution (with parameters a and b), the probabilities of the
% bounds [PLO, PUP], and the initial values for ab=[a0, b0]
% and returns the estimates of a and b (as well as XLO and XUP)
% by solving the non-linear functions in gampar(ab, XLO, XUP, PLO, PUP).\

%---------------------------- Gamma distribution ----------------------------------%
%--- p(x) = ( 1/(b^a Gamma(a) ) x^(a-1) exp(-x/b) for a>0, b>0, and x>=0.
%---    where a is the shape and b is the scale parameter.
%--- E(x) = a b;  var(x) = a b^2;
%--- Noninformative distribution: a,b -> 0.
%--- The density function is finite if a >= 1.
%-----------------------------------------------------------------------------------
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

if XLO >= XUP;
    error('the lower bound needs to be smaller than the upper bound')
elseif XLO<= 0;
    error('the support for Gamma distribution needs to be non-negative');
end;

if a0 <= 0 || b0 <= 0;
    error('the values for a0 and b0 need to be positive');
end;

disp(' ')
disp('*************** Convergence results for gamma density ***************')
options = optimset('Display', 'on','TolFun', 1.0e-10, 'TolX', 1.0e-10);
ab_values = fsolve('gampar', [a0, b0], options, XLO, XUP, PLO, PUP);
a = ab_values(1); b = ab_values(2);

% Alternatively, it is possible to constrain the search for the values of a
% and b in the positive range (with 0 being the explicit lower bound) by
% using the lsqnonlin function (see below) instead of the fsolve.  The
% tradeoff is that lsqnonlin is typically slower than fsolve.
% LB_a = 0; LB_b = 0; UB_a = Inf; UB_b = Inf;
% ab_values = lsqnonlin('gampar', [a0, b0], [LB_a, LB_b], [UB_a, UB_b],...
%     options, XLO, XUP, PLO, PUP);
% a = ab_values(1); b = ab_values(2);
