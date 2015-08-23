function P = invgamcdf(X, a, b);

% This function computes the cumulative density for the inverse gamma
% distribution with shape parameter a and scale parameter b.  It uses the
% following closed-form solution for the cdf of the inverse gamma
% distribution:
% ---------Cumulative density function for Inverse-Gamma distribution-------
% --- P(x; a, b) = G(b/x, a)/G(a), where G(b/x, a) is the upper incomplete
% gamma function and G(a) is the gamma function defined as
%    G(b/x, a) = int_{b/x}^{\infty} t^{a-1} e^{-t} dt.
% The Matlab definition is different, which is
%    gammainc(b/x, a) = (1.0/Gamma(a)) int_{0}^{b/x} t^{a-1} e^{-t} dt.
% Thus, we have
%    G(b/x, a)/G(a) = (1-gammainc(b./X, a)).

% First, check for input errors
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

if X<=0;
    error('the support for the inverse gamma function needs to be positive')
end;
if a<= 0 || b <=0;
    a = abs(a);
    b = abs(b);
          % abs() is used to allow continuous search for the fsolve function find_invgampar.m.
    %error('parameter requirements: a>0, b>0');
end;

P = (1-gammainc(b./X, a));
