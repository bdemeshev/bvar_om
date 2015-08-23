function logpdf_gam = fn_logpdf_gamma(x_gam, a_gam, b_gam)
%-----------------------------------------------------------------------------------
%---------------------------- Gamma distribution ----------------------------------%
%--- p(x) = ( b^a/Gamma(a) ) x^(a-1) exp(-bx) for a>0 and b>0.
%---    where a is shape and b is inverse scale (rate) parameter.
%--- E(x) = a/b;  var(x) = a/b^2;
%--- Noninformative distribution: a,b -> 0.
%--- The density function is finite if a >= 1.
%-----------------------------------------------------------------------------------
% Written by T. Zha; 12:20AM 02/18/2013
% Copyright (C) 2013 Tao Zha
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

if (a_gam == 1.0)
   logpdf_gam = a_gam*log(b_gam) - gammaln(a_gam) - b_gam*x_gam;
else
   logpdf_gam = a_gam*log(b_gam) - gammaln(a_gam) + (a_gam-1.0)*log(x_gam) - b_gam*x_gam;
end   