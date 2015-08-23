function logpdf_invgam = fn_logpdf_invgam(x_invgam, a_invgam, b_invgam)
%-----------------------------------------------------------------------------------
%------------------------ Inverse-Gamma distribution ------------------------------%
%--- p(x) = ( b^a/Gamma(a) ) x^(-a-1) exp(-b/x) for a>0 and b>0.
%---    where a is shape and b is scale parameter.
%--- E(x) = b/(a-1) for a>1;  var(x) = b^2/( (a-1)^2*(a-2) ) for a>2;
%--- Noninformative distribution: a,b -> 0.
%--- How to draw: (1) draw z from Gamma(a,b); (2) let x=1/z.
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

logpdf_invgam = a_invgam*log(b_invgam) - gammaln(a_invgam) - (a_invgam+1)*log(x_invgam) - b_invgam/x_invgam;

