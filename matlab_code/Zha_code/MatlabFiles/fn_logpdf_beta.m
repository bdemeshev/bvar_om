function logpdf_beta = fn_logpdf_beta(x_beta, a_beta, b_beta)
%-----------------------------------------------------------------------------------
%------------------------------- Beta distribution --------------------------------%
%--- p(x) = ( Gamma(a+b)/(Gamma(a)*Gamma(b)) ) x^(a-1) (1-x)^(b-1) for a>0 and b>0.
%--- E(x) = a/(a+b);  var(x) = a*b/( (a+b)^2*(a+b+1) );
%--- The density is finite if a,b>=1.
%--- Noninformative density: (1) a=b=1; (2) a=b=0.5; or (3) a=b=0.
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

logpdf_beta = log(betapdf(x_beta,a_beta,b_beta));

