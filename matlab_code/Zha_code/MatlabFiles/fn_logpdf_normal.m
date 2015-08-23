function logpdf_normal = fn_logpdf_normal(x_normal, a_normal, b_normal)
%-----------------------------------------------------------------------------------
%------------------------------ Normal distribution -------------------------------%
%--- p(x) = normalpdf(a,b^2), where a is E(x) and b^2=Var(x).
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

logpdf_normal = -0.5*log(2.0*pi*b_normal^2.0) - 0.5*(x_normal-a_normal)^2.0/b_normal^2.0;
   
