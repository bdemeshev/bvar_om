function y = mtpdf(x,xm,C,v,constIx)
% y = mtpdf(x,xm,C,v)
%    The pdf value for multivariate Student t distribution
%
% x:  p-by-draws matrix of values evaluated at where p=size(x,1) is # of variables
% xm: p-by-draws matrix of the mean of x
% C:  p-by-p Choleski square root of PDS S, which is the covariance matrix
%          in the normal case so that S = C*C'
% v (>0):  scaler -- the degrees of freedom
% constIx:  index for the constant.  1: constant (normalized); 0: no constant (unnormalized)
%----------
% y:  p-by-draws matrix of pdf's for multivariate Student t distribution with v degrees of freedom
%
%   Christian P. Robert, "The Bayesian Choice," Springer-Verlag, New York, 1994,
%      p. 382.
%
%   December 1998 by Tao Zha
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

[p,nx]=size(x);
z = C\(x-xm);


if constIx
   dSh=sum(log(diag(C)));   % (detSigma)^(1/2)
   %* Use gammaln function to avoid overflows.
   term = exp(gammaln((v + p) / 2) - gammaln(v/2) - dSh);
   y = (term*(v*pi)^(-p/2)) * (1 + sum(z.*z,1)/v).^(-(v + p)/2);
   y = y';
else
   y = (1 + sum(z.*z,1)/v).^(-(v + p)/2);
   y = y';
end
