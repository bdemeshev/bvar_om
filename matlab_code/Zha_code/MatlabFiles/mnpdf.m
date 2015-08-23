function y = mnpdf(x,xm,C,constIx)
% y = mnpdf(x,xm,C,constIx)
%    The pdf value for multivariate normal distribution
%
% x:  p-by-draws matrix of values evaluated at where p=size(x,1) is # of variables
% xm: p-by-draws matrix of the mean of x
% C:  p-by-p Choleski square root of PDS S -- the covariance matrix so that S = C*C'
% constIx: index for the constant.  1: constant (normalized); 0: no constant (unnormalized)
%----------
% y:  p-by-draws matrix of pdf's for multivariate normal distribution
%
%   Christian P. Robert, "The Bayesian Choice," Springer-Verlag, New York, 1994,
%      p. 381.
%
% November 1998 by Tao Zha
% rewritten by CAS 12/98 to take matrix x, return vector y
% Copyright (C) 1997-2012  Christopher A. Sims and Tao Zha 
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
   y = exp(-dSh-sum(z.*z,1)/2) / ((2*pi)^(p/2));
   y = y';
else
   y = exp(-sum(z.*z,1)/2);
   y = y';
end
