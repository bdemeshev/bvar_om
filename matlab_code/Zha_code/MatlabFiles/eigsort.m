function [vU2,dU2] = eigsort(U,desI)
%  E = EIGSORT(X,desI) is a vector containing the eigenvalues of a square
%    matrix X, where E is sorted descendingly if desI=1 or ascendingly otherwise
%
%  [V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
%    full matrix V whose columns are the corresponding eigenvectors so
%    that X*V = V*D.
%
% Written 9/6/98 by Tao Zha
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

[vU,dU]=eig(U);
dUvec = diag(dU);   % If dU is SPD, use "abs" in front of diag.  This is because,
           % even though dU is supposed to be positive, the numerical solution
			  % for a near-zero eigenvalue can be negative.
			  % *******************
			  %  Second thought: "abs" is not necessary, I think.  If negative, it
			  %      implies zero already -- thus the smallest number
[dUvec2,dUI2] = sort(dUvec);     % ascending
if desI==1      % if descending is required
	dUvec2 = flipud(dUvec2);
	dUI2 = flipud(dUI2);
end

if nargout==1
	vU2 = dUvec2;
	dU2 = NaN;
else
	vU2 = vU(:,dUI2);
	dU2 = diag(dUvec2);    % put it back to matrix form
end