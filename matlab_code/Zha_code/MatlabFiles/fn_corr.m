function corM = fn_corr(covM)
%  function corM = fn_corr(covM)
%
%  covM:  covariance matrix (input)
%  corM:  correlation matrix  (output)
%  6 September 1998
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


[nvar,jnk]=size(covM);
if nvar~=jnk
	error('The input matrix must be square!!')
end
corM=zeros(nvar);
for i=1:nvar
   for j=1:nvar
     corM(i,j)=covM(i,j) / sqrt( covM(i,i) * covM(j,j));
   end
end

