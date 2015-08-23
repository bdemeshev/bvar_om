function simtanzphi(X,k,nvar,U);
% simtanzphi(X,k,nvar,U)
%    Recursive function with global output argument ZPHI
%    Export ZPHI: PHI_k--the span of proj(R(k))(w(a(i))|i~=k), whose rank determines
%           the degree of simultaneity
%
% X: nvar-by-nvar matrix so that the 1st k columns are orthonormal
% k: <=nvar, orthonormal columns
% nvar: number of variables
% U: nvar-by-1 cell.  Each cell contains a set of orthonormal bases for the ith
%         column of idmat0s (restriction or R(i))
%-----------
% ZPHI:  global variable PHI_k--the span of proj(R(k))(w(a(i))|i~=k), whose rank determines
%           the degree of simultaneity
%
% Copyright (c) 1999 by D.F. Waggoner and T.A. Zha
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

global ZPHI

if k<nvar
   V=[X(:,1:k) U{k+1}];
   [Q,R,eindx]=qr(V,0);   % eindx: index vector
   Q1=Q;
   Q1(:,1:k) = Q(:,find(eindx<=k));  % keep the first k col's in V in Q1
   jnk = find(eindx>k);
   jnk1=length(Q(1,:))-k;
   Q1(:,k+1:length(Q(1,:))) = Q(:,jnk(1:jnk1));
   tmp = abs(diag(R));   % largest absoluate value in R
   jnk = find(tmp<max(tmp)*eps);   % index for the zero's in R
   if isempty(jnk)
      Y = Q1(:,k+1:length(tmp));
   else
      Y = Q1(:,k+1:min(jnk)-1);    % all meaningful orthonormal columns after Q(1:k)
   end
   %
   for ik=1:length(Y(1,:))
      X=[X(:,1:k) Y(:,ik)];
      simtanzphi(X,k+1,nvar,U);
   end
else
   Xn=X(:,nvar);
   %* Projection of Xn into R{n}
   Xnproj=U{nvar}*(Xn'*U{nvar})';
   ZPHI=[ZPHI Xnproj];
end
