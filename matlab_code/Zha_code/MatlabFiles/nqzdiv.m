function [A,B,P,R] = nqzdiv(stake,A,B,V,P)
%function [A,B,Q,Z] = nqzdiv(stake,A,B,V,P)
%
% Takes U.T. matrices A, B, V, and P and rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right
% corner, while the output A and B are diagonal, and P*A*R=G0, P*B*R=G1.
% Note: the input A and B are U.T.; the output A and B are diagonal.
%
% Modified by T.Zha, 5/26/97
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


[n jnk] = size(A);

ald = diag(diag(A));  % d: diagoanl
bed = diag(diag(B));
R = inv(V);

root = abs([diag(ald) diag(bed)]);
root(:,1) = root(:,1)-(root(:,1)<1.e-13).*(root(:,1)+root(:,2));
root(:,2) = root(:,2)./root(:,1);
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if (root(j,2) > stake | root(j,2) < -.1)
         m=j;
         break
      end
   end
   if (m==0)
      return
   elseif (m<i-1)
      pindx = 1:n;      % T.Z., 5/26/97
      pindx(m)=i-1;
      pindx(i-1)=m;
      %*** permutation begins
      ald = ald(pindx,pindx);
      bed = bed(pindx,pindx);
      P = P(:,pindx);   % column permutation, T.Z., 5/26/97
      R = R(pindx,:);   % row permutation
   end
end

%*** return the new Q, Z, A, B
A = diag(1./diag(bed));
B = diag(1./diag(ald));