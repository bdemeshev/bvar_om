function [z2,eu2]=gensys_z2(g0,g1,c,psi,gpi,div)
%[z2,eu2]=gensys_z2(g0,g1,c,psi,gpi,div)
%
% z2: n-by-k in general where k is the dimension of gpi or the number of expectational errors.
% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+gpi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.
% Space spanned by z2 is such that z2'*y_t = 0.
% If div is omitted from argument list, a div>1 is calculated.
%
%  eu2 = [0, 0] if the rows in Z2 corresponding to the fixed points are zeros (which happens only in the greater-than-2 regime case).
%    This value is set outside this function.
%  eu2 = [-2,-2] for coincident zeros.
% By Christopher A. Sims
% Corrected 10/28/96 by CAS
% Copyright (C) 1997-2012 Christopher A. Sims
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

eu2=NaN*ones(2,1);
realsmall=1e-6;
fixdiv=(nargin==6);
n=size(g0,1);
[a b q z v]=qz(g0,g1);
if ~fixdiv, div=1.01; end
nunstab=0;
zxz=0;
for i=1:n
% ------------------div calc------------
   if ~fixdiv
      if abs(a(i,i)) > 0
         divhat=abs(b(i,i))/abs(a(i,i));
	 % bug detected by Vasco Curdia and Daria Finocchiaro, 2/25/2004  A root of
	 % exactly 1.01 and no root between 1 and 1.02, led to div being stuck at 1.01
	 % and the 1.01 root being misclassified as stable.  Changing < to <= below fixes this.
         if 1+realsmall<divhat & divhat<=div
            div=.5*(1+divhat);
         end
      end
   end
% ----------------------------------------
   nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
   if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
      zxz=1;
   end
end

if ~zxz
   [a b q z]=qzdiv(div,a,b,q,z);
end
%gev=[diag(a) diag(b)];
if zxz
   disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
   eu2=[-2;-2];
   % correction added 7/29/2003.  Otherwise the failure to set output
   % arguments leads to an error message and no output (including eu2).
   G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];
   z2=[];
   return
end
%q1=q(1:n-nunstab,:);
%q2=q(n-nunstab+1:n,:);
%z1=z(:,1:n-nunstab)';

%if (nunstab>0)
%  z2=z(:,n-nunstab+1:n);
%else
%  z2=z(:,n-rank(gpi)+1:n);
%end

z2=z(:,n-rank(gpi)+1:n);
