function Q = SRestrictRWZalg(A0hatinv,Bhat,nvar,lags,irs)
% Rubio-Waggoner-Zha (RWZ) method of sign restrictions.  For related methods, see Canova, Faust, and Uhlig.
%      The detailed theoretical foundation of this algorithm can be found in Theorem 3 of Rubio, Waggoner, and Zha (RWZ)'s article "Regime Changes in the Euro Area."
%      Other M functions called by this function can be downloaded by clicking on Archived Matlab Library ZhaZippedCode on http://home.earthlink.net/~tzha02/programCode.html
% Strcutural VAR form:   Y*A0hat = X*Aphat + E, X where Y is T-by-nvar, A0hat is nvar-by-nvar, X is T-by-k (including the constant term and all other exogenous terms),
%   Aphat is k-by-nvar, and E is T-by-nvar.  Rows are in the order of 1st lag (with nvar variables) to lags (with nvar variables) plus the exogenous terms.
%   Note that columns of A0hat or Aphat correspond to equations.
% Inputs:
%   A0hatinv = inv(A0hat).
%   Bhat = Aphat*inv(A0hat).
%   nvar = number of endogenous variables.
%   lags = lag length.
%   irs = maximum number of periods in which sign restrictions are imposed.
% Outputs:
%   Q: orthogonal rotation matrix so that Q*A0hatinv or A0hat*Q' gives impulse responses that will satisfy
%      sign restrictions of Canova, Faust, and Uhlig.
%
% Modified Nov 2004 by T. Zha to
%     (1) correct the existing bugs;
%     (2) make the signs explicit to avoid the normalization problem when computing error bands;
%     (3) construct efficient way (i.e., earlier exit) to make all restrictions satisfied;
%     (4) construct efficient way to normalize.
% In this example, we have
%    Variables are in the following order: 1: y; 2: P; 3: R; 4: M3; 5: Exec (per $).
%    Shocks are in the following order: 1: AS; 2: AD; 3: MP; 4: MD; 5: Exec.
%
% Copyright (c) Rubio, Waggoner, and Zha, 2005.
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


nn = [nvar lags irs];
control=0;


Q=eye(nvar);

while control==0
   newmatrix=normrnd(0,1,nvar,nvar);
   [Q,R]=qr(newmatrix);
   for i=1:nvar;
       if R(i,i)<0
           Q(:,i)=-Q(:,i);
       end
   end
   imfhat = fn_impulse(Bhat,Q*A0hatinv,nn);
      %In the form that is congenial to RATS
   imf3hat=reshape(imfhat,size(imfhat,1),nvar,nvar);
      %imf3hat: row--steps, column--nvar responses, 3rd dimension--nvar shocks

   %=== Responses to a moneaty policy (MP) shock.
   % R>0, M<0, y<0, and P<0 for the irs periods.
   a = (imf3hat(1:irs,3,3) > 0) .* (imf3hat(1:irs,4,3) < 0) .* (imf3hat(1:irs,1,3) < 0) .* (imf3hat(1:irs,2,3) < 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock.
      am = (imf3hat(1:irs,3,3) < 0) .* (imf3hat(1:irs,4,3) > 0) .* (imf3hat(1:irs,1,3) > 0) .* (imf3hat(1:irs,2,3) > 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(3,:) = -Q(3,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end
   %--- R>0 and M<0 for the irs periods.
   %   a = (imf3hat(1:irs,3,3) > 0) .* (imf3hat(1:irs,4,3) < 0);
   %   if (max(a)==0)
   %      %--- Swiching the sign of the shock.
   %      am = (imf3hat(1:irs,3,3) < 0) .* (imf3hat(1:irs,4,3) > 0);
   %      if (min(am)==0)
   %         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
   %      else
   %         %--- Normalizing according to the switched sign.
   %         Q(3,:) = -Q(3,:);
   %      end
   %   elseif (min(a)==0)
   %      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   %   end

   %=== Responses to an money demand (MD) shock.
   % R>0 and M>0 for the irs periods.
   a = (imf3hat(1:irs,3,4) > 0) .* (imf3hat(1:irs,4,4) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (imf3hat(1:irs,3,4) < 0) .* (imf3hat(1:irs,4,4) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(4,:) = -Q(4,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to an aggregate demand (AD) shock.
   % P>0 and y>0 for the irs periods.
   a = (imf3hat(1:irs,1,2) > 0) .* (imf3hat(1:irs,2,2) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (imf3hat(1:irs,1,2) < 0) .* (imf3hat(1:irs,2,2) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(2,:) = -Q(2,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to an aggregate supply (AS) shock.
   % P>0 and y<0 for the irs periods.
   a = (imf3hat(1:irs,1,1) < 0) .* (imf3hat(1:irs,2,1) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (imf3hat(1:irs,1,1) > 0) .* (imf3hat(1:irs,2,1) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(1,:) = -Q(1,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to an exchange-rate shock (depreciation ==> exports >0 ==> y>0);
   % Ex>0 and y>0 for the irs periods.
   a= (imf3hat(1:irs,1,5) > 0) .* (imf3hat(1:irs,5,5) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (imf3hat(1:irs,1,5) < 0) .* (imf3hat(1:irs,5,5) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(5,:) = -Q(5,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %--- Terminating condition: all restrictions are satisfied.
   control=1;
end
