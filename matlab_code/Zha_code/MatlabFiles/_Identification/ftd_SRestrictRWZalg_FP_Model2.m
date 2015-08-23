function [Q, indx_fail] = ftd_SRestrictRWZalg_FP_Model2(A0hatinv,Bhat,nvar,lags,irs,ndraws_terminate)
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
%   ndraws_terminate = number of draws of Q before terminating the function.
% Outputs:
%   Q: orthogonal rotation matrix so that Q*A0hatinv or A0hat*Q' gives impulse responses that will satisfy
%      sign restrictions of Canova, Faust, and Uhlig.
%   indx_fail: 0: success; 1: fail.
%
% Copyright (c) Rubio, Waggoner, and Zha, 2005.
% Modified Nov 2004 by T. Zha to
%     (1) correct the existing bugs;
%     (2) make the signs explicit to avoid the normalization problem when computing error bands;
%     (3) construct efficient way (i.e., earlier exit) to make all restrictions satisfied;
%     (4) construct efficient way to normalize.
% In this example, we have
%    Variables are in the following order: 1: gov_exp; 2: gov_tax; 3: gdp.
%    Shocks are in the following order: 1: basic spending shock; 2: basic tax shock; 3: business shock.



nn = [nvar lags irs];
control=0;
cnt = 0;
indx_fail = 0;

Q=eye(nvar);

while control==0
   cnt = cnt + 1;
   if cnt > ndraws_terminate
      indx_fail = 1;
      break;
   end

   if (cnt == 1)
      Q = eye(nvar);
   else   
      newmatrix=normrnd(0,1,nvar,nvar);
      [Q,R]=qr(newmatrix);
      for i=1:nvar;
          if R(i,i)<0
              Q(:,i)=-Q(:,i);
          end
      end
   end
   imfhat = fn_impulse(Bhat,Q*A0hatinv,nn);
      %In the form that is congenial to RATS
   imf3hat=reshape(imfhat,size(imfhat,1),nvar,nvar);
      %imf3hat: row--steps, column--nvar responses, 3rd dimension--nvar shocks

   %=== Responses to a basic gov spending shock.
   % gov_exp>0 and gdp > 0 for the irs periods.  
   indx_both_g_gdp = 0;
   if (indx_both_g_gdp)
      a = (imf3hat(1:irs,1,1) > 0) .* (imf3hat(1:irs,3,1) > 0);
   else
      a = (imf3hat(1:irs,1,1) > 0);
   end      
   if (max(a)==0)
      %--- Swiching the sign of the shock.
      if (indx_both_g_gdp)
         am = (imf3hat(1:irs,1,1) <0) .* (imf3hat(1:irs,3,1) < 0);
      else
         am = (imf3hat(1:irs,1,1) <0);
      end   
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(1,:) = -Q(1,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to a basic gov tax shock.
   % gov_tax<0 for the irs periods.
   a = (imf3hat(1:irs,2,2) < 0) .* (imf3hat(1:irs,3,2) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (imf3hat(1:irs,2,2) > 0) .* (imf3hat(1:irs,3,2) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(2,:) = -Q(2,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to a normal business cycle shock.
   % gov_exp >0, gov_tax>0, and gdp>0 for the irs periods.
   %a = (imf3hat(1:irs,1,3) > 0) .* (imf3hat(1:irs,2,3) > 0) .* (imf3hat(1:irs,3,3) > 0);
   a = (imf3hat(1:irs,2,3) > 0) .* (imf3hat(1:irs,3,3) > 0);  %Mountford and Uhlig's article.
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      %am = (imf3hat(1:irs,1,3) < 0) .* (imf3hat(1:irs,2,3) < 0) .* (imf3hat(1:irs,3,3) < 0);
      am = (imf3hat(1:irs,2,3) < 0) .* (imf3hat(1:irs,3,3) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Q(3,:) = -Q(3,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end


   %--- Terminating condition: all restrictions are satisfied.
   control=1;
end
