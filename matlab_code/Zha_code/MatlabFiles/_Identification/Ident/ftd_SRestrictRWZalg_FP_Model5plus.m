function [Q, indx_fail] = ftd_SRestrictRWZalg_FP_Model5plus(A0hatinv,Bhat,nvar,lags,irs,ndraws_terminate)
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
%   irs = a vector of number of periods in which sign restrictions are imposed (vectors may correspond to different or same equations).
%         Note that the last element represent the maximum number of periods to be consider 
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
%    Shocks are in the following order: 1: unanticipated spending shock; 2: business shock; 3: anticipated spending shock.



nn = [nvar lags max(irs)];
control=0;           
%fail_shocks = zeros(nvar,1);  %0: no fail (success in satisfying sign restrictions); 1: fail (not satisfying sign restrictions)   
cnt = 0;
indx_fail = 0;

%Q=eye(nvar);

while (control==0)
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
      for i=1:nvar
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
   whichshock = 1;
   %--- gov_exp>0 and gdp > 0 for the irs(1) periods.
   amax = max([max(imf3hat(1:irs(1),1,whichshock) > 0); max(imf3hat(1:irs(1),3,whichshock) > 0)]);
   if (amax==0)   %The restrictions are satisfied.  
      %--- Swiching the sign of the shock and normalize.
      Q(whichshock,:) = -Q(whichshock,:);

      %%?????? Debugging.  The following lines should be got rid of later on.
      %disp('******** If no error pringing, we are fine. *******') 
      %am = (imf3hat(1:irs(1),1,whichshock) < 0) .* (imf3hat(1:irs(1),3,whichshock) < 0);
      %if (min(am)==0) 
      %   error('This canNOT happen.  Some bugs exist')
      %end
   else
      amin = min([min(imf3hat(1:irs(1),1,whichshock) > 0); min(imf3hat(1:irs(1),3,whichshock) > 0)]);       
      if (amin==0)
         continue;   %The restrictions are not satisfied.  Go to the beginning to redraw.
         %fail_shocks(1) = 1;  %The restrictions are not satisfied.  Go to the beginning to redraw.
      end   
   end
   
   %=== Responses to a normal business cycle shock.
   whichshock = 2;
   %%--- gov_exp >0, gov_tax>0, and gdp>0 for the irs periods.
   %a = (imf3hat(1:irs,1,3) > 0) .* (imf3hat(1:irs,2,3) > 0) .* (imf3hat(1:irs,3,3) > 0);
   %--- gov_tax>0 and gdp>0 for the irs(1) periods.
   amax = max([max(imf3hat(1:irs(1),2,whichshock) > 0); max(imf3hat(1:irs(1),3,whichshock) > 0)]);  %Mountford and Uhlig's article.
   if (amax==0)  %The restrictions are satisfied.
      %--- Swiching the sign of the shock and normalize.
      Q(whichshock,:) = -Q(whichshock,:);
   else
      amin = min([min(imf3hat(1:irs(1),2,whichshock) > 0); min(imf3hat(1:irs(1),3,whichshock) > 0)]);  %Mountford and Uhlig's article.
      if (amin==0)       
         continue;   %The restrictions are not satisfied.  Go to the beginning to redraw.
         %fail_shocks(whichshock) = 1;
      end
   end
   %continue;   %The restrictions are not satisfied.  Go to the beginning to redraw.
   %Qrow2 = Q(whichshock,:); 

   %=== Responses to an anticipated spending shock.
   whichshock = 3;
   %--- gov_exp <= 0 for irs(2) periods, gov_exp > 0 for additional irs(3) periods, and gdp>0 for irs(1) periods.
   nstps_tot = irs(2)+irs(3);
   amax = max([max(imf3hat(1:irs(2),1,whichshock) <= 0); max(imf3hat((irs(2)+1):nstps_tot,1,whichshock) > 0); max(imf3hat(1:irs(1),3,whichshock) > 0)]);  %Mountford and Uhlig's article.
   if (amax==0)  %The restrictions are satisfied.  
      %--- Swiching the sign of the shock and normalize.
      Q(whichshock,:) = -Q(whichshock,:); 
   else
      amin = min([min(imf3hat(1:irs(2),1,whichshock) <= 0); min(imf3hat((irs(2)+1):nstps_tot,1,whichshock) > 0); min(imf3hat(1:irs(1),3,whichshock) > 0)]);  %Mountford and Uhlig's article.   
      if (amin==0)       
         continue;   %The restrictions are not satisfied.  Go to the beginning to redraw.
         %fail_shocks(whichshock) = 1;
      end   
   end
   %continue;   %The restrictions are not satisfied.  Go to the beginning to redraw.
   %Qrow3 = Q(whichshock,:); 

   %%=== Responses to a basic gov tax shock.
   %% gov_tax<0 for the irs periods.
   %a = (imf3hat(1:irs,2,2) < 0) .* (imf3hat(1:irs,3,2) > 0);
   %if (max(a)==0)
   %   %--- Swiching the sign of the shock and normalize.
   %   am = (imf3hat(1:irs,2,2) > 0) .* (imf3hat(1:irs,3,2) < 0);
   %   if (min(am)==0)
   %      continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
   %   else
   %      %--- Normalizing according to the switched sign.
   %      Q(2,:) = -Q(2,:);
   %   end
   %elseif (min(a)==0)
   %   continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   %end
   
   %--- Terminating condition: all restrictions are satisfied.
   control = 1;
end
