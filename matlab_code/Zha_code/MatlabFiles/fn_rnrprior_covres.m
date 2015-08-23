function [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] ...
                      = fn_rnrprior_covres(nvar,q_m,lags,xdgel,mu,hpmsmd,indxmsmdeqn,nexo,asym0,asymp)
% [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] ...
%                      = fn_rnrprior_covres(nvar,q_m,lags,xdgel,mu,hpmsmd,indxmsmdeqn,nexo,asym0,asymp)
%
% More general than fn_rnrprior.m because when hpmsmd=0, fn_rnrprior_covres() is the same as fn_rnrprior().
%    Allows for prior covariances for the MS and MD equations to achieve liquidity effects.
%    Exports random Bayesian prior of Sims and Zha with asymmetric rior (but no linear restrictions yet)
%    See Waggoner and Zha's Gibbs sampling paper and TVBVAR NOTES p. 71k.0.
%
% nvar:  number of endogenous variables
% q_m:  quarter or month
% lags: the maximum length of lag
% xdgel: the general matrix of the original data (no manipulation involved)
%             with sample size including lags.  Used only to get variances of residuals for
%             the scaling purpose; NOT used for mu(5) and mu(6).
% mu: 6-by-1 vector of hyperparameters (the following numbers for Atlanta Fed's forecast), where
%          mu(5) and mu(6) are NOT used here.  See fn_dataxy.m for using mu(5) and mu(6).
%       mu(1): overall tightness and also for A0;  (0.57)
%       mu(2): relative tightness for A+;  (0.13)
%       mu(3): relative tightness for the constant term;  (0.1).  NOTE: for other
%               exogenous terms, the variance of each exogenous term must be taken into
%               acount to eliminate the scaling factor.
%       mu(4): tightness on lag decay;  (1)
%       mu(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
%       NOTE: for this subdirectory, mu(5) and mu(6) are not used.  See fn_dataxy.m for using mu(5) and mu(6).
% hpmsmd: 2-by-1 hyperparameters with -1<h1=hpmsmd(1)<=0 for the MS equation and 0<=h2=hpmsmd(2)<1 the MD equation.  Consider a1*R + a2*M.
%          The term h1*var(a1)*var(a2) is the prior covariance of a1 and a2 for MS, equivalent to penalizing the same sign of a1 and a2.
%          The term h2*var(a1)*var(a2) is the prior covariance of a1 and a2 for MD, equivalent to penalizing opposite signs of a1 and a2.
%          This will give us a liquidity effect.
% indxmsmdeqn: 4-by-1 index for the locations of the MS and MD equation and for the locations of M and R.
%               indxmsmdeqn(1) for MS and indxmsmdeqn(2) for MD.
%               indxmsmdeqn(3) for M and indxmsmdeqn(4) for R.
% nexo:  number of exogenous variables (if not specified, nexo=1 (constant) by default).
%         The constant term is always put to the last of all endogenous and exogenous variables.
% asym0: nvar-by-nvar asymmetric prior on A0.  Column -- equation.
%        If ones(nvar,nvar), symmetric prior;  if not, relative (asymmetric) tightness on A0.
% asymp: ncoef-1-by-nvar asymmetric prior on A+ bar constant.  Column -- equation.
%        If ones(ncoef-1,nvar), symmetric prior;  if not, relative (asymmetric) tightness on A+.
% --------------------
% Pi: ncoef-by-nvar matrix for the ith equation under random walk.  Same for all equations
% H0multi: nvar-by-nvar-by-nvar; H0 for different equations under asymmetric prior
% Hpmulti: ncoef-by-ncoef-by-nvar; H+ for different equations under asymmetric prior
% H0invmulti: nvar-by-nvar-by-nvar; inv(H0) for different equations under asymmetric prior
% Hpinvmulti: ncoef-by-ncoef-by-nvar; inv(H+) for different equations under asymmetric prior
%
% Tao Zha, February 2000.  Revised, September 2000, February 2003.
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



if nargin==7, nexo=1; end
ncoef = nvar*lags+nexo;  % number of coefficients in *each* equation, RHS coefficients only.

H0multi=zeros(nvar,nvar,nvar);  % H0 for different equations under asymmetric prior
Hpmulti=zeros(ncoef,ncoef,nvar);  % H+ for different equations under asymmetric prior
H0invmulti=zeros(nvar,nvar,nvar);  % inv(H0) for different equations under asymmetric prior
Hpinvmulti=zeros(ncoef,ncoef,nvar);  % inv(H+) for different equations under asymmetric prior

%*** Constructing Pi for the ith equation under the random walk assumption
Pi = zeros(ncoef,nvar);   % same for all equations
Pi(1:nvar,1:nvar) = eye(nvar);   % random walk

%
%@@@ Prepared for Bayesian prior
%
%
% ** monthly lag decay in order to match quarterly decay: a*exp(bl) where
% **  l is the monthly lag.  Suppose quarterly decay is 1/x where x=1,2,3,4.
% **  Let the decay of l1 (a*exp(b*l1)) match that of x1 (say, beginning: 1/1)
% **  and the decay of l2 (a*exp(b*l2)) match that of x2 (say, end: 1/5),
% **  we can solve for a and b which are
% **      b = (log_x1-log_x2)/(l1-l2), and a = x1*exp(-b*l1).
if q_m==12
   l1 = 1;   % 1st month == 1st quarter
   xx1 = 1;   % 1st quarter
   l2 = lags;   % last month
   xx2 = 1/((ceil(lags/3))^mu(4));   % last quarter
   %xx2 = 1/6;   % last quarter
   % 3rd quarter:  i.e., we intend to let decay of the 6th month match
   %    that of the 3rd quarter, so that the 6th month decays a little
   %    faster than the second quarter which is 1/2.
   if lags==1
      b = 0;
   else
      b = (log(xx1)-log(xx2))/(l1-l2);
   end
   a = xx1*exp(-b*l1);
end



%
% *** specify the prior for each equation separately, SZ method,
% ** get the residuals from univariate regressions.
%
sgh = zeros(nvar,1);        % square root
sgsh = sgh;              % square
nSample=size(xdgel,1);  % sample size-lags
yu = xdgel;
C = ones(nSample,1);
for k=1:nvar
   [Bk,ek,junk1,junk2,junk3,junk4] = sye([yu(:,k) C],lags);
   clear Bk junk1 junk2 junk3 junk4;
   sgsh(k) = ek'*ek/(nSample-lags);
   sgh(k) = sqrt(sgsh(k));
end
% ** prior variance for A0(:,1), same for all equations!!!
sg0bid = zeros(nvar,1);  % Sigma0_bar diagonal only for the ith equation
for j=1:nvar
   sg0bid(j) = 1/sgsh(j);    % sgsh = sigmai^2
end
% ** prior variance for lagged and exogeous variables, same for all equations
sgpbid = zeros(ncoef,1);     % Sigma_plus_bar, diagonal, for the ith equation
for i = 1:lags
   if (q_m==12)
      lagdecay = a*exp(b*i*mu(4));
   end
   %
   for j = 1:nvar
      if (q_m==12)
         % exponential decay to match quarterly decay
         sgpbid((i-1)*nvar+j) = lagdecay^2/sgsh(j);  % ith equation
      elseif (q_m==4)
         sgpbid((i-1)*nvar+j) = (1/i^mu(4))^2/sgsh(j);  % ith equation
      else
			error('Incompatibility with lags, check the possible errors!!!')
         %warning('Incompatibility with lags, check the possible errors!!!')
         %return
      end
   end
end
%


%=================================================
%   Computing the (prior) covariance matrix for the posterior of A0, no data yet
%=================================================
%
%
% ** set up the conditional prior variance sg0bi and sgpbi.
sg0bida = mu(1)^2*sg0bid;   % ith equation
sgpbida = mu(1)^2*mu(2)^2*sgpbid;
sgpbida(ncoef-nexo+1:ncoef) = mu(1)^2*mu(3)^2;
          %<<>> No scaling adjustment has been made for exogenous terms other than constant
sgppbd = sgpbida(nvar+1:ncoef);    % corresponding to A++, in a Sims-Zha paper

Hptd = zeros(ncoef);
Hptdi=Hptd;
Hptd(ncoef,ncoef)=sgppbd(ncoef-nvar);
Hptdinv(ncoef,ncoef)=1./sgppbd(ncoef-nvar);
             % condtional on A0i, H_plus_tilde


if nargin<9   % the default is no asymmetric information
   asym0 = ones(nvar,nvar);  % if not ones, then we have relative (asymmetric) tightness
   asymp = ones(ncoef-1,nvar);    % for A+.  Column -- equation
end

%**** Asymmetric Information
%asym0 = ones(nvar,nvar);  % if not ones, then we have relative (asymmetric) tightness
%asymp = ones(ncoef-1,nvar);    % pp: plus without constant.  Column -- equation
%>>>>>> B: asymmetric prior variance for asymp <<<<<<<<
%
%for i = 1:lags
%   rowif = (i-1)*nvar+1;
%   rowil = i*nvar;
%     idmatw0 = 0.5;   % weight assigned to idmat0 in the formation of asymp
%	if (i==1)
%     asymp(rowif:rowil,:)=(1-idmatw0)*ones(nvar)+idmatw0*idmat0;  % first lag
%		                 % note:  idmat1 is already transposed.  Column -- equation
%	else
%     %asymp(rowif:rowil,1:nvar) = (1-idmatw0)*ones(nvar)+idmatw0*idmat0;
%                % <<<<<<< toggle +
%                % Note: already transposed, since idmat0 is transposed.
%				     % Meaning: column implies equation
%     asymp(rowif:rowil,1:nvar) = ones(nvar);
%                % >>>>>>> toggle -
%	end
%end
%
%>>>>>> E: asymmetric prior variance for asymp <<<<<<<<


%=================================================
%   Computing the final covariance matrix (S1,...,Sm) for the prior of A0,
%      and final Gb=(G1,...,Gm) for A+ if asymmetric prior or for
%      B if symmetric prior for A+
%=================================================
%
for i = 1:nvar
   %------------------------------
   % Introduce prior information on which variables "belong" in various equations.
   % In this first trial, we just introduce this information here, in a model-specific way.
   % Eventually this info has to be passed parametricly.  In our first shot, we just damp down
   % all coefficients except those on the diagonal.

   %*** For A0
   factor0=asym0(:,i);
   sg0bd = sg0bida.*factor0;  %  Note, this only works for the prior variance Sg(i)
                      % of a0(i) being diagonal.  If the prior variance Sg(i) is not
                      % diagonal, we have to the inverse to get inv(Sg(i)).
   %sg0bdinv = 1./sg0bd;
   % *    unconditional variance on A0+
   H0td = diag(sg0bd);    % unconditional
   %=== Correlation in the MS equation to get a liquidity effect.
   if (i==indxmsmdeqn(1))
      H0td(indxmsmdeqn(3),indxmsmdeqn(4)) = hpmsmd(1)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
      H0td(indxmsmdeqn(4),indxmsmdeqn(3)) = hpmsmd(1)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
   else
      H0td(indxmsmdeqn(3),indxmsmdeqn(4)) = hpmsmd(2)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
      H0td(indxmsmdeqn(4),indxmsmdeqn(3)) = hpmsmd(2)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
   end
   H0tdinv = inv(H0td);
   %H0tdinv = diag(sg0bdinv);
   %
   H0multi(:,:,i)=H0td;
   H0invmulti(:,:,i)=H0tdinv;


   %*** For A+
   if ~(lags==0)  % For A1 to remain random walk properties
      factor1=asymp(1:nvar,i);
      sg1bd = sgpbida(1:nvar).*factor1;
      sg1bdinv = 1./sg1bd;
      %
      Hptd(1:nvar,1:nvar)=diag(sg1bd);
      Hptdinv(1:nvar,1:nvar)=diag(sg1bdinv);
      if lags>1
         factorpp=asymp(nvar+1:ncoef-1,i);
         sgpp_cbd = sgppbd(1:ncoef-nvar-1) .* factorpp;
         sgpp_cbdinv = 1./sgpp_cbd;
         Hptd(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbd);
         Hptdinv(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbdinv);
               % condtional on A0i, H_plus_tilde
      end
   end
   Hpmulti(:,:,i)=Hptd;
   Hpinvmulti(:,:,i)=Hptdinv;
end


