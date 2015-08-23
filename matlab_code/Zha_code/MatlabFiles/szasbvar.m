function [Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
             idmat0,idmatpp,H0invmulti,Hpinvmulti,xxhp] = szasbvar(idfile,q_m,lags,xdgel,mu)
% [Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
%      idmat0,idmatpp,H0invmulti,Hpinvmulti] = szasbvar(idfile,q_m,lags,xdgel,mu)
%
%    Estimating the Bayesian VAR of Sims and Zha with asymmetric prior (as)
%
% idfile:  Identification filename with rows corresponding to equations.
%              But all the output will be transposed to columns
%              corresponding to equations.
% q_m:  quarter or month
% lags: the maximum length of lag
% nhp:  number of haperparameters
% xdgel: the general matrix of the original data (no manipulation involved)
%                    with sample size including lags
% mu: 6-by-1 vector of hyperparameters (the following numbers for Atlanta Fed's forecast)
%       mu(1): overall tightness and also for A0;  (0.57)
%       mu(2): relative tightness for A+;  (0.13)
%       mu(3): relative tightness for the constant term;  (0.1)
%       mu(4): tightness on lag decay;  (1)
%       mu(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
% --------------------
% --------------------
% Gb: cell(nvar,1). Each cell, when postmultiplied by A0(:,i), is used to compute A+(:,i)
%      where A+ is k-by-m, where k--ncnoef, m--nvar, if asymmetric prior.
%      In particular, Aplus(:,i)=Gb{i}*A0(:,i), so Bh = Aplus/A0
% Sbd: cell(nvar,1). Sbd=diag(S(1), ..., S(m)).  Already divided by 'fss' for the
%        posterior of a0 or A0(:) in Waggoner and Zha when one has asymmetric prior.
%        Note,"bd" stands for block diagonal.
% Bh: reduced form B, k(mp+1)-by-m, column--equation. Bh=NaN if asymmetric prior.
%     In the form of y=X*B+U where B=[B1|B2| ...|Bp|C]
%       where Bi: m-by-m (i=1,..,p--lag length), C: 1-by-ml, constant.
% SpH:  divided by T, the final S in the Waggoner-Zha exponential part of p(A0|y)
%               SpH=NaN if asymmetric prior
% fss:  in-sample-size (for forecasting). Including dummy observations
% ndobs: number of dummy observations, ndobs=nvar+1
% phi:  X in the form of y = X*B+U. Row: nSmaple-lags+ndobs. Column: ncoef
% y:    y in the form y = X*B+U. Including the dummy observations too,
%         T-lags+ndobs-by-nvar.
% nvar: number of variables
% ncoef: number of coefficients in *each* equation. RHS coefficients only, nvar*lags+1
% xxhpc: chol(X'X+inv(H_p_tilde)): upper triangular but its transpose
%                                      is lower triagular Choleski
% a0indx: the index number for free parameters in A0.  Column meaning equation
% na0p:   number of free parameters in A0
% idmat0:  identification matrix for A0 with asymmetric prior;  column -- equation.
% idmatpp: asymmetric prior variance for A+ without constant; column -- equation.
% H0invmulti: nvar-by-nvar-by-nvar; inv(H0) for different equations under asymmetric prior
% Hpinvmulti: ncoef-by-ncoef-by-nvar; inv(H+) for different equations under asymmetric prior
% xxhp: ncoef-by-nocef, X'X+inv(H_p_tilde)
%
% Revisions by CAS 9/27/96:  If we transmit idmat, rather than a0indx, we can scale
%   the elements of idmat0 or idmatp so that it carries information about relative
%   prior variances.
% Revsions by TZ, 10/13/96:  Efficiency improvement by streamlining the previous
%   code according to the general setup in Sims and Zha's IER and according
%   to Zha's notes Forecast (0) (p.3) and Forecast (1) (p.9).
% Quick Revisions: May 2003.  See H1p_1 on lines 406-409.
%
% Copyright (c) December 1997 by C.A. Sims and T. Zha,
%
%  NOTE1: "nSample" is something I deleted as an input recently, so it may not
%    be compatible with old programs, 10/15/98 by TZ.
%  NOTE2: I put "mu" as an input arg and delete "nhp" as an input arg, so it may
%    not be compatible with old programs, 03/06/99 by TZ
%  NOTE3:  added three output arguments: H0invmulti and Hpinvmulti and xxhp.  9/27/99

% Copyright (C) 1997-2012 Christopher A. Sims and Tao Zha 
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


%
%@@@ Prepared for Bayesian VAR of Sims and Zha
%
%* Load identification and obtain idmat0
eval(['load ' idfile '.prn -ascii']);
eval(['idmat0=' idfile ';']);
idmat0=idmat0';   % so that each column corresponds to an equation
a0indx=find(idmat0);    % column meaning equation
na0p = length(a0indx);  % number of free parameters in A0
nhp = length(mu);    % total number of hyperparameters
nfp = na0p + nhp;     % total number of free parameters (including hyperparameters)
[nvar,neqn]=size(idmat0);
%
idmatpp = ones(nvar*lags,nvar);    % pp: plus without constant.  Column -- equation
%>>>>>> B: asymmetric prior variance for idmatpp <<<<<<<<
%
%for i = 1:lags
%   rowif = (i-1)*nvar+1;
%   rowil = i*nvar;
%		idmatw0 = 0.5;   % weight assigned to idmat0 in the formation of idmatpp
%	if (i==1)
%		idmatpp(rowif:rowil,:)=(1-idmatw0)*ones(nvar)+idmatw0*idmat0;  % first lag
%		                 % note:  idmat1 is already transposed.  Column -- equation
%	else
%   	%idmatpp(rowif:rowil,1:nvar) = (1-idmatw0)*ones(nvar)+idmatw0*idmat0;
%                % <<<<<<< toggle +
%                % Note: already transposed, since idmat0 is transposed.
%				     % Meaning: column implies equation
%   	idmatpp(rowif:rowil,1:nvar) = ones(nvar);
%                % >>>>>>> toggle -
%	end
%end
%
%>>>>>> E: asymmetric prior variance for idmatpp <<<<<<<<


%$$$ available computations
%
%@@@ original
% total number of the sample under study, including lags, etc. -- original sample size
nSample = size(xdgel,1); % the sample size (including lags, of course)
sb = lags+1;   % original beginning without dummies
sl = nSample;       % original last period without dummies
ssp = nSample - lags;  % original number of observations
%
ndobs=nvar+1;         % number of dummy observations
ncoef = nvar*lags+1;  % number of coefficients in *each* equation, RHS coefficients only.
%*** initializing for global variables
%%GlbAph=zeros(nvar,ncoef);
%Aplus=zeros(ncoef,nvar);
Gb = cell(nvar,1);      % Storing potential reduced-form Bh (k-by-m) or prepared
								%  for computing Aplus (k-by-m), where k--ncnoef, m--nvar
Sbd = cell(nvar,1);    % Storing for diag(S(1), ..., S(m)) for the posterior of a0
							  %   or vec(A0) when one has asymmetric prior.  Note,
							  %   "bd" stands for block diagonal.
SpH = ones(nvar,nvar);
%* flp:  forecast last period (with dummies);
%* fbp: forecast beginning period (with dummies).
flp = sl+ndobs;
%fbp = ndobs+lags+1;       % <<>> true begining by skipping dummies
fss = flp-lags;    % forecast sample size (with dummies).

% ** hyperparameters
%  mu = zeros(nhp,1);
%  mu(1) = 0.57;
%  mu(2) = 0.13;
%  mu(3) = 0.1;
%  mu(4) = 1;
%  mu(5) = 5;  %10;
%  mu(6) = 5;  %10;
%
% mu(1): overall tightness and also for A0;
% mu(2): relative tightness for A+;
% mu(3): relative tightness for the constant term;
% mu(4): tightness on lag decay;
% mu(5): weight on nvar sums of coeffs dummy observations (unit roots);
% mu(6): weight on single dummy initial observation including constant (cointegration, unit
%                                      roots, and stationarity).


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
% ** now run the VAR with
%
% **** nvar prior dummy observations with the sum of coefficients
% ** construct X for Y = X*B + U where phi = X: (T-lags)*k, Y: (T-lags)*nvar
% **    columns: k = # of [nvar for 1st lag, ..., nvar for last lag, const]
% **    Now, T=T+ndobs -- added with "ndobs" dummy observations
%
phi = zeros(fss,ncoef);
const = ones(fss,1);
const(1:nvar) = zeros(nvar,1);
phi(:,ncoef) = const;      % the first nvar periods: no or zero constant!
%
xdgelint = mean(xdgel(1:lags,:),1); % mean of the first lags initial conditions
%* Dummies
for k=1:nvar
   for m=1:lags
      phi(ndobs,nvar*(m-1)+k) = xdgelint(k);
      phi(k,nvar*(m-1)+k) = xdgelint(k);
      % <<>> multiply hyperparameter later
   end
end
%* True data
for k=1:lags
   phi(ndobs+1:fss,nvar*(k-1)+1:nvar*k) = xdgel(sb-k:sl-k,:);
   % row: T-lags; column: [nvar for 1st lag, ..., nvar for last lag, const]
   %                     Thus, # of columns is nvar*lags+1 = ncoef.
end
%
% ** Y with "ndobs" dummies added
y = zeros(fss,nvar);
%* Dummies
for k=1:nvar
   y(ndobs,k) = xdgelint(k);
   y(k,k) = xdgelint(k);
   % multiply hyperparameter later
end
%* True data
y(ndobs+1:fss,:) = xdgel(sb:sl,:);



%
% *** specify the prior for each equation separately, SZ method,
% ***
% *** obtaining the residuals for the univariate processes:
% ***                                   a total of "nvar" equations
% *** obtain the posterior peak of A0 which is Sims iden1 (Sims 1986)
% ** get the residuals from univariate regressions.
%
sgh = zeros(nvar,1);        % square root
sgsh = sgh;              % square
yu = xdgel;
C = ones(nSample,1);
for k=1:nvar
   [Bk,ek,junk1,junk2,junk3,junk4] = sye([yu(:,k) C],lags);
   clear Bk junk1 junk2 junk3 junk4;
   sgsh(k) = ek'*ek/ssp;
   %% sgh(k) = ek'*ek/(fss-7);  to match "sqrt(ess)" in RATS univariate regression
   sgh(k) = sqrt(sgsh(k));
end
% ** prior variance for alpha0, same for all equations!!!
al0b = zeros(nvar*neqn,1);    % prior mean for all A0
sg0bid = zeros(nvar,1);  % Sigma0_bar diagonal only
for j=1:nvar
   sg0bid(j) = 1/sgsh(j);    % sgsh = sigmai^2
end
% ** prior variance for alpha_plus, same for all equations
sgpbid = zeros(ncoef,1);     % Sigma_plus_bar, diagonal
for i = 1:lags
   if (q_m==12)
      lagdecay = a*exp(b*i);
   end
   %
   for j = 1:nvar
      if (q_m==12)
         % exponential decay to match quarterly decay
         sgpbid((i-1)*nvar+j) = lagdecay^2/sgsh(j);
         %%sgpbid((i-1)*nvar+j) = (1/i)^2/sgsh(j);
      elseif (q_m==4)
         sgpbid((i-1)*nvar+j) = (1/i^mu(4))^2/sgsh(j);
      else
			error('Incompatibility with lags, check the possible errors!!!')
         %warning('Incompatibility with lags, check the possible errors!!!')
         %return
      end
   end
end
%%A0bd = sqrt(sg0bid);
% commented out by T.A. Zha, 10/3/96, to avoid double-counting of scaling and the problem
%    of potential multiple peaks.
A0bd = zeros(nvar,1);
A0b = sparse(1:nvar,1:nvar,A0bd,nvar,nvar);
A0b = A0b';    % making a column in A0b correspond to an equation

%


%=================================================
%   Computing the (prior) covariance matrix for the posterior of A0, no data yet
%=================================================
%
%
% ** weight prior dummy observations
%phi(1:nvar,:) = (mu(5)^2/mu(1)^2)*phi(1:nvar,:);
%y(1:nvar,:) = (mu(5)^2/mu(1)^2)*y(1:nvar,:);
%phi(ndobs,:) = mu
% modified by CAS 8/6/96.  The dummy obs weights are naturally in std units, not var units.
%
phi(1:nvar,:) = 1*mu(5)*phi(1:nvar,:);    % standard Sims and Zha prior
y(1:nvar,:) = mu(5)*y(1:nvar,:);      % standard Sims and Zha prior
%----- The following prior designed for GLZ
%phi(1,:) = 1.00*mu(5)*phi(1,:);
%phi(2:nvar,:) = 1.02*mu(7)*phi(2:nvar,:);
%y(1,:) = mu(5)*y(1,:);
%y(2:nvar,:) = mu(7)*y(2:nvar,:);
%----------------------------
%
phi(nvar+1,:) = mu(6)*phi(nvar+1,:);
y(nvar+1,:) = mu(6)*y(nvar+1,:);


% ** set up the conditional prior variance sg0bi and sgpbi.
sg0bida = mu(1)^2*sg0bid;
sgpbida = mu(1)^2*mu(2)^2*sgpbid;
sgpbida(ncoef) = mu(1)^2*mu(3)^2;
sgppbd = sgpbida(nvar+1:ncoef);    % corresponding to A++, in a Sims-Zha paper
%%sgppbdi = 1./sgppbd;
%sgppb = diag(sgppbd);

Hptd = zeros(ncoef);
Hptdi=Hptd;
%%Hptd(nvar+1:ncoef,nvar+1:ncoef)=diag(sgppbd);
%%Hptdi(nvar+1:ncoef,nvar+1:ncoef)=diag(sgppbdi);
Hptd(ncoef,ncoef)=sgppbd(ncoef-nvar);
Hptdi(ncoef,ncoef)=1/sgppbd(ncoef-nvar);
             % condtional on A0i, H_plus_tilde

% **  some common terms
[u d v]=svd(phi,0); %trial
% xtx = phi'*phi; %trial
vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
xtx=vd*vd';
% M = eye(fss) - phi*(xtx\phi'); % not used except in line below %trial
yu = y'*u; %trial
cxyxx=yu*yu'; %trial
yty=y'*y;
ymy=yty-cxyxx; %trial
%ymy = y'*M*y;
%cxy = phi'*y; %trial
cxy = vd*yu'; %trial
cyx = cxy';
%cxpy = xtx\cxy; %not used further except in line below
%cxyxx = cxy'*cxpy;

% ** estimate A+_hat conditional on A0, ncoef*nvar, but not using full prior
%GlbAph = cxpy*A0h;


%=================================================
%   Computing the final covariance matrix (S1,...,Sm) for the posterior of A0,
%      with the data, and final Gb=(G1,...,Gm) for A+ if asymmetric prior or for
%      B if symmetric prior for A+
%=================================================
%
H0invmulti=zeros(nvar,nvar,nvar);  % inv(H0) for different equations under asymmetric prior
Hpinvmulti=zeros(ncoef,ncoef,nvar);  % inv(H+) for different equations under asymmetric prior
for i = 1:nvar
   stri = find(idmat0(:,i));
	         % CAS change 9/24/96, to avoid slim chance of an exact zero in x vector
            %     at some point during an optimization.
   %strm = length(stri);  % TZ, 11/27/97, no use any more

   %A0hb = A0h(stri,i)-A0b(stri,i);

   % ** set up the conditional prior variance sg0bi and sgpbi.
   % sg0bd =  zeros(stri,1); Commented out by CAS, 9/24/96.  Looks like this
	%   meant to have strm where it has stri, and in any case it is "overwritten" below.

   %------------------------------
   % Introduce prior information on which variables "belong" in various equations.
   % In this first trial, we just introduce this information here, in a model-specific way.
   % Eventually this info has to be passed parametricly.  In our first shot, we just damp down
   % all coefficients except those on the diagonal.
   factor0=idmat0(:,i);
   sg0bd =  sg0bida;  % added by TZ, 2/26/98.  I think at each equation, sg0bd must
                      % be refreshed.  Note, this only works for the prior variance Sg(i)
                      % of a0(i) being diagonal.  If the prior variance Sg(i) is not
                      % diagonal, we have to work on inv(Sg(i)) or sg0bdi directly.
   sg0bd(stri) = sg0bida(stri).*factor0(stri);
   sg0bdi = 1./sg0bd;
   %
   factor1=idmatpp(1:nvar,i);
   sg1bd = sgpbida(1:nvar).*factor1;
   sg1bdi = 1./sg1bd;
   %
   if lags>1
      factorpp=idmatpp(nvar+1:ncoef-1,i);
      sgpp_cbd = sgppbd(1:ncoef-nvar-1) .* factorpp;
      sgpp_cbdi = 1./sgpp_cbd;
   end

   % ** set up the unconditional prior variance on A0i and A+i
   %XX = zeros(nvar,strm);
   %XX(stri,:) = eye(strm);        % XX in Gelman, etel, on p479
   %  Commented out by TZ, 11/27/97. Streamline. No need for these any more.


   % * final conditional prior variance on A+i give that on A0i, and
   % *    unconditional variance on A0+
   H0td = diag(sg0bd);    % unconditional
   % H_~: td: tilde for ~
   % ** inverse and chol decomp
   H0tdi = diag(sg0bdi);
   H0invmulti(:,:,i)=H0tdi;

   %
   Hptd(1:nvar,1:nvar)=diag(sg1bd);
   Hptdi(1:nvar,1:nvar)=diag(sg1bdi);
   if lags>1
      Hptd(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbd);
      Hptdi(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbdi);
             % condtional on A0i, H_plus_tilde
   end
   Hpinvmulti(:,:,i)=Hptdi;

   % common terms
   xxhp = xtx+Hptdi;
	%Lxxhpc = chol(inv(xxhp))';   % lower triangular
	xxhpc = chol(xxhp);   % upper triangular but its transpose is lower triagular Choleski
	%A0hbH = A0hb'*H0tdi;
   %====== The following two lines seem incorrect.  The third line is supposed to correc the mistake.  May 2003.
   %  H1p_1 = zeros(ncoef,nvar);
   %  H1p_1(1:nvar,:) = Hptdi(1:nvar,1:nvar);
   H1p_1 = Hptdi(:,1:nvar);
   %%Hm = (cyx+H1p_1')*(xxhp\(cxy+H1p_1));
   Hm1 = (cyx+H1p_1')/xxhp;   % if symmetric prior, Bh=Hm1' -- reduced form k-by-m.
	Gb{i} = Hm1';      % k-by-m, where k--ncoef, m--nvar.
			  % If asymmetric prior, Gb is used to compute A+ where A+(i) = Gb(i)*a0(i)
			  %   where a0(i) is m-by-1 for the ith column of A0 or the ith equation.
   %Hm2 = cxy+H1p_1;
   Hm = Hm1*(cxy+H1p_1);
   %GlbAph(i,:) = A0h(stri,i)'*XX'*Hm1;     % alpha_plus_*_transpose
   %%alpMpyh = A0h(stri,i)'*XX'*(yty+Hptdi(1:nvar,1:nvar)-Hm)*XX;
   %Hss = XX'*(yty+Hptdi(1:nvar,1:nvar)-Hm)*XX;
	Hss = yty+Hptdi(1:nvar,1:nvar)-Hm;
   Sbd{i} = (H0tdi+Hss)./fss;       % m-by-m, where m--nvar.  If asymmetric prior, Scd is
				% usde to form diag(Sbd{1}, ..., Sbd{m}) for a0 or vec(A0) in the
            % posterior of A0.  Note, "bd" stands for block diagonal and divided by
            % "fss" already to make it compatible with SpH used in "a0lhfun" and make it
            % it easier according to the notations by Waggoner and Zha
end


%%=========================================
%% check if there is asymmetric prior on both A+ and A0
%%=========================================
%
spindx=1;    % A+ index for asymmetric prior
s0indx=1;    % A0 index for asymmetric prior
for i=2:nvar
	diffGb = Gb{i}-Gb{i-1};
	diffSbd = Sbd{i}-Sbd{i-1};
	if ~any(any(diffGb))
		spindx=spindx+1;
	end
   if ~any(any(diffSbd))
		s0indx=s0indx+1;
   end
end
%
if (spindx==nvar)
	Bh = Gb{1};            % reduced-form parameter. Note: does not depend on A0
else
	Bh = NaN;
end
%
if (s0indx==nvar)
	SpH = (H0tdi+Hss)/fss;   % the final S in the exponential part of p(A0|y)
		% divided by nobs in order to make Choleski decomp without optimization
		% or in the form conformable to Waggoner and Zha
else
	SpH=NaN;
end





%----------------------------------------------
%

% ***
% *** Form the inverse of covarinace to draw from t- or Gaussian distribution
% ***
%SpHs = SpH*fss;     % for the Wishart draw
%SpHsc = chol(SpHs);     % upper triangular where Sphsc' is
%	                        %  lower triangular Choleski, for Wishart draw
%SpHsic = chol(inv(SpHs))';     % inverse Choleski decomposition -- lower triangular
%A0hin = chol(SpH);    % upper triangular, inverse of A0h, each column
                         %   corresponds to an equation.


%@@@ The following can be used for other purposes such as forecasting
%
%swish = A0hin';       % each row corresponds to an equation
%A0h = inv(A0hin);     % each column corresponds to an equation
%xa0 = A0h(a0indx);

%%*** form Bh (ncoef*nvar)
%Aplus = Hm1t*A0h;      % estimate of Aplus -- the same as Astrar
%Bhp = Hm1t;
