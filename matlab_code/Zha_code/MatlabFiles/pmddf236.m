function of = pmddf236(x,idmat0,idmatp,fss,nvar,ncoef,phi, ...
                          y,A0b,sg0bid,sgpbid)
%function of = pmddf232(x,idmat,fss,nvar,ncoef,phi, ...
%                          y,A0b,sg0bid,sgpbid)
% Leeper-Sims-Zha BVAR setup
% general program to setup A0 matrix and compute the likelihood
% requires
% x (parameter vector)
% a0indx (matrix indicating the free parameters in A0)
% fss (forecast sample size)
% nvar (number of variables)
% ncoef (number of coefficients in a single equation in A+)
% phi (r.h.s. observations in the system for A+)
% y (l.h.s. observations for A0)
% A0b (initial prior mean on A0)
% sg0bid (diagonal of Sigma0_bar, unweighted prior vcv on the parameters in i-th equation in A0)
% sgpbid (diagonal of Sigma+_bar, unweighted prior vcv on the parameters in i-th equation in A+)
%
% Revisions by CAS 9/27/96:  If we transmit idmat, rather than a0indx, we can scale the elements of
% idmat so that it carries information about relative prior variances.
%
% Revsions by TZ, 10/13/96:  efficiency improvement by streamlining the previous code
%     according to the general setup in Sims and Zha "Bayesian Methods for ...".
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

global vaha0 vaMMa0 GlbAph FRESHFUNCTION
% FRESHFUNCTION added 8/20/90 by CAS.  It signals the pmddg23 whether the global variables
% it uses have been refreshed by a new pmddf231 call since the last gradient call.
a0indx=find(idmat0);
% pre-allocate GlbAph to be zeros(nvar,lags*nvar+1=ncoef)

nhp = 0;                 % <<>> 4 hyperparameters
na0p = length(a0indx);    % <<>> number of A0 parameters
nfp = na0p+nhp;

% *** initializing
vaha0 = zeros(na0p,1);
vaMMa0 = zeros(na0p,1);

%
A0h = zeros(nvar,nvar);
A0h(a0indx) = x(1:na0p);
%A0h = reshape(a0,nvar,nvar); CAS 9/24/96.  The reshape is not necessary if a0 starts as nvar x nvar
               %  restrictions on A0

[A0hl,A0hu] = lu(A0h);

% ** hyperparameters
%%mu = x(na0p+1:nfp);
mu = ones(5,1);
mu(1) = 1;
%mu(1) = 2;
mu(2) = 0.2;
mu(3) = 1;
%mu(4)=1;
mu(4)=10;
mu(5) =1;
%mu(5) = 40;
% results from ...\dummy1\foreh\pmdd6.out
%
% mu(1): overall tightness and also for A0;
% mu(2): relative tightness for A+;
% mu(3): relative tightness for the constant term;
% mu(4): weight on single dummy initial observation including constant;
% mu(5): weight on nvar sums of coeffs dummy observations.

% ** weight prior dummy observations
%phi(1:nvar,:) = (mu(5)^2/mu(1)^2)*phi(1:nvar,:);
%y(1:nvar,:) = (mu(5)^2/mu(1)^2)*y(1:nvar,:);
%phi(ndobs,:) = mu
% modified by CAS 8/6/96.  The dummy obs weights are naturally in std units, not var units.
phi(1:nvar,:) = mu(5)*phi(1:nvar,:);
y(1:nvar,:) = mu(5)*y(1:nvar,:);
phi(nvar+1,:) = mu(4)*phi(nvar+1,:);
y(nvar+1,:) = mu(4)*y(nvar+1,:);

% ** set up the conditional prior variance sg0bi and sgpbi.
sg0bida = mu(1)^2*sg0bid;
sgpbida = mu(1)^2*mu(2)^2*sgpbid;
sgpbida(ncoef) = mu(1)^2*mu(3)^2;
sgppbd = sgpbida(nvar+1:ncoef);    % corresponding to A++, in an SZ paper
%%sgppbdi = 1./sgppbd;
%sgppb = diag(sgppbd);

%
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


%%%%%%%
%%% build the objective function below
%%%%%%%
%
t1f = diag(abs(A0hu));
t1f = log(t1f);
% * without the prior |A0|^k
tt1 = (-1.0)*fss*sum(t1f);

tt3 = 0;

disp(sprintf('Starting loop (of): %g',toc))

stril = 0;
Hptd = zeros(ncoef);
Hptdi=Hptd;
%%Hptd(nvar+1:ncoef,nvar+1:ncoef)=diag(sgppbd);
%%Hptdi(nvar+1:ncoef,nvar+1:ncoef)=diag(sgppbdi);
Hptd(ncoef,ncoef)=sgppbd(ncoef-nvar);
Hptdi(ncoef,ncoef)=1/sgppbd(ncoef-nvar);
             % condtional on A0i, H_plus_tilde
for i = 1:nvar
   stri = find(idmat0(:,i)); %CAS change 9/24/96, to avoid slim chance of an exact zero in x vector
                              %at some point during an optimization.
   strm = length(stri);

   A0hb = A0h(stri,i)-A0b(stri,i);

   % ** set up the conditional prior variance sg0bi and sgpbi.
   % sg0bd =  zeros(stri,1); Commented out by CAS, 9/24/96.  Looks like this meant to have
   %                         strm where it has stri, and in any case it is "overwritten" below.
   %------------------------------
   % Introduce prior information on which variables "belong" in various equations.
   % In this first trial, we just introduce this information here, in a model-specific way.
   % Eventually this info has to be passed parametricly.  In our first shot, we just damp down
   % all coefficients except those on the diagonal.
   factor0=idmat0(:,i);
   sg0bd = sg0bida(stri).*factor0(stri);
   sg0bdi = 1./sg0bd;
   %sg0b = diag(sg0bd);
   %
   factor1=idmatp(1:nvar,i);
   sg1bd = sgpbida(1:nvar).*factor1;
   sg1bdi = 1./sg1bd;
   %sg1b = diag(sg1bd);
   %
   factorpp=idmatp(nvar+1:ncoef-1,i);
   sgpp_cbd = sgppbd(1:ncoef-nvar-1) .* factorpp;
   sgpp_cbdi = 1./sgpp_cbd;


   % ** set up the unconditional prior variance on A0i and A+i
   XX = zeros(nvar,strm);
   XX(stri,:) = eye(strm);        % XX in Gelman, etel, on p479
   %
   % * final conditional prior variance on A+i give that on A0i, and
   % *    unconditional variance on A0+
   H0td = diag(sg0bd);    % unconditional
   % H_~: td: tilde for ~
   % ** inverse and chol decomp
   H0tdi = diag(sg0bdi);

   %
   Hptd(1:nvar,1:nvar)=diag(sg1bd);
   Hptdi(1:nvar,1:nvar)=diag(sg1bdi);
   Hptd(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbd);
   Hptdi(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbdi);
             % condtional on A0i, H_plus_tilde

   % common terms
   xxhp = xtx+Hptdi;
   A0hbH = A0hb'*H0tdi;
   H1p_1 = zeros(ncoef,nvar);
   H1p_1(1:nvar,:) = Hptdi(1:nvar,1:nvar);
   %%Hm = (cyx+H1p_1')*(xxhp\(cxy+H1p_1));
   Hm1 = (cyx+H1p_1')/xxhp;
   Hm2 = cxy+H1p_1;
   GlbAph(i,:) = A0h(stri,i)'*XX'*Hm1;     % alpha_plus_*_transpose
   %%alpMpyh = A0h(stri,i)'*XX'*(yty+Hptdi(1:nvar,1:nvar)-Hm)*XX;
   alpMpyh = A0h(stri,i)'*XX'*(yty+Hptdi(1:nvar,1:nvar))*XX - GlbAph(i,:)*Hm2*XX;

   % ** 3rd bid term in i-th equation,
   % **       with all other 3rd big terms before the i-th equation
   tt3 = tt3 + A0hbH*A0hb + alpMpyh*A0h(stri,i);

   %%%%
   % *** passing the gradient to "pmdg6.m" (analytical gradient)
   %%%%
   %
   % ** daha0 = d((alpha0-alpha0_tilde)'*inv(H_0_tilde)*
   % **                   (alpha0-alpha0_tilde))/d(a0??)
   daha0 = 2.0*A0hbH;
   vaha0(stril+1:stril+strm) = daha0(:);
   %
   % ** daMMa0 = d(alpha0'*P'*(Y'Y+H_1^(-1)+H_m)*P*apha0) / d(a0?)
   daMMa0 = 2.0*alpMpyh;
   vaMMa0(stril+1:stril+strm) = daMMa0(:);

   stril = stril + strm;
end

%disp(sprintf('Loop %d end: %g', i, toc))
disp(sprintf('Loop end (of): %g', toc))

%e_tfat = toc

of = tt1 + 0.5*tt3;
FRESHFUNCTION=1;