% Applies to both linear and exclusion restrictions.
% (1) Marginal likelihood function p(Y) for constant structural VAR models, using Chib (1995)'s ``Marginal Likelihood from the Gibbs Output'' in JASA.
% (2) Conditional likelihood function f(Y|A0, A+) on the ML estimate for constant exclusion-identified models.
% See Forecast (II) pp.67-80.
%
% Tao Zha, September 1999.  Quick revisions, May 2003.  Final revision, September 2004.

msstart2    % start the program in which everyhting is initialized through msstart2.m
if ~indxEstima
   warning(' ')
   disp('You must set IxEstima=1 in msstart to run this program')
	disp('Press ctrl-c to abort now')
	pause
end

A0xhat = zeros(size(A0hat));
Apxhat = zeros(size(Aphat));
if (0)
   %Robustness check to see if the same result is obtained with the purterbation of the parameters.
   for k=1:nvar
      bk = Ui{k}'*A0hat(:,k);
      gk = Vi{k}'*Aphat(:,k);
      A0xhat(:,k) =  Ui{k}*(bk + 5.2*randn(size(bk)));   % Perturbing the posterior estimate.
      Apxhat(:,k) = Vi{k}*(gk + 5.2*randn(size(gk)));      % Perturbing the posterior estimate.
   end
else
   %At the posterior estimate.
   A0xhat = A0hat;   % ML estimate of A0
   Apxhat = Aphat;      % ML estimate of A+
end
%--- Rename variables.
YatYa = yty;
XatYa = xty;
ytx = xty';
YatXa = ytx;
XatXa = xtx;



%--------- The log value of p(A0,A+) at some point such as the peak ----------
vlog_a0p = 0;
Yexpt=0;  % exponential term for Y in p(Y|A0,A+) at some point such as the peak
Apexpt=0.0;  % 0.0 because we have chosen posterior estimate of A+ as A+*.  Exponential term for A+ conditional on A0 and Y
%======= Computing the log prior pdf of a0a+ and the exponential term for Y in p(Y|A0,A+).
for k=1:nvar
   a0k = A0xhat(:,k);  % meaningful parameters in the kth equation.
   apk = Apxhat(:,k);  % meaningful parameters in the kth equation.

   %--- Prior settings.
   S0bar = H0invtld{k};  %See Claim 2 on p.69b.
   Spbar = Hpinvtld{k};
   bk = Ui{k}'*a0k;  % free parameters in the kth equation.
   gk = Vi{k}'*apk;  % free parameters in the kth equation.
   gbark = Ptld{k}*bk;   % bar: prior

   %--- The exponential term for Y in p(Y|A0,A+)
   Yexpt = Yexpt - 0.5*(a0k'*YatYa*a0k - 2*apk'*XatYa*a0k + apk'*XatXa*apk);
   %--- The log prior pdf.
   vlog_a0p = vlog_a0p - 0.5*(size(Ui{k},2)+size(Vi{k},2))*log(2*pi) + 0.5*log(abs(det(S0bar))) + ...
         0.5*log(abs(det(Spbar))) - 0.5*(bk'*S0bar*bk+(gk-gbark)'*Spbar*(gk-gbark));
   %--- For p(A+|Y,a0) only.
   tmpd = gk - Pmat{k}*bk;
   Apexpt = Apexpt - 0.5*tmpd'*(Hpinv{k}*tmpd);
end
vlog_a0p

%--------- The log value of p(Y|A0,A+) at some point such as the peak. ----------
%--------- Note that logMarLHres is the same as vlog_Y_a, just to double check. ----------
vlog_Y_a = -0.5*nvar*fss*log(2*pi) + fss*log(abs(det(A0xhat))) + Yexpt
                 % a: given a0 and a+
logMarLHres = 0;   % Initialize log of the marginal likelihood (restricted or constant parameters).
for ki=1:fss   %ndobs+1:fss     % Forward recursion to get the marginal likelihood.  See F on p.19 and pp. 48-49.
   %----  Restricted log marginal likelihood function (constant parameters).
   [A0l,A0u] = lu(A0xhat);
   ada = sum(log(abs(diag(A0u))));   % log|A0|
   termexp = y(ki,:)*A0xhat - phi(ki,:)*Apxhat;   % 1-by-nvar
   logMarLHres = logMarLHres - (0.5*nvar)*log(2*pi) + ada - 0.5*termexp*termexp';  % log value
end
logMarLHres


%--------- The log value of p(A+|Y,A0) at some point such as the peak ----------
totparsp = 0.0;
tmpd = 0.0;
for k=1:nvar
   totparsp = totparsp + size(Vi{k},2);
   tmpd = tmpd + 0.5*log(abs(det(Hpinv{k})));
end
vlog_ap_Ya0 = -0.5*totparsp*log(2*pi) + tmpd + Apexpt;




%===================================
%  Compute p(a0,k|Y,ao) at some point such as the peak (in this situation, we simply
%   generate results from the original Gibbs sampler).  See FORECAST (2) pp.70-71
%===================================
%--- Global set up for Gibbs.
[Tinv,UT] = fn_gibbsrvar_setup(H0inv, Ui, Hpinv, Pmat, Vi, nvar, fss);
%
vlog_a0_Yao = zeros(nvar,1);
  % the log value of p(a0k|Y,ao) where ao: other a's at some point such as the peak of ONLY some a0's
vlog=zeros(ndraws2,1);
tic
for k=1:nvar
   bk = Ui{k}'*A0xhat(:,k);
   indx_ks=[k:nvar];  % the columns that exclude 1-(k-1)th columns
   A0gbs0 = A0hat;   % starting at some point such as the peak
   nk = n0(k);

   if k<nvar
      %--------- The 1st set of draws to be tossed away. ------------------
      for draws = 1:ndraws1
         if ~mod(draws,nbuffer)
            disp(' ')
            disp(sprintf('The %dth column or equation in A0 with %d 1st tossed-away draws in Gibbs',k,draws))
         end
         A0gbs1 = fn_gibbsrvar(A0gbs0,UT,nvar,fss,n0,indx_ks);
         A0gbs0=A0gbs1;    % repeat the Gibbs sampling
      end


      %--------- The 2nd set of draws to be used. ------------------
      for draws = 1:ndraws2
         if ~mod(draws,nbuffer)
            disp(' ')
            disp(sprintf('The %dth column or equation in A0 with %d usable draws in Gibbs',k,draws))
         end
         [A0gbs1, Wcell] = fn_gibbsrvar(A0gbs0,UT,nvar,fss,n0,indx_ks);
         %------ See p.71, Forecast (II).
         %------ Computing p(a0_k|Y,a_others) at some point such as the peak along the dimensions of indx_ks.
         Vk = Tinv{k}\Wcell{k};  %V_k on p.71 of Forecast (II).
         gbeta = Vk\bk;  % inv(V_k)*b_k on p.71 of Forecast (II) where alpha_k = b_k in our notation.
         [Vtq,Vtr]=qr(Vk',0);  %To get inv(V_k)'*inv(V_k) in (*) on p.71 of Forecast (II).
         %
         vlog(draws) = 0.5*(fss+nk)*log(fss)-log(abs(det(Vk)))-0.5*(nk-1)*log(2*pi)-...
                  0.5*(fss+1)*log(2)-gammaln(0.5*(fss+1))+fss*log(abs(gbeta(1)))-...
                  0.5*fss*bk'*(Vtr\(Vtr'\bk));

         A0gbs0=A0gbs1;    % repeat the Gibbs sampling
      end
      vlogm=max(vlog);
      qlog=vlog-vlogm;
      vlogxhat=vlogm-log(ndraws2)+log(sum(exp(qlog)));
      vlog_a0_Yao(k) = vlogxhat;
         % The log value of p(a0_k|Y,a_others) where a_others: other a's at some point such as the peak of ONLY some a0's
   else
      disp(' ')
      disp(sprintf('The last(6th) column or equation in A0 with no Gibbs draws'))
      [A0gbs1, Wcell] = fn_gibbsrvar(A0gbs0,UT,nvar,fss,n0,indx_ks)
      %------ See p.71, Forecast (II).
      %------ Computing p(a0_k|Y,a_others) at some point such as the peak along the dimensions of indx_ks.
      Vk = Tinv{k}\Wcell{k};  %V_k on p.71 of Forecast (II).
      gbeta = Vk\bk;  % inv(V_k)*b_k on p.71 of Forecast (II) where alpha_k = b_k in our notation.
      [Vtq,Vtr]=qr(Vk',0);  %To get inv(V_k)'*inv(V_k) in (*) on p.71 of Forecast (II).
      %
      vloglast = 0.5*(fss+nk)*log(fss)-log(abs(det(Vk)))-0.5*(nk-1)*log(2*pi)-...
               0.5*(fss+1)*log(2)-gammaln(0.5*(fss+1))+fss*log(abs(gbeta(1)))-...
               0.5*fss*bk'*(Vtr\(Vtr'\bk));
      vlog_a0_Yao(k) = vloglast;
   end
end
timimutes=toc/60
ndraws2

disp('Prior pdf -- log(p(a0hat, a+hat)):');
vlog_a0p
disp('LH pdf -- log(p(Y|a0hat, a+hat)):');
vlog_Y_a
disp('Posterior Kernal -- logp(ahat) + logp(Y|ahat):');
vlog_Y_a + vlog_a0p
disp('Posterior pdf -- log(p(a0_i_hat|a0_other_hat, Y)):');
vlog_a0_Yao
disp('Posterior pdf -- log(p(aphat|a0hat, Y)):');
vlog_ap_Ya0

%--------- The value of marginal density p(Y) ----------
disp(' ');
disp(' ');
disp('************ Marginal Likelihood of Y or Marginal Data Density: ************');
vlogY = vlog_a0p+vlog_Y_a-sum(vlog_a0_Yao)-vlog_ap_Ya0
