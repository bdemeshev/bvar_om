function [xdraw,timeminutes,nswitch] = a0onlysim(xinput)
% [xdraw,timeminutes,nswitch] = a0onlysim(xinput)
%        Export and plot the simulated pdfs of draws of only selected parameters;
%        Print and save Gelman's measures of B and W for all parameters;
%        Ref:  Waggoner and Zha "Does Normalization Matter for Inference"
%        See note Forecast (2)
%
% xinput{1}: nfp -- total number of free parameters
% xinput{2}: nvar -- number of variables
% xinput{3}: xhat -- ML estimate of free parameters in A0
% xinput{4}: hess1 -- Hessian of -logLH
% xinput{5}:Indxv -- index for selected variables of interest; normall first 2 are of our interest
%        to select variables, always check idmat0 to make sure.
%        When Indxv=[], xdraw is empty as well. When IndxGgraph=1, it plots
%         (1) pdf of 1st v for every buffer, (2) scattered plot of 1st and 2nd for every buffer,
%         (3) pdf of 1st v for all sequences; (4) scattered plot of 3rd and 4th for all sequences
%         (5) scattered plot of 1st and 2nd for al sequences.
% xinput{6}: IndxGraph - 1: plot graphs; 0: no graphs
% xinput{7}: idmat0 -- Index for non-zero elements in A0 with column to equation
% xinput{8}: nstarts -- # of starting points in Gibbs sampling
% xinput{9}: ndraws1 -- # of 1st loop to be discarded
% xinput{10}: ndraws2 -- # of 2nd loop for saving A0 draws
% xinput{11}: imndraws=nstarts*ndraws2
% xinput{12}: a0indx -- index number for non-zero elements in A0
% xinput{13}: tdf -- degrees of freedom for t-distribution
% xinput{14}: nbuffer -- interval for printing, plotting, and saving
% xinput{15}: Sbd -- nvar-by-nvar S{1}, ..., S{n} -- kind of covariance matrix for each simultaneous equation
%             Already divided by "fss."
% xinput{16}: nSample -- the original sample size including lags
% xinput{17}: IndxNmlr -- index for which normalization rule to choose
% xinput{18}: IndxGibbs -- index for WZ Gibbs; 1; Gibbs; 0: Metropolis
% xinput{19}: scf -- reduction scale factor for Metropolis jumping kernel
% xinput{20}: H_sr -- square root of the inverse of the covariance matrix
%             for free elements in A0 (nfp-by-nfp)
% xinput{21}: fss -- effective sample size == nSample-lags+# of dummy observations
%------------------
% xdraw: nstarts*ndraws2-by-length(Indxv) matrix of draws;
%        empty if Indxv=[]
% timeminutes:  minutes used for this simulation
%
% Written by Tao Zha 1999
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

nfp=xinput{1}; nvar=xinput{2}; xhat=xinput{3}; hess=xinput{4}; Indxv=xinput{5};
IndxGraph=xinput{6}; idmat0=xinput{7}; nstarts=xinput{8}; ndraws1=xinput{9}; ndraws2=xinput{10};
imndraws=xinput{11}; a0indx=xinput{12}; tdf=xinput{13}; nbuffer=xinput{14}; Sbd=xinput{15};
nSample = xinput{16}; IndxNmlr=xinput{17}; IndxGibbs=xinput{18}; scf=xinput{19}; H_sr=xinput{20};
fss=xinput{21};

if IndxGraph>length(Indxv)
   disp(' ')
   warning('when IndxGraph=1, Indxv must have at least 1 elements')
   disp('Press ctrl-c to abort')
   pause
end

Avhx_bs = zeros(nfp,1);
Avhx_bm = zeros(nfp,1);
Avhx_bj = zeros(nfp,1);
Avhx_cj = zeros(nfp,1);
Avhx_cm = zeros(nfp,1);
A0_h = zeros(nvar);
A0gbs = A0_h;    % drawn A0 in Gibbs sampling
Avhxm = zeros(nfp,1);
Avhxs = Avhxm;
A0xhat = zeros(nvar);
A0xhat(a0indx) = xhat;
%  A0hatw = zeros(nvar^2,nbuffer);

countJump = zeros(nstarts,1);


%===================================
%  Here begins with the big loop
%===================================
H1 = chol(hess);  % upper triangular so that H1' is a lower triangular decomp
baseW = H_sr;  %inv(H1);  %H_sr;   % covariance matrix without scaling
nswitch=0;  %<<>> total number of sign switches
A0inxhat = inv(A0xhat);   % inverse of ML estimate
a0indx0 = find(idmat0==0);    % index for all zero's in A0;
if isempty(Indxv)
   xdraw=[];
else
   xdraw=zeros(imndraws,length(Indxv));   %<<>> draws of selected variables
end


[cT,vR,kdf] = gibbsglb(Sbd,idmat0,nvar,fss);

tic
for starts = 1:nstarts
   starts
   if starts == 1
      A0gbs(a0indx) = xhat;   % from "load ..."
      if ~IndxGibbs   % Metropolist
         Avhx = xhat;
         hAvhx = a0asfun(Avhx,Sbd,fss,nvar,a0indx);
         hAvhx = -hAvhx;      % converted to logLH
      end
   else
      Avhx = baseW*randn(nfp,1);  %H_sr*randn(nfp,1);   % D: discarded sequence
      csq=randn(tdf,1);
      csq=sum(csq .* csq);
      Avhx = xhat+Avhx/sqrt(csq/tdf);
      %** Normalization by the choice of IndxNmlr
      A0gbs(a0indx) = Avhx;
      if ~IndxNmlr(5)
         [A0gbs,jnk] = nmlzvar(A0gbs,A0xhat,A0inxhat,IndxNmlr,nswitch,[]);
      else
         A0ingbs = inv(A0gbs);
         [A0gbs,jnk,jnk1] = nmlzvar(A0gbs,A0xhat,A0inxhat,IndxNmlr,nswitch,A0ingbs);
      end
      %
      if ~IndxGibbs   % Metropolist
         Avhx = A0gbs(a0indx);
         hAvhx = a0asfun(Avhx,Sbd,fss,nvar,a0indx);
         hAvhx = -hAvhx;      % converted to logLH
      end
   end
   %
   Avhxmm = zeros(nfp,1);
   Avhxss = zeros(nfp,1);
   cJump = 0;

   for draws = 1:ndraws1
      if IndxGibbs
         A0gbs = gibbsvar(A0gbs,cT,vR,nvar,fss,kdf);
      else     % Metropolis
         [Avhx,hAvhx,cJump] = smtplis(Avhx,hAvhx,tdf,cJump,scf,...
                       baseW,nfp,Sbd,fss,nvar,a0indx);
      end
   end

   wdraws=(starts-1)*ndraws2+0;
   for draws = 1:ndraws2
      drawsc = (starts-1)*ndraws2+draws;
      if IndxGibbs
         A0gbs = gibbsvar(A0gbs,cT,vR,nvar,fss,kdf);
         A0gbs(a0indx0) = 0; % set all zeros in A0gbs clean to avoid possible cumulative round-off errors
      else        % Metropolis
         [Avhx,hAvhx,cJump] = smtplis(Avhx,hAvhx,tdf,cJump,scf,...
                       baseW,nfp,Sbd,fss,nvar,a0indx);
         A0gbs(a0indx) = Avhx;
      end

      %*** call normalization so that A0_h is normalized
      if ~IndxNmlr(5)
         [A0_h,nswitch] = nmlzvar(A0gbs,A0xhat,A0inxhat,IndxNmlr,nswitch,[]);
      else
         A0ingbs = inv(A0gbs);
         [A0_h,nswitch,A0_hin] = nmlzvar(A0gbs,A0xhat,A0inxhat,IndxNmlr,nswitch,A0ingbs);
      end
      Avhx_norm = A0_h(a0indx);

      Avhxm = Avhxm + Avhx_norm;      % 1st step to overall mean of parameter
      Avhxs = Avhxs + Avhx_norm.^2;   % 1st step to overall 2nd moment of parameter

      %* compute the mean and 2nd moment
      %** Getting average of variances W and variance of means B/n -- B_n
      %*   see Gelman, p.331, my Shock(0), 12-13, and my Forecast (2), 28-31
      if (nstarts>1)
         Avhxmm = Avhxmm + Avhx_norm;      % n*(phi_.j)
         Avhxss = Avhxss + Avhx_norm.^2;   % 1st step to (phi_ij) for fixed j
      end

      if ~isempty(Indxv)
         xdraw((starts-1)*ndraws2+draws,:) = Avhx_norm(Indxv)';  %<<>>
      end
      %  A0hatw(:,drawsc-wdraws) = A0_h(:);
      if ~mod(draws,nbuffer)
         starts
         draws
         wdraws=drawsc
         %  fwriteid = fopen('outA0.bin','a');
         %  count = fwrite(fwriteid,A0hatw,'double');
         %  status = fclose('all');
         if IndxGraph
            %** 1-d pdf plot
            %  figure(2)
            %  histpdfg(xdraw((starts-1)*ndraws2+(1:draws),1),50,[],[],[]);
            %  pause(1)

            %** 2-d scatterplot
            figure(3)
            plot(xdraw((starts-1)*ndraws2+(1:draws),1),...
               xdraw((starts-1)*ndraws2+(1:draws),2),'.');
            drawnow
         end
      end
   end
   %
   if ~IndxGibbs
      countJump(starts,1) = cJump;
   end
   %
   %** Getting average of variances W and variance of means B/n -- B_n
   %**   see Gelman, p.331, my Shock(0), pp.12-13, and my Forecast (2), pp.28-31
   if (nstarts>1)
      Avhx_aj = Avhxmm/ndraws2;         %  (phi_.j)
      Avhx_bs = Avhx_bs + Avhx_aj.^2;  %  1st step to 2nd moment of (phi_.j)
      Avhx_bm = Avhx_bm + Avhx_aj;     % 1st step to (phi_..)
      %
      Avhx_bj = Avhxss/ndraws2;          % 2nd moment of (phi_ij) for fixed j
      Avhx_cj = Avhx_bj - Avhx_aj.^2;   %  ((n-1)/n)*S^2_j
      Avhx_cm = Avhx_cm + Avhx_cj;     %  sum of ((n-1)/n)*S^2_j across j
   end
end
timend = toc

if ~IndxGibbs
   countJump = countJump/ndraws2
end

Avhxm = Avhxm/(imndraws);
Avhxs = Avhxs/(imndraws);
Avhxv = Avhxs - Avhxm.^2;
Avhxs = sqrt(Avhxv);   % stardard deviation
A0hm = zeros(nvar);
A0hm(a0indx) = Avhxm   % mean
A0hv = zeros(nvar);
A0hv(a0indx) = Avhxv;  % varaince matrix
A0hs = zeros(nvar);
A0hs(a0indx) = Avhxs;   % standar deviation

%**** Getting Within-Sequence W and Between-Sequence B_n
%        see Gelman, p.331, my Shock(0), pp.12-13, and my Forecast (2), pp.28-31
if (nstarts>1)
   AvhxW = (ndraws2/(ndraws2-1))*Avhx_cm/nstarts;
                                 % W: average of j within-sequence variances
   AvhxB_n = (nstarts/(nstarts-1)) * ( Avhx_bs/nstarts - (Avhx_bm/nstarts).^2 );
                                 % B/n:  variance of J within-sequence means
   AvhxB = ndraws2*AvhxB_n;    % B
   %
   B_W1 = AvhxB ./ AvhxW;
   B1 = AvhxB;
   W1 = AvhxW;
   GR1 = sqrt((ndraws2-1)/ndraws2 + B_W1/ndraws2);
      % measure of Gelman reduction; need not be 1 to be accurate,
      %               contrary to what Gelman claims
   save outB_W B_W1 B1 W1 GR1 nstarts ndraws2 imndraws timend Avhxs ...
                  A0xhat A0hm A0hs A0hv IndxGibbs countJump nswitch

   titstr = ['J ' num2str(nstarts) ' n1 ' num2str(ndraws1) ...
            ' n2 ' num2str(ndraws2) ' timend minutes ' num2str(timend/60)];
   disp(' ')
   disp(titstr)
   disp('B/W sqrt(B) sqrt(W) Std(A0) GR')
   format short g
   [B_W1 sqrt(B1) sqrt(W1) Avhxs GR1]
else
   save outB_W nstarts ndraws2 imndraws timend Avhxs ...
                  A0xhat A0hm A0hs A0hv IndxGibbs countJump nswitch
end

disp(' ')
disp('nswitch nswitch/imndraws -- # of sign switches')
[nswitch nswitch/imndraws]
disp(' ')
disp('timend/60 minutes')
timeminutes = timend/60
nswitch=nswitch/imndraws;

if IndxGraph
   %
   %  figure(4)
   %  histpdfg(xdraw(:,1),50,[],[],[]);
   %  figure(5)
   %  plot(xdraw(:,3),xdraw(:,4),'.');
   %  figure(6)
   %  plot(xdraw(:,1),xdraw(:,2),'.');
end