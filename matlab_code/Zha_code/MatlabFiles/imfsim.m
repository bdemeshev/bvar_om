function imfsim(xinput)
% imfsim(xinput)
%        Save the simulated pdfs of impulse responses;
%        Print and save Gelman's measures of B and W for A0's,
%                    only when nstarts (# of starting points) >1.
%        Ref:  Waggoner and Zha "Does Normalization Matter for Inference"
%        See note Forecast (2)
%
% xinput{1}: nfp -- total number of free parameters
% xinput{2}: nvar -- number of variables
% xinput{3}: xhat -- ML estimate of free parameters in A0
% xinput{4}: hess -- Hessian of -logLH
% xinput{5}:Indxv -- index for selected variables of interest; normall first 2 are of our interest
%        to select variables, always check idmat0 to make sure
%        it plots: (1) pdf of 1st v for every buffer, (2) scattered plot of 1st and 2nd for every buffer,
%                  (2) pdf of 1st v for all sequences; (4) scattered plot of 3rd and 4th for all sequences
%                  (5) scattered plot of 1st and 2nd for al sequences.
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
% xinput{22}: idfile1 -- calls "iden6std."   Save stds. of both data and impulse responses in idfile1
% xinput{23}: xxhpc -- chol(X'X+inv(H_p_tilde)): upper triangular but its transpose
%                                      is lower triagular Choleski
% xinput{24}: ImfErr -- if 1, impulse response simulation; if 0, disable this simulation
% xinput{25}: ninv -- number of bins pre-specified to put each draw of impulse response
%                            into a proper bin (or small interval)
% xinput{26}: imstp -- # of steps for impulse responses
% xinput{27}: forep -- forecast periods (# of steps)
% xinput{28}: yact -- actual data (in log except R, U, etc.)
% xinput{29}: yactqg -- quarterly annualized growth in actual data
% xinput{30}: yactCalyg -- calendar annual growth in actual data
% xinput{31}: imfml -- imstp-by-nvar^2 ML impulse responses
% xinput{32}: forepq -- forecast periods for quarterly growth
% xinput{33}: forepy -- forecast periods for annual growth
% xinput{34}: ncoef -- k: # of coeffients per equation
% xinput{35}: Bhml -- ML reduced form parameter B (nvar-by-k)
% xinput{36}: lags -- # of lags
%------------------
% imfcntmulti: (ninv+2,imstp*nvar^2,nstarts) matrix
%           cnt: count; for impulse responses; multi (nstarts) sequences
% All output is saved in outB_W, including Range5, invc, ninv, and imfcntmulti
%
% Written by T. Zha 1999
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
fss=xinput{21}; idfile1=xinput{22}; xxhpc=xinput{23}; ImfErr=xinput{24}; ninv=xinput{25};
imstp=xinput{26}; forep=xinput{27}; yact=xinput{28}; yactqg=xinput{29}; yactCalyg=xinput{30};
imfml=xinput{31}; forepq=xinput{32}; forepy=xinput{33}; ncoef=xinput{34}; Bhml=xinput{35};
lags=xinput{36};

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

imfmean = zeros(imstp,nvar^2);
imfcntmulti = zeros(ninv+2,imstp*nvar^2,nstarts);
          % cnt: count; for impulse responses; multi (nstarts) sequences


%---------------------------------------------------
%  Specify the range for counting the empirical distribution
%
%** load the standard deviations of 6 variables, one for log(y), one for gq, and
%**   the third one for yg
eval(['load ' idfile1 '.prn -ascii']);
eval(['ABstd=' idfile1 ';']);
Range5 = cell(4,1);  % 4: log, qg, yg, and imf

%@@@ Tony's trick to expand the matrix
%
%** In order of log(y), qg, and yg for Range5{i} for i=1:3
Range5{1} =zeros(forep,nvar,2);  % 2: min and max
Range5{1}(:,:,1) = repmat(yact(length(yact(:,1)),:)-10*ABstd(1,:),[forep 1]);  % min, 30 std.
Range5{1}(:,:,2) = repmat(yact(length(yact(:,1)),:)+10*ABstd(1,:),[forep 1]);  % max, 30 std.
%
Range5{2} =zeros(forepq,nvar,2);  % 2: min and max
Range5{2}(:,:,1) = repmat(yactqg(length(yactqg(:,1)),:)-10*ABstd(2,:),[forepq 1]);  % min, 30 std.
Range5{2}(:,:,2) = repmat(yactqg(length(yactqg(:,1)),:)+10*ABstd(2,:),[forepq 1]);  % max, 30 std.
%
Range5{3} =zeros(forepy,nvar,2);  % 2: min and max
Range5{3}(:,:,1) = repmat(yactCalyg(length(yactCalyg(:,1)),:)-10*ABstd(3,:),[forepy 1]);  % min, 30 std.
Range5{3}(:,:,2) = repmat(yactCalyg(length(yactCalyg(:,1)),:)+10*ABstd(3,:),[forepy 1]);  % max, 30 std.
%
Range5{4} =zeros(imstp,nvar^2,2);  % 2: min and max
imfscale = repmat(ABstd(4,:),[1 nvar]);  % because nvar variables to 1, ..., nvar shocks
Range5{4}(:,:,1) = repmat(imfml(1,:)-5*imfscale,[imstp 1]);  % min, 5 std.
Range5{4}(:,:,2) = repmat(imfml(1,:)+5*imfscale,[imstp 1]);  % max, 5 std.
         % Range5(4)(:,:,1): imstp-by-nvar^2.  Column: nvar responses to 1st shock,
         %                        nvar responses to 2nd shock, ...
%**
invc = cell(4,1);     % interval length (used for counting later). 1st 3 cells have each
                      %    forep-by-nvar, and 4th has imstp-by-nvar^2.
for i=1:4
   invc{i} = Range5{i}(:,:,2) - Range5{i}(:,:,1);
end
hbin =  invc{4} ./ ninv;    % bin size for each point of impulse responses kkdf
imfloor = Range5{4}(:,:,1);



%===================================
%  Here begins with the big loop
%===================================
H1 = chol(hess);  % upper triangular so that H1' is a lower triangular decomp
baseW = H_sr;  %inv(H1);  %H_sr;   % covariance matrix without scaling
nswitch=0;  %<<>> total number of sign switches
A0inxhat = inv(A0xhat);   % inverse of ML estimate
a0indx0 = find(idmat0==0);    % index for all zero's in A0;
nn=[nvar lags imstp];

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
   imfcnt = zeros(ninv+2,imstp*nvar^2);   % cnt: count; for impulse responses

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
         A0_hin = inv(A0_h);
      else
         A0ingbs = inv(A0gbs);
         [A0_h,nswitch,A0_hin] = nmlzvar(A0gbs,A0xhat,A0inxhat,IndxNmlr,nswitch,A0ingbs);
      end
      Avhx_norm = A0_h(a0indx);

      %
      % *** normal draws for posterior Aplus conditional on A0h
      %
      %** the mean is Aplushm, and the covariance is inv(xxhp)=Lxxhpc*Lxxhpc'
      Apindm = randn(ncoef,nvar);
      %
      if ~all(all(finite(Bhml)))
         Aplushm=zeros(ncoef,nvar);
         for i=1:nvar
            Aplushm(:,i)=Gb{i}*A0_h(:,i);    % see Zha's forecast (1) p.9
                   % Here, Gb is used to compute A+ where A+(i) = Gb(i)*a0(i)
         end
         Bh_h = (Aplushm + xxhpc\Apindm)*A0_hin;
      else
         Bh_h = Bhml + (xxhpc\Apindm)*A0_hin;
      end

      if ImfErr
         swish_h = A0_hin';     % Switching back to the form A0*y(t) = e(t)
         imf_h = zimpulse(Bh_h,swish_h,nn);   % in the form that is congenial to RATS
         %  imf3_h=reshape(imf_h,size(imf_h,1),nvar,nvar);
         %      % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         imfmean = imfmean + imf_h;  % posterior mean
         imfcnt = empdfsort(imfcnt,imf_h,imfloor,hbin,ninv);
                       % sorted counts (prob.) in bins
      end

      Avhxm = Avhxm + Avhx_norm;      % 1st step to overall mean of parameter
      Avhxs = Avhxs + Avhx_norm.^2;   % 1st step to overall 2nd moment of parameter

      %* compute the mean and 2nd moment
      %** Getting average of variances W and variance of means B/n -- B_n
      %*   see Gelman, p.331, my Shock(0), 12-13, and my Forecast (2), 28-31
      if (nstarts>1)
         Avhxmm = Avhxmm + Avhx_norm;      % n*(phi_.j)
         Avhxss = Avhxss + Avhx_norm.^2;   % 1st step to (phi_ij) for fixed j
      end

      %  A0hatw(:,drawsc-wdraws) = A0_h(:);
      if ~mod(draws,nbuffer)
         starts
         draws
         wdraws=drawsc
         %  fwriteid = fopen('outA0.bin','a');
         %  count = fwrite(fwriteid,A0hatw,'double');
         %  status = fclose('all');
      end
   end
   %
   imfcntmulti(:,:,starts) = imfcnt;

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
timeminutes=timend/60

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

imfmean = imfmean/(imndraws);

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
   save outB_W B_W1 B1 W1 GR1 nstarts ndraws2 imndraws timeminutes Avhxs ...
                  A0xhat A0hm A0hs A0hv IndxGibbs countJump
   if ImfErr
      save outB_W nswitch imfml imfmean imfcntmulti Range5 ninv invc imstp nvar -append
   end

   titstr = ['J ' num2str(nstarts) ' n1 ' num2str(ndraws1) ...
            ' n2 ' num2str(ndraws2) ' timend minutes ' num2str(timend/60)];
   disp(' ')
   disp(titstr)
   disp('B/W sqrt(B) sqrt(W) Std(A0) GR')
   format short g
   [B_W1 sqrt(B1) sqrt(W1) Avhxs GR1]
else
   save outB_W nstarts ndraws1 ndraws2 imndraws timeminutes Avhxs ...
                  A0xhat A0hm A0hs A0hv IndxGibbs countJump
   if ImfErr
      save outB_W nswitch imfml imfmean imfcntmulti Range5 ninv invc imstp nvar -append
   end
end
