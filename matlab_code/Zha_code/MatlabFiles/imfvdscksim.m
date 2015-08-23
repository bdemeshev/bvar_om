function imfvdscksim(xinput,sa0indx,simfindx)
% imfvdscksim(xinput,sa0indx,simfindx)
%        Save the simulated pdfs of impulse responses, vd, shocks, and A0's;
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
% xinput{16}: nSample -- the original sample size including lags
% xinput{17}: IndxNmlr -- index for which normalization rule to choose
% xinput{18}: IndxGibbs -- index for WZ Gibbs; 1; Gibbs; 0: Metropolis
% xinput{19}: scf -- reduction scale factor for Metropolis jumping kernel
% xinput{20}: H_sr -- covariance matrix for free elements in A0 (nfp-by-nfp)
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
% xinput{37}: Psuedo -- 1: for Pseudo out-of-sample; 0: for in-sample (or real-time out-of-sample)
% xinput{38}: q_m = 4 or 12 -- quarterly (4) or monthly (12)
% xinput{39}: imf3ml -- ML impulse responses with row--steps, column--nvar responses,
%                 3rd dimension--nvar shocks
% xinput{40}: vlistlog -- sub index for log in vlist
% xinput{41}: vlistper -- sub index for percent in vlist
% xinput{42}: phi -- X in the form of y = X*B+U. Row: nSmaple-lags+ndobs. Column: ncoef
% xinput{43}: actup -- acutal periods for backing out structural shocks
% xinput{44}: A0ml -- ML A0; column-equation
% xinput{45}: Bhml -- ML Bh (k-by-nvar)
% xinput{46}: yrEnd -- the end year for the estimation period
% xinput{47}: qmEnd -- the end month or quarter for the estimation period
% xinput{48}: yrStart -- the start year for the estimation period
% xinput{49}: qmStart -- the start month or quarter for the estimation period
% sa0indx:  special a0 index.  1: use this; 0: disenable
% simfindx:   special impulse response index.  1: use it; 0: disenable
%------------------
% E.G.:   imfcntmulti: (ninv+2,imstp*nvar^2,nstarts) matrix
%           cnt: count; for impulse responses; multi (nstarts) sequences
% All output is saved in outB_W, including Range5, invc, ninv, imfcntmulti,
%      sckcorcntmulti, Avhxcntmulti, lzvdcntmulti
%
% Written by TAZ 1999
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
lags=xinput{36}; Psuedo=xinput{37}; q_m=xinput{38}; imf3ml=xinput{39}; vlistlog=xinput{40};
vlistper=xinput{41}; phi=xinput{42}; actup=xinput{43}; A0ml=xinput{44}; Bhml=xinput{45};
yrEnd=xinput{46}; qmEnd=xinput{47}; yrStart=xinput{48}; qmStart=xinput{49};


if Psuedo
   disp('Make sure (1) Psuedo=0 in msstart.m and (2) "actuap=nSample-lags" for strucutral shocks')
   disp('Press ctrl-c to abort now')
end

A0_h = zeros(nvar);
A0gbs = A0_h;    % drawn A0 in Gibbs sampling
Avhxml = xhat;   % ML estimate
Avhxmean = zeros(nfp,1);
Avhxs = Avhxmean;

A0xhat = zeros(nvar);
A0xhat(a0indx) = xhat;
%  A0hatw = zeros(nvar^2,nbuffer);

countJump = zeros(nstarts,1);
imstpyer = floor(imstp/q_m);   % yearly

imfmean = zeros(imstp,nvar^2);
imfs = imfmean;
imfyermean = zeros(imstpyer,nvar^2);
imfyers = imfyermean;
imfyer3ml = zeros(imstpyer,nvar,nvar);
imfyer3_h = zeros(imstpyer,nvar,nvar);


imfcntmulti = zeros(ninv+2,imstp*nvar^2,nstarts);
          % cnt: count; for impulse responses; multi (nstarts) sequences
Avhxcntmulti = zeros(ninv+2,nfp,nstarts);
          % cnt: count; for A0; multi (nstarts) sequences
imfyercntmulti = zeros(ninv+2,imstpyer*nvar^2,nstarts);
          % cnt: count; for impulse responses; multi (nstarts) sequences
lzvdcntmulti = zeros(ninv+2,imstp*nvar^2,nstarts);
          % cnt: count; for lz vd; multi (nstarts) sequences
lzvdyercntmulti = zeros(ninv+2,imstpyer*nvar^2,nstarts);
          % cnt: count; for lz vd; multi (nstarts) sequences
sckcorcntmulti = zeros(ninv+2,nvar^2,nstarts);

%---------------------------------------------------
%  Specify the range for counting the empirical distribution
%
%** load the standard deviations of 6 variables, one for log(y), one for gq, and
%**   the third one for yg
eval(['load ' idfile1 '.prn -ascii']);
eval(['ABstd=' idfile1 ';']);
Range5 = cell(9,1);  % 8: log, qg, yg, imf, Avhx, imfyer (yearly), lzvd,
                     % lzvdyer, sckcorv,

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
%
%*** for parameters A0's
Range5{5} =zeros(nfp,2);  % 2: min and max
Avhxscale = abs(Avhxml);
Range5{5}(:,1) = Avhxml-5*Avhxscale;  % min, 5 std.
Range5{5}(:,2) = Avhxml+5*Avhxscale;  % max, 5 std.

%
%*** for yearly impulse responses
yer3 = zeros(1+imstpyer,nvar,nvar);
for k=1:imstpyer
   yer3(k+1,:,:) = mean(imf3ml(1+q_m*(k-1):q_m*k,:,:),1);  % annual average
       % yer3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
end
imfyer3ml(:,vlistlog,:) = ( exp(yer3(2:1+imstpyer,vlistlog,:)-yer3(1:imstpyer,vlistlog,:)) - 1 ).*100;
imfyer3ml(:,vlistper,:) = yer3(2:1+imstpyer,vlistper,:) .* 100;
     % imfyer3ml: row--steps, column--nvar responses, 3rd dimension--nvar shocks
tmp = max(squeeze(max(imfyer3ml,[],1)),[],2);   % nvar-by-1
tmp = tmp';               % 1-by-nvar (variables)
imfyerml = reshape(imfyer3ml,imstpyer,nvar^2);
imfyermlmax = max(abs(imfyerml));     % for error bands (imfyercnt) later
%*** for annual impulse responses
Range5{6} =zeros(imstpyer,nvar^2,2);  % 2: min and max
tmp = repmat(tmp,[1 nvar]);  % because nvar variables to 1, ..., nvar shocks
imfyerscl = repmat(tmp,[imstpyer,1]);  % imstpyer-by-nvar^2
Range5{6}(:,:,1) = imfyerml-5*imfyerscl;  % min, 5 std.
Range5{6}(:,:,2) = imfyerml+5*imfyerscl;  % max, 5 std.
         % Range5(6)(:,:,1): imstpyer-by-nvar^2.  Column: nvar responses to 1st shock,
         %                        nvar responses to 2nd shock, ...

%*** lz variance decomposition (nunlike the traditional, non-cumulative).
tmp0=abs(imf3ml);
tmp1 = repmat(sum(tmp0,3),[1 1 nvar]);   % imstp-by-nvar^2
             % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
lzvd3ml = 100*(tmp0./tmp1);
lzvdml = reshape(lzvd3ml,imstp,nvar^2);
%** for lz vd (non-cumulative)
Range5{7} =zeros(imstp,nvar^2,2);  % 2: min and max
Range5{7}(:,:,1) = zeros(imstp,nvar^2);    % min, 5 std.
Range5{7}(:,:,2) = 100*ones(imstp,nvar^2);  % max, 5 std.
         % Range5(7)(:,:,1): imstp-by-nvar^2.  Column: nvar responses to 1st shock,
         %                        nvar responses to 2nd shock, ...

%*** lz annual variance decomposition (nunlike the traditional, non-cumulative).
tmp0=abs(imfyer3ml);
tmp1 = repmat(sum(tmp0,3),[1 1 nvar]);   % imstp-by-nvar^2
             % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
lzvdyer3ml = 100*(tmp0./tmp1);
lzvdyerml = reshape(lzvdyer3ml,imstpyer,nvar^2);
%** for lz annual vd
Range5{8} =zeros(imstpyer,nvar^2,2);  % 2: min and max
Range5{8}(:,:,1) = zeros(imstpyer,nvar^2);    % min, 5 std.
Range5{8}(:,:,2) = 100*ones(imstpyer,nvar^2);  % max, 5 std.
         % Range5(8)(:,:,1): imstpyer-by-nvar^2.  Column: nvar responses to 1st shock,
         %                                 nvar responses to 2nd shock, ...

%**** Correlations among structural shocks
philr = phi(size(phi,1)-actup+1,:);   % +1 is absolutely needed
yexa = yact(1:actup,:);
Estrexaml = fidcndexa(yexa,philr,A0ml,Bhml,nvar,lags,actup);  % actup-by-nvar
sckvarml = (Estrexaml'*Estrexaml)/actup;
sckcorml = corr(sckvarml);
%** for shock correlation sckcor
Range5{9} =zeros(nvar,nvar,2);  % 2: min and max
Range5{9}(:,:,1) = (-1)*ones(nvar,nvar);    % min, 5 std.
Range5{9}(:,:,2) = ones(nvar,nvar);  % max, 5 std.
         % Range5(9)(:,:,1): nvar^2.  Correlation among structural shocks




%**
invc = cell(9,1);     % interval length (used for counting later).
for i=[1:4 6:9]
   invc{i} = Range5{i}(:,:,2) - Range5{i}(:,:,1);
end
invc{5} = Range5{5}(:,2) - Range5{5}(:,1);
%
imfhbin =  invc{4} ./ ninv;  % bin size for each point of impulse responses
imfloor = Range5{4}(:,:,1);  % lowest point next to -infinity
Avhxhbin = invc{5} ./ ninv;  % bin size for each point of A0
Avhxfloor = Range5{5}(:,1);  % lowest point next to -infinity
imfyerhbin =  invc{6} ./ ninv;  % bin size for each point of annual impulse responses
imfyerfloor = Range5{6}(:,:,1);  % lowest point next to -infinity
lzvdhbin =  invc{7} ./ ninv;  % bin size for each point of variance decompositions
lzvdfloor = Range5{7}(:,:,1);  % lowest point next to -infinity
lzvdyerhbin =  invc{8} ./ ninv;  % bin size for each point of annul variance decompositions
lzvdyerfloor = Range5{8}(:,:,1);  % lowest point next to -infinity
sckcorhbin =  invc{9} ./ ninv;  % bin size for each point of shock correlations
sckcorfloor = Range5{9}(:,:,1);  % lowest point next to -infinity




%*** <<>> Specific requests
%*** compute prob(parameter>0) or joint prob. of sign matches in MP and MD
if sa0indx
   csix1a0 = [-1 -1 1 1 1 -1 1]';  % for 7 parameters in MP and MD
                                  % from 7th to 13th in Avhx_norm
   nspeca0 = 7+3;   % number of specific requests
         % additional 3: 1 for all MP parameters; 2 for all MD parameters;
         %               3 for all paramters in MP and MD
   cspeca0 = zeros(nspeca0,1);
end
%
if simfindx
   csix1imf = [-1 -1 1 -1 -1 1]';   % for 6 variables to MP shock
                                    % maybe at different horizons (esp. for inflation)
   csix1pimf = [1 1 1 24 36 24]';   % the periods for the 6 varialbes to MP shock
                                    % Pcm, M2, FFR, y, P, U.
   nspecimf = nvar+2;      % number of specific requests
               % additional 1: 1 for opposite signs of M2 and R
   cspecimf = zeros(nspecimf,1);
end




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
   cJump = 0;

   imfcnt = zeros(ninv+2,imstp*nvar^2);   % cnt: count; for impulse responses
   Avhxcnt = zeros(ninv+2,nfp);         % cnt: count; for A0's
   imfyercnt = zeros(ninv+2,imstpyer*nvar^2);   % cnt: count; for impulse responses
   lzvdcnt = zeros(ninv+2,imstp*nvar^2);   % cnt: count; for lz vd (non-cumulative)
   lzvdyercnt = zeros(ninv+2,imstpyer*nvar^2);   % cnt: count; for lz annual vd (non-cumulative)
   sckcorcnt = zeros(ninv+2,nvar^2);   % cnt: count; for shock correlations

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

      if sa0indx
         for k=1:7
            cspeca0(k) = cspeca0(k) + ((csix1a0(k)*Avhx_norm(6+k))>0);
         end
         %*** Joint tests
         j1=csix1a0;
         j2=Avhx_norm;
         %** MP equation
         mpall = ((j1(2)*j2(8))/(j1(3)*j2(9)))<0; % & ((j1(1)*j2(7))/(j1(3)*j2(9)))<0;

         cspeca0(8) = cspeca0(8) + mpall;
         %** MD equation
         mdall = ((j1(4)*j2(10))/(j1(5)*j2(11)))<0; %& ((j1(4)*j2(10))/(j1(6)*j2(12)))<0;
         cspeca0(9) = cspeca0(9) + mdall;
         mpdall = mpall & mdall;
         cspeca0(10) = cspeca0(10) + mpdall;
      end

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
         imf3_h=reshape(imf_h,size(imf_h,1),nvar,nvar);
             % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         imfmean = imfmean + imf_h;  % posterior mean
         imfs = imfs + imf_h.^2;   % posterior 2nd moment
         imfcnt = empdfsort(imfcnt,imf_h,imfloor,imfhbin,ninv);
                       % sorted counts (prob.) in bins

         %**** annula impulse responses
         for k=1:imstpyer
            yer3(k+1,:,:) = mean(imf3_h(1+q_m*(k-1):q_m*k,:,:),1);  % annual average
               % yer3:  initialized earlier already
               % yer3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         end
         imfyer3_h(:,vlistlog,:) = ( exp(yer3(2:1+imstpyer,vlistlog,:)-yer3(1:imstpyer,vlistlog,:)) - 1 ).*100;
         imfyer3_h(:,vlistper,:) = yer3(2:1+imstpyer,vlistper,:) .* 100;
               % imfyer3_h: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         imfyer_h = reshape(imfyer3_h,imstpyer,nvar^2);
         imfyermean = imfyermean + imfyer_h;  % posterior mean
         imfyers = imfyers + imfyer_h.^2;   % posterior 2nd moment
         imfyercnt = empdfsort(imfyercnt,imfyer_h,imfyerfloor,imfyerhbin,ninv);
                       % sorted counts (prob.) in bins

         %**** Leeper-Zha variance decomposition (non-cumulative)
         tmp0=abs(imf3_h);
         tmp = repmat(sum(tmp0,3),[1 1 nvar]);   % imstp-by-nvar^2
                      % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         lzvd3_h = 100*(tmp0./tmp);
         lzvd_h = reshape(lzvd3_h,imstp,nvar^2);
         lzvdcnt = empdfsort(lzvdcnt,lzvd_h,lzvdfloor,lzvdhbin,ninv);
                       % sorted counts (prob.) in bins

         %**** Leeper-Zha annual variance decomposition (non-cumulative)
         tmp0=abs(imfyer3_h);
         tmp = repmat(sum(tmp0,3),[1 1 nvar]);   % imstp-by-nvar^2
                      % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
         lzvdyer3_h = 100*(tmp0./tmp);
         lzvdyer_h = reshape(lzvdyer3_h,imstpyer,nvar^2);
         lzvdyercnt = empdfsort(lzvdyercnt,lzvdyer_h,lzvdyerfloor,lzvdyerhbin,ninv);
                       % sorted counts (prob.) in bins

         %**** Correlations among structural shocks
         Estrexa_h = fidcndexa(yexa,philr,A0_h,Bh_h,nvar,lags,actup);  % actup-by-nvar
         sckvar_h = (Estrexa_h'*Estrexa_h)/actup;
         sckcor_h = corr(sckvar_h);
         sckcorcnt = empdfsort(sckcorcnt,sckcor_h,sckcorfloor,sckcorhbin,ninv);

         if simfindx
            for k=1:6
               cspecimf(k) = cspecimf(k) + ((csix1a0(k)*imf_h(1,nvar+k))>0);
                                 % 1st dim in imf_h: periods
            end
            %*** Joint tests
            j1=csix1imf;
            j2=imf_h;
            j3=csix1pimf;
            %** M2 and R in MP equation
            mpall = ((j1(2)*j2(j3(2),nvar+2))/(j1(3)*j2(j3(3),nvar+3)))>0; % & ((j1(1)*j2(7))/(j1(3)*j2(9)))<0;
                          % 1st dim in j2: periods
            cspecimf(nvar+1) = cspecimf(nvar+1) + mpall;
            %** R and P in MP equation
            mpall = ((j1(3)*j2(j3(3),nvar+3))/(j1(5)*j2(j3(5),nvar+5)))>0; %
                          % 1st dim in j2: periods
            cspecimf(nvar+2) = cspecimf(nvar+2) + mpall;
            %
            %  %** MD equation
            %  mdall = ((j1(4)*j2(10))/(j1(5)*j2(11)))<0; %& ((j1(4)*j2(10))/(j1(6)*j2(12)))<0;
            %  cspeca0(9) = cspeca0(9) + mdall;
            %  mpdall = mpall & mdall;
            %  cspeca0(10) = cspeca0(10) + mpdall;
         end
      end

      Avhxmean = Avhxmean + Avhx_norm;      % 1st step to overall mean of parameter
      Avhxs = Avhxs + Avhx_norm.^2;   % 1st step to overall 2nd moment of parameter
      Avhxcnt = empdfsort(Avhxcnt,Avhx_norm,Avhxfloor,Avhxhbin,ninv);


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
   Avhxcntmulti(:,:,starts) = Avhxcnt;
   imfyercntmulti(:,:,starts) = imfyercnt;
   lzvdcntmulti(:,:,starts) = lzvdcnt;
   lzvdyercntmulti(:,:,starts) = lzvdyercnt;
   sckcorcntmulti(:,:,starts) = sckcorcnt;

   if ~IndxGibbs
      countJump(starts,1) = cJump;
   end
end
timend = toc
timeminutes=timend/60

if ~IndxGibbs
   countJump = countJump/ndraws2
end

Avhxmean = Avhxmean/(imndraws);
Avhxs = Avhxs/(imndraws);
Avhxs = sqrt(Avhxs - Avhxmean.^2);    % stardard deviation
A0hm = zeros(nvar);
A0hm(a0indx) = Avhxmean   % mean
A0hs = zeros(nvar);
A0hs(a0indx) = Avhxs;   % standar deviation

imfmean = imfmean/(imndraws);
imfs = imfs/imndraws;
imfs = sqrt(imfs - imfmean.^2);   % standard deviation

imfyermean = imfyermean/(imndraws);
imfyers = imfyers/imndraws;
imfyers = sqrt(imfyers - imfyermean.^2);   % standard deviation


save outB_W nstarts ndraws1 ndraws2 imndraws timeminutes Avhxml Avhxmean Avhxs ...
               Avhxcntmulti A0xhat A0hm A0hs IndxGibbs countJump nvar Range5 ...
               ninv invc nfp a0indx nswitch actup nSample lags yrEnd qmEnd ...
               yrStart qmStart q_m sckcorml sckcorcntmulti

if ImfErr
   if simfindx
      cspecimf = cspecimf/imndraws;
      save outB_W cspecimf -append
   end
   save outB_W imfml imfmean imfs imfcntmulti imstp cspecimf ...
         imfyerml imfyermean imfyers imfyercntmulti imstpyer ...
         lzvdml lzvdcntmulti lzvdyerml lzvdyercntmulti -append
end
%
if sa0indx
   cspeca0 = cspeca0/imndraws;
   save outB_W cspeca0 -append
end


