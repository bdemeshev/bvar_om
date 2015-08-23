%
% Get accurate B/W:  stage 1 in Metropolis algorithm
%
% November 1998

msstart    % start the program in which everyhting is initialized through msstart.m
if ~IxEstima
	warning('You must set IxEstima=1 in msstart to run this program')
	disp('Press ctrl-c to abort now')
	pause
end

%=================================================
% Generate error bands by first drawing from the
%  t-distribution for each column equation
%=================================================
%
seednumber = 472534;
randn('state',seednumber);
rand('state',seednumber);

nfp = length(a0indx);
Indxv = []; %[1:7];  %[3:7 15:19];
       % index for selected variables of interest; normall first 2 are of our interest
       % to select variables, always check idmat0 to make sure
       % it plots: (1) pdf of 1st v for every buffer, (2) scattered plot of 1st and 2nd for every buffer,
       %           (2) pdf of 1st v for all sequences; (4) scattered plot of 3rd and 4th for all sequences
       %           (5) scattered plot of 1st and 2nd for al sequences.
IndxGraph = 0;
IndxGibbs = 1;  % 1:  WZ Gibbs Sampler; 0: Straight Metropolis


%???????
%  H_sr = inv(Hsr1);   %eye(nfp); %  % upper triangular decomp

%*** generate draws of A0 from Gibbs or Metropolis
%xinput = cell(21,1);
%xinput{1}=nfp; xinput{2}=nvar; xinput{3}=xhat; xinput{4}=hess1; xinput{5}=Indxv;
%xinput{6}=IndxGraph; xinput{7}=idmat0; xinput{8}=nstarts; xinput{9}=ndraws1; xinput{10}=ndraws2;
%xinput{11}=imndraws; xinput{12}=a0indx; xinput{13}=tdf; xinput{14}=nbuffer; xinput{15}=Sbd;
%xinput{16}=nSample; xinput{17}=IndxNmlr; xinput{18}=IndxGibbs; xinput{19}=scf; xinput{20}=H_sr;
%xinput{21}=fss;
%%
%[xdraw,timeminutes,nswitch] = a0onlysim(xinput);
%save outxdraw xdraw nstarts ndraws2 ndraws1 imndraws IndxGibbs fss timeminutes Indxv ...
%                     IndxGibbs nswitch A0


%*** generate draws of A0 from importance sampling
%  xinput = cell(21,1);
%  xinput{1}=nfp; xinput{2}=nvar; xinput{3}=xhat; xinput{4}=hess1; xinput{5}=Indxv;
%  xinput{6}=imndraws; xinput{7}=a0indx; xinput{8}=tdf; xinput{9}=nbuffer;
%  xinput{10}=Sbd; xinput{11}=scf; xinput{12}=Hsr1; xinput{13}=fss;
%  %
%  [xdraw,mphv,sm,timeminutes] = a0impsmp(xinput);
%  mphvs=mphv/sm;
%  figure(2)
%  plot(mphvs)
%  max(mphv)
%  timeminutes
%  save outxdraw xdraw nstarts ndraws2 ndraws1 imndraws fss timeminutes Indxv ...
%           mphv sm mphvs

if 1    % Generate selected impulse responses from Gibbs sampling
   xinput = cell(35,1);
   xinput{1}=nfp; xinput{2}=nvar; xinput{3}=xhat; xinput{4}=hess1; xinput{5}=Indxv;
   xinput{6}=IndxGraph; xinput{7}=idmat0; xinput{8}=nstarts; xinput{9}=ndraws1; xinput{10}=ndraws2;
   xinput{11}=imndraws; xinput{12}=a0indx; xinput{13}=tdf; xinput{14}=nbuffer; xinput{15}=Sbd;
   xinput{16}=nSample; xinput{17}=IndxNmlr; xinput{18}=IndxGibbs; xinput{19}=scf; xinput{20}=H_sr;
   xinput{21}=fss;
   %------ Added June 2001 ------
   xinput{22}=Bhml; xinput{23}=xxhpc; xinput{24}=lags; xinput{25}=imstp; xinput{26}=ncoef;
   [xdraws,xdraws1,xdraws2,sel1,sel2,sel3,timeminutes,xdraw] = ftd_gibbsmp(xinput);
   xdrawsml = imf3ml(sel1,sel2,sel3);
   save outxdraws xdraws xdrawsml xdraws1 xdraws2 nstarts ndraws2 ndraws1 imndraws fss timeminutes ...
            sel1 sel2 sel3
else    % Generate selected impulse responses from importance sampling
   xinput = cell(36,1);
   xinput{1}=nfp; xinput{2}=nvar; xinput{3}=xhat; xinput{4}=hess1; xinput{5}=Indxv;
   xinput{6}=imndraws; xinput{7}=a0indx; xinput{8}=tdf; xinput{9}=nbuffer;
   xinput{10}=Sbd; xinput{11}=scf; xinput{12}=Hsr1; xinput{13}=fss;
   %------ Added June 2001 ------
   xinput{14}=Bhml; xinput{15}=xxhpc; xinput{16}=lags; xinput{17}=imstp; xinput{18}=ncoef;
   [xdraws,xdraws1,xdraws2,mphv,sm,sel1,sel2,sel3,timeminutes] = ftd_impsmp(xinput);
   mphvs=mphv/sm;
   plot(mphvs)
   save outxdraws xdraws xdrawsml xdraws1 xdraws2 nstarts ndraws2 ndraws1 imndraws fss timeminutes ...
            mphv sm mphvs sel1 sel2 sel3
end

%---------- Plot the selected impulse responses from Gibbs or importance sampling ---------
x1e = NaN*ones(size(xdraws1,1),2+size(xdraws1,2),1+2);
x1e(:,1,1)=[1:size(xdraws1,1)]';
x1e(:,2,1)=zeros(size(xdraws1,1),1);
x1e(:,3:end,1)=100*xdraws1;
x1e(:,3:end,2)=100*(xdraws1-sqrt(xdraws2));
x1e(:,3:end,3)=100*(xdraws1+sqrt(xdraws2));
set(gca,'ColorOrder',[0 0 0]); % turn the color off and set it to black
set(gca,'LineStyleOrder', '-|-|-|-.|--|:')
ftd_mseriesgraph(x1e,[1:size(x1e,2)-2],size(x1e,2)-2,1,12,[],[],[])



%  xinput = cell(36,1);
%  xinput{1}=nfp; xinput{2}=nvar; xinput{3}=xhat; xinput{4}=hess; xinput{5}=Indxv;
%  xinput{6}=IndxGraph; xinput{7}=idmat0; xinput{8}=nstarts; xinput{9}=ndraws1; xinput{10}=ndraws2;
%  xinput{11}=imndraws; xinput{12}=a0indx; xinput{13}=tdf; xinput{14}=nbuffer; xinput{15}=Sbd;
%  xinput{16}=nSample; xinput{17}=IndxNmlr; xinput{18}=IndxGibbs; xinput{19}=scf; xinput{20}=H_sr;
%  xinput{21}=fss; xinput{22}=idfile1; xinput{23}=xxhpc; xinput{24}=ImfErr; xinput{25}=ninv;
%  xinput{26}=imstp; xinput{27}=forep; xinput{28}=yact; xinput{29}=yactqg; xinput{30}=yactCalyg;
%  xinput{31}=imfml; xinput{32}=forepq; xinput{33}=forepy; xinput{34}=ncoef; xinput{35}=Bhml;
%  xinput{36}=lags;
%  %
%  imfsim(xinput);     % simulate draws of A0s and impusle responses.  All output is saved in outB_W
