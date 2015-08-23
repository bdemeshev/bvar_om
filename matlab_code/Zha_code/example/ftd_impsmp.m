function [xdraws,xdraws1,xdraws2,mphv,sm,sel1,sel2,sel3,timeminutes,xdraw] = ftd_impsmp(xinput)
% [xdraws,xdraws1,xdraws2,mphv,sm,sel1,sel2,sel3,timeminutes,xdraw] = ftd_impsmp(xinput)
%   Importance sampling for selected impulse responses.  Check <<>>1 lines for selected variables,
%      shocks, and steps.  We save all draws in xdraws.  6/2/01.
%
% xinput{1}: nfp -- total number of free parameters
% xinput{2}: nvar -- number of variables
% xinput{3}: xhat -- ML estimate of free parameters in A0
% xinput{4}: hess1 -- Hessian of -logLH for importance sampling
% xinput{5}: Indxv -- index for selected variables of interest; normall first 2 are of our interest
%        to select variables, always check idmat0 to make sure.
%        When Indxv=[], xdraw is empty as well. When IndxGgraph=1, it plots
%         (1) pdf of 1st v for every buffer, (2) scattered plot of 1st and 2nd for every buffer,
%         (3) pdf of 1st v for all sequences; (4) scattered plot of 3rd and 4th for all sequences
%         (5) scattered plot of 1st and 2nd for al sequences.
% xinput{6}: imndraws=nstarts*ndraws2
% xinput{7}: a0indx -- index number for non-zero elements in A0
% xinput{8}: tdf -- degrees of freedom for t-distribution
% xinput{9}: nbuffer -- interval for printing, plotting, and saving
% xinput{10}: Sbd -- nvar-by-nvar S{1}, ..., S{n} -- kind of covariance matrix for each simultaneous equation
%             Already divided by "fss."
% xinput{11}: scf -- reduction scale factor for Metropolis jumping kernel
% xinput{12}: Hsr1 -- square root of covariance matrix for free elements in A0 (nfp-by-nfp)
%             for importance sampling.
% xinput{13}: fss -- effective sample size == nSample-lags+# of dummy observations
%-----------------
% xdraws: length(sel1)-by-length(sel2)-by-length(sel3)-by-imndraws. Draws of selected impusle responses.
% xdraws1:  mean of the selected impulse responses.
% xdraws2:  variance of the selected impulse responses.
% mphv:  imndraws-by-1 vector of unscaled weights
% sm:   sum of the weights
% sel1:  Selected number of steps.  Check a <<>>1 line.
% sel2:  Selected number of variables.  Check a <<>>1 line.
% sel3:  Selected number of shocks.  Check a <<>>1 line.
% timeminutes:  minutes used for this simulation
% xdraw: imndraws-by-length(Indxv) matrix of draws; empty if Indxv is empty.
%
% Written by Tao Zha 1999

nfp=xinput{1}; nvar=xinput{2}; xhat=xinput{3}; hess1=xinput{4}; Indxv=xinput{5};
imndraws=xinput{6}; a0indx=xinput{7}; tdf=xinput{8}; nbuffer=xinput{9};
Sbd=xinput{10}; scf=xinput{11}; Hsr1=xinput{12}; fss=xinput{13};
Bhml=xinput{14}; xxhpc=xinput{15}; lags=xinput{16}; imstp=xinput{17}; ncoef=xinput{18};

sm = 0.0;
mphv=zeros(imndraws,1);
if isempty(Indxv)
   xdraw=[];
else
   xdraw=zeros(imndraws,length(Indxv));   %<<>> draws of selected variables
end

nn=[nvar lags imstp];
sel1 = [1:imstp];  % <<>>1 Selected number of steps.
sel2 = [2 3 4 5];  % <<>>1 Selected number of variables.
sel3 = [2];   % <<>>1 Selected number of shocks.
xdraws = zeros(length(sel1),length(sel2),length(sel3),imndraws);  % For impusle responses.
xdraws1 = zeros(length(sel1),length(sel2),length(sel3));  % 1st moment.
xdraws2 = xdraws1;  % 2nd moment.
tic
for draws=1:imndraws
   %
   if ~mod(draws,nbuffer)
      draws
      %  fwriteid = fopen('outA0.bin','a');
      %  count = fwrite(fwriteid,A0hatw,'double');
      %  status = fclose('all');
   end


   %** draw free elements Avh in A0 and hyperparameters from t-dist
   Avhz1 = Hsr1\randn(nfp,1);     % normal draws
   csq=randn(tdf,1);
   csq=sum(csq .* csq);
   Avhz = xhat+Avhz1/sqrt(csq/tdf);   % Robert, p.382

   % *** compute weights ***
   % * t-dist density, having taken log
   tdhz = -0.5*(length(Avhz)+tdf)*log(1+(Avhz-xhat)'*(hess1*(Avhz-xhat))/tdf);
   % * actual density, having taken log
   hAvhz = a0asfun(Avhz,Sbd,fss,nvar,a0indx);
   hAvhz = -hAvhz;      % converted to logLH

   % * mph: m prim hat in Kloek & van Dijk, p.5
   mph = exp(hAvhz-tdhz-scf);  % scf: scaling factor to prevent overflow
	mphv(draws) = mph;

   sm=sm+mph;
   % * matrix of draws for selected variables
   if ~isempty(Indxv)
      xdraw(draws,:) = Avhz(Indxv)';
   end

   %------- Selected impulse responses -------------
   A0_h = zeros(nvar);
   A0_h(a0indx)=Avhz;
   A0_hin = inv(A0_h);
   %
   Apindm = randn(ncoef,nvar);
   Bh_h = Bhml + (xxhpc\Apindm)*A0_hin;
   swish_h = A0_hin';     % Switching back to the form A0*y(t) = e(t)
   imf_h = zimpulse(Bh_h,swish_h,nn);   % in the form that is congenial to RATS
   imf3_h=reshape(imf_h,size(imf_h,1),nvar,nvar);
       % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
   xdraws(:,:,:,draws) = imf3_h(sel1,sel2,sel3);  % Selected impulse responses.
   xdraws1 = xdraws1 + imf3_h(sel1,sel2,sel3)*mph;  % 1st moment.
   xdraws2 = xdraws2 + (imf3_h(sel1,sel2,sel3).^2)*mph;  % 2st moment.
end
timeminutes = toc/60

xdraws1=xdraws1/sm;
xdraws2=xdraws2/sm;
xdraws2=xdraws2-xdraws1.^2;
