% 10/24/97
% Distance Method of Waggoner and Zha
% Modified from Sims and Zha's code
% Copyright (c) 1997 Tao Zha


% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')


%global xxhp Hm1t Hm1 Hm SpH FRESHFUNCTION
%
%* =================================================
%* ====== Beginning of the script ==================
%* =================================================
%
%* The available data considered
%
q_m = 12;   % quarters or months
yrBin=1959;   % beginning of the year
qmBin=1;    % begining of the quarter or month
yrFin=2001;   % final year
qmFin=1;    % final quarter
nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
       % total number of the available data -- this is all you have

%*** Load data and series
load xd29      % the default name for the variable is "xdd".
[nt,ndv]=size(xdd);
if nt~=nData
   disp(' ')
   warning(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   %disp(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   disp('Press ctrl-c to abort')
   return
end
%--------
%*** Atlanta model variables
%1  PCU   CPI-U: All Items (SA, 1982-84=100)
%2  PCUSLFE  CPI-U: All Items Less Food and Energy (SA, 1982-84=100)
%3  FFED  Federal Funds [effective] Rate (% p.a.)
%4  FTBS3 3-Month Treasury Bills, Secondary Market (% p.a.)
%5  FCM10 10-Year Treasury Bond Yield at Constant Maturity (% p.a.)
%6  FM2   Money Stock: M2 (SA, Bil.$)
%7 FM1   Money Stock: M1 (SA, Bil.$)
%8 CBM   Personal Consumption Expenditures (SAAR, Bil.$)
%9  CBHM  Personal Consumption Expenditures (SAAR, Bil.Chn.1996$)
%10 LR Civilian Unemployment Rate (SA, %)
%11 LE Civilian Employment: Sixteen Years & Over (SA, Thousands)
%12 LANAGRA  Employees on Nonfarm Payrolls (SA, Thousands)
%13 IP Industrial Production Index (SA, 1992=100)
%14 TRST  Retail Sales (SA, Mil.$)
%15 NAPMC Natl Assn of Purchasing Managers: Mfg: Composite Index, SA by BEA (%)
%16 SP3000   PPI: Finished Goods (SA, 1982=100)
%17 SP1000   PPI: Crude Materials for Further Processing (SA, 1982=100)
%18 SP1600   PPI: Crude Materials less Energy (SA, 1982=100)
%19 PZALL KR-CRB Spot Commodity Price Index: All Commodities (1967=100)
%20 PZRAW KR-CRB Spot Commodity Price Index: Raw Industrials (1967=100)
%21 PZTEXP   Spot Oil Price: West Texas Intermediate [Prior'82=Posted Price] ($/Barrel)
%22 PZRJOC   FIBER Industrial Materials Index: All Items (1990=100)
%23 SP500 Stock Price Index: Standard & Poor's 500 Composite  (1941-43=10)
%24 SP500E   Stock Price Index: Standard & Poor's 500 Composite (EOM, 1941-43=10)
%25 FARAT Adjusted Reserves of Depository Institutions (SA, Mil.$)
%26  rgdpmon  Real GDP (monthly, chain $92)
%27  dgdpmon  Deflator GDP (monthly interpolated, chain $92)
%28  dpcemon  PCE price index (monthly interpolated)
%29  imfcom   IMF commodity price index (all commodities excluding energy prices)

logindx = [1:2 6:9 11:14 16:29];
xdd(:,logindx) = log(xdd(:,logindx));
pctindx = [3:5 10 15];
xdd(:,pctindx)=.01*xdd(:,pctindx);
%
vlist = [20 6 3 26 1 10];    % regarding "xdd", Poil (or imfcom), M2, FFR, GDP, CPI (or PCE), U
vlistlog = [1 2 4 5];       % subset of "vlist"
vlistper = [3 6];           % subset of "vlist"

xlab = {'Inf'
        'MS'
        'MD'
        'y'
        'P'
        'U'};

ylab = {'Pcom'
        'M2'
        'FFR'
        'y'
        'P'
        'U'};

xdd_per = xdd(:,vlist);


%* A specific sample is considered for estimation
%*   Sample period 59:7-82:9, forecast period 82:10-84:9
yrStart=1959;
qmStart=1;
[yrEnd,qmEnd,forep,forepq,forepy,forelabel] = pararc(q_m);
nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmStart+1);
if qmEnd == q_m     % end of the year
   nSampleCal=nSample;            % Cal: calendar year
else
   nSampleCal=(yrEnd-1-yrStart)*q_m + (q_m-qmStart+1);   % Cal: calendar year
end

%* More script variables
%
lags = 13;        % <<>>
idfile='iden6s';
%  automatic decay code (monthly data), only two options: lags = 6 or 13
forepq = forep/3;      % quarterly
actup = 4*q_m;     % <<>> actual periods before forecasting (20 years)
%actup = 12*floor(nSample/12);     % <<>> actual periods before forecasting (8 years)
%actup = 48;     % <<>> actual periods before forecasting (4 years)
actupq = actup;  %actup/3;   % quarterly
actupy = actup/q_m;   % four years
imstp = 4*q_m;      % <<>>  impulse responses (4 years)
ninv = 1000;   % the number of intervals for counting impulse responses
nhp = 6;          % <<>> number of hyperparameters for estimation
%%scf = 2.4/sqrt(nvar);       % scf^2*Sigma (covaraince)
indxGimfml = 1;  % 1: graph ML impulse responses; 0: no graph
scf = 0.25;           % scf^2*Sigma (covaraince)
ndraws1=15;         % 1500, 1st part of Monte Carlo draws
ndraws2=2*ndraws1;         % 2nd part of Monte Carlo draws
ndraws=3*ndraws2         % a total number of Monte Carlo draws
nstarts=3;         % number of starting points
imndraws = nstarts*ndraws2;
tdf = 3;          % degrees of freedom for t-dist
ga = tdf/2;      % asymmetry parameter in Gamma
gb = 2/tdf;      % normalized parameter in Gamma
%
%* =================================================
%* ====== End of the script ========================
%* =================================================


if (q_m==12)
   nStart=(yrStart-yrBin)*12+qmStart-qmBin;  % positive number of months at the start
   nEnd=(yrEnd-yrFin)*12+qmEnd-qmFin;     % negative number of months towards end
elseif (q_m==4)
   nStart=(yrStart-yrBin)*4+qmStart-qmBin;  % positive number of months at the start
   nEnd=(yrEnd-yrFin)*4+qmEnd-qmFin;     % negative number of months towards end
else
   disp('Warning: this code is only good for monthly/quarterly data!!!')
   return
end
%
if nEnd>0 | nStart<0
   disp('Warning: this particular sample consider is out of bounds of the data!!!')
   return
end
%
xdgel=xdd(nStart+1:nData+nEnd,vlist);  % gel: general term for selected xdd
xdata=xdd(nStart+1:nData,vlist);

baddata = find(isnan(xdgel));
if ~isempty(baddata)
   warning('Some data are actually unavailable.')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end
%

[Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
             idmat0,idmatpp] = szasbvar(idfile,q_m,lags,nhp,xdgel);

%[A0hin,Hm1t,fss,ndobs,phi,y,nvar,ncoef,SpH,SpHsc,xxhpc,a0indx,na0p,...
%             idmat0] = szbvar(idfile,q_m,lags,nSample,nhp,xdgel)

%e = y - phi*Hm1t;      % from Y = XB + U, e: (T-lags)*nvar
%Sigu = e'*e/fss;


if all(all(finite(SpH)))
   asInforA0='No asymmetric prior on A0';
   % ***
   % *** Form the inverse of covarinace to draw from t- or Gaussian distribution
   % ***
   %if isempty( find(idmat0-triu(ones(nvar))) )
   %if isempty( find(idmat0-tril(ones(nvar))) )
	if 0
      SpHU=SpH;
      save SpHUout SpHU        % unrestricted case with idfile='iden6a'
   	SpHs = SpH*fss;     % for the Wishart draw
   	SpHsc = chol(SpHs);     % upper triangular where Sphsc' is
                             %  lower triangular Choleski, for Wishart draw
		A0in = chol(SpH);
		A0 = inv(A0in);
   	%SpHsic = chol(inv(SpHs))';     % inverse Choleski decomposition -- lower triangular
   	%A0hin = chol(SpH);    % upper triangular, inverse of A0h, each column
                         %   corresponds to an equation.
	else
   	%
   	% ** Get the ML estimate of A0
   	%
   	rand('state',fix(100*sum(clock)));
   	%rand('state',238465);
      x = 50*rand(na0p,1);
   	H0 = eye(na0p,na0p);
      crit = 1.0e-10;
   	nit = 10000;
   	%
   	tic
   	[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = ...
           csminwel('a0lhfun',x,H0,'a0lhgrad',crit,nit,SpH,fss,nvar,a0indx);
			                        % <<>> pmdf?
   	e_t = toc
		A0 = zeros(nvar);
		A0(a0indx)=xhat;
		A0in = inv(A0);
   	fhat
   	xhat
   	ghat
   	itct
   	fcount
   	retcodehat
		save outm e_t xhat ghat fhat idmat0 idmatpp itct itct fcount retcodehat
	end
else
   asInforA0='Yes, asymmetric prior on A0';
   %
   % ** Get the ML estimate of A0
   %
   rand('state',fix(100*sum(clock)));
   %rand('state',238465);
   x = 50*rand(na0p,1);
   H0 = eye(na0p,na0p);
   crit = 1.0e-10;
   nit = 10000;
   %
   tic
   [fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = ...
           csminwel('a0asfun',x,H0,'a0asgrad',crit,nit,Sbd,fss,nvar,a0indx);
			                        % <<>> pmdf?
   e_t = toc
	A0 = zeros(nvar);
	A0(a0indx)=xhat;
	A0in = inv(A0);
   fhat
   xhat
   ghat
   itct
   fcount
   retcodehat
	save outm e_t xhat ghat fhat idmat0 idmatpp itct itct fcount retcodehat
end

%==================
%  Impulse responses first
%==================
%A0(4,2) = -xhat(7);   % output in MD
%A0(5,2)=-xhat(7);   % price in MD

%<<<<<<<<<<< Beginning: positive on diag(A0) <<<<<<<<<<<<
%a0dpindx = find(diag(A0)<0);
%A0(:,a0dpindx) = -A0(:,a0dpindx);
%xhat = A0(a0indx);   % peak of posterior squeezing out all 0's or linear
										  % restrictions.
%A0in = inv(A0);
%<<<<<<<<<<< End: positive on diag(A0) <<<<<<<<<<<<


%>>>>>>>>>>>> Beginning: normalize on diagonal of inv(A0) >>>>>>>>
%
%a0dpindx = find(diag(A0in)<0);
%A0(:,a0dpindx) = -A0(:,a0dpindx);
%A0in(a0dpindx,:) = -A0in(a0dpindx,:);
% $$$ on the responses to MS shock, because we want R (2,3) to be positive
%
%  if (A0in(2,3)<0)
%     % making R positive
%     A0(:,2) = -A0(:,2);
%     A0in(2,:) = -A0in(2,:);
%  end
%
xhat = A0(a0indx);
%>>>>>>>>>>>> End: normalize on diagonal of inv(A0) >>>>>>>>

swish = A0in';       % each row corresponds to an equation

%--------------- Get Hessian at the peak ------------------
% * stepsize
%grdh = eps^(1/3) * (max([abs(xhat) ones(size(xhat,1),1)]'))' .* (abs(xhat) ./ xhat);
grdh = 1e-04;
%hess = hesscd('a0asfun',xhat,grdh,Sbd,fss,nvar,a0indx);
%check = chol(hess);

hess1 = gradcd('a0asgrad',xhat,grdh,Sbd,fss,nvar,a0indx);
Hsr1 = chol(hess1);     % note, upper triangular, for MC draws later


%=====================
% Now to get Bh (k-by-m) for impulse responses, where k=nocef, m=novar
%=====================
%
%Bh = Hm1t;


if ~all(all(finite(Bh)))
   asInforBh='Yes, asymmetric prior on lagged A+ and thus Bh';
   Aplus=zeros(ncoef,nvar);
   for i=1:nvar
      Aplus(:,i)=Gb{i}*A0(:,i)
   end
   Bh = Aplus/A0;
else
   asInforBh='No asymmetric prior on lagged A+ or Bh';
   Aplus = Bh*A0;
end

size(Bh)
size(idmatpp)
save idenml xhat idmat0 idmatpp A0 A0in Gb Sbd Aplus Bh asInforA0 asInforBh hess1 ...
            Hsr1
                      % <<>>  Most important save for later use!!!

% ** impulse responses
nn = [nvar lags imstp];
%imf = zimpulse(Bh,swish,nn);    % in the form that is congenial to RATS
[vd,str,imf] = errors(Bh,swish,nn);
vd3=reshape(vd,size(vd,1),nvar,nvar);
         % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
vd3s=permute(vd3,[1 3 2]);
         % imf3s: permuted so that row--steps, column--nvar shocks,
			%                                3rd dimension--nvar responses
			% Note: reshape(imf3s(1,:,:),nvar,nvar) = A0in  (columns -- equations)


scaleout = imcgraph(imf,nvar,imstp,xlab,ylab,indxGimfml)   % nvar-by-2 matrix


temp = max(abs(scaleout)')';
imfscale = repmat(temp,[1 2]);
imfscale(:,1) = -imfscale(:,1);
save iden6im.prn imfscale -ascii


return %<<>>

%----------------------------------------
% Tests for LR, HQ, Akaike, Schwarz
%----------------------------------------
%
load SpHUout

logLHU=-fss*sum(log(diag(chol(SpHU)))) -0.5*fss*nvar       % unrestricted logLH

logLHR=-fhat                                % restricted logLH
tra = reshape(SpHU,nvar*nvar,1)'*reshape(A0*A0',nvar*nvar,1)

S=2*(logLHU-logLHR)
SC = (nvar*(nvar+1)/2 - length(a0indx)) * log(fss)
df=(nvar*(nvar+1)/2 - length(a0indx))

norm(A0'*SpHU*A0-diag(diag(ones(6))))/6
norm(SpHU-A0in'*A0in)/6
SpHR=A0in'*A0in;


logLHU=-fss*sum(log(diag(chol(A0in'*A0in)))) -0.5*fss*nvar       % unrestricted logLH

corU=zeros(nvar);
for i=1:nvar
   for j=1:nvar
     corU(i,j)=SpHU(i,j) / sqrt( SpHU(i,i) * SpHU(j,j));
   end
end
corU

corR=zeros(nvar);
for i=1:nvar
   for j=1:nvar
     corR(i,j)=SpHR(i,j) / sqrt( SpHR(i,i) * SpHR(j,j));
   end
end
corR


%[vU,dU]=eig(SpHU)
%[vR,dR]=eig(A0in'*A0in)
