%
%   This .m file is called for point graphics or error bands and
%   It starts for both data and Bayesian estimation (when IxDataonly==1,
%          no Bayesian estimation but only data), which allows you to set
% (a) available data range,
% (b) sample range,
% (c) rearrangement of actual data such as mlog, qg, yg
% (d) conditions of shocks 'Cms',
% (c) conditions of variables 'nconstr',
% (e) soft conditions 'nbancon,'
% (f) produce point conditional forecast (at least conditions on variables).
%
% 10/24/97  -- begining of this enterprise
% ML Distance Method of Waggoner and Zha's "Normalization, Probability Distribution, ..."
% See Zha's Forecast (2) pp.28-58
%
% Copyright (c) October 1998 by Eric Leeper and Tao Zha
% Revision 10/19/98


% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')


format short g     % format

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
if (qmEnd==q_m)     % end of the year
   nSampleCal=nSample;            % Cal: calendar year
	qmSub=0;         % Sub: substraction or substitute
else
   nSampleCal=(yrEnd-1-yrStart)*q_m + (q_m-qmStart+1);   % Cal: calendar year
                        % the sample ends in December of the previous calendar year
	qmSub=qmEnd;     % Sub: substraction or substitute
end

%%%%----------------------------------------
%  Some settings
%
lags = 13;   % automatic decay code (monthly data), only two options: lags = 6 or 13
nindv = 2;    % number of individual variables
vindv = [1 2];  %[2 5];      % 2nd (M2) and 5th (inflation) variable
Aband = 1;    % 1: error bands with only A0 and A+ random.
Sband = 1;    % 1: error bands with only shocks e's random.
Gbook = 0;    % 1: for the Goldbook; 0: for academic analysis
indxGimfml = 1;  % 1: plot ML impulse responses; 0: no plot
IxEstima = 1;  %  applied below in this .m file.
       %1: Bayesian estimation; 0: no estimation and data only
lmqyIndx=[0 0 0 1];  % Index for which output: monthly level, mg, qg, and calendar yg
	                  % This index is designed solely for draws to save time.
		               % For ML, simply use [1 1 1 1]
IndxNmlr = [1 0 0 0 0];  % imported by nmlzvar.m
    % Index for which normalization rule to choose
    % Only one of the elments in IndxNmlr can be non-zero
    % IndxNmlr(1): ML A distance rule (supposed to be the best)
    % IndxNmlr(2): ML Ahat distance rule (to approximate IndxNmlr(1))
    % IndxNmlr(3): ML Euclidean distance rule (not invariant to scale)
    % IndxNmlr(4): Positive diagonal rule
    % IndxNmlr(5): Positive inv(A) diagonal rule (if ~IndxNmlr(5), no need for A0inu,
    %                                      so let A0inu=[])

%%%%----------------------------------------
% Hard conditions directly on variables
%
gDLSIx = 0;  % 1: graph point forecast on variables; 0: disable
nconstr1=forep;      % number of the 1st set of constraints
nconstr2=forepy;     % number of the 2nd set of constraints
nconstr=nconstr1+nconstr2;   % q: 4 years -- 4*12 months.
                         % When 0, no conditions directly on variables <<>>
%nconstr=forep
nconstr=0
eq_ms = 2; %[];      % location of MS equation; if [], all shocks
PorR = 3*ones(nconstr,1);   % the variable conditioned.  3: FFR; 5: CPI
%PorR = [3 5];   % the variable conditioned.  3: FFR; 5: CPI
%PorR = 3;

%%%%----------------------------------------
% Conditions directly on future shocks
%
Cms = 0     % 1: condition on ms shocks; 0: disable this and "fidcnderr.m" gives
             %   unconditional forecasts if nconstr = 0 as well;  <<>>
nCms = 4;   % number of the stance of policy; 0 if no tightening or loosening
eq_Cms = 2;  % location of MS shocks
TLindx = 1*ones(1,nCms);  % 1-by-nCms vector; 1: tightening; 0: loosen
TLnumber = [0.5 0.5 0 0];  %94:4 % [2 2 1.5 1.5]; %79:9  %[1.5 1.5 1 1]; 90:9
                          % 1-by-nCms vector; cut-off point for MS shocks
TLmean = zeros(1,nCms);
              % unconditional, i.e., 0 mean, for the final report in the paper
if Cms
	eq_ms = [];
	   % At least at this point, it makes no sense to have DLS type of eq_ms; 10/12/98
	if all(isfinite(TLnumber))
		for k=1:nCms
	 		TLmean(k) = lcnmean(TLnumber(k),TLindx(k));
			             % shock mean magnitude. 1: tight; 0: loose
							 % Never used for any subsequent computation but
							 %   simply used for the final report in the paper.
			%TLnumber(k) = fzero('lcutoff',0,[],[],TLmean(k))
                % get an idea about the cutoff point given TLmean instead

		end
	end
else
	nCms = 0;   % only for the use of the graph by msprobg.m
	TLnumber = NaN*ones(1,nCms);
	            % -infinity, only for the use of the graph by msprobg.m
end

%%%%----------------------------------------
% Soft conditions on variables
%
nbancon = 0  % # of band condtions; when 0, disable this option
  % Note (different from "fidencon") that each condition corres. to variable
banact = 1;    % 1: use infor on actual; 0:  preset without infor on actual
if nbancon
	banindx = cell(nbancon,1);  % index for each variable or conditon
	banstp = cell(nbancon,1);    % steps:  annual in general
	banvar = zeros(nbancon,1);    % varables:  annual in general
   banval = cell(nbancon,1);    % band value (each variable occupy a cell)
	badval{1} = zeros(length(banstp{1}),2);   % 2: lower or higher bound

	banstp{1} = 1:4;      % 3 or 4 years
	banvar(1) = 3;      % 3: FFR; 5: CPI
	if ~banact
   	for i=1:length(banstp{1})
      	banval{1}(i,:) = [5.0 10.0];
   	end
	end
end
%

pause(0.1)
disp(' ')
disp('For uncondtional forecasts, set nconstr = Cms = nbancon = 0')
pause(1)

%%%%----------------------------------------
% The rest of settings
%
ImfErr = 1;   % 1: generate error bands for impulse responses;
              % 0: no bands (i.e., only ML)
              % Note, error bands can be activated only when Aband = 1
				  % Set Aband to be zero when out-of-sample forecasts are generated.
				  % Set Cms=1 and nconstr=0 to improve efficiency (not yet efficient).
Psuedo = 0;     % 1: psuedo out-of-sample; 0: real time out-of-sample
Rform = 0;   % <<>> 1: reduced form: triangular, import "iden6a"
if Rform
   idfile='iden6r';
else
   idfile='iden6s';
end
Binwrite=0;   % 1 activate binary write; 0 disable it
%
if (Aband | Sband)
	Graphfore = 0;   % 1: graphics in "fore_gh.m"; 0: no graphics in the execution.
else
	Graphfore = 1;   % 1: graphics in "fore_gh.m"; 0: no graphics in the execution.
end
%actup = 12*floor(nSample/12);     % <<>> actual periods before forecasting (8 years)
actup = 7*q_m;     % In general, need 3 years to allow at least 1 calendar year growth
%actup = nSample-lags;     % for the whole sample, set Pseudo=0 and parac.m set to
								  % the end of the data for xinsample.mat and insampleg.m
%actup = 0;
if Psuedo
	Sexact={'Forecast' actup+forep};   % 'Back out exact shocks; # -- forecast steps
                          % If Sexact{2}=0, no shocks will be backed out
else
	Sexact={'History ' actup};   % Back out exact shocks; # -- actual steps
                          % If Sexact{2}=0, no shocks will be backed out
end
imstp = forep;      % <<>>  impulse responses (4 years)
ninv = 1000;   % the number of intervals for counting impulse responses
nhp = 6;          % <<>> number of hyperparameters for estimation
%scf = 2.4/sqrt(nvar);       % for Metropolis, scf^2*Sigma (covaraince)
%scf = 0.10; %0.40; 1.8;  %0.25;    % for Metropolis, scf^2*Sigma (inv(covaraince))
scf = 15340;%structural  4640;%choleski    % for importance sampling
uscf = 30;      % scale factor for starting points from uniform distributions

nbuffer = 1000;        % a block or buffer of draws (buffer) that
                      %       is saved to the disk (not memory)
ndraws1=0*3*nbuffer         % 1st part of Monte Carlo draws
ndraws2=50*nbuffer         % 2nd part of Monte Carlo draws
%ndraws=ndraws1+ndraws2;         % a total number of Monte Carlo draws
nstarts=1;         % number of starting points
imndraws = nstarts*ndraws2;   % total draws for impulse responses or forecasts
%n2draws=30;
%if nconstr
%	imndraws=1;
%	n2draws=4000;
%end
tdf = 9;          % degrees of freedom for t-dist
%
%* =================================================
%* ====== End of the script ========================
%* =================================================



%(1)--------------------------------------
% Further data analysis
%(1)--------------------------------------
%
%*** Month, quarter, and year of the actual series used for output
%      All these may include forecast horizon, depending on Psuedo
tem = q_m*yrEnd + qmEnd - actup;
yrIni = floor(tem/q_m);
qmIni = mod(tem,q_m) + 1;
              % monthly (for quarterly, need to modifiy quarterly stuff below)
actupq = floor((actup - mod(qmEnd,3)) / 3);  % quarters between qmIni and qmEnd
actupy = floor((actup - qmEnd*(qmEnd~=12)) / q_m);
                                         % calendar years between qmIni and qmEnd
nyears = floor(actup/q_m) + forepy*Psuedo+1;
     % +1 and others are more than the number of years we need, in order to be enough
nmonths = actup + forep*Psuedo;
years = yrIni + [0:nyears];
tem = ( kron(years,[1 zeros(1,11)])+kron(ones(1,nyears+1),[0 2:12]) )';
myears = tem(qmIni:qmIni+nmonths-1);   % to match output "yact"
%myears = myears(2:length(myears));    % to match output "yactmg"  XXXX Important
%
if (q_m == 12)     % monthly so as to convert to quarterly
	for i=2:4
		if ~mod(qmIni+i,3)
	   	qIni = (qmIni+i)/3;
		   nquarters = actupq + forepq*Psuedo;
		   tem = ( kron(years,[1 0 0 0])+kron(ones(1,nyears+1),[0 2 3 4]) )';
		   qyears = tem(qIni:qIni+nquarters-1);
		   qyears = qyears(5:length(qyears));    % to match output "yactqg"
		end
	end
end
%
yearsCal = yrEnd-actupy+(qmEnd==12):yrEnd+forepy*Psuedo-(qmEnd~=12);
yearsCal = yearsCal(2:length(yearsCal))';     % to match output "yactCalyg"
%
if Psuedo
	%----------------
	%  Monthly, quarterly, and calendar years for the forecast horizon
	%   This is only valid if Psuedo
	%----------------
	fbegm = qmEnd+1;   % +1 because the first month out of sample
	if fbegm > q_m    % the first forecast month into the next year
		fbegm = 1;     % reset
		fbegy = yrEnd+1;   % the first forecast year for the graph
	else
		fbegy = yrEnd;    % the first forecast year for the graph
	end
	fbegq = ceil(fbegm/3);   % turn to the first quarter out of sample
	ffiny = fbegy+forepy-1;   % the forecast final year

	%**** Monthly (I guess the same for log, level, or growth because of forecast)
	%       For acutal data, myears for log and level is +1 over mg
	ibegy = find(myears==fbegy);
	ifiny = find(myears==ffiny);
	myearsFore = myears(ibegy+fbegm-1:ibegy+fbegm+forep-2);

	%**** quarterly growth
	ibegy = find(qyears==fbegy);
	qyearsFore = qyears(ibegy+fbegq-1:ibegy+fbegq+forepq-2);

	%**** calendar annual growth
	yearsCalFore = [fbegy:fbegy+forepy-1]';
end

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
%***  Note, both xdgel and xdata have the same start with the specific sample
xdgel=xdd(nStart+1:nData+nEnd,vlist);
      % gel: general xdd within sample (nSample)
xdata=xdd(nStart+1:nData,vlist);
      % beyond sample into forecast horizon until the end of the data yrFin:qmFin
if ~(nSample==size(xdgel,1))
	warning('The sample size (including lags) and data are incompatible')
	disp('Check to make sure nSample and size(xdgel,1) are the same')
	return
end

baddata = find(isnan(xdgel));
if ~isempty(baddata)
   warning('Some data for this selected sample are actually unavailable.')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end

%-------------
% Rearranged actual data as mlog, mg, qg, and calendar yg.
%-------------
xinput = cell(15,1);
xinput{1}=xdata; xinput{2}=nSample; xinput{3}=nSampleCal;
xinput{4}=actup; xinput{5}=actupq; xinput{6}=actupy; xinput{7}=vlist;
xinput{8}=vlistlog; xinput{9}=vlistper; xinput{10}=q_m; xinput{11}=forep;
xinput{12}=forepq; xinput{13}=forepy; xinput{14}=Psuedo; xinput{15}=qmEnd;
[yact,yactmg,yactqg,yactCalyg,yactCal] = datactcon(xinput);
%============= Actual Data ============
actualmyear=[myears yact];
%actualqyear=[qyears yactqg];
actualCal=[yearsCal yactCalyg];


if IxEstima
	%(2)----------------------------------------------------------------------------
	% Estimation
	% ML forecast and impulse responses
	% Hard or soft conditions for conditional forecasts
	%(2)----------------------------------------------------------------------------
	%
	%*** Bayesian coefficient estimation with asymmetric prior
	[Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
	             idmat0,idmatpp] = szasbvar(idfile,q_m,lags,nhp,xdgel);
	neqn=nvar;

	if Rform
	   asInforA0='No asymmetric prior on A0';
	   % ***
	   % *** Form the inverse of covarinace to draw from t- or Gaussian distribution

	   SpHs = SpH*fss;     % for the Wishart draw
	   SpHsc = chol(SpHs);     % upper triangular where Sphsc' is
	                       %  lower triangular Choleski, for Wishart draw
	   A0in = chol(SpH);
	   A0 = inv(A0in);
	   %SpHsic = chol(inv(SpHs))';     % inverse Choleski decomposition -- lower triangular
	   %A0hin = chol(SpH);    % upper triangular, inverse of A0h, each column
	                      %   corresponds to an equation.

	   xhat = A0(a0indx);   % peak of posterior without 0's

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
	   %
	else
	   load idenml    % SpH, Bh, Gb, Sbd, etc.
	end

	%* prepare the inverse of the covariance for all the free elements in A0
	Sbdinc = cell(nvar,1);    % the inverse of the covariance
	Sp_h = sparse(na0p,na0p);    % sparse matrix of zeros

	Sp_cov = sparse(na0p,na0p);   % covariance sparse covariance
	H_sr = sparse(na0p,na0p);    % sparse matrix of zeros; sr -- square root
			               % square root of the inverse of covariance matrix
	stril = 0;
	for i = 1:nvar
		Sbdinc{i} = Sbd{i}*fss;         % the inverse of the covariance
	                         %    in the exponential part of p(A0|y)
		%
	   stri = find(idmat0(:,i));
	   strm = length(stri);
	   XX = zeros(nvar,strm);
	   XX(stri,:) = eye(strm);        % XX in Gelman, etel, on p479
	   Sp_h(stril+1:stril+strm,stril+1:stril+strm) = XX'*Sbdinc{i}*XX;
	   Sp_cov(stril+1:stril+strm,stril+1:stril+strm) = ...
	                   inv(Sp_h(stril+1:stril+strm,stril+1:stril+strm));
	   H_sr(stril+1:stril+strm,stril+1:stril+strm) = ...
	            chol(Sp_cov(stril+1:stril+strm,stril+1:stril+strm))';
	                    % lower triagular decomposition of covariance

	   stril = stril + strm;
	end


	%==================
	%  Impulse responses before forecast
	%==================
	%A0 = zeros(nvar);
	%A0(a0indx)=xhat;
	%A0(4,2) = -xhat(7);   % output in MD
	%A0(5,2)=-xhat(7);   % price in MD
	%A0in = inv(A0);

	swish = A0in';       % each row corresponds to an equation
	Bhml = Bh;      % from idenml.mat. Bh has been taken care of
	                % for asymmetric priors on lagged relationships.
	%Bh = Hm1t;    % no longer have Hm1t in this new szasbvar

	% ** impulse responses
	nn = [nvar lags imstp];
	imfml = zimpulse(Bhml,swish,nn);    % in the form that is congenial to RATS
	%[vd,str,imf] = errors(Bh,swish,nn);
	imf3ml=reshape(imfml,size(imfml,1),nvar,nvar);
	         % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
	imf3sml=permute(imf3ml,[1 3 2]);
	         % imf3s: permuted so that row--steps, column--nvar shocks,
				%                                3rd dimension--nvar responses
				% Note: reshape(imf3s(1,:,:),nvar,nvar) = A0in  (columns -- equations)
	%** M2 growth and Inflation calculation, inflation (5th variable) response to 2nd shock
	tmp1 = imf3ml(1:imstp-1,vindv,2);    % 2nd shock
	tmp1 = [zeros(1,nindv);tmp1];
	tmp2 = imf3ml(:,vindv,2);   % 2nd shock
	infrm2ml = ( (exp(tmp2-tmp1)).^12 -1)*100;    % monthly change (annualized)
	%figure
	%plot(1:imstp,infrm2ml)
	%figure
	%plot(1:imstp,tmp2)
	%** Annual inflation
	ntmp = floor(imstp/q_m);
	tmp = zeros(1+ntmp,nindv);
	for i=1:ntmp
		tmp(i+1,:) = mean(imf3ml(1+q_m*(i-1):q_m*i,vindv,2),1);  % annual average, 2nd shock
	end
	infra2ml = ( exp(tmp(2:1+ntmp,:)-tmp(1:ntmp,:)) - 1 ).*100;
	infra2max = max(abs(infra2ml));     % for error bands (infracnt) later
	%figure
	%plot(1:ntmp,infra2ml)

   scaleout = imcgraph(imfml,nvar,imstp,xlab,ylab,indxGimfml);
	disp('Click main Matlab screen to put the figure behind')
	pause(1)
	imfstd = max(abs(scaleout)');   % row: nvar (largest number); used for standard deviations

	%**** save stds. of both data and impulse responses in idfile1
	temp = [std(yact); std(yactqg); std(yactCalyg); imfstd];  %<<>>
	save iden6std.prn temp -ascii   % export to the file "iden6std.prn", 4-by-nvar
	%
	if (Aband | Sband)
		idfile1='iden6std';
	end



	%=====================================
	% Now, out-of-sample forecasts. Note: Hm1t does not change with A0.
	%=====================================
	%
	% * updating the last row of X (phi) with the current (last row of) y.
	phil = phi(size(phi,1),:);
	phil(nvar+1:ncoef-1) = phil(1:ncoef-1-nvar);
	enSample = size(y,1);    % effective number of samples
	phil(1:nvar) = y(enSample,:);
	%
	%*** ML unconditional point forecast
	nn = [nvar lags forep];
	yforeml = forecast(Bhml,phil,nn);    % forep-by-nvar, in log
	%*** Converted to mlevle, mg, qg, and calendar yg
   %
   %  [yforelml,yforemgml,yforeqgml,yforeCalygml] = fore_cal(yforeml,xdata,nvar,nSample,...
   %                    nSampleCal,forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,...
   %                    vlistper,[1 1 1 1]);
   %


	%-------------------------------------------------
	% Setup for point conditional forecast
	% ML Conditional Forecast
	%-------------------------------------------------
	%
	%% See Zha's note "Forecast (1)" p. 5, RATS manual (some errors in RATS), etc.
	%
	%% Some notations:  y(t+1) = y(t)B1 + e(t+1)inv(A0). e(t+1) is 1-by-n.
	%%    Let r(t+1)=e(t+1)inv(A0) + e(t+2)C + .... where inv(A0) is impulse
	%%          response at t=1, C at t=2, etc. The row of inv(A0) or C is
	%%          all responses to one shock.
	%%    Let r be q-by-1 (such as r(1) = r(t+1)
	%%                 = y(t+1) (constrained) - y(t+1) (forecast)).
	%%    Use impulse responses to find out R (k-by-q) where k=nvar*nsteps
	%%        where nsteps the largest constrained step.  The key of the program
	%%        is to creat R using impulse responses
	%%    Optimal solution for shock e where R'*e=r and e is k-by-1 is
	%%                 e = R*inv(R'*R)*r.
	%
	ylast = y(enSample,:);
	%
	%---------- B: last 12 months data (not calendar) ----------
	%
	%indx12 = size(y,1)-q_m+1:size(y,1);
	%ylast12 = y(indx12,:);        % last 12 months data
	%
	%---------- B: last 12 months in calendar year ----------
	%
	if (qmEnd == q_m)
		indx12 = enSample-q_m+1:enSample;
		ylast12Cal = y(indx12,:);
		yCal_1 = zeros(1,length(vlist));
	else
		indx12 = enSample-qmEnd-q_m+1:enSample-qmEnd;
		ylast12Cal = y(indx12,:);
		%
		tem = y(enSample-qmEnd+1:enSample,:);
	   yCal_1 = sum(tem,1);
	         % note, dimension "1" in sum(ytem,1) is necessary because when
	         %   qmEnd=1, sum(ytem,1) will give us a nomber, not a vector.  But this
				%   is not what we want.
	end
	%

	if Gbook                  % for Atlanta Fed's Goldbook
		if (nconstr > 0) & (length(PorR)>1)
			%*** initializing
			stepcon=cell(nconstr,1);  % initializing, value y conditioned
			valuecon=zeros(nconstr,1);  % initializing, value y conditioned
			varcon=zeros(nconstr,1);  % initializing, endogous variables conditioned
			varcon(:)=PorR;     % 3: FFR; 5: CPI
			%varcon(1:nconstr/2)=PorR(1);     % 2: M2
			%varcon(nconstr/2+1:nconstr)=PorR(2);     % 3: FFR

		   for i=1:5
		      stepcon{i}=1;
		   end
			%
			for i=6:9
				stepcon{i}=2;
			end
		   %
		   %k=0;
		   %for i=10:9+12
		   %   k=k+1;
		   %   stepcon{i}=2+k;    % first 2 months have been conditioned
		   %end
		   %stepcon{10}=[3:12]';   % average inflation or FFR
		   %stepcon{11}=[13:24]';   % average inflation or FFR
		   %stepcon{12}=[25:36]';   % average inflation or FFR
		   %stepcon{13}=[37:48]';   % average inflation or FFR
		   steps_ms = stepcon{1};

		   valuecon(1:9) = [
		     4.6062
		     8.3101
		     0.0556
		     5.0870
		     0.0470
		     8.3178
		     0.0551
		     5.0876
		     0.0460
		   		];
		   %valuecon(10)=0.055;    % average FFR
		   %valuecon(10)=(12*(mean(ylast12Cal(:,5))+log(1+1.3968/100))-5.087-5.0876)/10;   % average inflation
		   %valuecon(11)=mean(ylast12Cal(:,5))+log(1+1.3968/100)+log(1+1.7758/100);
		   %valuecon(12)=valuecon(11)+log(1+2.33/100);
		   %valuecon(13)=valuecon(12)+log(1+2.33/100);

		   %for i=1:nconstr/2
		   %   stepcon{i}=i;    % for M2
		   %end
		   %
		   if Psuedo
		      for i=1:nconstr/2
		         valuecon(i) = yact(actup+i,varcon(i));
		         %valuecon(i) = 0.060;      % 95:01
		         %valuecon(i) = (0.0475+0.055)/2;   % 94:10
		      end
		   end

			%----------------------------
			% Second constraint
			%----------------------------
		   %for i=(nconstr/2+1):nconstr
		   %   stepcon{i}=i-nconstr/2;    % for FFR
		   %end
		   %
		   if Psuedo
		      for i=(nconstr/2+1):nconstr
					k = i-nconstr/2;
		         valuecon(i) = yact(actup+k,varcon(i));
		         %valuecon(i) = 0.060;      % 95:01
		         %valuecon(i) = (0.0475+0.055)/2;   % 94:10
		      end
		   end
		end
	elseif (nconstr > 0)
		%*** initializing
		stepcon=cell(nconstr,1);  % initializing, value y conditioned
		valuecon=zeros(nconstr,1);  % initializing, value y conditioned
		varcon=zeros(nconstr,1);  % initializing, endogous variables conditioned
	   varcon(:)=PorR;     % 3: FFR; 5: CPI
	   %varcon(1:nconstr1)=PorR(1);     % 3: FFR
	   %varcon(nconstr1+1:nconstr1+nconstr2)=PorR(2);     % 5: CPI

	   for i=1:nconstr
	      stepcon{i}=i;      % FFR
	   end
		%
	   %bend=12;
	   %stepcon{1}=[1:bend]'; % average over
	   %stepcon{nconstr1+1}=[1:q_m-qmSub]';  % average over the remaing months in 1st forecast year
	   %stepcon{nconstr1+2}=[q_m-qmSub+1:q_m-qmSub+12]';    % average over 12 months next year
	   %stepcon{nconstr1+3}=[q_m-qmSub+13:q_m-qmSub+24]';    % average over 12 months. 3rd year
		%stepcon{nconstr1+4}=[q_m-qmSub+25:q_m-qmSub+36]';    % average over 12 months. 4th year

	   steps_ms=[1:nconstr];
	           % Only used for fidencond.m when Aband=Sband=0.  It is valid
	           %   only if length(eq_ms)=1.  But as Zha (1998) aruges,
	           %   letting length(eq_ms)=1 is seldom a good way to analyze policy
	           %   effects if steps_ms is too long.

	   if Psuedo
	      for i=1:nconstr
	         valuecon(i) = yact(actup+i,varcon(i));
	         %valuecon(i) = mean( yact(actup+1:actup+bend,varcon(i)) );
	         %valuecon(i) = 0.060;      % 95:01
	         %valuecon(i) = (0.0475+0.055)/2;   % 94:10
	      end
			%
	      %for i=nconstr1+1:nconstr1+nconstr2
	      %   i=1;
	      %   valuecon(nconstr1+i) = ( ( mean(ylast12Cal(:,varcon(nconstr1+i)),1) + ...
	      %        log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100) )*q_m - ...
	      %        yCal_1(:,varcon(nconstr1+i)) ) ./ length(stepcon{nconstr1+i});
	      %                    % the same as unconditional "yactCalyg" 1st calendar year
	      %   i=2;
	      %   valuecon(nconstr1+i) = mean(ylast12Cal(:,varcon(nconstr1+i))) +  ...
	      %        log(1+yactCalyg(yAg-yFg+1,varcon(nconstr1+i))/100) ...
	      %                       + log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
	      %                           % the same as actual "yactCalgy" 2nd calendar year
	      %   i=3;
	      %   valuecon(nconstr1+i) = valuecon(nconstr1+i-1) + ...
	      %                               log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
	      %                           % the same as actual "yactCalgy" 3rd calendar year
	      %   %i=4;
	      %   %valuecon(nconstr1+i) = valuecon(nconstr1+i-1) + ...
	      %   %                            log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
	      %                           % the same as actual "yactCalgy" 4th calendar year
	      %end
		else
	      %valuecon(i) = 0.060;      % 95:01
	   end
	else
		valuecon = [];
		stepcon = [];
		varcon = [];
	end

	nstepsm = 0;   % initializing, the maximum step in all constraints
	for i=1:nconstr
	   nstepsm = max([nstepsm max(stepcon{i})]);
	end

	if nbancon
		if Psuedo & banact
			for i=1:length(banstp{1})
				banval{1}(1:length(banstp{1}),1) = ...
			   	 yactCalyg(yAg-yFg+1:yAg-yFg+length(banstp{1}),banvar(1)) - 2;
				banval{1}(1:length(banstp{1}),2) = ...
				    yactCalyg(yAg-yFg+1:yAg-yFg+length(banstp{1}),banvar(1)) + 2;
			end
		end
	end


	%===================================================
	% ML conditional forecast
	%===================================================
	%[yhat,Estr,rcon] = fidencond(valuecon,stepcon,varcon,nconstr,nstepsm,nvar,lags,...
	%                                 yforeml,imf3sml,phil,Bhml,eq_ms,[]);
	                       %   The fidencond.m needs to be improved, 10/19/98
   %/*
   %  [yhat,Estr,rcon] = fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
   %           nconstr,eq_ms,nvar,lags,phil,0,0,yforeml,imf3sml,Bhml,...
   %           Cms,TLindx,TLnumber,nCms,eq_Cms);
   %  yhatml = yhat;
   %  %*** Converted to mlevle, mg, qg, and calendar yg
   %  [yhatlml,yhatmgml,yhatqgml,yhatCalygml] = fore_cal(yhat,xdata,nvar,nSample,...
   %                    nSampleCal,forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog,...
   %                    vlistper,[1 1 1 1]);
   %  yhatCalygml
   %  %
   %  %*** 1: Plot annual conditional forecast
   %  if gDLSIx & nconstr
   %     if eq_ms
   %        conlab = ['MS DLS ' ylab{PorR(1)} ' for ' num2str(nconstr) ' pers'];
   %     else
   %        conlab = ['Redu DLS ' ylab{PorR(1)} ' for ' num2str(nconstr) ' pers'];
   %     end
   %
   %     gfore(yrEnd,qmEnd,q_m,Psuedo,yearsCal,yactCalyg,forepy,...
   %                          yhatCalygml,xlab,ylab,forelabel,conlab)
   %  end
   %
   %
   %
   %  %---------------------------------------------------------------------
   %  %  Backed out exact shocks given parameter values and actual data
   %  %  In-sample shocks and, if Psuedo, out-of-sample shocks as well
   %  %---------------------------------------------------------------------
   %  if all(Sexact{1}=='Forecast') & Sexact{2}
   %     philr = phi(size(phi,1)-actup+1,:);   % +1 is absolutely needed
   %     yexa = yact(1:Sexact{2},:);
   %     Estrexa = fidcndexa(yexa,philr,A0,Bhml,nvar,lags,forep,Sexact);
   %  elseif Sexact{2}
   %     philr = phi(size(phi,1)-Sexact{2}+1,:);   % +1 is absolutely needed
   %     yexa = yact(actup-Sexact{2}+1:actup,:);
   %     Estrexa = fidcndexa(yexa,philr,A0,Bhml,nvar,lags,forep,Sexact);
   %  end
   %  ExactSmyear=[myears Estrexa];
   %  %
   %  save outXshock myears Estrexa yact qyears yactqg ...
   %           yearsCal yactCalyg A0 Bhml xdata
   %              % When parac.m is set right, it will be saved for actaul
   %              %   and exactly backed-out shocks for the entire history
   %
end
