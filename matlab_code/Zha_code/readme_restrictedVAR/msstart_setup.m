% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')


%===========================================
% Exordium I
%===========================================
format short g     % format
%
q_m = 4;   % quarters or months
yrBin=1970;   % beginning of the year
qmBin=1;    % begining of the quarter or month
yrFin=2005;   % final year
qmFin=4;    % final month or quarter
nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
       % total number of the available data -- this is all you have

%*** Load data and series
load data_dsge_jeske.prn     % the default name for the variable is "xdd".
xdd = data_dsge_jeske;
clear data_dsge_jeske;
[nt,ndv]=size(xdd);
if nt~=nData
   disp(' ')
   warning(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   %disp(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   disp('Press ctrl-c to abort')
   return
end
%--------
%1   energy price (p),
%2   durables investment (id),
%3   fixed investment (ik)
%4   output (y)

vlist = [1:4];    % 1: p; 2: id; 3: ik; 4: y.
varlist={'P', 'I_d','I_k','y'};
vlistlog = [ ];   % subset of "vlist.  Variables in log level so that differences are in **monthly** growth, unlike R and U which are in annual percent (divided by 100 already).
vlistper = [1:4];           % subset of "vlist"
idfile_const='ftd_cholesky';   %Only used by msstart2.m.
ylab = varlist;
xlab = varlist;

%----------------
nvar = length(vlist);   % number of endogenous variables
nlogeno = length(vlistlog)  % number of endogenous variables in vlistlog
npereno = length(vlistper)  % number of endogenous variables in vlistper
if (nvar~=(nlogeno+npereno))
   disp(' ')
   warning('Check xlab, nlogeno or npereno to make sure of endogenous variables in vlist')
   disp('Press ctrl-c to abort')
   return
elseif (nvar==length(vlist))
   nexo=1;    % only constants as an exogenous variable.  The default setting.
elseif (nvar<length(vlist))
   nexo=length(vlist)-nvar+1;
else
   disp(' ')
   warning('Make sure there are only nvar endogenous variables in vlist')
   disp('Press ctrl-c to abort')
   return
end



%------- A specific sample is considered for estimation -------
yrStart=1970;
qmStart=1;
yrEnd=2005;
qmEnd=4;
nfyr = 4;   % number of years for forecasting
if nfyr<1
   error('To be safe, the number of forecast years should be at least 1')
end
ystr=num2str(yrEnd);
forelabel = [ ystr(3:4) ':' num2str(qmEnd) ' Forecast'];

nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmStart+1);
if qmEnd==q_m
   E1yrqm = [yrEnd+1 1];  % first year and quarter (month) after the sample
else
   E1yrqm = [yrEnd qmEnd+1];  % first year and quarter (month) after the sample
end
E2yrqm = [yrEnd+nfyr qmEnd];   % end at the last month (quarter) of a calendar year after the sample
[fdates,nfqm]=fn_calyrqm(q_m,E1yrqm,E2yrqm);   % forecast dates and number of forecast dates
[sdates,nsqm] = fn_calyrqm(q_m,[yrStart qmStart],[yrEnd qmEnd]);
   % sdates: dates for the whole sample (including lags)
if nSample~=nsqm
   warning('Make sure that nSample is consistent with the size of sdates')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end
imstp = 4*q_m;    % <<>>  impulse responses (4 years)
nayr = 4; %nfyr;  % number of years before forecasting for plotting.


%------- Prior, etc. -------
lags = 4        % number of lags
indxC0Pres = 0;   % 1: cross-A0-and-A+ restrictions; 0: idfile_const is all we have
            % Example for indxOres==1: restrictions of the form P(t) = P(t-1).
Rform = 0;  % 1: contemporaneous recursive reduced form; 0: restricted (non-recursive) form
Pseudo = 0;  % 1: Pseudo forecasts; 0: real time forecasts
indxPrior = 1;  % 1: Standard Sims-Zha Bayesian prior; 0: no prior
indxDummy = indxPrior;  % 1: add dummy observations to the data; 0: no dummy added.
ndobs = 0;  % No dummy observations for xtx, phi, fss, xdatae, etc.  Dummy observations are used as an explicit prior in fn_rnrprior_covres_dobs.m.
%if indxDummy
%   ndobs=nvar+1;         % number of dummy observations
%else
%   ndobs=0;    % no dummy observations
%end
%=== The following mu is effective only if indxPrior==1.
mu = zeros(6,1);   % hyperparameters
mu(1) = 1;
mu(2) = 0.5;
mu(3) = 0.1;
mu(4) = 1.0;
mu(5) = 1.0;
mu(6) = 1.0;
%--- Default value.
%mu(1) = 1;
%mu(2) = 0.5;
%mu(3) = 0.1;
%mu(4) = 1.0;
%mu(5) = 1.0;
%mu(6) = 1.0;
%   mu(1): overall tightness and also for A0;
%   mu(2): relative tightness for A+;
%   mu(3): relative tightness for the constant term;
%   mu(4): tightness on lag decay;  (1)
%   mu(5): weight on nvar sums of coeffs dummy observations (unit roots);
%   mu(6): weight on single dummy initial observation including constant
%           (cointegration, unit roots, and stationarity);
%
%
hpmsmd = [0.0; 0.0];
indxmsmdeqn = [1; 2; 1; 2];


tdf = 3;          % degrees of freedom for t-dist for initial draw of the MC loop
nbuffer = 100;        % a block or buffer of draws (buffer) that is saved to the disk (not memory)
ndraws1=1*nbuffer;         % 1st part of Monte Carlo draws
ndraws2=10*ndraws1         % 2nd part of Monte Carlo draws
seednumber = 0; %7910;    %472534;   % if 0, random state at each clock time
           % good one 420 for [29 45], [29 54]
if seednumber
   randn('state',seednumber);
   rand('state',seednumber);
else
   randn('state',fix(100*sum(clock)));
   rand('state',fix(100*sum(clock)));
end
%  nstarts=1         % number of starting points
%  imndraws = nstarts*ndraws2;   % total draws for impulse responses or forecasts
%<<<<<<<<<<<<<<<<<<<





