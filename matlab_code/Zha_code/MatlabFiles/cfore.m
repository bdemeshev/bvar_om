% 10/24/97
% Distance Method of Waggoner and Zha
% Modified from Sims and Zha's code
% Copyright (c) 1997 Tao Zha
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
yrBin=59;   % beginning of the year
qmBin=1;    % begining of the quarter or month
yrFin=97;   % final year
qmFin=12;    % final quarter
%tnvar = 2;   % total number of variables
nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
% total number of the available data -- this is all you have
%
%* Load data and series
%
load xd24a      % the default name for the variable is "xdd".
[nt,ndv]=size(xdd);
if nt~=nData
   warning(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   %disp(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   return
end
%1  CPI-U
%2  FFR
%3  T-bill3
%4  Treasure note 10
%5  M2
%6  M1
%7  Nominal PCE
%8  real PCE
%9  Unemployment Rate
%10 IMF Commodity Price Index
%11 Civilians Employed: 16 & over
%12 Nonfarm Payroll Employment
%13 IP
%14 Retail Sales (Nominal)
%15 NAPM Composit Index
%16 Total Reserves
%17 PPI-finished goods
%18 PPI-Crude materials
%19 PPI-Crude materials less energy
%20 CRB Spot Commodity Index -- all commodities
%21 CRB Spot Commodity Index -- raw industrials
%22 PCE price index
%23 rgdpmon Real GDP (monthly, chain $92)
%24 dgdpmon Deflator GDP (monthly, chain $92)


%1  IMF CP
%2  M2
%3  FFR
%3  real GDP
%4  CPI-U
%5  U

%>>>>>>>>>>>>>>>>>>
logindx = [1 5:8 10:14 16:24];
xdd(:,logindx) = log(xdd(:,logindx));
pctindx = [2:4 9 15];
xdd(:,pctindx)=.01*xdd(:,pctindx);       % make it a general term for the following use
%
vlist = [10 5 2 23 1 9];    % regarding "xdd", IMF-CP, M2, FFR, GDP, CPI, U
vlistlog = [1 2 4 5];       % subset of "vlist"
vlistper = [3 6];           % subset of "vlist"
%<<<<<<<<<<<<<<<<<<<

idfile='iden6';

xlab = {'Inf'
        'MS'
        'FFR'
        'y'
        'P'
        'U'};

ylab = {'Pcm'
        'M2'
        'FFR'
        'y'
        'P'
        'U'};

xdd_per = xdd(:,vlist);

x1 = 'Pcm';
x2 = 'M2';
x3 = 'FFR';
x4 = 'GDP';
x5 = 'CPI';
x6 = 'U';
x7 = 'R10';


baddata = find(isnan(xdd_per));
if ~isempty(baddata)
   warning('Some data are actually unavailable.')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end
%

%* A specific sample is considered for estimation
%*   Sample period 59:7-82:9, forecast period 82:10-84:9
yrStart=59;
qmStart=1;
[yrEnd,qmEnd,forep,forepy,forelabel] = pararc;
nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmStart+1);
if qmEnd == q_m     % end of the year
   nSampleCal=nSample;            % Cal: calendar year
else
   nSampleCal=(yrEnd-1-yrStart)*q_m + (q_m-qmStart+1);   % Cal: calendar year
end

%* More script variables
%
lags = 13;        % <<>>
%  automatic decay code (monthly data), only two options: lags = 6 or 13
forepq = forep/3;      % quarterly
actup = 5*48;     % <<>> actual periods before forecasting (20 years)
%actup = 12*floor(nSample/12);     % <<>> actual periods before forecasting (8 years)
actup = 48;     % <<>> actual periods before forecasting (4 years)
actupq = actup/3;   % quarterly
actupy = actup/12;   % four years
imstp = 48;      % <<>>  impulse responses (4 years)
ninv = 1000;   % the number of intervals for counting impulse responses
nhp = 6;          % <<>> number of hyperparameters for estimation
%%scf = 2.4/sqrt(nvar);       % scf^2*Sigma (covaraince)
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

[Gb,Sbd,Bh,SpH,fss,ndobs,phi,y,nvar,ncoef,xxhpc,a0indx,na0p,...
             idmat0,idmatpp] = szasbvar(idfile,q_m,lags,nSample,nhp,xdgel);

% * the largest matrix in this file  <<>>
yforew = zeros(ndraws,forep*nvar);     % preallocating
yforeqgw = zeros(ndraws,forepq*nvar);     % preallocating
yforeCalygw = zeros(ndraws,forepy*nvar);     % preallocating
% * the largest matrix in this file  <<>>
%%imfcnt = zeros(ninv+2,imstp*nvar^2);   % cnt: count

load idenml    % xhat ghat fhat, etc.
%load outiden    % xhat ghat fhat

%==================
%  Impulse responses first
%==================
%A0 = zeros(nvar);
%A0(a0indx)=xhat;
%A0(4,2) = -xhat(7);   % output in MD
%A0(5,2)=-xhat(7);   % price in MD
A0in = inv(A0);
swish = A0in';       % each row corresponds to an equation
%Bh = Hm1t;    % no longer have Hm1t in this new szasbvar

% ** impulse responses
nn = [nvar lags imstp];
imf = zimpulse(Bh,swish,nn);    % in the form that is congenial to RATS
%[vd,str,imf] = errors(Bh,swish,nn);

scaleout = imcgraph(imf,nvar,imstp,xlab,ylab)


%%%%
%$$$ Out-of-sample forecasts. Note: Hm1t does not change with A0.
%%%%
%
% ** out-of-sample forecast, from 82:4 to 84:3 (flp+1:flp+forep)
% * updating the last row of X (phi) with the current (last row of) y.
phil = phi(size(phi,1),:);
phil(nvar+1:ncoef-1) = phil(1:ncoef-1-nvar);
phil(1:nvar) = y(size(y,1),:);
ylast = y(size(y,1),:);
indx12 = size(y,1)-q_m+1:size(y,1);
ylast12 = y(indx12,:);        % last 12 months data
nn = [nvar lags forep];
%
yfore = forecast(Bh,phil,nn);    % forep-by-nvar


%>>>>>>>>>>>>>>>
yforel=yfore;
yforel(:,vlistlog) = exp(yfore(:,vlistlog));
figure;
t2=1:forep;
for i = 1:nvar
   subplot(nvar/2,2,i)
   plot(t2,yforel(:,i),'--')
   %title('solid-actual, dotted-forecast');
   %title(eval(['forelabel']));
   %ylabel(eval(['x' int2str(i)]));
	ylabel(char(ylab(i)))
end
%<<<<<<<<<<<<<<<

%%%%%%%%
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
%%%%%%%%

nconstr=4;   % q: 4 years -- 4*12 months
eq_ms = 2;      % location of MS equation
%eq_ms = [];     % all shocks
%*** initializing
stepcon=cell(nconstr,1);  % initializing, value y conditioned
valuecon=zeros(nconstr,1);  % initializing, value y conditioned
varcon=zeros(nconstr,1);  % initializing, endogous variables conditioned
%
stepcon{1}=[1:12]';    % average over 12 months.
stepcon{2}=[13:24]';    % average over 12 months.
stepcon{3}=[25:36]';    % average over 12 months.
stepcon{4}=[37:48]';    % average over 12 months.
%
%for i=1:nconstr
%   stepcon{i}=i;
%end
%
chk1 = mean(yfore(stepcon{1},3))
chk2 = mean(yfore(stepcon{2},3))
chk3 = mean(yfore(stepcon{3},3))
chk4 = mean(yfore(stepcon{4},3))
Ro=[chk1 chk2 chk3 chk4];
%
chk1 = exp( (sum(yfore(stepcon{1},5))-sum(ylast12(:,5))) ./ q_m )
chk2 = exp( (sum(yfore(stepcon{2},5))-sum(yfore(stepcon{1},5))) ./ q_m )
chk3 = exp( (sum(yfore(stepcon{3},5))-sum(yfore(stepcon{2},5))) ./ q_m )
chk4 = exp( (sum(yfore(stepcon{4},5))-sum(yfore(stepcon{3},5))) ./ q_m )
%
%valuecon(:)=0.055;
%
%>>>>>>>>>>>>>>>>> E: Condition on funds rate path >>>>>>>>>>>>> Toggle
%delta=0.0010;
%valuecon(1) = mean(yfore(stepcon{1},3))+2*delta;
%valuecon(2) = mean(yfore(stepcon{2},3))+2*delta;
%valuecon(3) = mean(yfore(stepcon{3},3))-2*delta;
%valuecon(4) = mean(yfore(stepcon{4},3))-2*delta;
%valuecon(1) = 0.055;
%valuecon(2) = 0.050;
%valuecon(3) = 0.0475;
%valuecon(4) = 0.045;
%>>>>>>>>>>>>>>>>> E: Condition on funds rate path >>>>>>>>>>>>>
%
%<<<<<<<<<<<<<<<< B: Condition on inflation path <<<<<<<<<< Toggle
%delta=0.0010;
%valuecon(1)=mean(ylast12(:,5))+log(chk1-0*delta);
%valuecon(2)=valuecon(1)+log(chk2-2*delta);
%valuecon(3)=valuecon(2)+log(chk3-6*delta);
%valuecon(4)=valuecon(3)+log(chk4-12*delta);
%        % 5: CPI; 2.5%: annual inflation over 12 months, geometric means
%$$$ very good results -- following
valuecon(1)=mean(ylast12(:,5))+log(chk1);
valuecon(2)=valuecon(1)+log(1.020);
valuecon(3)=valuecon(2)+log(1.02);
valuecon(4)=valuecon(3)+log(1.02);
%        % 5: CPI; 2.5%: annual inflation over 12 months, geometric means
%<<<<<<<<<<<<<<<< E: Condition on inflation path <<<<<<<<<<

nstepsm = 0;   % initializing, the maximum step in all constraints
for i=1:nconstr
   nstepsm = max([nstepsm max(stepcon{i})]);
end
varcon(:)=5;     % 3: FFR; 5: CPI
%
imf3=reshape(imf,size(imf,1),nvar,nvar);
         % imf3: row-steps, column-nvar responses, 3rd dimension-nvar shocks
imf3s=permute(imf3,[1 3 2]);
         % imf3s: permuted so that row-steps, column-nvar shocks,
			%                                       3rd dimension-nvar responses

[yhat,Estr] = fidencond(valuecon,stepcon,varcon,nconstr,nstepsm,nvar,lags,...
                                 yfore,imf3s,phil,Bh,eq_ms);

chk1 = mean(yhat(stepcon{1},3))
chk2 = mean(yhat(stepcon{2},3))
chk3 = mean(yhat(stepcon{3},3))
chk4 = mean(yhat(stepcon{4},3))
Rh=[chk1 chk2 chk3 chk4];

chk1 = exp( (sum(yhat(stepcon{1},5))-sum(ylast12(:,5))) ./ q_m )
chk2 = exp( (sum(yhat(stepcon{2},5))-sum(yhat(stepcon{1},5))) ./ q_m )
chk3 = exp( (sum(yhat(stepcon{3},5))-sum(yhat(stepcon{2},5))) ./ q_m )
chk4 = exp( (sum(yhat(stepcon{4},5))-sum(yhat(stepcon{3},5))) ./ q_m )

%chk1 = mean(yhat(1:12,3))
%chk2 = mean(yhat(13:24,3))
%chk3 = mean(yhat(25:36,3))
%chk4 = mean(yhat(36:48,3))
%Rh=[chk1 chk2 chk3 chk4];

%chk1 = exp( (sum(yhat(1:12,5))-sum(ylast12(:,5))) ./ q_m )
%chk2 = exp( (sum(yhat(13:24,5))-sum(yhat(1:12,5))) ./ q_m )
%chk3 = exp( (sum(yhat(25:36,5))-sum(yhat(13:24,5))) ./ q_m )
%chk4 = exp( (sum(yhat(37:48,5))-sum(yhat(25:36,5))) ./ q_m )


idiff=mean(yfore(:,3))-mean(yhat(:,3))
mean(yfore(:,3))
mean(yhat(:,3))
figure
plot(1:48,yfore(:,3),1:48,yhat(:,3),':')
figure
plot(1:4,Ro,1:4,Rh,':')

%>>>>>>>>>>>>>>>
yhatl=yhat;
yhatl(:,vlistlog) = exp(yhat(:,vlistlog));
figure;
t2=1:forep;
for i = 1:nvar
   subplot(nvar/2,2,i)
   plot(t2,yhatl(:,i),'--')
   %title('solid-actual, dotted-forecast');
   %title(eval(['forelabel']));
   %ylabel(eval(['x' int2str(i)]));
	ylabel(char(ylab(i)))
end
%<<<<<<<<<<<<<<<


% inputs needed.
%yfore=yhat;

%===================================================
%%% Converting to calendar years and all at level
%===================================================
[yforeml,yforeqgml,yforeCalygml] = fore_cal(yhat,xdata,nvar,nSample,...
                  nSampleCal,forep,forepq,forepy,q_m,qmEnd,vlist,vlistlog);

%==================
%  Note
%=================
% ? 1--median; l--lower bound; h--upper bound: between the bound: 2/3 probability
% yfore?         % monthly, level
% yforeqg?       % quarterly, growth rate
% yforeCalyg?    % calendar year, growth rate

%save outw ndraws yfore1 yforeqg1 yforeCalyg1 yforel yforeqgl yforeh yforeqgh ...
%                        yforeCalygl yforeCalygh IUbeta

yforeCalygml


%----------------------------------------------------------------------
%=================  Graphics =====================
%----------------------------------------------------------------------
%

%[yactCalyg,yforeCalygml,yAg,yFg] = fore_gh(xdata,nvar,nSample,nSampleCal,...
%               yforeml,yforeqgml,yforeCalygml,actup,actupq,actupy,...
%		         vlist,vlistlog,vlistper,q_m,forep,ylab);

xinput = cell(16,1);
xinput{1}=xdata; xinput{2}=nvar; xinput{3}=nSample; xinput{4}=nSampleCal;
xinput{5}=yforeml; xinput{6}=yforeqgml; xinput{7}=yforeCalygml; xinput{8}=actup;
xinput{9}=actupq; xinput{10}=actupy; xinput{11}=vlist; xinput{12}=vlistlog;
xinput{13}=vlistper; xinput{14}=q_m; xinput{15}=forep; xinput{16}=ylab;
[yactCalyg,yforeCalygml,yAg,yFg] = fore_gh(xinput);


% Key Macroeconomic Variables: GDP, CPI, U
%******* From Goldbook July 1-2, 1997 FOMC
yforeBlue = [
  3.6  2.4  5.0
  2.5  2.6  5.0
      ];    % real GDP, CPI-U, U,
yforeMacro = [
  3.7  2.3  4.9
  2.1  2.4  4.9
      ];    % real GDP, CPI-U, U
yforeGold = [
  3.8  2.4  5.0
  2.5  2.5  4.8
  2.4  2.7  4.7
      ];    % real GDP, CPI-U, U
t3 = yAg+1:yAg+length(yforeBlue(:,1));
t4 = yAg+1:yAg+length(yforeGold(:,1));
%
keyindx = [4:nvar 3 2];      %  GdP, CPI, U, FFR, M2
count=0;
t1 = 1:yAg;
t2 = yAg:yAg+yFg;
for i = keyindx
   count = count+1;
   subplot(3,2,count)
   %plot(t1,yactqg(:,i),t2,yforeqg(:,i),':')
   if (i==3) | (i==2)
      plot(t1,yactCalyg(:,i),t2,[yactCalyg(length(t1),i);yforeCalygml(:,i)],'--')
      %title('solid-actual, dotted-forecast');
      %xlabel(eval(['forelabel']));
		%ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
   else
      plot(t1,yactCalyg(:,i),t2,[yactCalyg(length(t1),i);yforeCalygml(:,i)],'--',...
               t3,yforeBlue(:,count),'o',t3,yforeMacro(:,count),'d',...
               t4,yforeGold(:,count),'^')
               %title('solid-actual, dotted-forecast');
      %xlabel(eval(['forelabel']));
      %ylabel(eval(['x' int2str(i)]));
		ylabel(char(ylab(i)))
   end
end

actual=yactCalyg(:,keyindx)     %  GDP, CPI, U, M2
mode=yforeCalygml(:,keyindx)   %  GDP, CPI, U, M2
%low=yforeCalygl(:,keyindx)
%high=yforeCalygh(:,keyindx)
yforeBlue
yforeMacro
yforeGold