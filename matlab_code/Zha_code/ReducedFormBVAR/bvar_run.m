% This M file gives unconditional forecasts from a BVAR with the standard Sims and Zha prior.
%
% Written by Tao Zha, February 2004; Revised May 2005.



bvar_setup;
%(1)--------------------------------------
% Further data analysis
%(1)--------------------------------------
%
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
if ~(nSample==size(xdgel,1))
   warning('The sample size (including lags) and data are incompatible')
   disp('Check to make sure nSample and size(xdgel,1) are the same')
   return
end
%
baddata = find(isnan(xdgel));
if ~isempty(baddata)
   warning('Some data for this selected sample are actually unavailable.')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end
%
if qmBin==1
   yrB = yrBin; qmB = qmBin;
else
   yrB = yrBin+1; qmB = 1;
end
yrF = yrFin; qmF = qmFin;
[Mdate,tmp] = fn_calyrqm(q_m,[yrBin qmBin],[yrFin qmFin]);
xdatae=[Mdate xdd(1:nData,vlist)];
      % beyond sample into forecast horizon until the end of the data yrFin:qmFin
      % Note: may contain NaN data.  So must be careful about its use

%=========== Obtain prior-period, period-to-last period, and annual growth rates
[yactyrge,yactyre,yactqmyge,yactqmge,yactqme] = fn_datana(xdatae,q_m,vlistlog,vlistper,[yrB qmB],[yrF qmF]);
qdates = zeros(size(yactqmyge,1),1);
for ki=1:length(qdates)
   qdates(ki) = yactqmyge(1,1) + (yactqmyge(1,2)+ki-2)/q_m;
end
for ki=1:nvar
   figure
   plot(qdates, yactqmyge(:,2+ki)/100)
   xlabel(varlist{ki})
end
save outactqmygdata.prn yactqmyge -ascii



%===========  Write the output on the screen or to a file in an organized way ==============
%disp([sprintf('%4.0f %2.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',yactyrge')])
spstr1 = 'disp([sprintf(';
spstr2 = '%4.0f %2.0f';
yactyrget=yactyrge';
for ki=1:length(vlist)
   if ki==length(vlist)
      spstr2 = [spstr2 ' %8.3f\n'];
   else
      spstr2 = [spstr2 ' %8.3f'];
   end
end
spstr = [spstr1 'spstr2' ', yactyrget)])'];
eval(spstr)

%
fid = fopen('outyrqm.prn','w');
%fprintf(fid,'%4.0f %2.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',yactyrge');
fpstr1 = 'fprintf(fid,';
fpstr2 = '%4.0f %2.0f';
for ki=1:nvar
   if ki==nvar
      fpstr2 = [fpstr2 ' %8.3f\n'];
   else
      fpstr2 = [fpstr2 ' %8.3f'];
   end
end
fpstr = [fpstr1 'fpstr2' ', yactyrget);'];
eval(fpstr)
fclose(fid);



%(2)----------------------------------------------------------------------------
% Estimation
% ML forecast and impulse responses
% Hard or soft conditions for conditional forecasts
%(2)----------------------------------------------------------------------------
%
%* Arranged data information, WITHOUT dummy obs when 0 after mu is used.  See fn_rnrprior_covres_dobs.m for using the dummy
%    observations as part of an explicit prior.
[xtx,xty,yty,fss,phi,y,ncoef,xr,Bh] = fn_dataxy(nvar,lags,xdgel,mu,0,nexo);
if qmStart+lags-ndobs>0
   qmStartEsti = rem(qmStart+lags-ndobs,q_m);   % dummy observations are included in the sample.
   if (~qmStartEsti)
      qmStartEsti = q_m;
   end
   yrStartEsti = yrStart + floor((qmStart+lags-ndobs)/(q_m+0.01));
     % + 0.01 (or any number < 1)  is used so that qmStart+lags-ndobs==?*q_m doesn't give us an extra year forward.
else
   qmStartEsti = q_m + rem(qmStart+lags-ndobs,q_m);   % dummy observations are included in the sample.
   if (qmStart+lags-ndobs==0)
      yrStartEsti = yrStart - 1;   % one year back.
   else
      yrStartEsti = yrStart + floor((qmStart+lags-ndobs)/(q_m-0.01));
     % - 0.01 (or any number < 1)  is used so that qmStart+lags-ndobs==-?*q_m give us an extra year back.
   end
end
dateswd = fn_dataext([yrStartEsti qmStartEsti],[yrEnd qmEnd],xdatae(:,[1:2]));  % dates with dummies
phie = [dateswd phi];
ye = [dateswd y];


if Rform
   Ui=cell(nvar,1); Vi=cell(ncoef,1);
   for kj=1:nvar
      Ui{kj} = eye(nvar);  Vi{kj} = eye(ncoef);
   end
else
   %* Obtain linear restrictions
   eval(['[Ui,Vi,n0,np,ixmC0Pres] = ' idfile_const '(lags,nvar,nexo,0);'])
   if min(n0)==0
      disp(' ')
      warning('A0: restrictions in dlrprior.m give no free parameter in one of equations')
      disp('Press ctrl-c to abort')
      pause
   elseif min(np)==0
      disp(' ')
      warning('Ap: Restrictions in dlrprior.m give no free parameter in one of equations')
      disp('Press ctrl-c to abort')
      pause
   end
end

if indxPrior
   %*** Obtains asymmetric prior (with no linear restrictions) with dummy observations as part of an explicit prior (i.e,
   %      reflected in Hpmulti and Hpinvmulti).  See Forecast II, pp.69a-69b for details.
   [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior(nvar,q_m,lags,xdgel,mu);

   %*** Combines asymmetric prior with linear restrictions
   [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Ui,Vi,Pi,H0multi,Hpmulti,nvar);

   %*** Obtains the posterior matrices for estimation and inference
   [Pmat,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0invtld,Hpinvtld,Ui,Vi);

   if Rform
      %*** Obtain the ML estimate
      A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
      A0hat=inv(A0hatinv);
      a0indx = find(A0hat);
   else
      %*** Obtain the ML estimate
      %   load idenml
      x = 10*rand(sum(n0),1);
      H0 = eye(sum(n0));
      crit = 1.0e-9;
      nit = 10000;
      %
      tic
      [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = ...
            csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Ui,nvar,n0,fss,H0inv);
      endtime = toc

      A0hat = fn_tran_b2a(xhat,Ui,nvar,n0)
      A0hatinv = inv(A0hat);
      fhat
      xhat
      grad
      itct
      fcount
      retcodehat
      save outm endtime xhat A0hat A0hatinv grad fhat itct itct fcount retcodehat
   end
else
   %*** Obtain the posterior matrices for estimation and inference
   [Pmat,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Ui,Vi);

   if Rform
      %*** Obtain the ML estimate
      A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
      A0hat=inv(A0hatinv);
      a0indx = find(A0hat);
   else
      %*** Obtain the ML estimate
      %   load idenml
      x = 10*rand(sum(n0),1);
      H0 = eye(sum(n0));
      crit = 1.0e-9;
      nit = 10000;
      %
      tic
      [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = ...
            csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Ui,nvar,n0,fss,H0inv);
      endtime = toc

      A0hat = fn_tran_b2a(xhat,Ui,nvar,n0)
      A0hatinv = inv(A0hat);
      fhat
      xhat
      grad
      itct
      fcount
      retcodehat
      save outm endtime xhat A0hat A0hatinv grad fhat itct itct fcount retcodehat
   end
end

%**** impulse responses
swish = A0hatinv;       % each column corresponds to an equation
if Rform
   xhat = A0hat(a0indx);
   Bhat=Pmat{1};
   Fhat = Bhat*A0hat;
   ghat = NaN;
else
   xhat = fn_tran_a2b(A0hat,Ui,nvar,n0);
   [Fhat,ghat] = fn_gfmean(xhat,Pmat,Vi,nvar,ncoef,n0,np);
   Bhat = Fhat/A0hat;   % ncoef-by-nvar reduced form lagged parameters.
end
nn = [nvar lags imstp];
imfhat = fn_impulse(Bhat,swish,nn);    % in the form that is congenial to RATS
imf3hat=reshape(imfhat,size(imfhat,1),nvar,nvar);
         % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
imf3shat=permute(imf3hat,[1 3 2]);
         % imf3s: permuted so that row--steps, column--nvar shocks,
         %                                3rd dimension--nvar responses
         % Note: reshape(imf3s(1,:,:),nvar,nvar) = A0in  (columns -- equations)

%--- Graphing impulse responses.
figure
scaleout = fn_imcgraph(imfhat,nvar,imstp,xlab,ylab,1);
subtitle('Impulse responses')
imfstd = max(abs(scaleout)');   % row: nvar (largest number); used for standard deviations

%
%  %**** save stds. of both data and impulse responses in idfile1
%  temp = [std(yactqmyge(:,3:end)); std(yactyrge(:,3:end)); imfstd];  %<<>>
%  save idenyimstd.prn temp -ascii   % export forecast and impulse response to the file "idenyimstd.prn", 3-by-nvar
%  %
%  %**** save stds. of both data and impulse responses in idfile1
%  temp = [std(yactqmyge(:,3:end)); std(yactyrge(:,3:end)); imfstd];  %<<>>
%  save idenyimstd.prn temp -ascii   % export forecast and impulse response to the file "idenyimstd.prn", 3-by-nvar
%  if IndxParR
%     idfile1='idenyimstd';
%  end

%=====================================
% Now, out-of-sample forecasts. Note: Hm1t does not change with A0.
%=====================================
%
% * updating the last row of X (phi) with the current (last row of) y.
tcwx = nvar*lags;  % total coefficeint without exogenous variables
phil = phi(size(phi,1),:);
phil(nvar+1:tcwx) = phil(1:tcwx-nvar);
phil(1:nvar) = y(end,:);
%*** exogenous variables excluding constant terms
if (nexo>1)
   Xexoe = fn_dataext([yrEnd qmEnd],[yrEnd qmEnd],xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
   phil(1,tcwx+1:tcwx+nexo-1) = Xexoe(1,3:end);
end
%
%*** ML unconditional point forecast
nn = [nvar lags nfqm];
if nexo<2
   yforehat = fn_forecast(Bhat,phil,nn);    % nfqm-by-nvar, in log
else
   Xfexoe = fn_dataext(fdates(1,:),fdates(end,:),xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
   yforehat = fn_forecast(Bhat,phil,nn,nexo,Xfexoe(:,3:end));    % nfqm-by-nvar, in log
end
yforehate = [fdates yforehat];
%
yact1e = fn_dataext([yrEnd-nayr 1],[yrEnd qmEnd],xdatae(:,1:nvar+2));
if Pseudo
   %yact2e = fn_dataext([yrEnd-nayr 1],E2yrqm,xdatae);
   yact2e = fn_dataext([yrEnd-nayr 1],[fdates(end,1) q_m],xdatae(:,1:nvar+2));
else
   yact2e=yact1e;
end
yafhate = [yact1e; yforehate];  % actual and forecast
%
%===== Converted to mg, qg, and calendar yg
%
[yafyrghate,yafyrhate,yafqmyghate] = fn_datana(yafhate,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
               % actual and forecast growth rates
[yact2yrge,yact2yre,yact2qmyge] = fn_datana(yact2e,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
               % only actual growth rates
disp('--------- Unconditional out-of-sample forecasts:---------------')
yafyrghate


%---- Graphing the forecasts
keyindx = [1:nvar];
conlab=['unconditional'];

figure
yafyrghate(:,3:end) = yafyrghate(:,3:end)/100;
yact2yrge(:,3:end) = yact2yrge(:,3:end)/100;
fn_foregraph(yafyrghate,yact2yrge,keyindx,rnum,cnum,q_m,ylab,forelabel,conlab)



