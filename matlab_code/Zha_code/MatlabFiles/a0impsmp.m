function [xdraw,mphv,sm,timeminutes] = a0impsmp(xinput)
% [xdraw,mphv,sm,timeminutes] = a0impsmp(xinput)
%        Export the simulated pdfs of draws of only selected parameters
%           with importance sampling techniques
%
% xinput{1}: nfp -- total number of free parameters
% xinput{2}: nvar -- number of variables
% xinput{3}: xhat -- ML estimate of free parameters in A0
% xinput{4}: hess1 -- Hessian of -logLH for importance sampling
% xinput{5}:Indxv -- index for selected variables of interest; normall first 2 are of our interest
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
%------------------
% xdraw: imndraws-by-length(Indxv) matrix of draws;
%        empty if Indxv is empty.
% timeminutes:  minutes used for this simulation
% mphv:  imndraws-by-1 vector of unscaled weights
% sm:   sum of the weights
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

nfp=xinput{1}; nvar=xinput{2}; xhat=xinput{3}; hess1=xinput{4}; Indxv=xinput{5};
imndraws=xinput{6}; a0indx=xinput{7}; tdf=xinput{8}; nbuffer=xinput{9};
Sbd=xinput{10}; scf=xinput{11}; Hsr1=xinput{12}; fss=xinput{13};


sm = 0.0;
mphv=zeros(imndraws,1);
if isempty(Indxv)
   xdraw=[];
else
   xdraw=zeros(imndraws,length(Indxv));   %<<>> draws of selected variables
end

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
   %
end
timeminutes = toc/60
