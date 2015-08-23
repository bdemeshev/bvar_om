function [pabove,pbelow]=sbcontest(xinput)
%  [pabove,pbelow]=sbcontest(xinput)
%    Small Sample Bayesian Contour Test (sbcontest) for overidentified models
%
% xinput{1}: idenmlh -- handle for the file idenml.mat
% xinput{2}: SpHUouth -- handle for the file sphuout.mat
% xinput{3}: fss -- effective sample size == nSample-lags+# of dummy observations
% xinput{4}: imndraws=nstarts*ndraws2
% xinput{5}: nvar -- # of endogenous variables
% xinput{6}: nbuffer -- interval of whcih for printing, plotting, saving, etc.
%--------------
% pabove:  probablity of the reduced-form LH value above the overidentified LH peak value
% pbelow:  probablity of the reduced-form LH value below the overidentified LH peak value
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

idenmlh = xinput{1}; SpHUouth = xinput{2}; fss = xinput{3}; imndraws = xinput{4}; nvar = xinput{5};
nbuffer = xinput{6};

neqn = nvar;   %<<>> number of equations
eval(['load ' idenmlh]);
eval(['load ' SpHUouth]);

pabove=0;
pbelow=0;
load idenml    % fhat at ML for identified version
fhat=-fhat;
      % this value is strictly smaller than that at the peak of reduced-form LH
disp('If you havent run idenml with Rform=1, you must do so now')
disp('Press ctrl-c to abort now or any other key to continue')
disp(' ')
pause
load SpHUout    % this is the output file by running "idenml" with Rform=1
                % with output SpHU
SpHs = SpHU*fss;     % for the Wishart draw
SpHsc = chol(SpHs);     % upper triangular where Sphsc' is
                        %  lower triangular Choleski, for Wishart draw
tic
for draws=1:imndraws
   if ~mod(draws,nbuffer)
      draws
   end
   ranw = randn(neqn,fss-neqn-1);
   ranw = SpHsc\ranw;   % inv(SpHsc) is upper triagular Choleski of inv(SpHs)
   %ranw = SpHsic*ranw;
   % normal draws (T-nvar-1)*nvar, with variance inv(SpHs) (note, NOT inv(SpH))
   sinh = ranw*ranw';    % Wishart draws for A0h*A0h'
   A0_h = chol2(sinh);   % upper triangular Choleski
   a0indx = find(A0_h);

   xhat_h = A0_h(a0indx);
   jnk = a0lhfun(xhat_h,SpHU,fss,nvar,a0indx);
   fhat_h = -jnk;

   if fhat_h > fhat
      pabove = pabove+1;
   else
      pbelow = pbelow+1;
   end
end
timend=toc;
timeminutes=timend/60

pabove = pabove/imndraws
pbelow = pbelow/imndraws
