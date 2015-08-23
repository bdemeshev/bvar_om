function [hd,ehat] = history(uhat,Bh,A0,nn)
% Computing historical decompositions with
%                [hd,ehat] = history(uhat,Bh,A0,nn)
%   where hd: historical decompositions, dstp-by-nvar^2 matrix.
%                Column: 1st variable decompositions attributed to nvar 
%                  (cumulative) shocks, 2nd variable decompositions 
%                  attributed to nvar (cumulative) shocks, and so on.  
%                Row:  steps of decompositions "dstp".
%         ehat:  structural one-step forecast errors 
%                                    of each variable, dstp-by-nvar.
%         uhat: dstp-by-nvar reduced form one-step forecast errors 
%                                    begining with the decomposition period.
%         Bh: the estimated reduced form coefficient in the form 
%             Y(T*nvar) = XB + U, X: T*k, B: k*nvar.  The matrix
%             form or dimension is the same as "Bh" from the function "sye".
%         A0: comtemporaneous A in the structural model A(L)y(t) = e(t).
%         nn: # of inputs -- [nvar,lags,dstp (steps of decompositions)].
%    ===== Note: this function can be easily modified to get variance
%    =====       decompositions and impulse responses.
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

nvar = nn(1);
lags = nn(2);
dstp = nn(3);   % number of steps for decompositions
imstep = dstp;   % steps of impulse responses
tcwc = nvar*lags;     % total coefficients without constant
swish = inv(A0);


Ah = Bh';        
% Row: nvar equations 
% Column: 1st lag (with nvar variables) to 
%                      lags (with nvar variables) + const = k.
ehat = uhat*A0';
% dstp-by-nvar matrix:  structural one-step forecast errors.

hd = zeros(dstp,nvar*nvar);
% Column: 1st varaible to nvar shocks, 2nd variables to nvar shock, and so on.  
% Row:  steps of decompositions.
imf = zeros(dstp,nvar*nvar);          % impulse responses
% different format from "hd", i.e., nvar variables to 1st shock,
%      nvar variables to 2nd shock, etc.  To make the format the same
%      as "hd", just change Mtem to Mtem' right before stacking "imf".
%      See also "impulseo.m" in \toolbox\cstz
M = zeros(nvar*(lags+1),nvar);    % stacked matrices for impulse responses
% Stack M0;M1;M2;...;Mlags
MH = zeros(nvar*dstp,nvar);   
% same as M, but stacked in a different fashion after t > lags so that
%       all the cumulations are recorded for historical decompositions.
M(1:nvar,:) = swish;
Mtem = M(1:nvar,:);    % temporary M -- impulse responses.  
%
Hdc = M(1:nvar,:)*diag(ehat(1,:));  
% first period historical decompositions.
% * put in the form of "hd", as well as "imf"
Hdc = Hdc';
hd(1,:) = Hdc(:)';
imf(1,:) = Mtem(:)';   


%
% ** beginning with the second period.  Note, 1st period is t=0
%
t = 1;
ims1 = min([dstp-1 lags]);
while t <= ims1
   % ** impulse response functions
   Mtem = zeros(nvar,nvar);
   for k = 1:t
      Mtem = Ah(:,nvar*(k-1)+1:nvar*k)*M(nvar*(t-k)+1:nvar*(t-k+1),:) + Mtem;
      % Row: nvar equations, each for the nvar variables at tth lag
   end
   M(nvar*t+1:nvar*(t+1),:) = Mtem;
   imf(t+1,:) = Mtem(:)';   
   % stack imf with each step, Row: 6 var to 1st shock, 6 var to 2nd shock, etc.  
   % ** historical decompositions
   Hdc = M(1:nvar,:)*diag(ehat(t+1,:));  
   for k = 1:t
      Hdc = Hdc + M(nvar*k+1:nvar*(k+1),:)*diag(ehat(t+1-k,:));
      % Row: nvar variables; Column: nvar shocks
   end
   Hdc = Hdc';
   hd(t+1,:) = Hdc(:)';
   t= t+1;
end

MH(1:nvar*(lags+1),:) = M;
for t = lags+1:imstep-1
   % ** impulse response functions
   M(1:nvar*lags,:) = M(nvar+1:nvar*(lags+1),:);
   Mtem = zeros(nvar,nvar);
   for k = 1:lags
      Mtem = Ah(:,nvar*(k-1)+1:nvar*k)*M(nvar*(lags-k)+1:nvar*(lags-k+1),:) + Mtem;
      % Row: nvar equations, each for the nvar variables at tth lag
   end
   M(nvar*lags+1:nvar*(lags+1),:) = Mtem;
   imf(t+1,:) = Mtem(:)';   
   % stack imf with each step, Row: 6 var to 1st shock, 6 var to 2nd shock, etc.
   % ** historical decompositions
   MH(nvar*t+1:nvar*(t+1),:) = Mtem;
   Hdc = MH(1:nvar,:)*diag(ehat(t+1,:));  
   for k = 1:t
      Hdc = Hdc + MH(nvar*k+1:nvar*(k+1),:)*diag(ehat(t+1-k,:));
      % Row: nvar variables; Column: nvar shocks
   end
   Hdc = Hdc';
   hd(t+1,:) = Hdc(:)';
end
