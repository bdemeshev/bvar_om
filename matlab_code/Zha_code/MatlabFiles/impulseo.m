function imf = impulseo(Bh,swish,nn)
% impulse: computing impulse functions with
%                imf = impulseo(Bh,swish)
%  where imf is in a format that is conveninet to generate the RATS graphics;
%           Bh is the estimated reduced form coefficient in the form 
%               Y(T*nvar) = XB + U, X: T*k, B: k*nvar.  The matrix
%               form or dimension is the same as "Bh" from the function "sye";
%           swish is the inv(A0) in the structural model A(L)y(t) = e(t).
%           nn is the numbers of inputs [nvar,lags,# of impulse responses].
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
imstep = nn(3);   % number of steps for impulse responses

Ah = Bh';        
% Row: nvar equations 
% Column: 1st lag (with nvar variables) to lags (with nvar variables) + const = k.

imf = zeros(imstep,nvar*nvar);        
% Column: 1st variable response to nvar shocks, 2nd variable response to 
%         nvar shocks, and so on.  
% Row:  steps of impulse responses. 
M = zeros(nvar*(lags+1),nvar);
% Stack M0;M1;M2:...;Mlags
M(1:nvar,:) = swish;
Mtem = M(1:nvar,:);    % temporary M.  
% first (initial) responses to 1 standard deviation shock.  Row: responses; Column: shocks
% * put in the form of "imf"
Mtem = Mtem';
imf(1,:) = Mtem(:)';

t = 1;
ims1 = min([imstep-1 lags]);
while t <= ims1
   Mtem = zeros(nvar,nvar);
   for k = 1:t
      Mtem = Ah(:,nvar*(k-1)+1:nvar*k)*M(nvar*(t-k)+1:nvar*(t-k+1),:) + Mtem;
      % Row: nvar equations, each for the nvar variables at tth lag
   end
   M(nvar*t+1:nvar*(t+1),:) = Mtem;
   Mtem = Mtem';
   imf(t+1,:) = Mtem(:)';   
   % stack imf with each step, Row: 1st variable to nvar shocks, 
   %                                2nd variable to nvar shocks, etc.
   t= t+1;
end

for t = lags+1:imstep-1
   M(1:nvar*lags,:) = M(nvar+1:nvar*(lags+1),:);
   Mtem = zeros(nvar,nvar);
   for k = 1:lags
      Mtem = Ah(:,nvar*(k-1)+1:nvar*k)*M(nvar*(lags-k)+1:nvar*(lags-k+1),:) + Mtem;
      % Row: nvar equations, each for the nvar variables at tth lag
   end
   M(nvar*lags+1:nvar*(lags+1),:) = Mtem;
   Mtem = Mtem';
   imf(t+1,:) = Mtem(:)';   
   % stack imf with each step, Row: 1st variable to nvar shocks, 
   %                                2nd variable to nvar shocks, etc.
end
