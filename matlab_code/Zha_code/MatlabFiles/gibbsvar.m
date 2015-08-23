function A0gbs = gibbsvar(A0gbs,cT,vR,nvar,fss,kdf,Indxcol)
% A0gbs = gibbsvar(A0gbs,cT,vR,nvar,fss,kdf)
%    One-step Gibbs sampler for structural VARs -- simultaneous equations approach
%    Ref.:  D.F. Waggoner and T. Zha: "Does Normalization Matter for Inference?"
%    See Note Forecast (2) pp. 44-51
%
% A0gbs:  the last draw of A0 matrix
% cT{i}: nvar-by-nvar where T'*T=Sbd{i} which is kind of covariance martrix
%          divided by fss already
% vR{i}: nvar-by-q{i} -- orthonormral basis for T*R, which is obtained through
%          single value decomposition of Q*inv(T).  See gibbsglb.m
% nvar:  rank of A0 or # of variables
% fss:  effective sample size == nSample (T)-lags+# of dummy observations
% kdf:  polynomial power in the Gamma or 1d Wishart distribution
% Indxcol: a vector of random columns this Gibbs draws.
%           When this input is not supplied, the Gibbs draws all columns
%------------------
% A0bgs:  new draw of A0 matrix in this Gibbs step
%
% Written by Tao Zha; Copyright (c) 1999 by Waggoner and Zha
% NOTE: Added Indxcol on 2/13/00 so that previous programs may not be compatible.
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

if (nargin==6), Indxcol=[1:nvar]; end

%---------------- Local loop for Gibbs given last A0gbs ----------
%* uR{i}: nvar-by-q{i} -- orthonormal with first q(i)-1 vectors lies in the
%          span(T*a(j)|j~=i)
%*** Constructing u(1),...,u(q{i}) at each Gibbs step
%
sw0 = zeros(nvar,1);
for k=Indxcol            % given last A0gbs and general new A0bgs
   X = cT{k}*A0gbs;    % given the latest updated A0gbs
   X(:,k) = 0;    % want to find non-zero sw s.t., X'*sw=0
   [jL,Ux] = lu(X');
   jIx0 = min(find(abs(diag(Ux))<eps)); % if isempty(jIx0), then something is wrong here
   %
   sw0(jIx0+1:end) = 0;
   sw0(jIx0) = 1;
   jA = Ux(1:jIx0-1,1:jIx0-1);
   jb = Ux(1:jIx0-1,jIx0);
   jy = -jA\jb;
   sw0(1:jIx0-1) = jy;
   sw = sw0/sqrt(sum(sw0.^2));
   %
   lenk = length(vR{k}(1,:));
   gkb = zeros(lenk,1);  % greek beta's
   uR = zeros(nvar,lenk);
   sx = zeros(nvar,lenk);
   sx(:,1) = vR{k}(:,1);
   for ki = 1:lenk-1
      wxv = [sw'*sx(:,ki);sw'*vR{k}(:,ki+1)];  % w'*x and w'*v(+1)
      dwxv = sqrt(sum(wxv.^2));
      if (dwxv<eps)
         uR(:,ki)=sx(:,ki); sx(:,ki+1)=vR{k}(:,ki+1);
      else
         wxv = wxv/dwxv;
         uR(:,ki) = wxv(1)*vR{k}(:,ki+1) - wxv(2)*sx(:,ki);
         sx(:,ki+1) = wxv(2)*vR{k}(:,ki+1) + wxv(1)*sx(:,ki);
      end
   end
   uR(:,lenk) = sx(:,lenk);  % uR now constructed
   %
   %--------- Gibbs loop ----------
   %*** draw independently beta's that combine uR to form a's (columns of A0)
   jcon = sqrt(1/fss);
   gkb(1:lenk-1) = jcon*randn(lenk-1,1);
   %* gamma or 1-d Wishart draw
   jnk = jcon*randn(kdf+1,1);
   if rand(1)<0.5
      gkb(lenk) = sqrt(jnk'*jnk);
   else
      gkb(lenk) = -sqrt(jnk'*jnk);
   end
   %
   %*** form new a(i) - ith column of A0
   A0gbs(:,k) = (cT{k}\uR)*gkb;
end
