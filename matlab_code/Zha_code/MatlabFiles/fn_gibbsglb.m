function [cT,vR,kdf] = fn_gibbsglb(Sbd,idmat0,nvar,fss)
% [cT,vR,kdf] = gibbsglb(Sbd,idmat0,nvar,fss)
%    Global setup outside the Gibbs loop -- c.f. gibbsvar
%    Ref.:  D.F. Waggoner and T.A. Zha: "Does Normalization Matter for Inference?"
%    See Note Forecast (2) pp. 44-51
%
% Sbd: cell(nvar,1). Sbd=diag(S(1), ..., S(m)).  Already divided by 'fss' for the
%        posterior of a0 or A0(:) in Waggoner and Zha when one has asymmetric prior.
%        Note,"bd" stands for block diagonal.
% idmat0:  identification matrix for A0 with asymmetric prior;  column -- equation.
% nvar:  rank of A0 or # of variables
% fss: effective sample size (in the exponential term) --
%            # of observations + # of dummy observations
%                                   (or nSample - lags + # of dummy observations)
%-------------
% cT{i}: nvar-by-nvar where T'*T=Sbd{i} which is kind of covariance martrix
%          divided by fss already
% vR{i}: nvar-by-q{i} -- orthonormral basis for T*R, which is obtained through
%          single value decomposition of Q*inv(T)
% kdf:  the polynomial power in the Gamma or 1d Wishart distribution, used in
%          "gibbsvar.m"
%
% Written by Tao Zha; Copyright (c) 1999 by Waggoner and Zha
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


kdf = fss;  %2;     %fss;

cT = cell(nvar,1);
for k=1:nvar
   cT{k} = chol(Sbd{k});   % upper triangular but lower triangular Choleski
end

Q = cell(nvar,1);
       % Q{i}: nvar-by-nvar with rank nvar-q(i) with ith equation a and
       %       q(i) which is # of non-zero elements in a.  Note: Q*a=0.

vR = cell(nvar,1);
for k=1:nvar
   Q{k} = diag(idmat0(:,k)==0);   % constructing Q{k}
   %
   [jnk1,d1,v1] = svd(Q{k}/cT{k});
   d1max=max(diag(d1));
   if d1max==0
      Idxk=1:nvar;
   else
      Idxk = find(diag(d1)<eps*d1max);
   end
   lenk = length(find(idmat0(:,k)));
   if ( length(Idxk)<lenk )
      warning('Dimension of non-zero A0(:,k) is different from svd(Q*inv(T))')
      disp('Press ctrl-c to abort')
      pause
   else
      jnk1 = v1(:,Idxk);
      vR{k} = jnk1(:,1:lenk);   % orthonormal basis for T*R
   end
end
