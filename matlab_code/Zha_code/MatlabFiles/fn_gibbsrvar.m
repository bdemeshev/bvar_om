function [A0gbs, Wcell] = fn_gibbsrvar(A0gbs,UT,nvar,fss,n0,Indxcol)
% [A0gbs, Wcell] = fn_gibbsrvar(A0gbs,UT,nvar,fss,n0,Indxcol)
%    One-step Gibbs sampler for restricted VARs in the structural form (including over-identified cases).
%    Reference: "A Gibbs sampler for structural VARs" by D.F. Waggoner and T. Zha,  ``
%               Journal of Economic Dynamics & Control (JEDC) 28 (2003) 349-366.
%    See Note Forecast (2) pp. 44-51, 70-71, and Theorem 1 and Section 3.1 in the WZ JEDC paper.
%
% A0gbs:  the last draw of A0 matrix
% UT: cell(nvar,1) -- U_i*T_i in the proof of Theorem 1 where
%            (1) a_i = U_i*b_i with b_i being a vector of free parameters
%            (2) T_i (q_i-by-q_i) is from T_i*T_i'= S_i = inv(H0inv{i}/T) where H0inv is the inverse of
%                the covariance martrix NOT divided by fss and S_i is defined in (14) on p.355 of the WZ JEDC paper.
% nvar:  rank of A0 or # of variables
% fss:  effective sample size == nSample (T)-lags+# of dummy observations
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
% Indxcol: a row vector indicating random columns this Gibbs draws.
%           When this input is not supplied, the Gibbs draws all columns
%------------------
% A0bgs: nvar-by-nvar.  New draw of A0 matrix in this Gibbs step
% Wcell: cell(nvar,1).  In each cell, columns of Wcell{i} form an orthonormal basis w_1,...,w_qi in Section 3.1 in the WZ paper.  Added 9/04.
%
% Written by Tao Zha, August 2000. Revised, September 2004.
% Copyright (c) by Waggoner and Zha
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


if (nargin==5), Indxcol=[1:nvar]; end

%---------------- Local loop for Gibbs given last A0gbs ----------
Wcell = cell(length(Indxcol),1);
w = zeros(nvar,1);
for ki=Indxcol     % given last A0gbs and generate new A0bgs
   X = A0gbs;    % WZ's Section 4.3
   X(:,ki) = 0;   % want to find non-zero sw s.t., X'*w=0


   %*** Solving for w and getting an orthonormal basis.  See Theorem 1 and also p.48 in Forecast II
   [jL,Ux] = lu(X');
   jIx0 = min(find(abs(diag(Ux))<eps)); % if isempty(jIx0), then something is wrong here
   w(jIx0+1:end) = 0;  % if jIx0+1>end, no effect on w0
   w(jIx0) = 1;
   jA = Ux(1:jIx0-1,1:jIx0-1);
   jb = Ux(1:jIx0-1,jIx0);
   jy = -jA\jb;
   w(1:jIx0-1) = jy;
      % Note: if jIx0=1 (which almost never happens for numerical stability reasons), no effect on w.

   %*** Constructing orthonormal basis w_1, ..., w_qi at each Gibbs step
   w0 = UT{ki}'*w;
   w1 = w0/sqrt(sum(w0.^2));
   [W,jnk] = qr(w1);   % Columns of W form an orthonormal basis w_1,...,w_qi in Section 3.1 in the WZ paper

   %*** Draw beta's according to Theorem 1 in the WZ JEDC paper.
   gkbeta = zeros(n0(ki),1);  % qi-by-1: greak beta's
   jstd = sqrt(1/fss);
   gkbeta(2:end) = jstd*randn(n0(ki)-1,1);  % for beta_2, ..., beta_qi
   %--- Unnormalized (i.e., not normalized) gamma or 1-d Wishart draw of beta_1 in Theorem 1 of the WZ JEDC paper.
   jr = jstd*randn(fss+1,1);
   if rand(1)<0.5
      gkbeta(1) = sqrt(jr'*jr);
   else
      gkbeta(1) = -sqrt(jr'*jr);
   end

   %*** Getting a new ki_th column in A0
   A0gbs(:,ki) = UT{ki}*(W*gkbeta);   % n-by-1: U_i*T_i*W*beta;
   Wcell{ki} = W;  %q_i-by-1.
end
