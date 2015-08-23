function [g,badg] = fn_a0sfreegrad(b,Uistar,Uibar,nvar,nStates,n0,Tkave,Sgm0tldinvave)
% [g,badg] = fn_a0sfreegrad(b,Uistar,Uibar,nvar,nStates,n0,Tkave,Sgm0tldinvave)
%   Analytical gradient for fn_a0sfreefun.m when using csminwel.m.
%    The case of no asymmetric prior and no lag restrictions.
%    Note: (1) columns correspond to equations; s stands for state.
%    See TBVAR NOTE p.34a.
%
% b: sum(n0)*nStates-by-1 vector of free A0 parameters, vectorized from the sum(n0)-by-nStates matrix.
% Uistar: cell(nvar,1).  In each cell, nvar*nStates-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.  See p.33.
% Uibar: cell(nvar,1).  In each cell, we have nvar*nStates-by-qi*nStates, rearranged
%           from Uistar or Ui.  See p.33.
% nvar:  Number of endogeous variables.
% nStates:  NUmber of states.
% n0: nvar-by-1.  n0(i)=qi where ith element represents the number of free A0 parameters in ith equation for each state.
% Tkave: nStates-by-1 of sample sizes (excluding lags but including dummies if provided) for different states k,
%           averaged (ave) over E-step draws.  For T_k.  See p.38.
% Sgm0tldinvave:  nvar*nStates-by-nvar*nStates.  The matrix inv(\Sigma~_0) averaged (ave)
%         over E-step draws. Same for all equations because of no asymmetric prior and no lag
%         restrictions.  Resembles old SpH in the exponent term in posterior of A0,
%         but NOT divided by fss (T) yet.  See p.38.
%----------------
% g: sum(n0)*nStates-by-1 analytical gradient for fn_a0sfreefun.m.
%       Vectorized frmo the sum(n0)-by-nStates matrix.
% badg: 0, the value that is used in csminwel.m.
%
% Tao Zha, March 2001
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



n0=n0(:);
A0=zeros(nvar,nvar,nStates);
n0cum = [0;cumsum(n0)];
b=reshape(b,n0cum(end),nStates);

g = zeros(n0cum(end),nStates);   % which will be vectorized later.
badg = 0;


%*** The derivative of the exponential term w.r.t. each free paramater
for kj = 1:nvar
   lenbjs=length(n0cum(kj)+1:n0cum(kj+1));
   bj = zeros(nStates*lenbjs,1);
   for si=1:nStates
      bj((si-1)*lenbjs+1:si*lenbjs) = b(n0cum(kj)+1:n0cum(kj+1),si);  % bj(si).  See p.34a.
      A0(:,kj,si) = Uistar{kj}((si-1)*nvar+1:si*nvar,:)*b(n0cum(kj)+1:n0cum(kj+1),si);
   end
   gj = (Uibar{kj}'*Sgm0tldinvave*Uibar{kj})*bj;
   for si=1:nStates
      g(n0cum(kj)+1:n0cum(kj+1),si) = gj((si-1)*lenbjs+1:si*lenbjs);
   end
end

%*** Add the derivative of -T_klog|A0(k)| w.r.t. each free paramater
for si=1:nStates     % See p.34a.
   B=inv(A0(:,:,si)');
   for ki = 1:n0cum(end)   % from 1 to sum(q_i)
      n = max(find( (ki-n0cum)>0 ));  % note, 1<=n<=nvar equations.
      g(ki,si) = g(ki,si) - Tkave(si)*B(:,n)'*Uistar{n}((si-1)*nvar+1:si*nvar,ki-n0cum(n));  % See p.34a.
   end
end
g = g(:);  % vectorized the same way as b is vectorized.
