function of = fn_a0sfreefun2(b,Uistar,Uibar,nvar,nStates,n0,n0cumsum,tvstate,tvtot,constot,...
                   tvinx,constinx,indxTV,indxConst,Tkave,Sgm0tldinvaveConst,Sgm0tldinvaveTV)
%  Negative logPosterior function for regime-switching vectorized free A0 parameters,
%    which are b's in the WZ notation. The case of no asymmetric prior and no lag restrictions.
%    It improves fn_a0sfreefun2.m in several aspects:
%       (a) allows some equations to have constant parameters.
%    It differs from fn_a0cfreefun_tv.m in several aspects:
%       (a) does not deal with equations with only structural variances time varying;
%       (b) cannot deal with lag restrictions;
%       (c) only deal with the marginal distribution of A0_s.
%    Note: (1) columns correspond to equations; (2) s stands for state.
%    See Time-Varying BVAR NOTE pp.34-34c,40.
%
% b: sum(n0)*nStates-by-1 vector of free A0 parameters, vectorized from the sum(n0)-by-nStates matrix.
% Uistar: cell(nvar,1).  In each cell, nvar*nSates-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.  See p.33.
% Uibar: cell(nvar,1).  In each cell, we have nvar*nStates-by-qi*nStates, rearranged
%           from Uistar or Ui.  See p.33.
% nvar:  Number of endogeous variables.
% nStates:  NUmber of states.
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation for each state.
% n0cumsum: Equal to [0;cumsum(n0)];
% tvstate:  Number of time-varying parameters for all equations for each state.
% tvtot: Total number of time-varying parameters.
% constot: Total number of constant parameters.
% tvinx: Vectorized index to include time-varying parameters for all equations under *each* state.
% constinx:  Vectorized index to include constant parameters for all equations.
% indxTV:  Index of acending order for equations with time-varying parameters.
% indxConst: Index of ascending order for equations with constant parameters.  When [], all equations
%            are time varying; when [1:nvar], all equations have constant parameters.
% Tkave: nStates-by-1 of sample sizes (excluding lags but including dummies if provided) for different states k,
%           averaged (ave) over E-step draws.  For T_k.  See p.40.
% Sgm0tldinvaveConst
% Sgm0tldinvaveTV
% Sgm0tldinvave:  nvar*nStates-by-nvar*nStates.  The matrix inv(\Sigma~_0) averaged (ave)
%         over E-step draws. Same for all equations because of no asymmetric prior and no lag
%         restrictions.  Resembles old SpH in the exponent term in posterior of A0,
%         but NOT divided by fss (T) yet.  See p.40.
%----------------
% of:  Objective function (negative logPosterior).
%
% This function is called by tveml*.m which is an old program compared with szeml*.m.  Thus,
%   use fn_a0cfreefun_tv.m if possible.
% Tao Zha, August 2001.
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
b_shat = zeros(n0cumsum(end),nStates);
b_shat(tvinx,:) = reshape(b(1:tvtot),tvstate,nStates);  % Time-varying parameter matrix.
b_shat(constinx,:) = repmat(b(tvtot+1:end), [1 nStates]);  % Constant parameter matrix.

tra = 0.0;
for kj = 1:nvar
   lenbjs=length(n0cumsum(kj)+1:n0cumsum(kj+1));
   bj = zeros(nStates*lenbjs,1);
   a0j = zeros(nStates*nvar,1);  % See p.33.
   for si=1:nStates
      bj((si-1)*lenbjs+1:si*lenbjs) = b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);  % bj(si).  See p.34a.
      A0(:,kj,si) = Uistar{kj}((si-1)*nvar+1:si*nvar,:)*b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);
      a0j((si-1)*nvar+1:si*nvar) = A0(:,kj,si);   % a0j(si).  See p.33.
   end
   %tra = tra + 0.5*a0j'*Sgm0tldinvave*a0j;   % negative exponential term
   if find(kj==indxTV)   % For time-varying equations.
      tra = tra+0.5*bj'*(Uibar{kj}'*Sgm0tldinvaveTV*Uibar{kj})*bj;
   else     % % For constant parameter equations.
      tra = tra+0.5*bj'*(Uibar{kj}'*Sgm0tldinvaveConst*Uibar{kj})*bj;
   end
end

ada=0.0;
for si=1:nStates     % See p.34a.
   [A0l,A0u] = lu(A0(:,:,si));
   ada = ada - Tkave(si)*sum(log(abs(diag(A0u))));    % negative log determinant of A0 raised to power T
end

of = ada + tra;
