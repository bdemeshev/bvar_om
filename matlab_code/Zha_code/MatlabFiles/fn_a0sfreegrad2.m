function [g,badg] = fn_a0sfreegrad2(b,Uistar,Uibar,nvar,nStates,n0,n0cumsum,tvstate,tvtot,constot,...
                     tvinx,constinx,indxTV,indxConst,Tkave,Sgm0tldinvaveConst,Sgm0tldinvaveTV)
%   Analytical gradient for both time-varying and constant parameters for fn_a0sfreefun2.m when using csminwel.m.
%    The case of no asymmetric prior and no lag restrictions.
%    Note: (1) columns correspond to equations; s stands for state.
%    See Time-Varying BVAR NOTE pp.34a-34c.
%
% b: tvtot+constot-by-1 vector of free A0 parameters, vectorized from b_shat for both time-varying and constant parameter matrices.
% Uistar: cell(nvar,1).  In each cell, nvar*nStates-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters
%           where in each state, nvar-by-qi is the same. With this transformation, we have
%           ai = Ui*bi or Ui'*ai = bi where ai is a vector of total original parameters and
%           bi is a vector of free parameters. When no restrictions are imposed, we have Ui = I.
%           There must be at least one free parameter left for the ith equation.  See p.33.
% Uibar: cell(nvar,1).  In each cell, we have nvar*nStates-by-qi*nStates, rearranged
%           from Uistar or Ui.  See p.33.
% nvar:  Number of endogenous variables.
% nStates:  NUmber of states.
% n0: nvar-by-1.  n0(i)=qi where ith element represents the number of free A0 parameters in ith equation for each state.
% n0cumsum=[0;cumsum(n0)]: nvar+1-by-1.  Cumulative sums of n0 with zero added.
% tvstate=length(tvinx): Number of time-varying parameters for all equations under *each* state.
% constot=length(constinx): Total number of constant parameters for all equations.
% tvtot = tvstate*nStates: Total number of time-varying parameters for all equations under all states.
% tvinx: Vectorized index to include time-varying parameters for all equations under *each* state.
% constinx: Vectorized index to include constant parameters for all equations.
% indxTV: Index vector of ascending order for equations with time-varying parameters.
% indxConst: Index vector of ascending order for equations with constant parameters.  When [], all equations
%       are time varying; when [1:nvar], all equations have constant parameters.
% Tkave: nStates-by-1 of sample sizes (excluding lags but including dummies if provided) for different states k,
%           averaged (ave) over E-step draws.  For T_k.  See p.40.
% Sgm0tldinvaveConst
% Sgm0tldinvaveTV
% Sgm0tldinvave:  nvar*nStates-by-nvar*nStates.  The matrix inv(\Sigma~_0) averaged (ave)
%         over E-step draws. Same for all equations because of no asymmetric prior and no lag
%         restrictions.  Resembles old SpH in the exponent term in posterior of A0,
%         but NOT divided by fss (T) yet.  See p.40.
%----------------
% g: tvtot+constot-by-1 analytical gradient for fn_a0sfreefun2.m.  Vectorized first from the
%     time-varying parameter matrix and then from constant parameter and constant parameter matrix.
% badg: Always 0 --- the value that is used in csminwel.m.
%
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
B0=A0;
b_shat = zeros(n0cumsum(end),nStates);
b_shat(tvinx,:) = reshape(b(1:tvtot),tvstate,nStates);  % Time-varying parameter matrix.
b_shat(constinx,:) = repmat(b(tvtot+1:end), [1 nStates]);  % Constant parameter matrix.
%
gmix = zeros(n0cumsum(end),nStates);
      %  Gradient matrix mixing time-varying and constant parameters.
      %  Rows of gmix, indexed by constinx, will not be used but instead treated
      %    different under "for kj=indxConst."
gtv = zeros(tvstate,nStates);   % Gradient wrt time-varying parameters, which will be vectorized later.
gconst1 = zeros(constot,1);    % The first gradient wrt constant parameters for the exponetial term.
gconst2=zeros(constot,1);   % The second gradient wrt constant parameters for -T_klog|A0(k)|.


%---------------------------------------------------
%  Time-varying situation.  See p. 34a.
%---------------------------------------------------
%*** The derivative of the exponential term w.r.t. each free time-varying paramater.
for kj = 1:nvar    % For all equations
   lenbjs=length(n0cumsum(kj)+1:n0cumsum(kj+1));
   bj = zeros(nStates*lenbjs,1);
   for si=1:nStates
      bj((si-1)*lenbjs+1:si*lenbjs) = b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);  % bj(si).  See p.34a.
      A0(:,kj,si) = Uistar{kj}((si-1)*nvar+1:si*nvar,:)*b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);
   end
   gj = (Uibar{kj}'*Sgm0tldinvaveTV*Uibar{kj})*bj;
   if find(kj==indxTV)   % For time-varying equations.
      for si=1:nStates
         gmix(n0cumsum(kj)+1:n0cumsum(kj+1),si) = gj((si-1)*lenbjs+1:si*lenbjs);
      end
   end
end

%*** Add the derivative of -T_klog|A0(k)| w.r.t. each free time-varying paramater.
for si=1:nStates     % See p.34a.
   B0(:,:,si)=inv(A0(:,:,si)');
   for ki = tvinx   % Those from 1 to sum(q_i) that pertain to time-varying parameters.
      n = max(find( (ki-n0cumsum)>0 ));  % note, 1<=n<=nvar equations.
      gmix(ki,si) = gmix(ki,si) - Tkave(si)*B0(:,n,si)'*Uistar{n}((si-1)*nvar+1:si*nvar,ki-n0cumsum(n));  % See p.34a.
   end
end

%---------------------------------------------------
%  Constant-parameter situation.  See pp. 34b-34c.
%---------------------------------------------------
kcount = 0;   % Counts.
for kj = indxConst    % Equations
   klocat=kcount+n0(kj);   % Losition for the last constant element in the jth equation.
   Bbar = repmat(eye(n0(kj)),[1 nStates]);   % See p.34b.
   lenbjs=length(n0cumsum(kj)+1:n0cumsum(kj+1));
   bj = zeros(nStates*lenbjs,1);
   for si=1:nStates
      bj((si-1)*lenbjs+1:si*lenbjs) = b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);  % bj(si).  See p.34a.
      A0(:,kj,si) = Uistar{kj}((si-1)*nvar+1:si*nvar,:)*b_shat(n0cumsum(kj)+1:n0cumsum(kj+1),si);
   end
   gconst1(kcount+1:klocat) = Bbar*( (Uibar{kj}'*Sgm0tldinvaveConst*Uibar{kj})*bj );
   kcount=klocat;
end
%*** Add the derivative of -T_klog|A0(k)| w.r.t. each free time-varying paramater.
kcount = 0;   % Counts.
for ki = constinx   % Those from 1 to sum(q_i) that pertain to constant parameters.
   kcount=kcount+1;
   n = max(find( (ki-n0cumsum)>0 ));  % note, 1<=n<=nvar equations.
   gi=0;   % Initializes a gradient summed across the states.  See p. 34c.
   for si=1:nStates     % See p.34c.
      gi = gi - Tkave(si)*B0(:,n,si)'*Uistar{n}((si-1)*nvar+1:si*nvar,ki-n0cumsum(n));  % See p. 34a.
   end
   gconst2(kcount)=gi;
end


%---------------------------------------------------
%  Conclusion.
%---------------------------------------------------
gtv = gmix(tvinx,:);   % Extract the gradient wrt time-varying parameters.
gconst = gconst1+gconst2;  % Gradient wrt to constant parameters.
%
g = [gtv(:);gconst];  % Vectorized the same way as b is vectorized.
badg = 0;    % This term is used for Sims's csminwel program.
