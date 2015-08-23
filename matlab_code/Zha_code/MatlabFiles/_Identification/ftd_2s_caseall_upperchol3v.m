function [Ui,Vi,n0,np] = ftd_2s_caseall_upperchol3v(lags,nvar,nStates,indxEqnTv_m,nexo)
% Case 2&3: Policy a0j and a+j have only time-varying structural variances -- Case 2.
%           All policy and nonpolicy a0j's and the corresponding constant terms are completely time-varying and only the scale
%             of each variable in d+j,1** (excluding the constant term) is time-varying -- Case 3.
%
% Variables: Pcom, M2, FFR, y, P, U.  Equations: information, policy, money demand, y, P, U.
% Exporting orthonormal matrices for the deterministic linear restrictions
%     (equation by equation) with time-varying A0 and D+** equations.
% See Model II.3 on pp.71k-71r in TVBVAR NOTES (and Waggoner and Zha's Gibbs sampling paper and TVBVAR NOTES p.58).
%
% lags: Maximum length of lag.
% nvar:  Number of endogeous variables.
% nStates:  Number of states.
% indxEqnTv_m: nvar-by-2. Stores equation characteristics.
%                   1st column: labels of equations [1:nvar]'.
%                   2nd column: labels of time-varying features with
%                        1: indxConst -- all coefficients are constant,
%                        2: indxStv -- only shocks are time-varying,
%                        3: indxTva0pv -- a0 are freely time-varying and each variable i for d+ is time-varying only by the scale lambda_i(s_t).
%                        4: indxTva0ps -- a0 are freely time-varying and only the scale for the whole of d+ is time-varying.
%                        5: indxTv -- time-varying for all coeffficients (a0 and a+) where the lag length for a+ may be shorter.
% nexo:  number of exogenous variables.  If nexo is not supplied, nexo=1 as default for a constant.
%          So far this function is written to handle one exogenous variable, which is a constant term.
%-----------------
% Ui: nvar-by-1 cell.  In each cell, nvar*nStates-by-qi*si orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters
%           within the state and si is the number of free states.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation in the order of [a_i for 1st state, ..., a_i for last state].
% Vi: nvar-by-1 cell.  In each cell, k*nStates-by-ri*si orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters within the state and si is the number of free states.
%           With this transformation, we have fi = Vi*gi or Vi'*fi = gi where fi is a vector of total original
%           parameters and gi is a vector of free parameters.  The ith equation is in the order of [nvar variables
%           for 1st lag and 1st state, ..., nvar variables for last lag and 1st state, const for 1st state, nvar
%           variables for 1st lag and 2nd state, nvar variables for last lag and 2nd state, const for 2nd state, and so on].
% n0: nvar-by-1, whose ith element represents the number of free A0 parameters in ith equation in *all states*.
% np: nvar-by-1, whose ith element represents the number of free D+ parameters in ith equation in *all states*.
%
% Tao Zha, February 2003



Ui = cell(nvar,1);  % initializing for contemporaneous endogenous variables
Vi = cell(nvar,1);  % initializing for lagged and exogenous variables
n0 = zeros(nvar,1); % ith element represents the number of free A0 parameters in ith equation in all states.
np = zeros(nvar,1); % ith element represents the number of free D+ parameters in ith equation in all states.

if (nargin==3)
   nexo = 1;  % 1: constant as default where nexo must be a nonnegative integer
end


n = nvar*nStates;
kvar=lags*nvar+nexo;  % Maximum number of lagged and exogenous variables in each equation under each state.
k = kvar*nStates;  % Maximum number of lagged and exogenous variables in each equation in all states.

Qi = zeros(n,n,nvar);   % 3rd dim: nvar contemporaneous equations.
Ri = zeros(k,k,nvar);    % 1st and 2nd dims: lagged and exogenous equations.
   % Row corresponds to equation with nvar variables for state 1, ..., nvar variables for state nState.
   %        0 means no restriction.
   %        1 and -1 or any other number means the linear combination of the corresponding parameters is restricted to 0.
   %        1 (only 1) means that the corresponding parameter is restricted to 0.

%nfvar = 6;   % number of foreign (Granger causing) variables
%nhvar = nvar-nfvar;  % number of home (affected) variables.


%-------------------------------------------------------------
%  Beginning the manual input of the restrictions one quation at a time for A0_s.
%-------------------------------------------------------------
%

%======== The first equation ===========
eqninx = 1;
nreseqn = 2;  % Number of linear restrictions for A0(:,eqninx) for each state.
if (indxEqnTv_m(eqninx, 2)<=2)
   %**** For constant A0_s.    In the order of [a0j(1),...,a0j(nStates)] for the 2nd dim of Qi.
   Qi(1:(nStates-1)*nvar+nreseqn,:,eqninx) = [
      1  0  0      -1  0  0
      0  1  0       0 -1  0
      0  0  1       0  0 -1

      0 0 0       0 1 0
      0 0 0       0 0 1
                         ];
   %**** For constant D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   for si=1:nStates-1
      for ki=1:kvar
         Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
      end
   end
else    % Time-varying equations at least for A0_s.  For D+_s, constant-parameter equations in general.
   %**** For time-varying A0_s.    In the order of [a0j(1),...,a0j(nStates)] for the 2nd dim of Qi.
   Qi(1:nreseqn*nStates,:,eqninx) = [
      0 1 0       0 0 0
      0 0 1       0 0 0

      0 0 0       0 1 0
      0 0 0       0 0 1
                         ];

   %**** For D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   if (indxEqnTv_m(eqninx, 2)==3)    % For constant D+** except the constant term.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar-1   % -1: no restrictions on the constant term, which is freely time-varying.
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   elseif (indxEqnTv_m(eqninx, 2)==4)    % For constant D+**.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   else
      error('.../ftd_2s_caseall_*.m:  Have not got time to deal with the simple case indxEqnTv_m(eqninx, 2)=5.')
   end
end


%======== The second equation ===========
eqninx = 2;
nreseqn = 1;  % Number of linear restrictions for A0(:,eqninx) for each state.
if (indxEqnTv_m(eqninx, 2)<=2)
   %**** For constant A0_s.    In the order of [a0j(1),...,a0j(nStates)] for the 2nd dim of Qi.
   Qi(1:(nStates-1)*nvar+nreseqn,:,eqninx) = [
      1  0  0     -1  0  0
      0  1  0      0 -1  0
      0  0  1      0  0 -1

      0 0 0      0 0 1
                         ];
   %**** For constant D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   for si=1:nStates-1
      for ki=1:kvar
         Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
      end
   end
else    % Time-varying equations at least for A0_s.  For D+_s, constant-parameter equations in general.
   %**** For time-varying A0_s.    In the order of [a0j(1),...,a0j(nStates)] for the 2nd dim of Qi.
   Qi(1:nreseqn*nStates,:,eqninx) = [
      0 0 1       0 0 0

      0 0 0      0 0 1
                         ];

   %**** For D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   if (indxEqnTv_m(eqninx, 2)==3)    % For constant D+** except the constant term.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar-1   % -1: no restrictions on the constant term, which is freely time-varying.
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   elseif (indxEqnTv_m(eqninx, 2)==4)    % For constant D+**.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   else
      error('.../ftd_3s_case3a.m:  Have not got time to deal with the simple case indxEqnTv_m(eqninx, 2)=5.')
   end

   %==== For freely time-varying A+ for only the first 6 lags.
   %====       Lagged restrictions: zeros on all lags except the first 6 lags in the MS equation.
   %  nlagsno0 = 6;   % Number of lags to be nonzero.
   %  for si=1:nStates
   %     for ki = 1:lags-nlagsno0
   %        for kj=1:nvar
   %           Ri(kvar*(si-1)+nlagsno0*nvar+nvar*(ki-1)+kj,kvar*(si-1)+nlagsno0*nvar+nvar*(ki-1)+kj,2) = 1;
   %        end
   %     end
   %  end
   %**** For constant D+_s except the first two lags and the constant term.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   %  for si=1:nStates-1
   %     for ki=[2*nvar+1:kvar-1]
   %        Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
   %     end
   %  end
end


%======== The third equation (money demand) ===========
eqninx = 3;
nreseqn = 0;  % Number of linear restrictions for the equation for each state.
if (indxEqnTv_m(eqninx, 2)<=2)
   %**** For constant A0_s.    In the order of [a0j(1),...,a0j(nStates)] for the 2nd dim of Qi.
   Qi(1:(nStates-1)*nvar+nreseqn,:,eqninx) = [
      1  0  0      -1  0  0
      0  1  0       0 -1  0
      0  0  1       0  0 -1
                         ];
   %**** For constant D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   for si=1:nStates-1
      for ki=1:kvar
         Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
      end
   end
else    % Time-varying equations at least for A0_s.  For D+_s, constant-parameter equations in general.
   %**** For D+_s.  In the order of [aj+(1),...,aj+(nStates)] for the 2nd dim of Ri.
   if (indxEqnTv_m(eqninx, 2)==3)    % For constant D+** except the constant term.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar-1   % -1: no restrictions on the constant term, which is freely time-varying.
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   elseif (indxEqnTv_m(eqninx, 2)==4)    % For constant D+**.  In the order of [dj**(1),...,dj**(nStates)] for the 2nd dim of Ri.
      for si=1:nStates-1
         for ki=1:kvar
            Ri(kvar*(si-1)+ki,[kvar*(si-1)+ki si*kvar+ki],eqninx) = [1 -1];
         end
      end
   else
      error('.../ftd_2s_caseall_*.m:  Have not got time to deal with the simple case indxEqnTv_m(eqninx, 2)=5.')
   end
end



for ki=1:nvar   %  initializing loop for each equation
   Ui{ki} = null(Qi(:,:,ki));
   Vi{ki} = null(Ri(:,:,ki));
   n0(ki) = size(Ui{ki},2);
   np(ki) = size(Vi{ki},2);
end
