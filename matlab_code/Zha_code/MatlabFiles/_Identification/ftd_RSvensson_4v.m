function [Ui,Vi,n0,np,ixmC0Pres] = ftd_reac_function_4v(lags,nvar,nexo,indxC0Pres)
%  vlist = [ff+ch fh dpgdp ffr)
%
%    Exporting orthonormal matrices for the deterministic linear restrictions (equation by equation)
%    See Waggoner and Zha's Gibbs sampling paper.
%
%   HERE FIRTS 3 EQUATIONS ARE AR2 AND THE LAST EQUATION IS AN UNRESTRICTED
%   REACTION FUNCTION 2 lags
%
% q_m:  quarter or month
% lags: the maximum length of lag
% nvar:  number of endogeous variables
% nexo:  number of exogenous variables.  If nexo is not supplied, nexo=1 as default for a constant
% indxC0Pres: index for cross-A0-A+ restrictions.  if 1: cross-A0-and-A+ restrictions; 0: idfile is all we have
%                Example for indxOres==1: restrictions of the form P(t) = P(t-1).
%                These restrictions have to be manually and carefully keyed in.
%-----------------
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
% Vi: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
% n0: nvar-by-1, ith element represents the number of free A0 parameters in ith equation
% np: nvar-by-1, ith element represents the number of free A+ parameters in ith equation
% ixmC0Pres:  neq_cres-by-1 cell.  Effective only if indxC0Pres=1, otherwise equals NaN.
%             neq_cres is the number of equations in which cross-A0-A+ restrictions occur.
%             In the jth cell representing equation, we have 4 columns:
%               1st: the jth column (equation) of A+ or A0: f_j or a_j
%               2nd: the ith element f_j(i) -- the ith element in the jth column of A+
%               3rd: the hth element a_j(h) -- the hth element in the jth column of A0
%               4th: the number s such that f_j(i) = s * a_j(h) holds.
%
% Tao Zha, May 2000
% BN

nvar=4;
lags=4;
nexo=1;

Ui = cell(nvar,1);  % initializing for contemporaneous endogenous variables
Vi = cell(nvar,1);  % initializing for lagged and exogenous variables
n0 = zeros(nvar,1); % ith element represents the number of free A0 parameters in ith equation
np = zeros(nvar,1); % ith element represents the number of free A+ parameters in ith equation

if (nargin==2)
   nexo = 1;  % 1: constant as default where nexo must be a nonnegative integer
elseif (nargin==3)
   indxC0Pres = 0;  % default is no cross-A0-and-A+ restrictions.
end

k = lags*nvar+nexo;  % maximum number of lagged and exogenous variables in each equation

Qi = zeros(nvar,nvar,nvar);   % for nvar contemporaneous equations
Ri = zeros(k,k,nvar);    % for nvar lagged and exogenous equations
  % Row corresponds to equation. 0 means no restriction.
  %                              1 means exclusion restriction such that the corresponding parameter is restricted to 0.

%nfvar = 6;   % number of foreign (Granger causing) variables
%nhvar = nvar-nfvar;  % number of home (affected) variables.


%-------------------------------------------------------------
%  Beginning the manual input of the restrictions one quation at a time
%-------------------------------------------------------------
%
%======== The first equation ===========
Qi(1:3,:,1) = [
    0 1 0 0
    0 0 1 0
    0 0 0 1
        ];

%======== The second equation ===========
Qi(1:2,:,2) = [
    0 0 1 0
    0 0 0 1
        ];

%======== The third equation =========== NOTE THAT WE FORBID A
%CONTEMPORANEOUS IMPACT OF OUTPUTON PRICES TO AVOID A CONSTRAINT THAT
%INVOLVE A0 and Aplus
Qi(1:3,:,3) = [
    1 0 0 0
    0 1 0 0
    0 0 0 1
        ];

%======== The fourth equation ===========


% Restrictions on the A+ in order to focus strictly on the reaction fucntion

% indicates free parameterers X i
%	Ap = [
%      X  X    X  X
%	   X  X    X  X
%     -a1 -b1  X  X
%      a1 b1   0  X  (1st lag)
%      X  X    X  X
%	   X  X    X  X
%     -a2 -b2  X  X
%      b2  b2  0  X  (2nd lag)
%      X   0   X  X
%	   X  X    X  X
%     -a3 -b3  X  X
%      a3  a3  0  X  (3rd lag)
%      X  X    X  X
%	   X  X    X  X
%     -a4 -b4  X  X
%      a4  b4  0  X  (4th lag)
%      X  X    X  X  (constant terms)
%			  ];

k=nvar*lags+nexo;
Ri = zeros(k,k,nvar);
% constraints on IS curve /conso+corporate investment
for nv=1:2
for ll=1:lags
Ri(ll,3+lags*(ll-1),nv)=1;
Ri(ll,4+lags*(ll-1),nv)=1;
end
end

% constraints on IS curve /conso+corporate investment only on the long run
% impact
% for nv=1:2
% for ll=1:lags
% Ri(1,3+lags*(ll-1),nv)=1;
% Ri(1,4+lags*(ll-1),nv)=1;
% end
% end


% constraints on Ph curve / inflation does not react to interest rates
for ll=1:lags
Ri(ll,4+lags*(ll-1),3)=1;
end


for n=1:nvar   %  initializing loop for each equation
   Ui{n} = null(Qi(:,:,n));
   Vi{n} = null(Ri(:,:,n));
   n0(n) = size(Ui{n},2);
   np(n) = size(Vi{n},2);
end



%(2)-------------------------------------------------------------
%  Cross-A0-and-A+ rerestrictions one quation at a time
%    i.e., the first, second, ..., kjth, ..., equation
%(2)-------------------------------------------------------------
%
if indxC0Pres
   neq_cres = 3;   % the number of equations in which cross-A0-A+ restrictions occur.
   ixmC0Pres = cell(neq_cres,1);  % in each cell representing equation, we have 4 columns:
           % 1st: the jth column (equation) of A+ or A0: f_j or a_j
           % 2nd: the ith element f_j(i) -- the ith element in the jth column of A+
           % 3rd: the hth element a_j(h) -- the hth element in the jth column of A0
           % 4th: the number s such that f_j(i) = s * a_j(h) holds.
   %** 1st equation
   ixmC0Pres{1} = [1 2 2 1
                   1 7 1 1];
   %** 2nd equation
   ixmC0Pres{2} = [2 2 2 2];
   %** 3rd equation
   ixmC0Pres{3} = [3 7 1 1
                   3 2 2 1];


%         % 4 columns.
%   ncres = 5;  % manually key in the number of cross-A0-A+ restrictions

%           % 1st: the jth column (equation) of A+ or A0: f_j or a_j
%           % 2nd: the ith element f_j(i) -- the ith element in the jth column of A+
%           % 3rd: the hth element a_j(h) -- the hth element in the jth column of A0
%           % 4th: the number s such that f_j(i) = s * a_j(h) holds.
else
   ixmC0Pres = NaN;
end

