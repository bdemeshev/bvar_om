function [Ui,Vi,n0,np,ixmC0Pres] = ftd_upperchol3v(lags,nvar,nexo,indxC0Pres)
%vlist = [20 6 3 44 1 10];    % regarding "xdd", Pcom (Poil or imfcom), M2, FFR, GDP, CPI (or PCE), and U.
%
%    Exporting orthonormal matrices for the deterministic linear restrictions (equation by equation)
%    See Waggoner and Zha's Gibbs sampling paper.
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
Qi(1:2,:,1) = [
    0 1 0
    0 0 1
        ];

%======== The second equation ===========
Qi(1:1,:,2) = [
    0 0 1
        ];


%======== The third equation ===========



%===== Lagged restrictions in foreign (Granger causing) block
%nfbres = lags*(nvar-nfvar);  % number of block restrictions in each foreign equation
%bfor = zeros(nfbres,k);  % each foreign equation
%cnt=0;
%for ki = 1:lags
%   for kj=1:nvar-nfvar
%      cnt=cnt+1;
%      bfor(cnt,nvar*(ki-1)+nfvar+kj) = 1;
%   end
%end
%%
%if cnt~=nfbres
%   error('Check lagged restrictions in foreign equations!')
%end
%%
%for kj=1:nfvar
%   Ri(1:nfbres,:,kj) = bfor;
%end


%===== Lagged restrictions in home (affected) block
%
%~~~~~ selected domestic equations
%dlrindx = nfvar+1; %[nfvar+1 nfvar+2];  % index for relevant home equations
%rfvindx = []; %[6];  %[1 2 3 5];  % index for restricted foreign variables (e.g., Poil, M2, FFR, P).
%%nf2hvar = nfvar-length(rfvindx);  % number of free parameters -- foreign variables entering the home sector
%nhbres = lags*length(rfvindx);  % number of block restrictions in each home equation
%bhom = zeros(nhbres,k);  % each home equation
%cnt=0;
%for ki = 1:lags
%   for kj=1:length(rfvindx)
%      cnt=cnt+1;
%      bhom(cnt,nvar*(ki-1)+rfvindx(kj)) = 1;
%   end
%end
%%
%if cnt~=nhbres
%   error('Check lagged restrictions in domestic equations!')
%end
%%
%for kj=dlrindx
%   Ri(1:nhbres,:,kj) = bhom;
%end


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

