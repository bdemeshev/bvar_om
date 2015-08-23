function [Ui,Vi,n0,np,ixmC0Pres] = ftd_upperchol5v_LagRes_1st2v(lags,nvar,nexo,indxC0Pres)
%  vlist = [127 124 93 141 21];    % 1: GDP; 2: GDP deflator 124 (consumption deflator 79); 3: R; 4: M3 141 (M2 140); 5: exchange rate 21.
%  varlist={'y', 'P', 'R', 'M3', 'Ex'};
%
%    Exporting orthonormal matrices for the deterministic linear restrictions (equation by equation)
%    See Waggoner and Zha's Gibbs sampling paper.
%
% q_m:  quarter or month
% lags: the maximum length of lag
% nvar:  number of endogeous variables
% nexo:  number of exogenous variables.  If nexo is not supplied, nexo=1 as default for a constant.
%        Our prior at this point allows only nexo=1 (constant term).  
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

%-------------------------------------------------------------
%  Beginning the manual input of the restrictions one quation at a time
%-------------------------------------------------------------
%
%======== The first equation for Ramer's spending dummy ===========
Qi(1:4,:,1) = [
    0 1 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 1
        ];

%======== The second equation for Romer's tax dummy ===========
Qi(1:3,:,2) = [
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 1
        ];

%======== The third equation ===========
Qi(1:2,:,3) = [
    0 0 0 1 0
    0 0 0 0 1
        ];


%======== The fourth equation ===========
Qi(1:1,:,4) = [
    0 0 0 0 1
         ];


%======== The fifth equation ===========




%****************************************************************
%*** Block lagged restrictions in exogenous (Granger causing) block    
%****************************************************************
indx_lag_res_exoblock = 1; %1: allowing for lagged resrictions in the exogenous block; 0: not allowing for lagged restrictions.
if (indx_lag_res_exoblock)                            
   indx_order_exo = 1; %1: ordering exogenous block first; 0: ordering endogenous block first. 
   n_exo_nvar = 2;   % number of exgoenous (Granger causing) variables -- in our case, 2 dummy variables.
      
   %~~~~~~ The following is automated. ~~~~~~
   n_eng_nvar = nvar-n_exo_nvar;  % number of endogenous (affected) variables -- in our case, 3 variables: spending, tax, and gdp.
   if (indx_order_exo)
      loc_exovars = [1:n_exo_nvar];     % locations for relevant endogenous equations
   else
      loc_exovars = [n_eng_nvar+1:nvar];     % locations for relevant endogenous equations   
   end
   %
   n_exo_bres = lags*n_eng_nvar;  % number of block restrictions in each exogenous equation
   b_exo_res = zeros(n_exo_bres,k);  % each exogenous equation
   cnt=0;
   for ki = 1:lags
      for kj= 1:n_eng_nvar 
         cnt=cnt+1;     
         if (indx_order_exo) %Exgoneous block of variables ordered first (restrictions on endogenous block).
            b_exo_res(cnt,nvar*(ki-1)+n_exo_nvar+kj) = 1;
         else  %Endognoues block of variables ordered first (restrictions on endogenous block).
            b_exo_res(cnt,nvar*(ki-1)+kj) = 1;            
         end
      end
   end
   %
   if cnt~=n_exo_bres
      error('Check lagged restrictions in exogenous equations!')
   end
   %
   for kj=loc_exovars
      Ri(1:n_exo_bres,:,kj) = b_exo_res;
   end
end


%****************************************************************
%*** Lagged restrictions in endogenous (affected) block    
%****************************************************************
indx_lag_res_engblock = 0; %1: allowing for lagged resrictions in the endogenous block; 0: not allowing for lagged restrictions.
if (indx_lag_res_engblock)
   %~~~~~~ The following is automated. ~~~~~~
   if (indx_order_exo)
      loc_engvars = [n_exo_nvar+1:nvar]; %[n_exo_nvar+1 n_exo_nvar+2];  % locations for relevant endogenous equations
   else
      loc_engvars = [1:n_eng_nvar];    % locations for relevant endogenous equations
   end  
   %
   n_eng_bres = lags*n_exo_nvar;  % number of block restrictions in each endogenous equation
   b_eng_res = zeros(n_eng_bres,k);  % each endogenous equation
   cnt=0;
   for ki = 1:lags
      for kj=1:n_exo_nvar
         cnt=cnt+1;
         if (indx_order_exo) %Exgoneous block of variables ordered first (restrictions on exogenous block).
            b_eng_res(cnt,nvar*(ki-1)+kj) = 1;
         else  %Endognoues block of variables ordered first (restrictions on exogenous block).
            b_eng_res(cnt,nvar*(ki-1)+n_eng_nvar+kj) = 1;
         end            
      end
   end
   %
   if cnt~=n_eng_bres
      error('Check lagged restrictions in endogenous equations!')
   end
   %
   for kj=loc_engvars
      Ri(1:n_eng_bres,:,kj) = b_eng_res;
   end
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

