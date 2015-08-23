function of = fn_a0cfreefun_tv(b,nvar,nStates,n0cumsum,Ui,Tkave,Del00invcell_ave,dpDelp0cell_ave)
% of = fn_a0cfreefun_tv(b,nvar,nStates,n0cumsum,Ui,Tkave,Del00invcell_ave,dpDelp0cell_ave)
%   Negative logPosterior function for squeesed free A0 parameters (bar tv structural variances) in
%     a new SZ time-varying model, which are vectorized as b's in the WZ notation.  In other words,
%     A0 parameters can be fully time varying, but the time varying structural variances will be
%     excluded.  It differs from fn_a0sfreefun*.m in several aspects:
%        (a) can deal with equations with only structural variances time varying;
%        (b) take care of lag restrictions;
%        (c) conditional on the values of all other parameters including A+.
% Note: (1) columns correspond to equations; (2) c stands for constant or for an exclusion of
%           time varying structural variances.
% See TBVAR NOTE p. 61.
%
% b: n0cumsum(end)-by-1 vector of free constant A0 parameters, vectorized from b_ihatcell.
% nvar:  Number of endogeous variables.
% nStates:  Number of states.
% n0cumsum: [0;cumsum(n0)] where n0 is nvar-by-1 and its ith element represents the number of
%           free constant A0 parameters in ith equation for *all states*.
% Ui: nvar-by-1 cell.  In each cell, nvar*nStates-by-(qi+si) orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters
%           within the state and si is the number of free parameters across the states.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
% Tkave: nStates-by-1 of sample sizes (excluding lags but including dummies if provided) for different states k,
%           averaged (ave) over E-step draws.  For T_k.  See p.61.  In the Gibbs-Metropolis step, there is
%           no need to make Tk the same as Tkave (See Chib and Jeliazkov 2001).
% Del00invcell_ave:  Quardratic term for b0_j in each cell.  See p.61. When fn_a0cfreegrad_tv.m is used,
%           always make this systemetric by using (Del00invcell_ave'+Del00invcell_ave)/2.
%           Also note that in the Gibbs-Metropolis step, there is no need to make Del00invcell the same
%           as Del00invcell_ave (See Chib and Jeliazkov 2001).
% dpDelp0cell_ave:  Cross d+_j and b0_j term in each cell.  See p.61. In the Gibbs-Metropolis step, there is no need
%           to make dpDelp0cell the same as dpDelp0cell_ave (See Chib and Jeliazkov 2001).
%----------------
% of:  Objective function (negative logPosterior).
%
% This function is called by szeml*.m which is a bettern program than tveml*.m.  Thus,
%   use this function in place of fn_a0sfreefun2.m if possible.
% See fn_a0cfreegrad_tv.m for analytical gradient for this function.
%
% Tao Zha, September 2001
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



A0_shat=zeros(nvar,nvar,nStates);    %  Arranged A0 matrix across states.

tra = 0.0;
for kj = 1:nvar
   bj = b(n0cumsum(kj)+1:n0cumsum(kj+1));
   A0_shat(:,kj,:) = reshape(Ui{kj}*bj,nvar,nStates);
   tra = tra+0.5*( (bj'*Del00invcell_ave{kj})*bj - 2*dpDelp0cell_ave{kj}*bj );  % Negative exponential term
end

ada=0.0;
for si=1:nStates     % See p.34a.
   [A0l,A0u] = lu(A0_shat(:,:,si));
   ada = ada - Tkave(si)*sum(log(abs(diag(A0u))));    % Negative log determinant of A0 raised to power T
end

of = ada + tra;
