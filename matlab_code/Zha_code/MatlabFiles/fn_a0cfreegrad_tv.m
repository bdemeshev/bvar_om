function [g,badg] = fn_a0cfreegrad_tv(b,nvar,nStates,n0cumsum,Ui,Tkave,Del00invcell_ave,dpDelp0cell_ave)
% [g,badg] = fn_a0cfreegrad_tv(b,nvar,nStates,n0cumsum,Ui,Tkave,Del00invcell_ave,dpDelp0cell_ave)
%   Analytical gradient for fn_a0cfreegrad_tv.m when using csminwel.m.
%   Note: (1) columns correspond to equations; (2) c stands for constant.
%   See TBVAR NOTE pp. 61-62, 34a-34c.
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
%           averaged (ave) over E-step draws.  For T_k.  See p.61.
% Del00invcell_ave:  Quardratic term for b0_j in each cell.  See p.61.
% dpDelp0cell_ave:  Cross d+_j and b0_j term in each cell.  See p.61.
%----------------
% g: n0cumsum(end)-by-1 analytical gradient for fn_a0cfreegrad_tv.m.
% badg: 0, the value that is used in csminwel.m.
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
B_shat=zeros(nvar,nvar,nStates);   % inv(A0_shat);

badg = 0;
g = zeros(size(b(:)));

%**** The derivative of the exponential term w.r.t. each free constanta A0 paramater. See pp. 34a and 62.
for kj = 1:nvar
   indxn0j = [n0cumsum(kj)+1:n0cumsum(kj+1)];  % Index for the parameters in the ith equation.
   bj = b(indxn0j);
   g(indxn0j) = Del00invcell_ave{kj}*bj - dpDelp0cell_ave{kj}';
   A0_shat(:,kj,:) = reshape(Ui{kj}*bj,nvar,nStates);
end


%for si=1:nStates     % See p.34a.

%end
%**** Add the derivative of -Tk_ave*log|A0(k)| w.r.t. each free constanta A0 paramater.  See pp. 62 and 34a.
for kj = 1:nvar
   indxn0j = [n0cumsum(kj)+1:n0cumsum(kj+1)];  % Index for the parameters in the ith equation.
   for si=1:nStates     % See p.34a.
      B_shat(:,:,si)=inv(A0_shat(:,:,si));   % See p.62.
      g(indxn0j) = g(indxn0j) - Tkave(si)*( B_shat(kj,:,si)*Ui{kj}((si-1)*nvar+1:si*nvar,:) )';  % See p.62.
   end
end

