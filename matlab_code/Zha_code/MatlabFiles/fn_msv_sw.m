function [G1cell_sw, err] = fn_msv_sw(G1cell_0, A11cell, A12cell, A21cell, A22cell, Hcell, P)
%G1cell_sw = fn_msv_sw(G1cell_0, A11cell, A12cell, A21cell, A22cell, Hcell,P)
%
% Model:  X_{t+1} = A11_{t+1} X_t + A12_{t+1} x_t + C_{t+1} episilon_{t+1}
%         E_t(H_{t+1} x_{t+1}) = A21_t X_t + A22_t x_t
% Solution: x_t = G1_t X_t where G1_t is nx-by-nX.
%
%Output:
% G1cell_sw: solution.
% err: if the value > sqrt(eps), then not converged ==> no solution.
%Inputs:
% G1cell_0: initial guess (starting point) for the Svensson-William algorithm.
%           If G1cell_0 is empty, G1_0 will be randomly selected.
% A11cell, A12cell, A21cell, A22cell, and Hcell are all ns-by-1 cells -- the model coefficients depending on regime.
%   In each cell, the dimensions must be
%      A11=zeros(nX,nX); A12=zeros(nX,nx);
%      A21=zeros(nx,nX); A22=zeros(nx,nx);
%      H=diag(ones(nx,1));
% P: ns-by-ns transition matrix with each column summing up to 1.
%
% See ZhaNotes on FWZ MSV solution, p. AA.1.  For the SW solution, see (7) on p. AA.1.
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

%--- Convergence criteria
tol = 1.0e-09;
maxiter = 1e+04;

%--- Dimensions.
ns = length(A11cell);   %number of states (regimes)
nX = size(A11cell{1},1);
nx = size(A22cell{1},1);

%--- Initialization.
if isempty(G1cell_0)
   G1cell_0 = cell(ns,1);
   for si=1:ns
      G1cell_0{si} = 1000.0*randn(nx,nX);
   end
end

%====== Iterative procedure. ======
Gerr = zeros(ns,1);
err = 1.0;
iter = 1;
while ( (iter<maxiter) && (err>tol) )
   for j=1:ns
      M1 = zeros(nx,nx);
      M2 = zeros(nx,nX);
      for k=1:ns
         HGk = Hcell{k}*G1cell_0{k};
         M1 = M1 + P(k,j)*HGk*A12cell{k};
         M2 = M2 + P(k,j)*HGk*A11cell{k};
      end  % k loop
      Mden = A22cell{j} - M1;
      Mnum = M2 - A21cell{j};
      G1cell_sw{j} = Mden\Mnum;
      Gerr(j)=norm(G1cell_sw{j} - G1cell_0{j});
      G1cell_0{j} = G1cell_sw{j};  % Updated for the next iteration.
   end  %j loop
   err=sum(Gerr);
   iter=iter+1;
end % err loop

