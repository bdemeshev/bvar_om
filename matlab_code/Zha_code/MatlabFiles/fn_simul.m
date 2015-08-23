function simul = fn_simul(G1, impact, nsim, shocks, x0);
% Inputs:
%  G1: n-by-n;
%  impact: n-by-r;
%  nsim: length of time series for simulation.
%  shocks: r-by-nsim, exogenous driving processes
%  x0: n-by-1, initial values for the variables of interest
%---
% Outputs:
%  simul: nsim-by-n.  
%
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

[n,r] = size(impact);

[n1,n2] = size(G1);
if (n1 ~= n2) || (n1 ~= n)
   error('fn_simul.m: make sure that (1) G1 is square and (2) size(G1,1) = size(impact,1)');
end

simul_tmp = zeros(n, nsim);

simul_tmp(:, 1) = x0;
for ti = 2:nsim;
    simul_tmp(:, ti) = G1*simul_tmp(:, ti-1) + impact*shocks(:, ti);
end;
simul = simul_tmp';
