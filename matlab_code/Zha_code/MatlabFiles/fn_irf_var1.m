function irfs = fn_irf_var1(G1, impact, nsteps);
%irfs = fn_irf_var1(G1, impact,nsteps);
% Inputs:
%  G1: n-by-n;
%  impact: n-by-r;
%  nsteps: number of steps for impulse responses.
%---
% Outputs:
%  irfs: nsteps-by-n-by-r.  For each shock of r shocks, impulse responses are nsteps-by-n.
%
% See fn_vds.m and fn_uncondfcst_var1.m.
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
   error('fn_irf_var1.m: make sure that (1) G1 is square and (2) size(G1,1) = size(impact,1)');
end

irfs = zeros(nsteps,n,r);


%---- Impuse responses at the first step.
M = impact;
for ri=1:r
   irfs(1,:,ri) = M(:,ri)';
end

for ti = 2:nsteps
   M = G1*M;
   for ri=1:r
      irfs(ti,:,ri) = M(:,ri)';
   end
end
