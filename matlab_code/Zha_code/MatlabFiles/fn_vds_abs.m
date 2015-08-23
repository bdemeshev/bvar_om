function vds = fn_vds_abs(irfs);
%vds = fn_vds(irfs)
%
% Inputs:
%  irfs: nsteps-by-n-by-r.  For each shock of r shocks, impulse responses are nsteps-by-n.
%---
% Outputs:
%  vds: nsteps-by-n-by-r.  For each shock of r shocks, variance decompositions (%) are nsteps-by-n.
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

[nsteps, n,r] = size(irfs);

vds = zeros(nsteps,n,r);

vds_square_cum = cumsum(abs(irfs),1);
vds_square_sum = repmat(sum(vds_square_cum,3),[1 1 r]);
vds = (vds_square_cum ./ vds_square_sum) .* 100;

