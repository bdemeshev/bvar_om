function [uniform_cell, uniform_matrix3D, names_fields] = fn_numstruct2numcell_nummatrix(uniform_ps)
%[uniform_cell, uniform_matrix3D, names_fields] = fn_numstruct2numcell_nummatrix(uniform_ps)
%   
% Inputs:
%    uniform_ps: a structure where each field has the same nrows-by-ncols matrix.  This function does NOT work
%                  if each field has different data types.
% Outputs:
%    uniform_cell:      nfields cells where each cell has a nrows-by-ncols matrix.  
%    uniform_matrix3D: 3-D array: nrows-by-ncols-by-nfields
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

names_fields = fieldnames(uniform_ps);
nfields = length(names_fields); 
uniform_cell = struct2cell(uniform_ps);
nrows = size(uniform_cell{1},1);
ncols = size(uniform_cell{1},2);

uniform_matrix3D = zeros(nrows,ncols,nfields);
for (ni=1:nfields)
   uniform_matrix3D(:,:,ni) = uniform_cell{ni};
end