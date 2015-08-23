function fn_fprintvector(fid, vec, ncols, indxFloat)
% Prints a row vector to an ascii file indexed by fid without any spacing.
%
% Inputs:
%   fid:  Ascii file id.  Example:  fid = fopen('outdatainp_3s_stv_tvms6lags.prn','a');
%   vec:  A row vector to be written to the file.
%   ncols:  Number of columns of vec.
%   indxFloat:  1 if double;
%               2 if single;
%               0 if integer.
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

%
if size(vec,1)~=1
   error('fn_fprintvector(): The vector must be a row vector');
end
if ncols~=size(vec,2)
   error('fn_fprintvector(): The column number supplied match that of the vector');
end
for kj=1:ncols
   if (indxFloat == 1)
      fprintf(fid,' %.16e ',vec(kj));
   elseif (indxFloat == 2)
      fprintf(fid,' %.8e ',vec(kj));
   else
      fprintf(fid,' %d ',vec(kj));
   end
   if (kj==ncols)
      fprintf(fid,'\n');
   end
end
