function fn_printmatrix4tex(M, nrows, ncols, indxFloat)
% Prints the matrix to a screen for a tex file.
%
% Inputs:
%   M:  The matrix to be written to the file.
%   nrows:  Number of rows of M.
%   ncols:  Number of columns of M.
%   indxFloat:  1 if double;
%               2 if single;
%               3 if only 2 significant digits
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
if nrows~=size(M,1)
   error('fn_fprintmatrix(): Make sure the row number supplied match that of the matrix');
end
if ncols~=size(M,2)
   error('fn_fprintmatrix(): Make sure the column number supplied match that of the matrix');
end
for ki=1:nrows
   for kj=1:ncols
      fprintf(1,' & ');
      if (indxFloat == 1)
         fprintf(1,' %.16e ',M((kj-1)*nrows+ki));
      elseif (indxFloat == 2)
         fprintf(1,' %.8e ',M((kj-1)*nrows+ki));
      elseif (indxFloat == 3)
         fprintf(1,' %5.2f ',M((kj-1)*nrows+ki));
      else
         fprintf(1,' %d ',M((kj-1)*nrows+ki));
      end
      if (kj==ncols)
         fprintf(1,'\\\\\n');
      end
   end
   if (ki==nrows)
      fprintf(1,'\n\n');
   end
end
