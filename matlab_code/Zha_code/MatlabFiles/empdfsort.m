function imfcnt = empdfsort(imfcnt,imf_h,imfloor,hbin,ninv)
% imfcnt = empdfsort(imfcnt,imf_h,imfloor,hbin,ninv)
%    Variablewise (but multivariate) empirical probability distribution with counts
%       sorted into bins variablewise
%    Note that since a particular draw (imf_h) can be below "imfloor" and above "imceiling"
%       (=(imceiling-imfloor)*hbin), this function allows ninv+2 bins for each variable
%
% imfcnt:  initial ninv+2-by-k matrix that records counts in each bin given a column in imfcnt
%          if k==1, then only one variable is considered.
% imf_h:  particular draw, needed not to be a 1-by-k row vector
% imfloor: the low value of imf but imf_h can be below "imfloor" and above "imceiling"
% hbin:  bin size = (imceilling-imfloor)/ninv
% ninv:  number of bins between "imfloor" and "imceiling"
%------------------
% imfcnt:  updated ninv+2-by-k matrix of counts (probability)
%
% January 1999 TAZ
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


%** find the bin locations and arrange them to the order of 1, 2, ...
countInt = floor( (imf_h-imfloor) ./ hbin ); % imstp-by-nvar^2
       % bin locations from <0, 0, 1,..., ninv-1, >=ninv, a total of ninv+2 bins
countInt = 2+countInt(:)';  % row vector, 1-by-imstp*nvar^2, see my shock (1), pp.1
                    % move everyting by 2 so as to take account of <0

countInt(find(countInt<2)) = 1;     % set <0 or -infinity at 1 to start
countInt(find(countInt>=ninv+2)) = ninv+2;  % set >=ninv+2 or +infinity at ninv+2 to end
countInt = countInt + (0:length(countInt)-1)*(ninv+2);  % index to fill the matrix with prob. (counts)
             % The term after "+": with every count, skip ninv+2 to keep
             %    each column in "imfcnt" with only one element (which is probability)
imfcnt(countInt) = imfcnt(countInt) + 1;
