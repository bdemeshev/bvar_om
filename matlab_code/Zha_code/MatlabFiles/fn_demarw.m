function [imfl,imfh,imfl1,imfh1, imfmedian] = demarw(xpo,xprob)
%  [imfl,imfh,imfl1,imfh1] = fn_demarw(xpo,xprob)
%    Demarcate the .68 and .90 error bands given positions and properly scaled probabilities (weights).
%    Calls fn_histwpdfg() first.
%
% xpo:  ndraws-by-n.  Centered position
% xprob:  ndraws-by-n.  Must be properly scaled probabilities corresponding to xpo.
%-----------------
% imfl:  lower .90 bound, n-by-1
% imfh:  higher .90 bound, n-by-1
% imfl1: lower .68 bound, n-by-1
% imfh1: higher .68 bound, n-by-1
% imfmedian: .50 estimate, n-by-1.
%
% Note: fn_histwpdfg() must be called first to get xpo and xprob.
% August 1999 Tao Zha; Revised, August 2000
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


[ndraws,n]=size(xpo);

%$$$ .68 and .95 probability bands
imfl = zeros(n,1);     % preallocating
imfh = zeros(n,1);     % preallocating
imfl1 = zeros(n,1);     % preallocating
imfh1 = zeros(n,1);     % preallocating
imfmedian = zeros(n,1);     % preallocating
imfpo = zeros(4,n);    % 4 positions: l,h,l1,h1.
%
%tic
tem = cumsum(xprob);    % cumulative probability
tem = tem .* 100;
%
%@@@ the following operations are valid only because tem are increasing!
for k = 1:n
   %
   %@@@ the following operations are valid only because tem are increasing!
   %** 5% low tail
   if isempty(max(find(tem(:,k)<5)))
      imfl(k) = xpo(1,k);
   else
      imfl(k) = xpo(max(find(tem(:,k)<5)),k);
   end
   %** 16% low tail
   if isempty(max(find(tem(:,k)<16)))
      imfl1(k) = xpo(1,k);
   else
      imfl1(k) = xpo(max(find(tem(:,k)<16)),k);
   end
   %** 50% estimate
   if isempty(max(find(tem(:,k)<50)))
      imfmedian(k) = xpo(1,k);
   else
      imfmedian(k) = xpo(max(find(tem(:,k)<50)),k);
   end
   %** 5% high tail
   imfh(k) = xpo(min(find(tem(:,k)>95)),k);
   %** 16% low tail
   imfh1(k) = xpo(min(find(tem(:,k)>84)),k);
end
