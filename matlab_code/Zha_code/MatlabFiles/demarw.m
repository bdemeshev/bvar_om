function [imfl,imfh,imfl1,imfh1] = demarw(xpo,xprob)
%  [imfl,imfh,imfl1,imfh1] = demarw(xpo,xprob)
%    Demarcate the .68 and .95 error bands given prob. (weights) and positions
%
% xpo:  ndraws-by-nvar.  Centered position
% xprob:  ndraws-by-nvar.  Properly scaled probabilities corresponding to xpo.
%-----------------
% imfl:  lower .95 bound, nvar-by-1
% imfh:  higher .95 bound, nvar-by-1
% imfl1: lower .68 bound, nvar-by-1
% imfh1: higher .68 bound, nvar-by-1
%
% Copyright (c) August 1999 Tao Zha
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


[ndraws,nvar]=size(xpo);

%$$$ .68 and .95 probability bands
imfl = zeros(nvar,1);     % preallocating
imfh = zeros(nvar,1);     % preallocating
imfl1 = zeros(nvar,1);     % preallocating
imfh1 = zeros(nvar,1);     % preallocating
imfpo = zeros(4,nvar);    % 4 positions: l,h,l1,h1.
%
%tic
tem = cumsum(xprob);    % cumulative probability
tem = tem .* 100;
%
%@@@ the following operations are valid only because tem are increasing!
for k = 1:nvar
   %
   %@@@ the following operations are valid only because tem are increasing!
   %** 2.5% low tail
   if isempty(max(find(tem(:,k)<2.5)))
      imfl(k) = xpo(1,k);
   else
      imfl(k) = xpo(max(find(tem(:,k)<2.5)),k);
   end
   %** 16% low tail
   if isempty(max(find(tem(:,k)<16)))
      imfl1(k) = xpo(1,k);
   else
      imfl1(k) = xpo(max(find(tem(:,k)<16)),k);
   end
   %** 2.5% high tail
   imfh(k) = xpo(min(find(tem(:,k)>97.5)),k);
   %** 16% low tail
   imfh1(k) = xpo(min(find(tem(:,k)>84)),k);
end
