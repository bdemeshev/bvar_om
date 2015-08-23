function [Myrqm,nMyrqm] = fn_calyrqm(q_m,Byrqm,Eyrqm)
% [Myrqm,nMyrqm] = fn_calyrqm(q_m,Byrqm,Eyrqm)
%
%    Given the beginning and end years and quarters (months), export a matrix of all years and
%       quarters (months) for these years and in between
%
% q_m:  4 if quarterly and 12 if monthly
% Byrqm:  [year quarter(month)] -- all integers, the begining year and quarter (month)
% Eyrqm:  [year quarter(month)] -- all integers, the end year and quarter (month)
%-------------------
% Myrqm:  matrix of all years and quarters (months) between and incl. Byrqm and Eyrqm
% nMyrqm:  number of data points incl. Byrqm and Eyrqm
%
% Tao Zha, April 2000
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


if ~isempty(find(Byrqm-round(Byrqm))) | (q_m-round(q_m)) | ~isempty(find(Byrqm-round(Byrqm)))
   error('argin qm, Byrqm, or Eyrqm must of integer')
elseif Byrqm(1)>Eyrqm(1)
   error('Eyrqm(1) must be equal to or greater than Byrqm(1)')
elseif Byrqm(1)==Eyrqm(1)
   if Byrqm(2)>Eyrqm(2)
      error('Eyrqm(2) must be equal to or greater than Byrqm(2) because of the same year')
   end
end


Yr = Byrqm(1)+[0:Eyrqm(1)-Byrqm(1)]';

if length(Yr)>=2   %  there are years and quarters (months) between Byrqm and Eyrqm
   n=length(Yr)-2;
   C=zeros(n*q_m,2);
   C(:,1) = kron(Yr(2:end-1),ones(q_m,1));
   C(:,2) = kron(ones(n,1),[1:q_m]');

   %* initialize a matrix of years and quarters (months) including Byrqm and Eyrqm
   Myrqm = zeros((q_m-Byrqm(2)+1)+Eyrqm(2)+n*q_m,2);

   %* Years in between
   n1=q_m-Byrqm(2)+1;
   n2=Eyrqm(2);
   Myrqm(n1+1:end-n2,:) = C;
   %* Beginning year
   for k=1:n1
      Myrqm(k,:) = [Byrqm(1) Byrqm(2)+k-1];
   end
   %* End year
   for k=1:n2
      Myrqm(end-Eyrqm(2)+k,:) = [Eyrqm(1) k];
   end
else     %* all the data are in the same calendar year
   n1=Eyrqm(2)-Byrqm(2)+1;
   Myrqm = zeros(n1,2);
   for k=1:n1
      Myrqm(k,:) = [Byrqm(1) Byrqm(2)+k-1];
   end
end

nMyrqm = size(Myrqm,1);
