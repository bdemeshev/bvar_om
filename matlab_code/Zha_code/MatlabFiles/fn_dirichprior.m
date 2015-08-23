function [Galpha,stdthe] = fn_dirichprior(Glamda,std_maxv)
% [Galpha,stdthe] = fn_dirichprior(Glamda,std_maxv)
%   Setup of the Dirichlet prior on the transition matrix P.  Given the mode of each random
%     variable theta_j, solve for the (hyper)parameter \alpha where one of \alpha, say, \alpha_j
%     is free for controlling the overall variance of \alpha consisdent with the given modes.
%
% Glamda: h-by-1 vector of prior modes for the random variables \theta.
%         For each column of P (transitional matrix), Glamda'll be different with a large weight on a diaganol.
% std_maxv:  prior standard deviation for the selected random variable
%            (here the one with the maximum (largest) mode).
%----------
% Galpha: h-by-1 vector of parameters' values for the Dirichlet distribtuion.  If one of the
%    elements in Galpha is Inf, we mean a degenerate case where stdthe is set to zeros.
% stdthe:  h-by-1 vector of standard deviations for all random variables \theta.
%
% See Gelman, et al p. 476 and TVBvar Note pp. 24-31
% Tao Zha, March 2001
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


Glamda = Glamda(:);   % Column vector!  Prior modes for the random variables \theta.
[maxlam,indxlam] = max(Glamda);  % index for the selected hyperparameter that is freely set to match the variance of the corresponding random variable.
ndiri = length(Glamda);   % the number of random variables in the Dirichlet distribution.
Galpha=zeros(ndiri,1);  % \alpha's: hyperparameters that need to be choosen. Must be all greater than 1.
stdthe = zeros(ndiri,1);   % Standard deviations for each \theta
if indxlam>ndiri
   warning('Make sure indxlam is less than the total number of parameters in the Dirichlet')
   error(' ')
end
if (maxlam==1)   % We treat this as a degenerate case for convenient coding.
   Galpha(indxlam)=Inf;     % Note the use of Inf, NOT NaN.  Because
         %  0(or any finite number)/Inf=0 while 0/NaN gives a nonsense (i.e., NaN again).
         %  See tveml.m and pp. 27a-27b for more information.
   return   % The function stops here.
end


%*** Given std_maxv, find the corresponding hyperparameter Galpha_j by soloving
%      a polynomial of order 3.  See Note (*) page 30.
vj=std_maxv^2;  % variance of the selected random variable \theta.
Gsi = (1-Glamda(indxlam))/Glamda(indxlam);  % unchanged value, given Glamda.
coevec = zeros(4,1);  % 4 because we're going to consider a polynomial of order 3 (+ constant).
coevec(1) = vj*(1+Gsi)^3;
coevec(2) = vj*(1+Gsi)^2*(3*(ndiri-Gsi)-2) - Gsi;
coevec(3) = (ndiri-Gsi-1)*(vj*(1+Gsi)*(3*(ndiri-Gsi)-1)-1);
coevec(4) = vj*(ndiri-Gsi)*(ndiri-Gsi-1)^2;

Galphaj = max(roots(coevec)); % j: jth element, set to match the variance.  Must be greater than 1.
if (Galphaj<=1)
   disp(' ')
   warning('Error!  Make sure std_maxv is small enough to have a large Galphaj > 1')
   error(' ')
end

%*** Solve for the rest of \alpha's.  See Note page 27.
A = eye(ndiri) - repmat(Glamda,[1 ndiri]);
C = ones(ndiri,1) + (Galphaj-ndiri)*Glamda;
A1=A; C1=C;
A1(indxlam,:)=[]; A1(:,indxlam)=[]; C1(indxlam)=[];
indxrest=1:ndiri; indxrest(indxlam)=[];
Galpha(indxrest) = A1\C1;
Galpha(indxlam) = Galphaj;

%========= Debugging or checking the properties ==========
%disp(' ')
%disp('2nd-to-first-or-rest ratio: it should not change with the last element of \alpha')
%(Galpha(2)-1)/(Galpha(1)-1)
%(Galpha(2)-1)/(sum(Galpha([1 3]))-2)
%disp(' ')
%disp('Solution for \alpha')
%Galpha

for k=1:ndiri
   tmp=Galpha;
   tmp(k)=[];
   stdthe(k)=sqrt( Galpha(k)*sum(tmp)/(sum(Galpha)^2*(sum(Galpha)+1)) );
end
%disp(' ')
%disp('Standard deviations for each \theta')
%stdthe
