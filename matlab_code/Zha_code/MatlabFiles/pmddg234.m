function [g,badg] = pmddg234(x,idmat0,idmat1,fss,nvar, ...
                               ncoef,phi,y,A0b,sg0bid,sgpbid);
% Leeper-Sims-Zha BVAR setup, analytical gradient for "csminwel.m"
% general program to setup A0 matrix and compute the likelihood
% requires x (parameter vector), a0indx (matrix indicating the free
% parameters in A0), fss (forecast sample size), nvar (number of variables),
% ncoef (number of coefficients in a single equation in A+), phi (r.h.s.
% observations in the system for A+), y (l.h.s. observations for A0),
% ymy (Y'*M*Y), xd0 (X(0) for dummy observations y*(1)), ys1 (dummy initial
% observations y*(1)), xtx (X'X or phi'*phi), A0b (initial prior on A0),
% sg0bid (diagonal of Sigma0_bar on the parameters in i-th equation in A0,
% initial prior covariace matrix), sgpbid (diagonal of Sigma+_bar on the
% parameters in i-th equation in A0, initial prior covariance matrix).
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

global vaya0 vaha0 vaxah vahad vbiga FRESHFUNCTION
% CAS added FRESHFUNCTION 8/20/96 to allow use of pmddg23 in a numerical Hessian calculation, where
% it is invoked repeatedly without new function evaluations.
if ~FRESHFUNCTION
   fjnk = pmddf234(x,idmat0,idmat1,fss,nvar,ncoef,phi,y,A0b,sg0bid,sgpbid);
end
FRESHFUNCTION=0;
badg = 0;
%
a0indx=find(idmat0);
nhp = 0;                 % <<>> 4 hyperparameters
na0p = length(a0indx);    % <<>> number of A0 parameters
nfp = na0p+nhp;
g = zeros(nfp,1);
A0h= zeros(nvar);
A0h(a0indx) = x(1:na0p);
               %  restrictions on A0

%disp(sprintf('Starting loop (gradient): %g',toc))


%%%%%%%
%%% analytical gradient for A0 below
%%%%%%%
%
% ** dla0 = dlog|A0|/d(a0??)
vdla0 = zeros(na0p,1);
dla0 = inv(A0h)';
dla0 = dla0(:);
vdla0 = dla0(a0indx);

% without the prior |A0|^k
g = -fss*vdla0 + 0.5*vaya0 + 0.5*vaha0 + 0.5*vaxah + 0.5*vahad - 0.5*vbiga;

%disp(sprintf('Loop end (gradient): %g', toc))