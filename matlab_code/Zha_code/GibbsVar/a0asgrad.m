function [g,badg] = a0asgrad(x,s,nobs,nvar,a0indx);
% Computes analytical gradient for use in csminwel.m
%      function [g,badg] = a0asgrad(x,s,nobs,nvar,a0indx);
%
%  x (parameter vector),
%  s (diag(S1,...,Sm)): note, as in "a0lhfun", already divided by "nobs"
%  nobs (no of obs),
%  nvar (no of variables),
%  a0indx (matrix indicating the free parameters in A0, and each column in A0 corresponds
%                    to an equation)


a0 = zeros(nvar);
%%g = zeros(nvar*nvar,1);     % 4/27/97: not necessary.
badg = 0;

a0(a0indx) = x;

b1=-nobs*inv(a0');
b2 = zeros(nvar,nvar);
for i=1:nvar
   b2(:,i) = nobs*s{i}*a0(:,i);
end

%b = -nobs*inv(a') + nobs*s*a;
b = b1(:)+b2(:);
g = b(a0indx);