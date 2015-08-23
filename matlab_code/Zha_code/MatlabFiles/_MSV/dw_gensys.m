function [G1,impact,gev,eu]=dw_gensys(g0,g1,psi,pi,div)
% function [G1,impact,gev,eu]=gensys(g0,g1,c,psi,pi,div)
% System given as
%        g0*y(t)=g1*y(t-1)+psi*e(t)+pi*eta(t),
% with e an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+impact*e(t).
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Daniel Waggoner

realsmall=1e-8;

G1=[];
impact=[];
gev=[];

n=size(pi,1);

if nargin == 4
    div=1.0;
end

[a b q z]=qz(g0,g1);

for i=1:n
  if (abs(a(i,i)) < realsmall) & (abs(b(i,i)) < realsmall)
    disp('Coincident zeros.')
    eu=[-2;-2];
    gev=[diag(a) diag(b)];
    return
  end
end

[a b q z]=qzsort(a,b,q,z);

suppress=0;
i=n;
while (i > 0) & (abs(b(i,i)) > abs(a(i,i)*div))
    suppress=suppress+1;
    i=i-1;
end

[G1,impact,eu]=qzcomputesolution(a,b,q,z,psi,pi,suppress);

