function [G1,impact,gev,eu]=msv_simple(g0,g1,psi,pi)
% function [G1,C,impact,gev,eu]=msv_one(g0,g1,psi,pi)
% System given as
%        g0*y(t)=g1*y(t-1)+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+impact*z(t).
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  
% eu=[-2,-2] for coincident zeros.
% By Daniel Waggoner -- Based on code by Christopher A. Sims

realsmall=1e-10;

G1=[];
impact=[];
gev=[];

n=size(pi,1);
s=size(pi,2);

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

% determine number unstable roots
nunstable=0;
i=n;
while i > 0
    if (abs(b(i,i)) > abs(a(i,i)))
        nunstable=nunstable+1;
        i=i-1;
    else
        break;
    end
end

eu=[0;0];
gev=[diag(a) diag(b)];

% No bounded MSV solutions?
if nunstable > s
    disp('No bounded MSV solutions')
    return;
end

% Multiple or unique bounded MSV solutions?
if nunstable < s
    disp('Multiple bounded MSV solutions should exist');
else
    disp('Unique bounded MVS solution');
end

[G1,impact,eu]=qzcomputesolution(a,b,q,z,psi,pi,s);
