function [G1,impact,gev,eu]=msv_one(g0,g1,psi,pi)
% function [G1,C,impact,gev,eu]=msv_one(g0,g1,psi,pi)
% System given as
%        g0*y(t)=g1*y(t-1)+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+impact*z(t).
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Daniel Waggoner -- Based on code by Christopher A. Sims

realsmall=1e-10;

n=size(pi,1);
s=size(pi,2);

[a b q z]=qz(g0,g1);

for i=1:n
  if (abs(a(i,i)) < realsmall) & (abs(b(i,i)) < realsmall)
    disp('Coincident zeros.')
    eu=[-2;-2];
    gev=[diag(a) diag(b)];
    G1=[]
    impact=[]
    return
  end
end

[a b q z]=qzsort(a,b,q,z);

% determine complex conjugate pairs
pairs=zeros(n,1);
i=1;
while i < n
  if abs(b(i+1,i+1)*conj(a(i,i)) - a(i+1,i+1)*conj(b(i,i))) < realsmall
    pairs(i)=1;
    pairs(i+1)=-1;
    i=i+2;
  else
    i=i+1;
  end
end

% determine number unstable roots
nunstable=0;
i=1;
while i <= n
  if pairs(i) == 0
    if (abs(b(i,i)) > abs(a(i,i)))
      nunstable=n-i+1;
      break;
    else
      i=i+1;
    end
  else
    if (abs(b(i,i)) > abs(a(i,i))) | (abs(b(i+1,i+1)) > abs(a(i+1,i+1)))
      nunstable=n-i+1;
      break;
    else
      i=i+2;
    end
  end
end

% No bounded MSV solutions?
if nunstable > s
    disp('No bounded MSV solutions')
    eu=[0;0];
    gev=[diag(a) diag(b)];
    G1=[]
    impact=[]
    return;
end

% Multiple bounded MSV solutions?
if nunstable < s
    disp('Multiple bounded MSV solutions should exist');

    % Get roots to suppress
    idx=zeros(n,1);
    idx_size=zeros(n,1);
    k=1;
    idx(1)=n-nunstable;
    cont=1;
    while cont > 0
        if (pairs(idx(k)) == -1)
            idx_size(k)=2;
        else
            idx_size(k)=1;
        end

        d=nunstable;
        for i=1:k
            d=d+idx_size(i);
        end

        if (d == s)
            cont=cont+1;
            [an bn qn zn]=qzmoveindex(a,b,q,z,pairs,idx,k,n-nunstable);
            [G1,impact,eu]=qzcomputesolution(an,bn,qn,zn,psi,pi,s);
            if (eu(1) == 1) & (eu(2) == 1)
                gev=[diag(an) diag(bn)];
                return;
            else
                disp('MSV solution does not exist for the ordering:');
                abs(diag(bn)./diag(an))
                eu
            end
        end

        if (d >= s)
            if idx(k) > idx_size(k)
                idx(k)=idx(k)-idx_size(k);
            else
                if k > 1
                    k=k-1;
                    idx(k)=idx(k)-idx_size(k);
                else
                    disp('No MSV solutions exist');
                    return;
                end
            end
        else
            if idx(k) > idx_size(k)
                idx(k+1)=idx(k)-idx_size(k);
                k=k+1;
            else
                if k > 1
                    k=k-1;
                    idx(k)=idx(k)-idx_size(k);
                else
                    disp('No MSV solutions exist');
                    return;
                end
            end
        end
    end
end

% Unique MSV solution?
disp('Unique bounded MVS solution');
[G1,impact,eu]=qzcomputesolution(a,b,q,z,psi,pi,s);
gev=[diag(a) diag(b)];

