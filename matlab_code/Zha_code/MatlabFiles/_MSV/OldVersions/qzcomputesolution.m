function [G1,impact,eu]=qzcomputesolution(a,b,q,z,psi,pi,suppress)
% function [G1,C,impact,eu]=qzcomputesolution(a,b,q,z,c,psi,pi,suppress)
% System given as
%        (q'*a*z')*y(t)=(q'*b*z')*y(t-1)+psi*e(t)+pi*eta(t),
% with e an exogenous variable process and eta an endogenously determined
% one-step-ahead expectational error.  q and z are unitary matrices and a
% and b are upper triangular matrices that are computed via a generalized
% Schur decomposition.  It is assumed that the ordering of the system
% is such that complex conjugate pairs are together.  If a complex
% conjugate pair is split, then G1 and impact will be complex.
% Returned system is
%        y(t)=G1*y(t-1)+impact*z(t).
% eu(1)=1 for existence, eu(2)=1 for uniqueness.
%
% The last suppress roots are suppressed.
%
% By Daniel Waggoner -- Based on code by Christopher A. Sims

realsmall=sqrt(eps);

eu=[0;0];
G1=[];
impact=[];

n=size(pi,1);
s=size(pi,2);

q1=q(1:n-suppress,:);
q2=q(n-suppress+1:n,:);

[u1 d1 v1]=svd(q2*pi,'econ');
if suppress > 0
    m=sum(diag(d1) > realsmall*d1(1,1));
else
    m=0;
end
u1=u1(:,1:m);
d1=d1(1:m,1:m);
v1=v1(:,1:m);

if (m == s)
    eu(2)=1;
else
    eu(2)=0;
end

[u2 d2 v2]=svd(q2*psi,'econ');
if suppress > 0
    m=sum(diag(d2) > realsmall*d2(1,1));
else
    m=0;
end
u2=u2(:,1:m);
d2=d2(1:m,1:m);
v2=v2(:,1:m);

if norm((eye(suppress) - u1*u1')*u2) < realsmall
    eu(1)=1;
else
    eu(1)=0;
    return;
end


% Compute impact and G0
a11_inv=inv(a(1:n-suppress,1:n-suppress));
x=a11_inv*q1*pi*v1*diag(1./diag(d1))*u1';
impact=[a11_inv*q1*psi - x*q2*psi; zeros(suppress,size(psi,2))];

G1=zeros(n,n);
G1(1:n-suppress,1:n-suppress)=a11_inv*b(1:n-suppress,1:n-suppress);

% this is to make the answer agree with gensys - this piece does not effect
% the answer, at least after the initial period.
G1(1:n-suppress,n-suppress+1:n)=a11_inv*b(1:n-suppress,n-suppress+1:n) - x*b(n-suppress+1:n,n-suppress+1:n);

% Convert back to y
impact=z*impact;
G1=z*G1*z';

% Is the solution real?
if (suppress == 0) | (suppress == n) | (abs(b(n-suppress+1,n-suppress+1)*conj(a(n-suppress,n-suppress)) - a(n-suppress+1,n-suppress+1)*conj(b(n-suppress,n-suppress))) > realsmall)
    impact=real(impact);
    G1=real(G1);
end

%--------------------------------------------------------------------------
% The code below check that a solution is obtained with stable manifold
% z1.  Should be commented out for production runs.
%--------------------------------------------------------------------------
% z2=z(:,n-suppress+1:n);
% z1=null(z2');
% A=q'*a*z';
% B=q'*b*z';
%
% n=norm(A*(G1*z1) - B*z1);
% disp('norm(A*(g1*z1) - B*z1) -- should be zero');
% disp(n);
%
% eta=-v1*diag(1./diag(d1))*u1'*q2*psi;
% n=norm(A*impact - psi - pi*eta);
% disp('norm(A*impact - psi - pi*eta) -- should be zero');
% disp(n);
%
% n=norm(z2'*G1*z1);
% disp('norm(z2''*G1*z1) -- should be zero');
% disp(n);
%
% n=norm(z2'*G1);
% disp('norm(z2''G1) -- best if zero');
% disp(n);
%
% n=norm(z2'*impact);
% disp('norm(z2''*impact) -- should be zero');
% disp(n);
%
% n=norm(G1*z2);
% disp('norm(G1*z2)');
% disp(n);
%
% disp('press any key to continue');
% pause




