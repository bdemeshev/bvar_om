function [G1,impact,eu]=qzcomputesolution(a,b,q,z,psi,pi,suppress)
% System given as
%
%     (q'*a*z')*y(t) = (q'*b*z')*y(t-1) + psi*epsilon(t) + pi*eta(t),
%
% with epsilon(t) an exogenous process and eta(t) an endogenously
% determined one-step-ahead expectational error.  The matrices q and z are
% unitary and a and b are upper triangular.  Furthermore, it is assumed
% that if a11 is the upper (n-suppress) x (n-suppress) block of a, then a11
% is non-singular.
%
% The span of the returned solution is equal to the span of the first
% (n-suppress) columns of z.
%
% If t
% It is assumed that the ordering of
% the system is such that complex conjugate pairs are together.  If a
% complex conjugate pair is split, then G1 and impact will be complex.
%
% Returned solution is
%
%        y(t) = G1*y(t-1) + impact*z(t).
%
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

% find the dimension of the column space of q2*pi
if suppress > 0
    % all singular values small - q2*pi is the zero matrix
    if d1(1,1) < realsmall
        m=0;
    else
        % exclude the small singular values relative to largest
        m=sum(diag(d1) > realsmall*d1(1,1));
    end
else
    m=0;
end
u1=u1(:,1:m);
d1=d1(1:m,1:m);
v1=v1(:,1:m);

if (m == suppress)
    eu(2)=1;
end

% find othogonal basis for span of q2*psi
[u2 d2 v2]=svd(q2*psi,'econ');
if suppress > 0
    % all singular values small - q2*psi is the zero matrix
    if d2(1,1) < realsmall
        m=0;
    else
        % exclude the small singular values relative to largest
        m=sum(diag(d2) > realsmall*d2(1,1));
    end
else
    m=0;
end
u2=u2(:,1:m);
d2=d2(1:m,1:m);
v2=v2(:,1:m);

% is the projection of the column space of q2*psi on to the column space of
% the column space of q2*pi the identity mapping?
if norm((eye(suppress) - u1*u1')*u2) > realsmall
    return
end

eu(1)=1;

% Compute impact and G0
a11_inv=inv(a(1:n-suppress,1:n-suppress));
x=a11_inv*q1*pi*v1*diag(1./diag(d1))*u1';
impact=[a11_inv*q1*psi - x*q2*psi; zeros(suppress,size(psi,2))];

G1=zeros(n,n);

G1(1:n-suppress,1:n-suppress)=a11_inv*b(1:n-suppress,1:n-suppress);

% Uncomment the line below to make the answer agree with gensys - this
% piece does not effect the answer, at least after the initial period.
G1(1:n-suppress,n-suppress+1:n)=a11_inv*b(1:n-suppress,n-suppress+1:n) - x*b(n-suppress+1:n,n-suppress+1:n);

% Convert back to y
impact=z*impact;
G1=z*G1*z';

% Is the solution real?
% if (suppress == 0) | (suppress == n) | (abs(b(n-suppress+1,n-suppress+1)*conj(a(n-suppress,n-suppress)) - a(n-suppress+1,n-suppress+1)*conj(b(n-suppress,n-suppress))) > realsmall)
%     impact=real(impact);
%     G1=real(G1);
% end

%--------------------------------------------------------------------------
% The code below check that a solution is obtained with stable manifold
% z1.  Should be commented out for production runs.
%--------------------------------------------------------------------------
z2=z(:,n-suppress+1:n);
z1=null(z2');
%if check_solution_AR(G1,impact,q'*a*z',q'*b*z',psi,pi,z1) == 1
%    eu(1)=3;
%end
