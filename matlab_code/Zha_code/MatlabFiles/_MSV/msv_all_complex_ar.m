function [G1,Impact,gev,z2,err]=msv_all_complex_ar(A,B,Psi,Pi)
% System given as
%
%        A*y(t) = B*y(t-1) + psi*epsilon(t) + pi*eta(t),
%
% with epsilon(t) an exogenous process and eta(t) endogenously determined
% one-step-ahead expectational errors.  Returned solutions are
%
%        y(t)=G1{i}*y(t-1)+Impact{i}*epsilon(t).
%
% The span of the each of the solutions returned will be n-s and the
% solution will be uniquely determined by its span.  gev returns the
% generalized eigenvalues, of which the last s were suppressed.  z2 is the
% subspace that is perpendicular to the span of the solution.
%
% err is the error return.
%  err = 0 -- success, at least one msv solution found.
%  err = 1 -- no msv solutions
%  err = 2 -- degenerate system
%
% Solutions are obtained via calls to qzcomputesolution().
%
% By Daniel Waggoner

realsmall=sqrt(eps);

n=size(Pi,1);
s=size(Pi,2);

[a,b,q,z]=qz(A,B,'complex');

G1=cell(0,1);
Impact=cell(1,1);
gev=cell(0,1);
gev{1}=[diag(a) diag(b)];
z2=cell(0,1);
err=1;

% determine if the system is degenerate and count the number of diagonal
% elements of a that are non-zero.
kk=0;
for i=1:n
  if abs(a(i,i)) < realsmall
      if abs(b(i,i)) < realsmall
          disp('Coincident zeros.')
          err=2;
          return
      end
  else
      kk=kk+1;
  end
end

% too many roots must be suppressed. no msv type solutions.
if (n-kk > s)
    return
end

% sort the generalized eigenvalues
[d,id]=sort(abs(diag(b)./diag(a)),'descend');
[id,idx]=sort(id);
[a0,b0,q0,z0]=ordqz(a,b,q,z,idx);
a=a0; b=b0; q=q0; z=z0;

idx=ones(n,1);
for i=0:s-1
    idx(n-i)=0;
end

ii=0;
cont=1;
cntmax = 500;
cntnumber = 1;
while (cont == 1) && (cntnumber <=cntmax)
    % compute solution
    [g1,impact,eu]=qzcomputesolution(a,b,q,z,Psi,Pi,s);

    % save solution if unique
    if (eu == [1;1])
        ii=ii+1;
        G1{ii,1}=g1;
        Impact{ii,1}=impact;
        gev{ii,1}=[diag(a) diag(b)];
        z2{ii,1}=z(:,n-s+1:n);
        err=0;
    end

    % increment idx
    cont=0;
    jj=1;
    while (jj < kk) & (idx(jj) == 0)
        idx(jj)=1;
        jj=jj+1;
    end
    for j=jj+1:kk
        if idx(j) == 0
            idx(j)=1;
            cont=1;
            break;
        end
    end
    for i=1:jj
        idx(j-i)=0;
    end

    % reorder roots
    [a,b,q,z]=ordqz(a0,b0,q0,z0,idx);

    %--- Ad hoc addtion.
    cntnumber = cntnumber+1;
end

