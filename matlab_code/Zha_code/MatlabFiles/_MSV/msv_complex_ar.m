function [g1,impact,gev,z2,err]=msv_complex_ar(A,B,Psi,Pi)
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
%  err = [2;2] -- degenerate system
%  err = [1;1] -- success, unique msv solution found.
%  err = [1;0] -- multiple msv solutions
%  err = [3;?] -- msv solutions are unbounded
%
% Solutions are obtained via calls to msv_all_complex_ar().
%
% By Daniel Waggoner

realsmall=sqrt(eps);

% find all msv type solutions
[G1,Impact,Gev,Z2,err]=msv_all_complex_ar(A,B,Psi,Pi);

if err ~= 0
    % no solution exists
    g1=[];
    impact=[];
    gev=Gev{1,1};
    z2=[];
    if err == 2
        err=[-2;-2];
    else
        err=[0;0];
    end
else
    % solution exists
    g1=G1{1};
    impact=Impact{1,1};
    gev=Gev{1,1};
    z2=Z2{1,1};

    % uniqueness
    if size(G1,1) == 1
        err=[1,1];
    else
        n=size(Pi,1);
        m=n-size(Pi,2)+1;
        x1=min(abs(gev(m:n,2)./gev(m:n,1)));
        x2=min(abs(Gev{2,1}(m:n,2)./Gev{2,1}(m:n,1)));
        if abs(x1 - x2) > realsmall
            err=[1,1];
        else
            err=[1,0];
        end
    end

    % boundedness
    m=size(Pi,1)-size(Pi,2);
    if max(abs(gev(1:m,2)./gev(1:m,1))) > 1+1.0e-9
        err(1)=3
    end
end
