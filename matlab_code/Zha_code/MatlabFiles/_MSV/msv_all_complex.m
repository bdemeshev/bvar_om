function [Gamma,Omega]=msv_all_complex(A1,A2,A3,B1,B2,B3,C1)
% System given as
%
%        E_t[y(t+1)] = A1*u(t) + A2*z(t) + A3*y(t)
%
%             z(t+1) = B1*u(t) + B2*z(t) + B3*y(t)
%
%             u(t+1) = C1*u(t) + epsilon(t+1)
%
% for t >= 0.  The dimension of u(t) is r, the dimension of z(t) is k, and
% the dimension of y(t) is s.  E_t[epsilon(t+1)]=0.
%
% Solution technique:
%
% Assume a solution of the form y(t) = Gamma*u(t) + Omega*z(t).  This
% implies that 
%
%  Gamma*C1*u(t) + Omega*(B1*u(t) + B2*z(t) + B3*(Gamma*u(t) + Omega*z(t)))
%     = A1*u(t) + A2*z(t) + A3*(Gamma*u(t) + Omega*z(t))
%
% Gathering terms, we have that 
%
%  (1) (Omega*B2 + Omega*B3*Omega - A2 - A3*Omega)*z(t) = 0
%
%  (2) (Gamma*C1 + Omega*B1 + Omega*B3*Gamma - A1 - A3*Gamma)*u(t) = 0
%
% Let 
%
%    X = [ -B2 -B3
%          -A2 -A3 ]
%
% and let X = U*T*U' be a Schur decomposition of X, where U is unitary, T
% is upper triangular, and ' denotes the conjugate transpose.  If
%
%    U = [ U11  U12
%          U21  U22 ] 
%
% and U11 is invertiable, then 
%
%  (3)  Omega = U21*inv(U11) 
%
% will be a solution of (1) for any value of z(t) in k-dimensional
% euclidean space.  On the other hand, if Omega is a solution of (1) for
% all z(t) in k-dimensional euclidean space, then Omega is of the form of
% (3) for some Schur decomposition of X.
% 
% Given Omega, the solution of the second equation is
%
% vec(Gamma) = inv(kron(C1',eye(s)) + kron(eye(r),Omega*B3 - A3))
%                                                       *vec(A1 - Omega*B1)
%
% The solution Omega of (1) depends only on the column space
% 
%    [ U11
%      U21 ]
%
% which we denote by V.  It is well known that the eigenvalues of X appear
% along the diagonal of T.  If the eigenvalues of X are distinct, then the
% ordering of the eigenvalues uniquely determines U.  However, when the
% dimension of the eigenspace of some eigenvalue is greater than one, this
% is no longer true and reordering the Schur decomposition is no longer
% sufficient to capture all possible solutions.
%                                                        
% By Daniel Waggoner 

realsmall=sqrt(eps);

r=size(A1,2);
k=size(A2,2);
s=size(A3,2);
kk=k+s;

Gamma=cell(0,1);
Omega=cell(0,1);

if k == 0
    Y=kron(C1',eye(s)) + kron(eye(r),-A3);
    if rank(Y) == size(Y,1)
        Omega{1,1}=zeros(s,k);
        Gamma{1,1}=reshape(inv(Y)*reshape(A1,s*r,1),s,r);

        %------------------------------------------------------------------
        % Check if solution - comment out for production runs
        % E_t[y(t+1)] = E_t[Gamma*u(t+1) + Omega*z(t+1)]
        %             = Gamma*C1*u(t) + Omega*(B1*u(t) + B2*z(t) 
        %                                   + B3*(Gamma*u(t) + Omega*z(t)))
        %             = A1*u(t) + A2*z(t) + A3*(Gamma*u(t) + Omega*z(t))
        %------------------------------------------------------------------
        n1=norm(Gamma*C1 - (A1 + A3*Gamma));
        if n1 > realsmall
            disp('Should be zero - press any key to continue');
            disp('Gamma*C1 - (A1 + A3*Gamma)');
            disp(n1);
            pause;
        end
        %------------------------------------------------------------------

    end
else
    % get and order Schur decomposition
    X=[ -B2 -B3
        -A2 -A3
        ];

    [U,T]=schur(X,'complex');
    [d,id]=sort(abs(diag(T)),'descend');
    [id,idx]=sort(id);
    [U,T]=ordschur(U,T,idx);
    u=U; t=T;
    
    idx=ones(s+k,1);
    for i=0:s-1
        idx(s+k-i)=0;
    end
    
    ii=1;
    cont=1;
    while cont == 1
        
        u11=u(1:k,1:k);
        u21=u(k+1:k+s,1:k);

        if rank(u11) == k
            omega=u21/u11;
            Y=kron(C1',eye(s)) + kron(eye(r),omega*B3 - A3);
            if rank(Y) == size(Y,1)
                Omega{ii,1}=omega;
                Gamma{ii,1}=reshape(Y\reshape(A1 - omega*B1,s*r,1),s,r);
                ii=ii+1;

                %----------------------------------------------------------
                % Check if solution - comment out for production runs
                % E_t[y(t+1)] = E_t[Gamma*u(t+1) + Omega*z(t+1)]
                %             = Gamma*C1*u(t) + Omega*(B1*u(t) + B2*z(t) 
                %                           + B3*(Gamma*u(t) + Omega*z(t)))
                %             = A1*u(t) + A2*z(t) + A3*(Gamma*u(t) 
                %                                             + Omega*z(t))
                %----------------------------------------------------------
                n1=norm(Gamma{ii-1,1}*C1 + Omega{ii-1,1}*(B1 + B3*Gamma{ii-1,1}) - (A1 + A3*Gamma{ii-1,1}));
                n2=norm(Omega{ii-1,1}*(B2 + B3*Omega{ii-1,1}) - (A2 + A3*Omega{ii-1,1}));
                if (n1 > realsmall) || (n2 > realsmall)
                    i
                    idx
                    disp('Both should be zero - press any key to continue');
                    disp('Gamma*C1 + Omega*(B1 + B3*Gamma) - (A1 + A3*Gamma)');
                    disp(n1);
                    disp('Omega*(B2 + B3*Omega) - (A2 + A3*Omega)');
                    disp(n2);
                    pause
                end
                %----------------------------------------------------------             
            end
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
        
        [u,t]=ordschur(U,T,idx);
    end
end


