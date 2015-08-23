function [ F1 F2 G1 G2 V ] = MSV_Newton(P, A, B, Psi, s, x)
%
% Computes MSV solution of 
%
%   A(s(t) x(t) = B(s(t) x(t-1) + Psi(s(t)) epsilon(t) + Pi eta(t)
%
% using Newton's Method.  Assumes that Pi' = [zeros(s,n-s)  eye(s)]
%
% Assumes that A{i} is invertable for 1 <= i <= h

max_count=1000;
tol=1e-5;

h=size(P,1);
n=size(A{1},1);
r=size(Psi{1},2);

if nargin <= 5
    x=zeros(h*s*(n-s),1);
end
f=zeros(h*s*(n-s),1);

Is=eye(s);
I=eye(n-s);

% The following code would be used if Pi were passed instead of the
% assumption that Pi' = [zeros(s,n-s)  eye(s)].
% for i=1:h
%     [Q,R] = qr(Pi{i});
%     W=[zeros(n-s,s)  I; inv(R(1:s,1:s)); zeros(s,n-s)]*Q';
%     A{i}=W*A{i};
%     B{i}=W*B{i};
%     Psi{i}=W*Psi{i};
% end

U=cell(h,1);
for i=1:h
    U{i}=inv(A{i});
end
C=cell(h,h);
for i=1:h
    for j=1:h
        C{i,j}=P(i,j)*B{j}*U{i};
    end
end

cont=true;
count=1;
D=zeros(h*s*(n-s),h*s*(n-s));
X=cell(h,1);
for i=1:h
    X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
end
while cont
    for i=1:h
        for j=1:h
            W1=C{i,j}*[I; -X{i}];
            W2=W1(1:n-s,:);
            D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = kron(W2',Is);
            if i == j
                W1=zeros(s,n);
                for k=1:h
                    W1=W1+[X{k}  Is]*C{i,k};
                end
              W2=-W1(:,n-s+1:end);
              D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) = D((i-1)*s*(n-s)+1:i*s*(n-s),(j-1)*s*(n-s)+1:j*s*(n-s)) + kron(I,W2);
            end
        end
    end
    
    for i=1:h
        mf=zeros(s,n-s);
        for j=1:h
            mf=mf+[X{j} Is]*C{i,j}*[I; -X{i}];
        end
        f((i-1)*s*(n-s)+1:i*s*(n-s))=reshape(mf,s*(n-s),1);
    end
    
    y=D\f;
    x=x - y;
    
    if (count > max_count) || (norm(y) < tol)
        cont=false;
    end
    
    count=count+1;
    for i=1:h
        X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
    end
end

F1=cell(h,1);
F2=cell(h,1);
G1=cell(h,1);
G2=cell(h,1);
V=cell(h,1);
pi=[zeros(n-s,s); Is]; 
for i=1:h
    X=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
    V{i}=U{i}*[I; -X];
    W=[A{i}*V{i} pi];
    F=W\B{i};
    F1{i}=F(1:n-s,:);
    F2{i}=F(n-s+1:end,:);
    G=W\Psi{i};
    G1{i}=G(1:n-s,:);
    G2{i}=G(n-s+1:end,:);
end
  