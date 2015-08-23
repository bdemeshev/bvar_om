function [F1, F2, G1, G2, V, err] = fwz_msv_msre(P, A, B, Psi, s, x, max_count, tol)
%[ F1 F2 G1 G2 V err ] = fwz_msv_msre(P, A, B, Psi, s, x, max_count, tol)
% Computes MSV solution of 
%
%   A(s(t)) x(t) = B(s(t)) x(t-1) + Psi(s(t)) epsilon(t) + Pi eta(t)
%
% using Newton's Method.  Assumes that Pi' = [zeros(s,n-s)  eye(s)] and
% that A{i} is invertible.  P is the transition matrix and P(i,j) is the
% probability that s(t+1)=j given that s(t)=i.  Note that the rows of P
% must sum to one.  x is the initial value and if not passed is set to
% zero.  max_count is the maximum number of iterations on Newton's method
% before failing and tol is the convergence criterion.
%
% The solution is of the form
%
%   x(t) = V{s(t)}*F1{s(t)}*x(t-1) + V{s(t)}*G1{s(t)}*epsilon(t)
%
%  eta(t) = F2{s(t)}*x(t-1) + G2{s(t)}*epsilont(t)
%
% A positive value of err is the number of iterations needed to obtain 
% convergence and indicates success.  A negitive value of err is the number of 
% iterations before the method terminated without convergence and indicates
% failure.
%
% Copyright (C) 1997-2013 Roger Farmer, Daniel F. Waggoner, and Tao Zha
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


h=size(P,1);
n=size(A{1},1);
r=size(Psi{1},2);

if (nargin <= 7) || (tol <= 0)
    tol=1e-5;
end

if (nargin <= 6) || (max_count <= 0)
    max_count=1000;
end

if nargin <= 5
    x=zeros(h*s*(n-s),1);
end
f=zeros(h*s*(n-s),1);

Is=eye(s);
I=eye(n-s);

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
    
    if (count > max_count) || (norm(f) < tol)
        cont=false;
    else
        count=count+1;
        for i=1:h
            X{i}=reshape(x((i-1)*s*(n-s)+1:i*s*(n-s)),s,n-s);
        end
    end
end

if (norm(f) < tol)
    err=count;
else
    err=-count;
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
  