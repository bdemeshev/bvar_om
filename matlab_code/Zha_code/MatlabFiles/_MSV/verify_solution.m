function diff = verify_solution(P,A,B,Psi,s,F1,F2,G1,G2,V)
%
% Verifies that 
%
%     x(t) = V(s(t))(F1(s(t))x(t-1) + G1(s(t))epsilon(t)
%   eta(t) = F2(s(t))x(t-1) + G2(s(t))epsilon(t)
%
% is a solution of 
% 
%    A(s(t)) x(t) = B(s(t)) x(t-1) + Psi(s(t) epsilon(t) + Pi eta(t)
%
% where Pi' = [zeros(s,n-s)  eye(s)] by verifying that
%
%  [A(j)*V(j)  Pi] [ F1(j)   
%                    F2(j) ] = B(j)
%
%  [A(j)*V(j)  Pi] [ G1(j)   
%                    G2(j) ] = Psi(j)
%
%  (P(i,1)*F2(1) + ... + P(i,h)*F2(h))*V(i) = 0
%

h=size(P,1);
n=size(A{1},1);
r=size(Psi{1},2);

diff=0;
Pi=[zeros(n-s,s); eye(s)];
for j=1:h
    X=inv([A{j}*V{j} Pi]);
    
    tmp=norm(X*B{j} - [F1{j}; F2{j}]);
    if tmp > diff
        diff=tmp;
    end
    
    tmp=norm(X*Psi{j} - [G1{j}; G2{j}]);
    if (tmp > diff)
        diff=tmp;
    end
end

for i=1:h
    X=zeros(s,n);
    for j=1:h
        X=X+P(i,j)*F2{j};
    end
    tmp=norm(X*V{i});
    if tmp > diff
        diff=tmp;
    end;
end
