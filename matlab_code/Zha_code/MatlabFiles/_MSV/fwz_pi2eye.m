function [A_new, B_new, Psi_new] = fwz_pi2eye(A, B, Psi, Pi)
% [A, B, Psi] = fwz_pi2eye(A, B, Psi, Pi)
% Converts the representation 
%
%   A(s(t) x(t) = B(s(t) x(t-1) + Psi(s(t)) epsilon(t) + Pi{s(t)} eta(t)
%
% to a form compatible with fwz_msv_msre().  This is the form
%
%   A(s(t) x(t) = B(s(t) x(t-1) + Psi(s(t)) epsilon(t) + Pi eta(t)
%
% where Pi' = [zeros(s,n-s)  eye(s)]
%

for i=1:h
    [Q,R] = qr(Pi{i});
    W=[zeros(n-s,s)  I; inv(R(1:s,1:s)); zeros(s,n-s)]*Q';
    A_new{i}=W*A{i};
    B_new{i}=W*B{i};
    Psi_new{i}=W*Psi{i};
end