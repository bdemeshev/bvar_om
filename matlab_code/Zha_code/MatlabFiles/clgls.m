function [bhat,Vq,ACt,uqhat,q2] = clgls(a,yq,Xq,C,qm,mT,qT)
%function [bhat,Vq,ACt,uqhat,q2] = clgls(a,yq,Xq,C,qm,mT,qT)
%   bhat = inv[Xq' * inv(Vq) * Xq] * Xq' * inv(Vq) * yq
%   Vq = C*A*C'
%   ACt = A*C'
%   uqhat = yq - Xq * bhat
%
% Written by E.M. Leeper
% Modified by T. Zha, 5/6/97
% Copyright (C) 1997-2012 Eric Leeper and Tao Zha
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


%*** The following creation of Ahat may be inefficient, T.Z., 5/7/97
Ahat = eye(mT,mT);
i = 1;
arow = zeros(1,mT-1);
while i <= mT-1
 arow(1,i) = a^i;
 i = i + 1;
end
i = 1;
while i <= mT-1
 Ahat(i,i+1:mT) = arow(1,1:mT-i);
 i = i + 1;
end
Ahat = Ahat + Ahat' - eye(mT,mT);



%*** GLS to estimate bhat
ACt = Ahat*C';
Vq = C*ACt;
Xqt = Xq';
%CACinv = inv(CAhatC); commented out by T.Z., not necessary
%XqCACinv = Xq' * CACinv; commented out by T.Z., not necessary
bhat = ((Xqt/Vq)*Xq)\(Xqt*(Vq\yq));
uqhat = yq - Xq * bhat;

%*** compute the new quarterly coefficeint "q2"
uqlag = zeros(qT-1,1);
uqlag = uqhat(1:qT-1,1);
uqc = uqhat(2:qT,1);
[u d v]=svd(uqlag,0); %trial
dinv = 1./diag(d);    % inv(diag(d))
vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
xtxinv = vdinv*vdinv';
uy = u'*uqc; %trial
xty = vd*uy; %trial
%Bh = xtx\xty;        %inv(X'X)*(X'Y), k*m (ncoe*nvar).
q2 = xtxinv*xty;      % initial q