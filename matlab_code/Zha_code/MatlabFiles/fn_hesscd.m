function grdd = fn_hesscd(fcn,x0,grdh,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
                   P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20)
% Computing numerical hessian using a central difference with
%                    function grdd = hesscd(fcn,x0,grdh,Passed variables1)
%
%   fcn: a string naming the objective function.
%   x0: a column vector n*1, at which point the hessian is evaluated.
%   grdh: step size.
%   grdd: hessian matrix (second derivative), n*n.
%
% Written by Tao Zha, May 1998. Revised June 1002.
% Copyright (C) 1997-2012 Tao Zha
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


stps = eps^(1/3);
% eps: floating point relative accuracy or machine precision: 2.22e-16
% stps: step size recommended by Dennis and Schnabel: 6.006e-6

x0 = x0(:);
tailstr = ')';
for i=nargin-3:-1:1
   tailstr=[ ',P' num2str(i)  tailstr];
end
f0 = eval([fcn '(x0' tailstr]);

% ** initializations
k = length(x0);
grdd = zeros(k);

% ** Computation of stepsize (dh)
if all(grdh)
    dh = grdh;
else
    ax0 = abs(x0);
    if all(x0)
        dax0 = x0 ./ ax0;
    else
        dax0 = 1;
    end
    dh = stps * (max([ax0 (1e-2)*ones(k,1)]'))' .* dax0;
end

xdh = x0 + dh;
dh = xdh - x0;    % This increases precision slightly
dhm = dh(:,ones(k,1));
ee = eye(k) .* dhm;

i = 1;
while i <= k
    j = i;
    while j <= k

		  fune1 = eval([fcn '(x0 + ee(:,i) + ee(:,j)' tailstr]);
		  fune2 = eval([fcn '(x0 - ee(:,i) + ee(:,j)' tailstr]);
		  fune3 = eval([fcn '(x0 + ee(:,i) - ee(:,j)' tailstr]);
		  fune4 = eval([fcn '(x0 - ee(:,i) - ee(:,j)' tailstr]);
        grdd(i,j) = (fune1 - fune2 - fune3 + fune4)  / (4 * dh(i) * dh(j));

        if i ~= j
            grdd(j,i) = grdd(i,j);
        end

    j = j+1;
    end
    i = i+1;
end

