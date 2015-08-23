function [grdd, badg] = numgradcd(fcn,x0,varargin)
%function grdd = fn_gradcd(fcn,x0,varargin)
% computes numerical gradient of a single-valued function or Jacobian
%   matrix of a vector-valued function using a central difference with
%                    function grdd = gradcd(fcn,x0,grdh)
%
%   fcn: a string naming a vector-valued function (f:n-by-1 -> k-by-1).
%   x0: a column vector n-by-1, at which point the hessian is evaluated.
%   grdh: step size, n*1. Set as follows
%              step = eps^(1/3);
%              %step = 1e-04;
%              grdh = step * (max([abs(x) ones(length(x),1)]'))' .* (abs(x) ./ x);
%--------------------
%   grdd: n-by-k Jacobian matrix (gradients).
%
% Written by Tao Zha, 2002.
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
%stps = 1.0e-02;
% eps: floating point relative accuracy or machine precision: 2.22e-16
% stps: step size recommended by Dennis and Schnabel: 6.006e-6

x0 = x0(:);
%tailstr = ')';
%for i=nargin-3:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%end
if ischar(fcn)
    f0 = eval([fcn '(x0,varargin{:})']);
else
    f0 = fcn(x0,varargin{:});
end
%f0 = eval([fcn '(x0' tailstr]);

% ** initializations
n = length(x0);
k = length(f0);   % dimension of "fcn"

% ** Computation of stepsize (dh)
if (0)
    dh = 1.0e-06; %grdh;
else
    ax0 = abs(x0);
    if all(x0)
        dax0 = x0 ./ ax0;
    else
        dax0 = 1;
    end
    dh = stps * (max([ax0 ones(n,1)]'))' .* dax0;
end

xdh = x0 + dh;
dh = xdh - x0;    % This increases precision slightly
%
argplus = x0(:,ones(n,1));
argminus = argplus;
dnum = 1:n+1:n^2;    % positions of diagonals in vec(argplus).
argplus(dnum) = xdh;     % replace the diagonals of "argplus" by "xdh".
argminus(dnum) = x0-dh;    % replace the diagonals of "argplus" by "xdh".

grdd = zeros(k,n);   % preallocate to speed the loop.
badg=0;
i = 0;
while i ~= n
   i = i+1;
   if ischar(fcn)
       
    fp = eval([fcn '(argplus(:,i),varargin{:})']);
    fm = eval([fcn '(argminus(:,i),varargin{:})']);
   else
    
    fp = fcn(argplus(:,i),varargin{:});
    fm = fcn(argminus(:,i),varargin{:});   
   end
%    fp = eval([fcn '(argplus(:,i)' tailstr]);
%    fm = eval([fcn '(argminus(:,i)' tailstr]);
   g0  = fp - fm;

   if abs(g0)< 1e15
      grdd(:,i)=g0;
      % disp('good gradient')
   else
      disp('bad gradient ------------------------')
      % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see above
      grdd(:,i)=0;
      badg=1;
      % return
      % can return here to save time if the gradient will never be
      % used when badg returns as true.
   end

end
dhm = dh(:,ones(k,1));
dhm = dhm';      % k*n
grdd = grdd ./ (2*dhm);   % k-by-n.
grdd = grdd';  % n-by-k.

