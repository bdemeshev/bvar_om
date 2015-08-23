function Q = adapt(f,a,b,tol,trace,varargin)
%ADAPT  Numerically evaluate integral using adaptive Simpson rule.
%
%   Q = ADAPT('F',A,B) approximates the integral of F(X) from A to B
%   to machine precision.  'F' is a string containing the name of the
%   function.  Function F must return a vector of output values if given
%   a vector of input values.
%
%   Q = ADAPT('F',A,B,TOL) integrates to a relative error of TOL.
%
%   Q = ADAPT('F',A,B,TOL,TRACE) displays the left end point of the
%   current interval, the interval length and the partial integral.
%
%   Q = ADAPT('F',A,B,TOL,TRACE,P1,P2,...) allows coefficients P1, ...
%   to be passed directly to function F:  G = F(X,P1,P2,...).
%   To use default values for TOL or TRACE, you may pass in the empty
%   matrix ([]).
%
%   See also QUAD, QUAD8, DBLQUAD.

%  % Copyright (C) 1997-2012 Tao Zha
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

if (nargin < 4), tol = []; end;
if (nargin < 5), trace = []; end;
if (isempty(tol)), tol = 10*eps; end;
if (isempty(trace)), trace = 0; end;

x = [a (a+b)/2 b];
y = feval(f, x, varargin{:});
fa = y(1); fm = y(2); fb = y(3);
is = (b - a)/6 * (fa + 4*fm + fb);
s = sign(is); if (s == 0), s = 1; end;
is = s*(abs(is) + b - a)/2*tol/eps;
Q = adaptstp(f, a, b, fa, fm, fb, is, trace, varargin{:});