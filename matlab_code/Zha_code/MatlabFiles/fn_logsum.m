function [log_sum, log_max] = fn_logsum(log_sum, log_max, log_new)
%  [log_sum, log_max] = fn_logsum(log_sum, log_max, log_new)
%
%Outputs:
% log_sum:  updated log(sum of x_1, ..., x_{N+1})
% log_max:  updated max of log(x_1), ..., log(x_{N+1}).
%--------------
%Inputs:
% log_sum:  log(sum of x_1, ..., x_N)
% log_max:  max of log(x_1), ..., log(x_N).
% log_new:  log(x_{N+1}).
%
%Written by T. Zha; 12:20PM 06/28/2005
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


%=== Updates log_sumpdf with an additional pdf value.  See TVBVAR Notes p.81a.
if (log_max>=log_new)
   log_sum = log(  exp(log_sum-log_max) + exp(log_new-log_max)  ) + log_max;
else
   log_sum = log(  exp(log_sum-log_new) + 1.0  ) + log_new;
   log_max = log_new;
end
