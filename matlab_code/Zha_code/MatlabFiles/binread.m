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
nbuffer = 10;
ndraws=3*nbuffer;

n=3;
yhatw = zeros(n^2,ndraws);
seldraw = [3 7];

[fid,message]=fopen('outyhat.bin')
%[xd,count]=fread(fid,[n^2,length(seldraw)],'double');   % (1) working
[xd,count]=fread(fid,inf,'single');    % (2) working.  with 'single' or 'double'
fid
message
count
%xdd=reshape(xd,n^2,length(seldraw))   % (1) working
xdd=reshape(xd,n^2,ndraws)  % (2) working