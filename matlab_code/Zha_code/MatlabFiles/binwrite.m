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
yhatw = zeros(n^2,nbuffer);

wdraws=0;
for draws=1:ndraws
	tmp = draws*ones(n);
	yhatw(:,draws-wdraws) = tmp(:);

	if ~mod(draws,nbuffer)
      wdraws=draws

   	fwriteid = fopen('outyhat.bin','a');
      count = fwrite(fwriteid,yhatw,'single');   % 'single' or 'double'
   	status = fclose('all');
	end
end

