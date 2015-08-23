function Estrexa = fidcndexa(yexa,phil,A0_h,Bh_h,nvar,lags)
% Estrexa = fidcndexa(yexa,phil,A0_h,Bh_h,nvar,lags)
%   Exact structural shocks e(t) (backed out by conditioning on the path "yexa")
%
% yexa:  actup-by-nvar.  Actual data (all log except R, U, etc.) begin
%         at nSample-actup+1 and ends at nSample, where nSample is the period
%         for actual data (excluding dummy observations)
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags, prior to nSample-actup+1)
% A0_h:   A0, column-equation; in the form of y(t)A0=X(t)A+ + e(t)
% Bh_h:  reduced-form parameter matrix: k-by-nvar, y(t) = X(t)*Bh+u(t)
%                    where X(t) is k-by-nvar and y(t) is 1-by-nvar
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% ------
% Estrexa:  actup-by-nvar -- backed out exact shocks given parameter values
%           where actup = length(yexa(:,1)).
%
%% See Zha's note "Forecast (2)" p.11
%
%% Some notations:  y(t+1) = y(t)B1 + e(t+1)inv(A0). e(t+1) is 1-by-n.
%%    Let r(t+1)=e(t+1)inv(A0) + e(t+2)C + .... where inv(A0) is impulse
%%          response at t=1, C at t=2, etc. The row of inv(A0) or C is
%%          all responses to one shock.
%%    Let r be q-by-1 (such as r(1) = r(t+1)
%%                 = y(t+1) (constrained) - y(t+1) (forecast)).
%%    Use impulse responses to find out R (k-by-q) where k=nvar*nsteps
%%        where nsteps the largest constrained step.  The key of the program
%%        is to creat R using impulse responses
%%    Optimal solution for shock e where R'*e=r and e is k-by-1 is
%%                 e = R*inv(R'*R)*r.
%
% Copyright (c) April 1998 by Tao Zha
% Change the input arguments so that old programs may not be compatible, 03/19/98.
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

tcwc = nvar*lags;     % total coefficients without constant
actup = length(yexa(:,1));   % periods considered for structural shocks.
phis=phil;     % for exact backed out shocks

Estrsexa=zeros(actup,nvar);    % backed out exact shocks

for k=1:actup
	Estrexa(k,:) = (yexa(k,:)-phis*Bh_h) * A0_h;
	phis(nvar+1:tcwc) = phis(1:tcwc-nvar);
	phis(1:nvar) = yexa(k,:);
end
