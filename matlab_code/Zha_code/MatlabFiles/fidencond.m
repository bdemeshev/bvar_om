function [yhat,Estr,rcon,Ome,Rmat,u,v,d] = fidencond(valuecon,stepcon,varcon,nconstr,...
                         nstepsm,nvar,lags,yfore,imf3s,phil,Bh,eq_iden,steps_iden)
% Estimating conditional forecasting in the identified model
%function [yhat,Estr,rcon,Ome,Rmat,u,v,d] = fidencond(valuecon,stepcon,varcon,nconstr,...
%                         nstepsm,nvar,lags,yfore,imf3s,phil,Bh,eq_iden,steps_iden)
%
% valuecon:  vector of values conditioned
% stepcon:   sequence (cell) of steps conditioned; if length(stepcon{i}) > 1, the condition
%               is then an arithmetic average of log(y) over the stepcon{i} period.
% varcon:    vector of variables conditioned
% nconstr:   number of constraints
% nstepsm:   maximum number of steps in all constraints
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% yfore:  uncondtional forecasts: forep-by-nvar
% imf3s: 3-dimensional impulse responses matrix:
%                               impsteps-by-nvar shocks-by-nvar responses
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% Bh:  reduced-form parameter matrix: k-by-nvar, y(t) = X(t)*Bh+e(t)
%                    where X(t) is k-by-nvar and y(t) is 1-by-nvar
% eq_iden: identified equation or shock (in terms of number).
%          If eq_iden=[], then'fidencond' is, similar to RATS, to compute
%             forecasts with *all* shocks.
%          If length(eq_iden)=1, compute forecasts with only "MS" shocks.
%          If eq_iden = [a b c], a is # of constraints for all shocks (<=nconstr),
%             b is location of MS, and c is max steps for all shocks.
%             My recall is that this option has never been tested 9/18/98
% steps_iden:  a vector (set) of steps for identified shocks.  Valid only
%              if length(eq_iden)=1. Note, length(steps_iden) must nconstr.
% ------
% yhat:  conditional forecasts: forep-by-nvar
% Estr:  backed-out structural shocks (from N(0,1))
% rcon:  vector - the difference between valuecon and log(yfore) (unconditional forecasts)
% Rcon:  k-by-q (q constranits and k=nvar*max(nsteps)) so that
%                        Rcon'*e = rcon where e is k-by-1, where k=nvar*nstepm
% Rmat:  nstepsm-by-nvar (shocks), for only one constraint at a time. See Zha's
%                  Forecast (1), pp.5-6
% Ome:  k-by-k: covariance matrix of vectorized structural shocks vec(Estr)
%
% [u,d,v]:  svd(Rcon,0)
%
%% See Zha's note "Forecast (1)" pp.5-7, RATS manual (some errors in RATS), etc.
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
% Copyright (c) February 1998 by Tao Zha
% NOTE: the code needs to be improved, 10/19/98 (use lzpaper/fcstidcnd.m for the time being).
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


%IdenShock = ~isempty(eq_iden);   % if not empty, the shock is identified
forep=size(yfore,1);
impsteps=size(imf3s,1);
if (forep<nstepsm) | (impsteps<nstepsm)
	disp('Increase # of forecast or impulse steps!!')
   disp('Or decrease # of constraints (nconstr) or constrained steps (stepcon(i))!!')
	error('Maximum of conditional steps > # of forecast or impulse steps!!')
	%warning
	%return
end
if max(steps_iden) < nconstr
	disp('Increase # of identified steps or decrease # of constraints')
	error('length(steps_iden) > nconstr')
end

kfs = nvar*nstepsm;   % k -- fs: free shocks
%*** initializing
Rcon = zeros(kfs,nconstr);   % R: k-by-q
Econ = zeros(kfs,1);      % E: k-by-1
rcon = zeros(nconstr,1);   % r: q-by-1
%rcon=valuecon-diag(yfore(stepcon,varcon));  % another way is to use "loop" below.
A0in = reshape(imf3s(1,:,:),nvar,nvar);  % nvar shocks-by-nvar responses


for i=1:nconstr
   rcon(i)=valuecon(i)-mean(yfore(stepcon{i},varcon(i)));
	%rcon(i)=valuecon(i)-sum(yfore(stepcon{i},varcon(i)));
   Rmat = zeros(nstepsm,nvar);
           % Rmat: row--nstepsm, column--nvar shocks (here all shocks except
	        %     the identified one are set to zero) for a particular
           %     endogenous variable 'varcon(i)'.  See Zha Forcast (1), pp.6-7
	r2mat = zeros(nstepsm,1);   % simply one identified equation
           % Must be here inside the loop because it's matrix of one column of Rcon
   for j=1:length(stepcon{i})
      if (length(eq_iden)==1)
			r2mat(1:stepcon{i}(j)) = r2mat(1:stepcon{i}(j)) + ...
                                      imf3s(stepcon{i}(j):-1:1,eq_iden,varcon(i));
      elseif (length(eq_iden)>1)
         if (i<=eq_iden(1))
            Rmat(1:stepcon{i}(j),:) = Rmat(1:stepcon{i}(j),:) + ...
			                             imf3s(stepcon{i}(j):-1:1,:,varcon(i));
         else
            Rmat(1:eq_iden(3),:) = Rmat(1:eq_iden(3),:) + ...
                     imf3s(stepcon{i}(j):-1:stepcon{i}(j)-eq_iden(3)+1,:,varcon(i));

            Rmat(eq_iden(3)+1:stepcon{i}(j),eq_iden(2)) = ...
                         Rmat(eq_iden(3)+1:stepcon{i}(j),eq_iden(2)) + ...
                         imf3s(stepcon{i}(j)-eq_iden(3):-1:1,eq_iden(2),varcon(i));
         end
		else
			Rmat(1:stepcon{i}(j),:) = Rmat(1:stepcon{i}(j),:) + ...
			                             imf3s(stepcon{i}(j):-1:1,:,varcon(i));
	                % Rmat: row--nstepsm, column--nvar shocks (here all shocks are
						 %     *not* set to zero) for a particular endogenous
                   %     variable 'varcon(i)'.  See Zha Forcast (1), pp.6-7
		end
   end
	%
	if (length(eq_iden)==1)
      Rmat(steps_iden,eq_iden) = r2mat(steps_iden);   % for only one constraint at a time
	end
	Rmat=Rmat/length(stepcon{i});    % <<>>
	Rmatt = Rmat';   % Now, nvar-by-nstepsm. I think here is where RATS has an error
							  % i.e. "OVERR" is not transposed when overlaid to "CAPR"
	Rcon(:,i)=Rmatt(:);      % Rcon: k-by-q where q=nconstr
end

if nconstr
   [u d v]=svd(Rcon,0); %trial.   Rcon: k-by-q; u: k-by-q
	% rtr = Rcon'*Rcon; %trial
	% rtrinv = inv(Rcon'*Rcon); %trial
	vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
	dinv = 1./diag(d);    % inv(diag(d))
	vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
	rtr=vd*vd';       % R'*R
	rtrinv = vdinv*vdinv';   % inv(R'*R)
   Ome = eye(kfs) - u*u';       % note: I-u*u' = I - R*inv(R'*R)*R'; k-by-k

	Econ = Rcon*rtrinv*rcon;    % E = R*inv(R'R)*r; mean
else
	Econ = zeros(kfs,1);
	Rcon = NaN;
	rcon = NaN;
	u = NaN;
	d = NaN;
	v = NaN;
end


Estr = reshape(Econ,nvar,nstepsm);
Estr = Estr';   % transpose so that
          % Estr: structural shocks. Row--steps, Column--n shocks
Ures = Estr*A0in;     % nstepsm-by-nvar
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks
%Ures = Estr*A0in;

% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
% **       where phi = x(t+h-1) with last column being constant
%
tcwc = nvar*lags;     % total coefficients without constant
phi=phil;
%
yhat = zeros(forep,nvar);
for k=1:forep
	if (k<=nstepsm)
   	epsl = Ures(k,:);
   	yhat(k,:) = phi*Bh + epsl;
		%yhat(k,:) = phi*Bh;
   	phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
   	phi(1,1:nvar) = yhat(k,:);
	else
   	yhat(k,:) = phi*Bh;
   	phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
   	phi(1,1:nvar) = yhat(k,:);
	end
end
