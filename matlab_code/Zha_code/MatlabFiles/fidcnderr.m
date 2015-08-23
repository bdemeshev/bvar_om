function [yhat,Estr,rcon,Rcon,u,v,d] = fidcnderr(valuecon,stepcon,varcon,nconstr,...
               nstepsm,nvar,lags,phil,eq_iden,Aband,Sband,yforeml,imf3sml,Bhml,...
               yfore_h,imf3s_h,Bh_h,Cms,TLindx,TLnumber,nCms,eq_Cms)
%function [yhat,Estr,rcon,Rcon,u,v,d] = fidcnderr(valuecon,stepcon,varcon,nconstr,...
%               nstepsm,nvar,lags,phil,eq_iden,Aband,Sband,yforeml,imf3sml,Bhml,...
%               yfore_h,imf3s_h,Bh_h,Cms,TLindx,TLnumber,nCms,eq_Cms)
% To be done (8/25/98): (1) yforeml, imf3sml, and Bhml should not be there.
%                        (2) MS stuff need to be generalized.
%                        (3) perhaps, imf3s_h needs to be re-examined.
% Conditional forecasting in the identified model with one draw for error bands
%
% valuecon:  vector of values conditioned
% stepcon:   sequence (cell) of steps conditioned; if length(stepcon{i}) > 1, the condition
%               is then an arithmetic average of log(y) over the stepcon{i} period.
% varcon:    vector of variables conditioned
% nconstr:   number of constraints
% nstepsm:   maximum number of steps in all constraints
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% eq_iden: identified equation or shock (in terms of number). If equ_iden=[], then
%          'fidencond' is, similar to RATS, to compute forecasts with *all* shocks.
% Sband:  1: generate error bands from random shocks E; 0: no random shocks
% yfore:  uncondtional forecasts: forep-by-nvar
% imf3s: 3-dimensional impulse responses matrix:
%                               impsteps-by-nvar shocks-by-nvar responses
% Bh:  reduced-form parameter matrix: k-by-nvar, y(t) = X(t)*Bh+e(t)
%                    where X(t) is k-by-nvar and y(t) is 1-by-nvar
% ------
% yhat:  conditional forecasts: forep-by-nvar
% Estr:  backed-out structural shocks (from N(0,1))
% rcon:  vector - the difference between valuecon and log(yfore) (unconditional forecasts)
% Rcon:  k-by-q (q constranits and k=nvar*max(nsteps)) so that
%                        Rcon'*e = rcon where e is k-by-1
% [u,d,v]:  svd(Rcon,0)
%
%% See Zha's note "Forecast (1)" p. 5, RATS manual (some errors in RATS), etc.
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
% Copyright (c) March 1998 by Tao Zha
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



IdenShock = ~isempty(eq_iden);   % if not empty, the shock is identified

forep=size(yfore_h,1);
impsteps=size(imf3s_h,1);
if (forep<nstepsm) | (impsteps<nstepsm)
	disp('Increase # of forecast or impulse steps!!')
   disp('Or decrease # of constraints (nconstr) or constrained steps (stepcon(i))!!')
	error('Maximum of conditional steps > # of forecast or impulse steps!!')
end
kfs = nvar*nstepsm;   % k -- fs: free shocks
%*** initializing
Rcon = zeros(kfs,nconstr);   % R: k-by-q
Econ = zeros(kfs,1);      % E: k-by-1
rcon = zeros(nconstr,1);   % r: q-by-1
%rcon=valuecon-diag(yfore(stepcon,varcon));  % another way is to use "loop" below.
A0in = reshape(imf3s_h(1,:,:),nvar,nvar);  % nvar shocks-by-nvar responses


for i=1:nconstr
	if IdenShock
   	rcon(i)=length(stepcon{i})*valuecon(i) - ...
		                     sum(yforeml(stepcon{i},varcon(i)),1);  % <<>>
	else
		rcon(i)=length(stepcon{i})*valuecon(i) - ...
		                     sum(yfore_h(stepcon{i},varcon(i)),1);  %<<>>
	end
   Rmat = zeros(nstepsm,nvar);
	r2mat = zeros(nstepsm,1);   % simply one identified equation
           % Must be here inside the loop because it's matrix of one column of Rcon
   for j=1:length(stepcon{i})
		if IdenShock     % Rcon only at ML; assuming the Fed can't see all other shocks
      	Rmat(1:stepcon{i}(j),eq_iden) = Rmat(1:stepcon{i}(j),eq_iden) + ...
                 imf3sml(stepcon{i}(j):-1:1,eq_iden,varcon(i));
	                % Rmat: row--nstepsm, column--nvar shocks (here all shocks except
						 %     the identified one are set to zero) for a particular
                   %     endogenous variable 'varcon(i)'.  See Zha Forcast (1), pp.6-7
		else             % Rcon random with (A0,A+)
			Rmat(1:stepcon{i}(j),:) = Rmat(1:stepcon{i}(j),:) + ...
			        imf3s_h(stepcon{i}(j):-1:1,:,varcon(i));
	                % Rmat: row--nstepsm, column--nvar shocks (here all shocks are
						 %     *not* set to zero) for a particular endogenous
                   %     variable 'varcon(i)'.  See Zha Forcast (1), pp.6-7
		end
   end
	Rmatt = Rmat';   % Now, nvar-by-nstepsm. I think here is where RATS has an error
							  % i.e. "OVERR" is not transposed when overlaid to "CAPR"
	Rcon(:,i)=Rmatt(:);      % Rcon: k-by-q where q=nconstr
end

if nconstr
	[u d v]=svd(Rcon,0); %trial
	% rtr = Rcon'*Rcon; %trial
	% rtrinv = inv(Rcon'*Rcon); %trial
	vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
	dinv = 1./diag(d);    % inv(diag(d))
	vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
	rtr=vd*vd';       % R'*R
	rtrinv = vdinv*vdinv';   % inv(R'*R)

	Econ = Rcon*rtrinv*rcon;    % E = R*inv(R'R)*r; mean
else
	Econ = zeros(kfs,1);
	Rcon = NaN;
	rcon = NaN;
	u = NaN;
	d = NaN;
	v = NaN;
end


tcwc = nvar*lags;     % total coefficients without constant
phi=phil;
phis=phil;     % for exact backed out shocks
%
if (Sband==0)     % no random or random only from (A0,A+)
	Estr = reshape(Econ,nvar,nstepsm);
	Estr = Estr';   % transpose so that
	          % Estr: structural shocks. Row--steps, Column--n shocks
	Estr = [Estr;randn(forep-nstepsm,nvar)];
				 % Now, forep-by-nvar -- ready for forecasts
	Ures = Estr*A0in;     % nstepsm-by-nvar
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks

	% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
	% **       where phi = x(t+h-1) with last column being constant
	%
	yhat = zeros(forep,nvar);
	for k=1:forep
   	yhat(k,:) = phi*Bh_h + Ures(k,:);
		phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
		phi(1,1:nvar) = yhat(k,:);
      %
	end
else     % random from shocks E and possibly (A0,A+) depending on if imf3s is random
	if nconstr
		if IdenShock    % other shocks are indepedent of the eq_iden shock
			Osk = randn(kfs,1);    % other shocks
			for j=1:nstepsm
            Osk(nvar*(j-1)+eq_iden)=0;     % no shock to the MS or identified equation
			end
			Estr = Econ + Osk;   % Econ is non zero only at position
			                     %  eq_iden*j where j=1:nstepsm
		else
    		Ome = eye(kfs) - u*u';        % note, I-u*u' = I - R*inv(R'*R)*R'
    		%[u1 d1 v1] = svd(Ome);  % too slow
    		[u1 d1] = eig(Ome);
    		Stdcon = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of conditional variance
    							% see Zha's forecast (1), p.17
			if Cms
				if TLindx      % tight
    				Estr1 = Econ + Stdcon*randn(kfs,1);
					Estr2 = reshape(Estr1,nvar,nstepsm);
					Estr2 = Estr2';   % transpose so that
	                % Estr2: structural shocks. Row--nstepsm, Column--n shocks
					Estr = [Estr2;randn(forep-nstepsm,nvar)];
				       % Now, forep-by-nvar -- ready for forecasts
					Sindx = find(Estr(1:nCms,eq_Cms)<TLnumber);
					while ~isempty(Sindx)
		    			Estr1 = Econ + Stdcon*randn(kfs,1);
						Estr2 = reshape(Estr1,nvar,nstepsm);
						Estr2 = Estr2';   % transpose so that
	              		 % Estr2: structural shocks. Row--nstepsm, Column--n shocks
						Estr = [Estr2;randn(forep-nstepsm,nvar)];
				     		% Now, forep-by-nvar -- ready for forecasts
						Sindx = find(Estr(1:nCms,eq_Cms)<TLnumber);
					end
				else                     % loose
    				Estr1 = Econ + Stdcon*randn(kfs,1);
					Estr2 = reshape(Estr1,nvar,nstepsm);
					Estr2 = Estr2';   % transpose so that
	                % Estr2: structural shocks. Row--nstepsm, Column--n shocks
					Estr = [Estr2;randn(forep-nstepsm,nvar)];
				       % Now, forep-by-nvar -- ready for forecasts
					Sindx = find(Estr(1:nCms,eq_Cms)>TLnumber);
					while ~isempty(Sindx)
		    			Estr1 = Econ + Stdcon*randn(kfs,1);
						Estr2 = reshape(Estr1,nvar,nstepsm);
						Estr2 = Estr2';   % transpose so that
	              		 % Estr2: structural shocks. Row--nstepsm, Column--n shocks
						Estr = [Estr2;randn(forep-nstepsm,nvar)];
				     		% Now, forep-by-nvar -- ready for forecasts
						Sindx = find(Estr(1:nCms,eq_Cms)>TLnumber);
					end
				end
			else
  				Estr1 = Econ + Stdcon*randn(kfs,1);
				Estr2 = reshape(Estr1,nvar,nstepsm);
				Estr2 = Estr2';   % transpose so that
                % Estr2: structural shocks. Row--nstepsm, Column--n shocks
				Estr = [Estr2;randn(forep-nstepsm,nvar)];
			       % Now, forep-by-nvar -- ready for forecasts
			end
		end
	else
		if Cms
			if TLindx      % tight
				Estr = randn(forep,nvar);
				    % Now, forep-by-nvar -- ready for forecasts
				Sindx = find(Estr(1:nCms,eq_Cms)<TLnumber);
				while ~isempty(Sindx)
					Estr = randn(forep,nvar);
				   	% Now, forep-by-nvar -- ready for forecasts
					Sindx = find(Estr(1:nCms,eq_Cms)<TLnumber);
				end
			else                     % loose
				Estr = randn(forep,nvar);
				    % Now, forep-by-nvar -- ready for forecasts
				Sindx = find(Estr(1:nCms,eq_Cms)>TLnumber);
				while ~isempty(Sindx)
					Estr = randn(forep,nvar);
				   	% Now, forep-by-nvar -- ready for forecasts
					Sindx = find(Estr(1:nCms,eq_Cms)>TLnumber);
				end
			end
		else
			Estr = randn(forep,nvar);
			    % Now, forep-by-nvar -- ready for forecasts
		end
	end
	%
	%disp('HERE')

	Ures = Estr*A0in;     % nstepsm-by-nvar
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks

	% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
	% **       where phi = x(t+h-1) with last column being constant
	%
	yhat = zeros(forep,nvar);
	for k=1:forep
   	yhat(k,:) = phi*Bh_h + Ures(k,:);
		phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
		phi(1,1:nvar) = yhat(k,:);
	end
end