function [yhat,Estr,rcon,Rcon,u,v,d] = fn_fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
            nconstr,eq_ms,nvar,lags,phil,Aband,Sband,yfore_h,imf3s_h,A0_h,Bh_h,...
            forep,TLindx,TLnumber,nCms,eq_Cms)
% [yhat,Estr,rcon,Rcon,u,v,d] = fn_fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
%            nconstr,eq_ms,nvar,lags,phil,Aband,Sband,yfore_h,imf3s_h,A0_h,Bh_h,...
%            forep,TLindx,TLnumber,nCms,eq_Cms)
%
%   Note that the case Aband=1 is not finished yet.
%
%   Conditional forecasting in the identified model with or without error bands
%   It handles conditions on average values as well, so "valuecon" must be
%      expressed at average (NOT sum) level.
%   Aband is used only once when nconstr>0 and Aband=1, where Gibbs sampler may be used.  NOT yet finished.
%   Unconditional forecast when imf3s_h, etc is fixed and nconstr=0, where yfore_h must set to [].
%
% valuecon:  vector of values conditioned
% stepcon:   sequence (cell) of steps conditioned; if length(stepcon{i}) > 1, the condition
%               is then an arithmetic average of log(y) over the stepcon{i} period.
% varcon:    vector of variables conditioned
% nconstr:   number of DLS constraints.  If zero, it leads to unconditional forecasts.
% nstepsm:   maximum number of steps in all DLS constraints
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% Aband:  1: draws from A0 and Bh; 0: no draws.  The case Aband=1 and nconstr>0 has NOT been finished yet.
% Sband:  1: draws from random shocks E; 0: no random shocks
% yfore_h:  uncondtional forecasts: forep-by-nvar.  Never used when nconstr=0.
%            In this case, set it to [];
% imf3s_h: 3-dimensional impulse responses matrix: impsteps-by-nvar shocks-by-nvar responses
%            Never used when nconstr=0.  In this case, set it to [];
% A0_h:  A0 contemporaneous parameter matrix
% Bh_h:  reduced-form parameter matrix: k-by-nvar, y(t) = X(t)*Bh+e(t)
%                    where X(t) is k-by-nvar and y(t) is 1-by-nvar
% forep:  # of forecast periods (e.g., monthly for a monthly model)
% TLindx: 1-by-nCms vector of 1's and 0's, indicating tight or loose; 1: tighter, 0: looser
%       Used only when /* (MS draws) is activated.  Right now, MS shocks are deterministic.
% TLnumber: 1-by-nCms vector, lower bound for tight and upper bound for loose
% nCms: # of LZ conditions
% eq_Cms:  equation location of MS shocks
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
%%                 e = R*inv(R'*R)*r and k>=q
%
% Copyright (c) March 1998 by Tao Zha. Revised November 1998;
% 3/20/99 Disenabled draws of MS shcoks.  To enable it, activate /* part
% 3/20/99 Added A0_h and forep and deleted Cms as input argument.  Previous
%              programs may not be compatible.
% 3/15/2004 There are some BUG problems when calling fn_fcstcnd.m().
%
% Comparing with the version fn_fcstidcnd2.m, this is a better version.
%   See the note in fn_fcstidcnd2.m for explanation.  April, 2009.  TZ.
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

DLSIdShock = ~isempty(eq_ms);   % if not empty, the MS shock is identified as in DLS

impsteps=size(imf3s_h,1);
if (forep<nstepsm) | (impsteps<nstepsm)
	disp('Increase # of forecast or impulse steps!!')
   disp('Or decrease # of constraints (nconstr) or constrained steps (stepcon(i))!!')
	error('Maximum of conditional steps > # of forecast or impulse steps!!')
end
kts = nvar*nstepsm;   % k -- ts: total shocks some of which are restricted and others
							  %  are free.
%*** initializing
Rcon = zeros(kts,nconstr);   % R: k-by-q
Econ = zeros(kts,1);      % E: k-by-1
rcon = zeros(nconstr,1);   % r: q-by-1
%rcon=valuecon-diag(yfore(stepcon,varcon));  % another way is to use "loop" below.
tcwc = nvar*lags;     % total coefficients without constant
phi=phil;



%----------------------------------------------------
%  Form rcon, Rcon, and Econ (the mean of structural shocks)
%----------------------------------------------------
if nconstr
   A0in = reshape(imf3s_h(1,:,:),nvar,nvar);  % nvar shocks-by-nvar responses
	for i=1:nconstr
		rcon(i)=length(stepcon{i})*valuecon(i) - ...
		                     sum(yfore_h(stepcon{i},varcon(i)),1);  %<<>>
	   Rmat = zeros(nstepsm,nvar);
		r2mat = zeros(nstepsm,1);   % simply one identified equation
	         % Must be here inside the loop because it's matrix of one column of Rcon
	   for j=1:length(stepcon{i})
			if DLSIdShock   % Assuming the Fed can't see all other shocks within a month
	      	Rmat(1:stepcon{i}(j),eq_ms) = Rmat(1:stepcon{i}(j),eq_ms) + ...
	                 imf3s_h(stepcon{i}(j):-1:1,eq_ms,varcon(i));
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
		Rmatt = Rmat';  % Now, nvar-by-nstepsm. I think here is where RATS has an error
							 % i.e. "OVERR" is not transposed when overlaid to "CAPR"
		Rcon(:,i)=Rmatt(:);      % Rcon: k-by-q where q=nconstr
	end

	[u d v]=svd(Rcon,0); %trial
	               %???? Can we reduce the time by computing inv(R'*R) directly?
	% rtr = Rcon'*Rcon; %trial
	% rtrinv = inv(Rcon'*Rcon); %trial
	vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
	dinv = 1./diag(d);    % inv(diag(d))
	vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
	rtr=vd*vd';       % R'*R
	rtrinv = vdinv*vdinv';   % inv(R'*R)

	Econ = Rcon*rtrinv*rcon;    % E = R*inv(R'R)*r; the mean of structural shocks
else
	Econ = zeros(kts,1);  % the mean of shocks is zero under no variable condition
	Rcon = NaN;
	rcon = NaN;
	u = NaN;
	d = NaN;
	v = NaN;
end



%---------------------------------------
%  No uncertainty at all or only random (A0,A+)
%  In other words, no future shocks
%---------------------------------------
if (~Sband) %| (nconstr & (length(eq_ms)==1))
         % length(eq_ms)==1 implies one-one mapping between MS shocks and, say, FFR
         %  if nstepsm==nconstr.  If this condition does not hold, this procedure
         %  is incorrect.  I don't have time to fix it now (3/20/99).  So I use
         %  this as a proximation
	Estr = reshape(Econ,nvar,nstepsm);
	Estr = Estr';   % transpose so that
	          % Estr: structural shocks. Row--steps, Column--n shocks
   Estr = [Estr;zeros(forep-nstepsm,nvar)];
				 % Now, forep-by-nvar -- ready for forecasts
   Estr(1:nCms,eq_Cms) = TLnumber(:);
   Ures = Estr/A0_h;     % nstepsm-by-nvar
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

%---------------------------------------
%  With random future shocks and possibly (A0,A+) depending
%           on if imf3s_h is random
%---------------------------------------
else
	%--------------
	% Condition on variables and A random
	%--------------
	if nconstr & Aband
      warning(' ')
		disp('This situation (both E and A random) is still under construction')
      disp('It is closely related to Waggoner and Zha ReStat Gibbs sampling method')
		disp('Please press ctrl-c to abort')
		pause
	elseif nconstr
		%--------------
		% Condition on variables and DLS MS shock, no A random but S random
		%--------------
		if DLSIdShock    % other shocks are indepedent of the eq_ms shock
           % 3/20/99 The following may be problematic because Osk should depend
           %  on u (A0_h and Bh_h) in general.  I haven't worked out any good version
         %/*
         %  Osk = randn(kts,1);    % other shocks
         %  for j=1:nstepsm
         %     Osk(nvar*(j-1)+eq_ms)=0;     % no shock to the MS or identified equation
         %  end
         %  Estr = Econ + Osk;   % Econ is non zero only at position
         %                       %  eq_ms*j where j=1:nstepsm
         %  Estr = reshape(Estr,nvar,nstepsm);
         %  Estr = Estr';   % transpose so that
         %           % Estr: structural shocks. Row--steps, Column--n shocks
         %  Estr = [Estr;randn(forep-nstepsm,nvar)];
         %     % Now, forep-by-nvar -- ready for forecasts
         %
         disp('DLS')
         Ome = eye(kts) - u*u';        % note, I-u*u' = I - R*inv(R'*R)*R'
         %[u1 d1 v1] = svd(Ome);  % too slow
         [u1 d1] = eig(Ome);
         Stdcon = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of conditional variance
                        % see Zha's forecast (1), p.17
         tmp1=zeros(nvar,nstepsm);
         tmp1(eq_ms,:)=randn(1,nstepsm);
         tmp2=tmp1(:);
         %Estr1 = Econ + Stdcon*randn(kts,1);
         %jnk = reshape(Stdcon*tmp2,nvar,nstepsm)
         Estr1 = Econ + Stdcon*tmp2;
         Estr2 = reshape(Estr1,nvar,nstepsm);
         Estr2 = Estr2';   % transpose so that
            % Estr2: structural shocks. Row--nstepsm, Column--n shocks
         Estr = [Estr2;randn(forep-nstepsm,nvar)];
            % Now, forep-by-nvar -- ready for forecasts
		else
    		Ome = eye(kts) - u*u';        % note, I-u*u' = I - R*inv(R'*R)*R'
    		%[u1 d1 v1] = svd(Ome);  % too slow
    		[u1 d1] = eig(Ome);
    		Stdcon = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of conditional variance
    							% see Zha's forecast (1), p.17
			%--------------
			% Condition on variables and LZ MS shock, no A random but S random
			%   This section has not be tested yet, 10/14/98
			%--------------
         if nCms
 				Estr1 = Econ + Stdcon*randn(kts,1);
				Estr2 = reshape(Estr1,nvar,nstepsm);
				Estr2 = Estr2';   % transpose so that
                % Estr2: structural shocks. Row--nstepsm, Column--n shocks
				Estr = [Estr2;randn(forep-nstepsm,nvar)];
			       % Now, forep-by-nvar -- ready for forecasts
            Estr(1:nCms,eq_Cms) = TLnumber(:);

            %/* draw MS shocks
            %  for k=1:nCms
            %     if TLindx(k)     % tighter
            %        while (Estr(k,eq_Cms)<TLnumber(k))
            %           Estr(k,eq_Cms) = randn(1,1);
            %        end
            %     else        % looser
            %        while (Estr(k,eq_Cms)>TLnumber(k))
            %           Estr(k,eq_Cms) = randn(1,1);
            %        end
            %     end
            %  end
			%--------------
			% Condition on variables only, no A random but S random
			%--------------
			else
  				Estr1 = Econ + Stdcon*randn(kts,1);
				Estr2 = reshape(Estr1,nvar,nstepsm);
				Estr2 = Estr2';   % transpose so that
                % Estr2: structural shocks. Row--nstepsm, Column--n shocks
				Estr = [Estr2;randn(forep-nstepsm,nvar)];
			       % Now, forep-by-nvar -- ready for forecasts
			end
		end
  	%--------------
  	% Condition on LZ MS shocks only, S random and possibly A random depending on
   %                     if A0_h and Bh_h are random
  	%--------------
	else
      if nCms
			Estr = randn(forep,nvar);
				    % Now, forep-by-nvar -- ready for forecasts
         Estr(1:nCms,eq_Cms) = TLnumber(:);

         %/* draw MS shocks
         %  for k=1:nCms
         %     if TLindx(k)     % tighter
         %        while (Estr(k,eq_Cms)<TLnumber(k))
         %           Estr(k,eq_Cms) = randn(1,1);
         %        end
         %     else        % looser
         %        while (Estr(k,eq_Cms)>TLnumber(k))
         %           Estr(k,eq_Cms) = randn(1,1);
         %        end
         %     end
         %  end
		else
			Estr = randn(forep,nvar);    % Unconditional forecast
			    % Now, forep-by-nvar -- ready for forecasts
		end
	end
	%


   Ures = Estr/A0_h;     % nstepsm-by-nvar
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
