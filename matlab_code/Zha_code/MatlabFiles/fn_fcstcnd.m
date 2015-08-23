function [yhat,Estr,rcon,Rcon,u,v,d] = fn_fcstcnd(valuecon,stepcon,varcon,nstepsm,...
                      nconstr,eq_ms,nvar,lags,phil,Sband,yfore_h,imf3s_h,Bh_h,forep)
% ******* There are some BUG problems when calling this fucntion. ******* 3/15/2004.
% [yhat,Estr,rcon,Rcon,u,v,d] = fn_fcstcnd(valuecon,stepcon,varcon,nstepsm,...
%                      nconstr,eq_ms,nvar,lags,phil,Sband,yfore_h,imf3s_h,Bh_h,forep)
%
%   Conditional forecasting in the identified model with or without error bands
%   It handles conditions on average values as well, so "valuecon" must be
%      expressed at average (NOT sum) levels (i.e., arithmetic averages of log(y) over
%      the period of stepcon{i}.  5/22/01.
%   Unconditional forecast when imf3s_h, etc is fixed and nconstr=0.
%   Note that length(eq_ms)==1 implies one-one mapping between MS shocks and, say, FFR
%     if nstepsm==nconstr.  If this condition does not hold, this procedure is incorrect.
%     I don't have time to fix it now (3/20/99).  Meantime, consult or use the old code
%     fidencond.m.
%
% valuecon:  vector of values conditioned
% stepcon:   sequence (cell) of steps conditioned; if length(stepcon{i}) > 1, the condition
%               is then an arithmetic average of log(y) over the stepcon{i} period.
% varcon:    vector of variables conditioned
% nstepsm:   maximum number of steps in all DLS constraints
% nconstr:   number of DLS constraints
% eq_ms:  Equation location of MS shocks.  If [], all shocks.
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% Sband:  1: draws from random shocks E; 0: no random shocks  For now (4/27/01), no option
%    for Aband because I don't think it works best to do both Aband and Sband in one function.
% yfore_h:  uncondtional forecasts: forep-by-nvar.  Never used when nconstr=0.
%            In this case, may set it to [];
% imf3s_h: 3-dimensional impulse responses matrix: impsteps-by-nvar shocks-by-nvar responses
%            Never used when nconstr=0.  In this case, may set it to [];
% Bh_h:  reduced-form parameter matrix: k-by-nvar, y(t) = X(t)*Bh+e(t)
%                    where X(t) is k-by-nvar and y(t) is 1-by-nvar
% forep:  # of forecast periods (e.g., monthly for a monthly model)
% eq_Cms:  equation location of MS shocks
% ------
% yhat:  conditional forecasts: forep-by-nvar
% Estr:  backed-out structural shocks (from constrained Gaussians)
% rcon:  vector - the difference between valuecon and log(yfore) (unconditional forecasts)
% Rcon:  k-by-q (q constranits and k=nvar*max(nsteps)) so that
%                        Rcon'*e = rcon where e is k-by-1
% [u,d,v]:  svd(Rcon,0)
%
%% See Zha's note "Forecast (1)" p. 5, RATS manual (some errors in RATS), etc.
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
% See the old code fidencond.m.  I wond't use fn_fcstidcnd?.m, 5/22/01.
% Copyright (c) March 1998 by Tao Zha. Revised November 1998, May 2001 (Delete A0_h as
%    input arg so that previous programs may not be compatible).

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
A0in = reshape(imf3s_h(1,:,:),nvar,nvar);  % <<>>1 nvar shocks-by-nvar responses
if nconstr   % Conditional forecasts.
	for i=1:nconstr
      rcon(i)=length(stepcon{i})*valuecon(i) - sum(yfore_h(stepcon{i},varcon(i)),1);
                                   %<<>>2 Automatically taking care of average conditions.
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
   %
	[u d v]=svd(Rcon,0); %trial
	vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
	dinv = 1./diag(d);    % inv(diag(d))
	vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
	rtr=vd*vd';       % R'*R
	rtrinv = vdinv*vdinv';   % inv(R'*R)

	Econ = Rcon*rtrinv*rcon;    % E = R*inv(R'R)*r; the mean of structural shocks
else   % Unconditional forecasts.
	Econ = zeros(kts,1);  % the mean of shocks is zero under no variable condition
	Rcon = NaN;
	rcon = NaN;
	u = NaN;
	d = NaN;
	v = NaN;
end



%---------------------------------------
%  No uncertainty at all.  In other words, no future shocks.
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
   Ures = Estr*A0in;      % <<>>1 nstepsm-by-nvar
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks
          % Note: We don't use /A0_h so as to eliminate small discrepancies to be
          %   completely compatible with the computation of Rmat and Estr, which uses A0in.

	% ** reconstruct x(t) for y(t+h) = x(t+h-1)*B
	% **       where phi = x(t+h-1) with last column being constant
	%
	yhat = zeros(forep,nvar);
	for k=1:forep
   	yhat(k,:) = phi*Bh_h + Ures(k,:);
		phi(1,nvar+1:tcwc) = phi(1,1:tcwc-nvar);
		phi(1,1:nvar) = yhat(k,:);
	end
%---------------------------------------
%  With random future shocks.
%---------------------------------------
else
   if nconstr     % Conditional forecasts.
		%--------------
      % Condition on variables with all shocks backed out. Straight DLS forecasts.  No A random but S random.
		%--------------
      if ~DLSIdShock
    		Ome = eye(kts) - u*u';        % note, I-u*u' = I - R*inv(R'*R)*R'
    		%[u1 d1 v1] = svd(Ome);  % too slow
    		[u1 d1] = eig(Ome);
    		Stdcon = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of conditional variance
    							% see Zha's forecast (1), p.17
         Estr1 = Econ + Stdcon*randn(kts,1);   % Draws of constrained (conditioned) shocks.
         Estr2 = reshape(Estr1,nvar,nstepsm);
         Estr2 = Estr2';   % transpose so that
             % Estr2: structural shocks. Row--nstepsm, Column--n shocks
         Estr = [Estr2;randn(forep-nstepsm,nvar)];  % Second concatenated part: draws of free shocks.
             % Now, forep-by-nvar -- ready for forecasts
		%--------------
      % Condition on variables with identified MS shocks backed out, no A random but S random.
		%--------------
      else     % other shocks are indepedent of the eq_ms shock
           % 3/20/99 The following may be problematic because Osk should depend
           %  on u (A0_h and Bh_h) in general.  I have NOT worked out any good version.
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
      end
   else   % Unconditional forecasts.
     Estr = randn(forep,nvar);    % Unconditional draws.
	end

   Ures = Estr*A0in;     % <<>>1 nstepsm-by-nvar
			 % Ures: reduced-form residuals.  Row--steps; Column--n shocks
          % Note: We don't use /A0_h so as to eliminate small discrepancies to be
          %   completely compatible with the computation of Rmat and Estr, which uses A0in.

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
