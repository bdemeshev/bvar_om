function [yhat,Estr,rcon,Rcon,u,v,d] = fcstidcnd2(valuecon,stepcon,varcon,nstepsm,...
            nconstr,eq_ms,nvar,lags,phil,Aband,Sband,yfore_h,imf3s_h,A0_h,Bh_h,...
            forep,TLindx,TLnumber,nCms,eq_Cms,nconadd,eq_oth,radd,Radd)
% [yhat,Estr,rcon,Rcon,u,v,d] = fcstidcnd2(valuecon,stepcon,varcon,nstepsm,...
%            nconstr,eq_ms,nvar,lags,phil,Aband,Sband,yfore_h,imf3s_h,A0_h,Bh_h,...
%            forep,TLindx,TLnumber,nCms,eq_Cms,nconadd,eq_oth,radd,Radd)
%
% 5/2/99: This one is much, much slower than fidcnderr.m when eq_ms=[] and
%    nconstr>0 and length(eq_ms)*nstepsm>nconstr.  Seems to me that fcstidcnd2.m
%    is designed for the situation which eq_ms=2.  The slow speed is due to
%    "Radd" added.  When I have time, I need to check to see if I can speed up
%    this M file.
%
%   Conditional forecasting in the identified model with or without error bands
%   It handles conditions on average values as well, so "valuecon" must be
%      expressed at average (NOT sum) level.  Including unconditional forecasts
%      when nconstr=nCms=0.  Maybe slower than "fcstidcnd.m" when length(eq_ms)*nstepsm=nconstr
%      because of Radd added.  But "fsctidcnd.m" is incorrect when length(eq_ms)
%     *nstepsm>nconstr (situations like the average annual FFR constraint).  3/22/99
%   Ref.:  Zha Forecast (I), pp. 17-17c.
%
% valuecon:  vector of values conditioned
% stepcon:   sequence (cell) of steps conditioned; if length(stepcon{i}) > 1, the condition
%               is then an arithmetic average of log(y) over the stepcon{i} period.
% varcon:    vector of variables conditioned
% nconstr:   number of DLS constraints
% nstepsm:   maximum number of steps in all DLS constraints
% nvar:   number of variables in the BVAR model
% lags:   number of lags in the BVAR model
% phil:  the 1-by-(nvar*lags+1) data matrix where k=nvar*lags+1
%                 (last period plus lags before the beginning of forecast)
% eq_ms: identified equation or shock (in terms of number). If equ_ms=[], then "fidencond" or
%          'fcstidcond' is, similar to RATS, to compute forecasts with *all* shocks.
% Aband:  1: draws from A0 and Bh; 0: no draws
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
% nconadd:  number of DLS added constraints DIRECTLY on structural shocks
%            for identified version where eq_ms is not empty.  Note
%            nconadd=0 when eq_ms is empty.
% eq_oth:  index for other structural shocks than those in eq_ms.
%            If eq_ms=[], eq_oth=[].  Note length(eq_oth)+length(eq_ms)=nvar
% radd:  sparse nconadd-by-1;  nconadd values added to rcon in the text later
% Radd:  sparce nvar*nstepsm-by-nconadd; added to Rcon in the text later
% ------
% yhat:  conditional forecasts: forep-by-nvar
% Estr:  backed-out structural shocks (from N(0,1))
% rcon:  vector - the difference between valuecon and log(yfore) (unconditional forecasts)
% Rcon:  k-by-q (q constranits and k=nvar*max(nsteps)) so that
%                        Rcon'*e = rcon where e is k-by-1
% [u,d,v]:  svd(Rcon,0)
%
% Copyright (c) March 1998 by Tao Zha. Revised November 1998;
% 3/20/99 Disenabled draws of MS shcoks.  To enable it, activate /* part
% 3/20/99 Added A0_h and forep and deleted Cms as input argument.
% 3/21/99 Added nconadd and eq_oth.
%   Previous programs may not be compatible.
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



%
%% See Zha's note "Forecast (1)" p. 5, RATS manual (some errors in RATS), etc.
%% See also Zha Forecast (I), pp. 17-17c.
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
Rcon = zeros(kts,nconstr);   % R: k-by-q where q include possible additional
               % constraints directly on structural shocks for identified models.
               % These addtional constraints will be added later.
Econ = sparse(zeros(kts,1));      % E: k-by-1
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
	         % Must be here inside the loop because it's matrix of one column of Rcon
	   for j=1:length(stepcon{i})
			if DLSIdShock   % Assuming the Fed can't see all other shocks within a month
            Rmat(1:stepcon{i}(j),:) = Rmat(1:stepcon{i}(j),:) + ...
                    imf3s_h(stepcon{i}(j):-1:1,:,varcon(i));
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

   rcon = [rcon;radd];   % added nconadd constrained values for identified model
               % sparse because radd is sparse
   Rcon = [Rcon Radd];   % added constraints on shocks for identified version
               % sparse because Radd is sparse
   Rcont=Rcon';
   rinvrtr=Rcon/(Rcont*Rcon);
   Econ = rinvrtr*rcon;


   %/*  Too slow, I believe, when q is large.  3/21/99
   %  [u d v]=svd(Rcon,0); %trial
   %                 %???? Can we reduce the time by computing inv(R'*R) directly?
   %  % rtr = Rcon'*Rcon; %trial
   %  % rtrinv = inv(Rcon'*Rcon); %trial
   %  vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
   %  dinv = 1./diag(d);    % inv(diag(d))
   %  vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
   %  rtr=vd*vd';       % R'*R
   %  rtrinv = vdinv*vdinv';   % inv(R'*R)
   %
   %  Econ = Rcon*rtrinv*rcon;    % E = R*inv(R'R)*r; the mean of structural shocks
else
   Econ = sparse(zeros(kts,1));  % the mean of shocks is zero under no variable condition
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
if (~Sband) | (nconstr & length(eq_ms)*nstepsm==nconstr)
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
      disp('and related to "if DLSIdShock" in the following')
		disp('Please press ctrl-c to abort')
		pause
	elseif nconstr
		%--------------
		% Condition on variables and DLS MS shock, no A random but S random
		%--------------
		if DLSIdShock    % other shocks are indepedent of the eq_ms shock
         Ome=speye(kts)-rinvrtr*Rcont;
         Ome(find(abs(Ome)<1e-10))=0;  % tighter the tolerance is, the longer it
               % takes for svds to compute.  E.g., size(Ome,1)*eps is perhaps very
               % tight.  After all, when the number in Ome is very small relative to
               % 1 (which is the variance of structural shocks Estr), we can treat it
               % as zero.
         [u1 d1 v1] = svds(Ome,kts-nconstr-nconadd);
                 % We have exactly nconstr+nconadd zero singular values
         Stdcon = u1*diag(sqrt(diag(abs(d1))));  % find a square root of Ome
                 % see Zha's forecast (1), p.17
         Stdcon = [Stdcon sparse(zeros(kts,kts-size(d1,1)))];
         Estr1 = Econ + Stdcon*randn(kts,1);  % We must have rand(kts,1).  Assigning
               % some of randn(kts,1) to be zero according to Radd is WRONG.
               % For this argument, see Zha Forecast (I) p.17a
         Estr2 = reshape(Estr1,nvar,nstepsm);
                 % nvar-by-nstepsm;  Needed to be transposed later
         Estr = [Estr2';randn(forep-nstepsm,nvar)];  % transpose so that
            % Estr2: structural shocks. Row--nstepsm, Column--n shocks
            % Now, forep-by-nvar -- ready for forecasts
		else
         Ome=speye(kts)-rinvrtr*Rcont;
         Ome(find(abs(Ome)<1e-10))=0;  % tighter the tolerance is, the longer it
               % takes for svds to compute.  E.g., size(Ome,1)*eps is perhaps very
               % tight.  After all, when the number in Ome is very small relative to
               % 1 (which is the variance of structural shocks Estr), we can treat it
               % as zero.
         [u1 d1 v1] = svds(Ome,kts-nconstr-nconadd);
                 % We have exactly nconstr+nconadd zero singular values
         Stdcon = u1*diag(sqrt(diag(abs(d1))));  % find a square root of Ome
                 % see Zha's forecast (1), p.17
         Stdcon = [Stdcon sparse(zeros(kts,nconstr+nconadd))];

         %*** The following very inefficient
         %  Ome = eye(kts) - u*u';        % note, I-u*u' = I - R*inv(R'*R)*R'
         %  %[u1 d1 v1] = svd(Ome);  % too slow
         %  [u1 d1] = eig(Ome);
         %  Stdcon = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of conditional variance
         %                 % see Zha's forecast (1), p.17

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
         Estr(:,eq_Cms)=0;
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