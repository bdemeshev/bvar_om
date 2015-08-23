function [Avhx,AvhxD,hAvhx,hAvhxD,cJump,hAvhh] = smtplis2seq(Avhx,AvhxD,hAvhx,hAvhxD,...
                      A0xhat,tdf,cJump,scf,H_sr,nfp,Sbd,fss,nvar,a0indx,hAvhh)
% [Avhx,AvhxD,hAvhx,hAvhxD,cJump,hAvhh] = smtplis2(Avhx,AvhxD,hAvhx,hAvhxD,...
%                      A0xhat,tdf,cJump,scf,H_sr,nfp,Sbd,fss,nvar,a0indx,hAvhh)
%      Straight Metropolis Algorithm for identified VARs with 2 parallel sequences
%
% Avhx: previous draw of parameter x in 1st (kept) sequence
% AvhxD: previous draw of parameter x in 2nd (discarded) sequence
% hAvhx: lh value of previous draw in 1st (kept) sequence
% hAvhxD: lh vlaue of previous draw in 2nd (discarded) sequence
% A0xhat: ML estimate of A0
% tdf: degrees of freedom of the jump t-distribution
% cJump: old count for times of jumping
% scf: scale down factor for stand. dev. -- square root of covariance matrix H
% H_sr: square root of covariance matrix H
% nfp: number of free parameters
% Sbd: S in block diagonal covariance matrix in asymmetric prior case
% fss: effective sample size
% nvar: number of variables
% a0indx: index of locations of free elements in A0
% hAvhh: highest point of lh
%--------------
% Avhx: new draw of parameter x in 1st (kept) sequence
% AvhxD: new draw of parameter x in 2nd (discarded) sequence
% hAvhx: lh value of new draw in 1st (kept) sequence
% AvhxD: lh vlaue of new draw in 2nd (discarded) sequence
% cJump: new count for times of jumping
% hAvhh: highest point of lh
%
% November 1998 by Tao Zha
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

if ~(fix(tdf)-tdf==0)
	warning('tdf in msstart.m must be integer for drawing chi^2 from normal')
	disp('press ctrl-c to abort')
	pause
end

for k=1:2    % parallel sequences; 2nd to be discarded
	%** draw free elements Avh in A0 and hyperparameters from t-dist
	%gsg = gamrnd(ga,gb);     % G-sigma from Gamma draw
	%Avhz = (1/sqrt(gsg))*(scf*H_sr*randn(nfp,1));
   	% scf is used to control the acceptance ratio -- countJump
	Avhz1 = scf*H_sr*randn(nfp,1);     % normal draws
	%Avhz1 = 2*randn(nfp,1);     % normal draws
	csq=randn(tdf,1);
	csq=sum(csq .* csq);
	Avhz = Avhz1/sqrt(csq/tdf);

	if (k==1)
		Avhy = Avhx + Avhz;      % random walk chain -- Metropolis
	else
		Avhy = AvhxD + Avhz;     % random walk chain -- Metropolis
	end

   % ** actual density, having taken log
   hAvhy = a0asfun(Avhy,Sbd,fss,nvar,a0indx);
   hAvhy = -hAvhy;      % converted to logLH

	if (k==1)
   	mphxy = exp(hAvhy-hAvhx);
	else
		mphxy = exp(hAvhy-hAvhxD);
	end

	%** draw u from uniform(0,1)
	u = rand(1);
	Jump = min([mphxy 1]);
	if u <= Jump
      %** Normalization: get the point that is nearest to axhat or A0xhat
		Atem=zeros(nvar);
	   Atem(a0indx) = Avhy;
		Adiff = (Atem - A0xhat).^2;
                      % distance that may be far from axhat or A0xhat
      Adiffn = (-Atem - A0xhat).^2;
		                      % distance by chaning the sign of Atem
		cAdiff = sum(Adiff);    % each column summed up
		cAdiffn = sum(Adiffn);  % each column summed up
		cAindx = find(cAdiffn<cAdiff); % index for shorter distance
		Atemn = Atem;
		Atemn(:,cAindx) = -Atem(:,cAindx);
                      % find the shortest or nearest distance
		%** get the value of logPoster or logLH
		Avhy = Atemn(a0indx);
		%

		if (k==1)
   		Avhx = Avhy;
   		hAvhx = hAvhy;
   		if hAvhy > hAvhh
      		hAvhh = hAvhy;
      		Avhh = Avhy;
   		end
   		cJump = cJump+1;
		else
			AvhxD = Avhy;
			hAvhxD = hAvhy;
		end
	end
end