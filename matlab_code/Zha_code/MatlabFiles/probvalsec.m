function [probelow,valow,yhatmean] = probvalsec(yhatprob,yhatpo,ninv,forep,nvar,findex,...
                     sindex,valix,prolix)
% [probelow,valow,yhatmean] = probvalsec(yhatprob,yhatpo,ninv,forep,nvar,findex,...
%                     sindex,valix,prolix)
%
%   Probabilites (below) and values for selected variables and levels.
%   This program takes outputs from "histpdfcnt.m" which must be run first.
%   Compute (1) the probability of the x-value (such as FFR) that is below a level
%      prespecified by valix (e.g., probability of R below 5%); (2) the x-value
%     (such as median or .60 lower-tail value) below which the probability is at a level
%     prespecified by "prolix" (e.g., median when 0.50 is prespecified); (3) the mean
%     of the x-values (e.g., the mean of Pcm, M2, FFR, etc.).
%
% yhatprob:  2+ninv-by-forep*shockp(=1 here)*nvar.  Probability (NOT density) at each bin
% yhatpo:  2+ninv-by-forep*shockp(=1 here)*nvar.  Bin position (x-axis) in relation to yhat
% ninv:   the number of bins which are small interior intervals on the x-axis
% forep:  forecast periods -- 1st dim  (must be compatible with findex)
% nvar:   number of shocks or variables -- 2nd dim  (must be compatible with sindex)
% findex:  index for sected forecast periods 1st dim (c.f., forep)
% sindex:  index for selected shocks or variables, 2nd dim (c.f., nvar)
% valix:  length(findex)-by-length(sindex).  Selected values on the x-axis
% prolix:  length(findex)-by-length(sindex).  Selected (below) probabilites
%-----------
% probelow: length(findex)-by-length(sindex).  The probability of the x-value
%            (such as FFR) that is below a level prespecified by valix (e.g.,
%            probability of R below 5%).
% valow: length(findex)-by-length(sindex).  The x-value (such as median or .60
%            lower-tail value) below which the probability is at a level
%            prespecified by "prolix" (e.g., median when 0.50 is prespecified);
% yhatmean: forep-by-nvar.  The mean of forecasts of Pcm, M2, FFR, etc.
%
% 3/25/99 Tao A. Zha
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


yhatprob2=reshape(yhatprob,2+ninv,forep,nvar);
yhatpo2=reshape(yhatpo,2+ninv,forep,nvar);
%
%** Mean
yhatmean0 = sum(yhatprob.*yhatpo,1);
yhatmean = reshape(yhatmean0,forep,nvar);
%
probelow=zeros(length(findex),length(sindex));
            % probablity of the x-value that is below a specified level
            % E.g., probability of R below 5%.
valow=probelow;
            % x-value for the probability that is below a specified level
            % E.g., median when 0.50 is specified.
cumprob=cumsum(yhatprob2,1);   % cumulative probabilities
%
count1=0;
for k1=sindex     % forecast variables (or shocks)
   count1=count1+1;
   count2=0;
   for k2=findex    % forecast periods
      count2=count2+1;
      posix = max(find(yhatpo2(:,k2,k1)<valix(count2,count1)));
            % index for the position corr. to the x-value respeicied.
      if isempty(posix)
         probelow(count2,count1) = cumprob(1,k2,k1);
      else
         probelow(count2,count1) = cumprob(posix,k2,k1);
      end
      posix1 = max(find(cumprob(:,k2,k1)<prolix(count2,count1)));
            % index for the position corr. probability (below) specified.
      if isempty(posix1)
         valow(count2,count1) = yhatpo2(1,k2,k1);
      else
         valow(count2,count1) = yhatpo2(posix1,k2,k1);
      end
    end
end



