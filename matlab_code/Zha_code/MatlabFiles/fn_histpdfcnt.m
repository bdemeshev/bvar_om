function [imfpdf,imfpo,imfprob] = fn_histpdfcnt(imfcnt,ndrawscnt,imfloor,hbin,ninv,...
                          gIndx,ncom,kIndx,lval,hval)
% [imfpdf,imfpo,imfprob] = fn_histpdfcnt(imfcnt,ndrawscnt,imfloor,hbin,ninv,gIndx,ncom,kIndx,lval,hval)
%       From already ordered and binned (cnt: count) series (not pdf yet).
%       Produce (1) the dataset for generating p.d.f.;
%               (2) the dataset for probability (NOT density) at each bin;
%               (3) with option gIndx=1, graph density functions.
%
% imfcnt:  2+ninv-by-k.  Counted structural parameters, qmygcnt, ygcnt, or impulse responses.
% ndrawscnt:  a total number of draws used in imfcnt
% imfloor: k-element vector of low values of imf but imf_h can be below "imfloor" and above "imceiling"
% hbin:  k-element vector of bin lengths = (imceilling-imfloor)/ninv. Need not be a 1-by-k row vector.
% ninv:   the number of bins (small interior intervals) between ``imfloor'' and ``imfceiling''
% gIndx: 1: plot graphs of pdf's; 0: no plot.
%         If gIndx=0, ncom, kIndx, lval, and hval are irrelevant and no densities will be plotted.
% ncom:  number of combined intervals.  Large ncom implies large size of the bin.
%            Must set ncom so that ninv can be divided
% kIndx: index for selected variables.  Length(kIndx)<=k.
%         If kIndx=[], lval and hval are irrelevant and no densities will be plotted.
% lval: length(kIndx)-element vector of the lowest values on the axis;
%        The vector corresponds to kIndx.  This option is disenabled by setting it to [].
% hval: length(kIndx)-element vector of the highest values on the axis;
%        The vector corresponds to kIndx.  This option is disenabled by setting it to [].
%-----------------
% imfpdf:  2+ninv-by-k.  Density (NOT probability)
% imfpo:  2+ninv-by-k.  Bin position (x-axis) in relation to imfs3
% imfprob:  2+ninv-by-k.  Probability (NOT density) at each bin
%
% 27 August 1998 Tao Zha
% Revised, October 1998, March 1999, August 2000.
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


if (nargin==5), gIndx=0; end
if (nargin==8), lval=[]; hval=[]; end

imfloor=imfloor(:); hbin=hbin(:);
k=size(imfcnt,2);
invlengthM = repmat(hbin',[2+ninv,1]);
%invlengthM([1 2+ninv],:) = 1;  % August 2000.  I comment this out because it leads to an extorted picture
                          % when the rest of invlengthM is much smaller or larger than 1.
        % For the first interval (-inf, low bound) and last interval (high bound, +inf),
        % the length is set at 1.  Of course, theoretically, the interval length should be set
        % to infinity.  Adjust low and high bounds if 1 is not large enough compared with hbin.
imfprob = imfcnt ./ ndrawscnt;    % probability (NOT density)
imfpdf = imfprob ./ invlengthM;    % density

imfpo = [1:2+ninv]';    % positions for each forecast
imfpo = imfpo - 2;  % the first step to put back to original form
imfpo = repmat(imfpo,[1,k]);
imfloorM = repmat(imfloor',[2+ninv,1]);
imfpo = (imfpo .* invlengthM) + imfloorM;  % 2+ninv-by-k
     % the final step to put back to original form of impulse responses

if gIndx
   if mod(ninv,ncom)
      warning('Set ncom so that ninv can be divided')
      return
   end
   %
   ninv2=ninv/ncom;
   imfpdfn=zeros(2+ninv2,size(imfcnt,2));  %<<>>
   imfpon=imfpdfn;
   %
   for ik=1:ninv2   % from 2 to 1+ninv2 for imfpdfn and imfpon
      imfpdfn(1+ik,:) = sum(imfpdf(1+(ik-1)*ncom+1:1+ik*ncom,:))/ncom;
      imfpon(1+ik,:) = mean(imfpo(1+(ik-1)*ncom+1:1+ik*ncom,:));
   end

   ck1=0;
   for k1=kIndx
      ck1=ck1+1;
      figure
      plot(imfpon(2:1+ninv2,k1),imfpdfn(2:1+ninv2,k1))   % do not plot the values at -infty and +infty
      if ~isempty(lval) & ~isempty(hval)
         set(gca,'XLim',[lval(ck1) hval(ck1)])
      end
      grid
   end
end
