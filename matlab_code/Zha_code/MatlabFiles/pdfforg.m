function [imfpdf,imfpo,imfprob] = pdfforg(imfcnt,imndraws,forep,nvar,ninv,invc,Am)
% [imfpdf,imfpo,imfprob] = pdfforg(imfcnt,imndraws,forep,nvar,ninv,invc,Am)
%
%       Produce the dataset for generating p.d.f. and export the dataset for
%              probability (NOT density) at each bin
%
% imfcnt:  2+ninv-by-forep*nvar.  Counted logcnt, yhatqgcnt, or yhatCalygcnt.
%         In the case of impulse responses, forep=imstp and nvar=nvar^2
% ninv:   the number of small interior intervals for sorting.
% invc:  1-by-forep*nvar.  Whole inverval lenghth from min to max for one of
%                 (yhat, yhatqg, or yhatCalyg)
% Am:  1-by-forep*nvar.  Range5{i}(:,:,1) is lowest range for for one of
%                (yhat, yhatqg, or yhatCalyg)
%-----------------
% imfpdf:  2+ninv-by-forep*nvar.  Density
% imfpo:  2+ninv-by-forep*nvar.  Bin position (x-axis) in relation to imfs3
% imfprob:  2+ninv-by-forep*nvar.  Probability (NOT density) at each bin
%
% 27 August 1998 Tao Zha
% Revised, October 1998
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

invlength = invc ./ ninv;
invlengthM = repmat(invlength,[2+ninv,1]);
invlengthM([1 2+ninv],:) = 1;
        % first (-inf, low bound) and last (high bound, +inf), the interval is set
        % to be 1.  Of course, theoretically, the interval length should be set
		  % to infinity, but 1 is large enough compared with invlength.
imfprob = imfcnt ./ imndraws;    % probability (NOT density)
imfpdf = imfprob ./ invlengthM;    % density

imfpo = [1:2+ninv]';    % positions for each forecast
imfpo = imfpo - 2;  % the first step to put back to original form
imfpo = repmat(imfpo,[1,forep*nvar]);
invcM = repmat(invc,[2+ninv,1]);
AmM = repmat(Am,[2+ninv,1]);
imfpo = ((imfpo .* invcM) ./ ninv) + AmM;  % 2+ninv-by-forep*nvar
     % the final step to put back to original form of impulse responses