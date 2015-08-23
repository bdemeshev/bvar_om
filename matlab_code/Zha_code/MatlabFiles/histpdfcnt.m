function [imfpdf,imfpo,imfprob] = histpdfcnt(imfcnt,ndrawscnt,forep,shockp,nvar,ninv,invc,Am,...
                        ixdeng,findex,sindex,vindex,idxuniscl,lval,hval,ncom,gIndx)
% [imfpdf,imfpo,imfprob] = histpdfcnt(imfcnt,ndrawscnt,forep,shockp,nvar,ninv,invc,Am,...
%                   ixdeng,findex,sindex,vindex,idxuniscl,lval,hval,ncom,gIndx)
%       From already ordered and binned (cnt: count) series (not pdf yet).
%       Produce (1) the dataset for generating p.d.f., (2) the dataset for
%              probability (NOT density) at each bin, (3) with option ixdeng=1,
%              graphs of density functions
%
% imfcnt:  2+ninv-by-forep*shockp*nvar.  Counted logcnt, yhatqgcnt, or yhatCalygcnt.
%         In the case of impulse responses, forep=imstp and shockp=nvar.
%         In case of of A0(Avhx), forep=1, shockp=1, nvar=nfp
%         In case of forecasts, shockp=1
% ndrawscnt:  a total number of draws used in imfcnt
% forep:  forecast periods -- 1st dim  (must be compatible with findex)
% shockp:   number of shocks -- 2nd dim  (must be compatible with sindex)
% nvar:   number of responses (variables) -- 3rd dim  (must be compatible with vindex)
% ninv:   the number of small interior intervals for sorting.
% invc:  1-by-forep*shockp*nvar.  Whole inverval lenghth from min to max for one of
%                 (yhat, yhatqg, or yhatCalyg)
% Am:  1-by-forep*shockp*nvar.  Range5{i}(:,:,1) is lowest range for for one of
%                (yhat, yhatqg, or yhatCalyg)
% ixdeng:  index for density graphs.  1: enable; 0: disenable
% findex:  index for selected forecast periods, 1st dim (c.f., forep)
% sindex:  index for selected shocks, 2nd dim (c.f., shockp)
% vindex:  variable index, 3rd dim (c.f., nvar)
% idxuniscl:  1: scale on each graph by using lval and hval; 0: disenable this
% lval: (length(findex),length(sindex),length(vindex)); lowest point on the axis;
%        The number matches a total of findex, sindex, and vindex
% hval:  (length(findex),length(sindex),length(vindex)); highest point on the axis;
%        The number matches a total of findex, sindex, and vindex
% ncom:  number of combined intervals.  Large ncom implies large size of the bin.
%            Must set ncom so that ninv can be divided
% gIndx: 1: plot graphs of pdf's; 0: no plot
%-----------------
% imfpdf:  2+ninv-by-forep*shockp*nvar.  Density
% imfpo:  2+ninv-by-forep*shockp*nvar.  Bin position (x-axis) in relation to imfs3
% imfprob:  2+ninv-by-forep*shockp*nvar.  Probability (NOT density) at each bin
%
% 27 August 1998 Tao Zha
% Revised, October 1998, March 1999
% 3/24/99, added gIndx so that the previous programs may not be compatible.
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
imfprob = imfcnt ./ ndrawscnt;    % probability (NOT density)
imfpdf = imfprob ./ invlengthM;    % density

imfpo = [1:2+ninv]';    % positions for each forecast
imfpo = imfpo - 2;  % the first step to put back to original form
imfpo = repmat(imfpo,[1,forep*shockp*nvar]);
invcM = repmat(invc,[2+ninv,1]);
AmM = repmat(Am,[2+ninv,1]);
imfpo = ((imfpo .* invcM) ./ ninv) + AmM;  % 2+ninv-by-forep*shockp*nvar
     % the final step to put back to original form of impulse responses


if mod(ninv,ncom)
   warning('Set ncom so that ninv can be divided')
   return
end
%
ninv2=ninv/ncom;
imfpdfn=zeros(2+ninv2,forep*shockp*nvar);  %<<>>
imfpon=imfpdfn;
%
for ik=1:ninv2   % from 2 to 1+ninv2 for imfpdfn and imfpon
   imfpdfn(1+ik,:) = sum(imfpdf(1+(ik-1)*ncom+1:1+ik*ncom,:))/ncom;
   imfpon(1+ik,:) = mean(imfpo(1+(ik-1)*ncom+1:1+ik*ncom,:));
end

imfpdf4 = reshape(imfpdfn,2+ninv2,forep,nvar,shockp);
% imfpdf3: row--bin numbers, column--steps, 3rd dim--variables(responses), 4rd dime--shocks
imfpdf4s = permute(imfpdf4,[1 2 4 3]);
      % imf3pdfs: permuted so that
      %     row--bin numbers, column--steps, 3rd dim--shocks, 4rd dime--variables (responses)
imfpo4 = reshape(imfpon,2+ninv2,forep,nvar,shockp);
imfpo4s = permute(imfpo4,[1 2 4 3]);


if gIndx
   ck1=0;
   for k1=findex
      ck1=ck1+1;
      ck2=0;
      for k2=sindex
         ck2=ck2+1;
         ck3=0;
         for k3=vindex
            ck3=ck3+1;
            %figure
            if idxuniscl
               lpos = min(find(imfpo4s(2:1+ninv2,k1,k2,k3)>lval(ck1,ck2,ck3)));   % low position
               hpos = max(find(imfpo4s(2:1+ninv2,k1,k2,k3)<hval(ck1,ck2,ck3)));    % high position
               plot(imfpo4s(lpos+1:hpos+1,k1,k2,k3),imfpdf4s(lpos+1:hpos+1,k1,k2,k3))   %
                  % push everything forward by 1 because lpos and hpos start at 2
               set(gca,'XLim',[lval(ck1,ck2,ck3) hval(ck1,ck2,ck3)])
               grid
            else
               plot(imfpo4s(2:1+ninv2,k1,k2,k3),imfpdf4s(2:1+ninv2,k1,k2,k3))   %
               grid
            end
         end
      end
   end
end

%xlabel('The 1984 U Forecast')   %The 1982 GDP Growth Forecast')
%set(gca,'XLim',[lval hval])
%xact = [7.508 7.508];  %[-2.13 -2.13];
%yact = [0 0.5];  %[0 0.5];
%set(line(xact,yact),'Linestyle','-')
%title('Figure 8')
%line([-2 -2],[0 400])