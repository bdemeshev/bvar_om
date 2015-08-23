function [imfl,imfh,imfl1,imfh1] = demarcate(imfcnt,imndraws,forep,nvar,...
                          ninv,invc,Am)
%  [imfl,imfh,imfl1,imfh1] = demarcate(imfcnt,imndraws,forep,nvar,...
%                          ninv,invc,Am)
% imfcnt:  2+ninv-by-forep*nvar.  Counted logcnt, yhatqgcnt, or yhatCalygcnt.
%         In the case of impulse responses, forep=imstp and nvar=nvar^2
% imndraws: total number of draws allocated to "imfcnt"
% forep:  the number of forecast periods (for both impulse responses and forecasting)
% nvar:  the number of variables.
% ninv:   the number of small interior intervals for sorting.
% invc:  1-by-forep*nvar.  Whole inverval lenghth from min to max for one of
%                 (yhat, yhatqg, or yhatCalyg)
% Am:  1-by-forep*nvar.  Range5{i}(:,:,1) is lowest range for for one of
%                (yhat, yhatqg, or yhatCalyg)
%-----------------
% imfl:  lower .95 bound, forep-by-nvar
% imfh:  higher .95 bound, forep-by-nvar
% imfl1: lower .68 bound, forep-by-nvar
% imfh1: higher .68 bound, forep-by-nvar
%
% Copyright (c) March 1998 Tao Zha
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




%$$$ .68 and .95 probability bands
imfl = zeros(forep*nvar,1);     % preallocating
imfh = zeros(forep*nvar,1);     % preallocating
imfl1 = zeros(forep*nvar,1);     % preallocating
imfh1 = zeros(forep*nvar,1);     % preallocating
imfpo = zeros(4,forep*nvar);    % 4 positions: l,h,l1,h1.
%
%tic
tem = cumsum(imfcnt) ./ imndraws;    % cumulative probability
tem = tem .* 100;
clear imfcnt

uptail1=5;   % 2.5, interval
lowtail1=95;  % 97.5

%t_cum = toc
%
%tic
%
%@@@ the following operations are valid only because tem are increasing!
for k = 1:forep*nvar
   %
   %@@@ the following operations are valid only because tem are increasing!
   %** 2.5% low tail
   if isempty(max(find(tem(:,k)<uptail1)))
      imfpo(1,k) = 1;
   else
      imfpo(1,k) = max(find(tem(:,k)<uptail1));
   end
   %** 16% low tail
   if isempty(max(find(tem(:,k)<16)))
      imfpo(2,k) = 1;
   else
      imfpo(2,k) = max(find(tem(:,k)<16));
   end
   %** 2.5% high tail
   imfpo(3,k) = min(find(tem(:,k)>lowtail1));
   %** 16% low tail
   imfpo(4,k) = min(find(tem(:,k)>84));
end

tem = imfpo';
%save outprob.txt tem -ascii
%

imfpo = imfpo - 2;  % the first step to put back to original form
% * 2.5% low tail
imfs = ((imfpo(1,:) .* invc) ./ ninv) + Am;
    % the final step to put back to original form of impulse responses
imfl = imfs';
% * 16% low tail
imfs = ((imfpo(2,:) .* invc) ./ ninv) + Am;
    % the final step to put back to original form of impulse responses
imfl1 = imfs';
% * 2.5% high tail
imfs = ((imfpo(3,:) .* invc) ./ ninv) + Am;
    % the final step to put back to original form of impulse responses
imfh = imfs';
% * 16% high tail
imfs = ((imfpo(4,:) .* invc) ./ ninv) + Am;
    % the final step to put back to original form of impulse responses
imfh1 = imfs';
%
%e_t = toc


% *** write out final results
clear tem
imfl = reshape(imfl,forep,nvar);
imfh = reshape(imfh,forep,nvar);
imfl1 = reshape(imfl1,forep,nvar);
imfh1 = reshape(imfh1,forep,nvar);
%
%