function [A0n,a0dpindx,nswitch,A0inn] = fn_nmlzvar(A0u,A0xhat,A0inxhat,IndxNmlr,nswitch,A0inu)
% [A0n,nswitch,A0inn] = nmlzvar(A0u,A0xhat,A0inxhat,IndxNmlr,nswitch,A0inu)
%    Export normalized new draw A0 (and inv(A0) only if IndxNmlr(5)) and the number of sign switches
%    Ref.:  D.F. Waggoner and T.A. Zha: "Does Normalization Matter for Inference?"
%    See Note Forecast (2) pp. 52-53
%
% A0u:  unnormalized A0; column--equation
% A0xhat:  ML estimate or posterior mode of A0
% A0inxhat:  inv(A0xhat)
% IndxNmlr: index for which normalization rule to choose
%     Only one of the elments in IndxNmlr can be non-zero
%     IndxNmlr(1): ML A distance rule (supposed to be the best)
%     IndxNmlr(2): ML Ahat distance rule (to approximate IndxNmlr(1))
%     IndxNmlr(3): ML Euclidean distance rule (not invariant to scale)
%     IndxNmlr(4): Positive diagonal rule
%     IndxNmlr(5): Positive inv(A) diagonal rule (if ~IndxNmlr(5), no need for A0inu,
%                                      so we set A0inn=[])
%     IndxNmlr(6): Assigned postive rule (such as off-diagonal elements).  Added 1/3/00
% nswitch:  # of sign switches
% A0inu:   unnormalized inv(A0); used only if IndxNmlr(5)
%-----------------
% A0n:  normalized new A0; column--equation
% a0dpindx:  Index of the columns in A0 or A+ whose signs need to be switched.
% nswitch:  updated # of sign switches
% A0inn:   normalized inv(A0); used only if IndxNmlr(5)
%
% Written by Tao Zha
% 1/3/00: added IndxNmlr(6) so that previous programs may not be compatible.
% 10/11/01:  added a0dpindx in front of nswitch as an output argument so that previous programs may not be compatible.

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


A0inn = [];    % no output for normalized A0in unless IndxNmlr(5)

if (length(find(IndxNmlr))>1)
   warning('You cannot choose more than one normalization rule at a time')
   disp('Press ctrl-c to abort')
   pause
elseif isempty(find(IndxNmlr))      % no normalization
   A0n=A0u; nswitch=0;
elseif IndxNmlr(1)
   a0dpindx = find(diag(A0u\A0xhat)<0);
   A0n = A0u;
   if ~isempty(a0dpindx)
      A0n(:,a0dpindx) = -A0u(:,a0dpindx);   % normalized new draws
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
elseif IndxNmlr(2)
   a0dpindx = find(diag(A0inxhat*A0u)<0);
   A0n = A0u;
   if ~isempty(a0dpindx)
      A0n(:,a0dpindx) = -A0u(:,a0dpindx);   % normalized new draws
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
elseif IndxNmlr(3)
   Adiff = (A0u - A0xhat).^2;  % distance that may be far from axhat or A0xhat
   Adiffn = (-A0u - A0xhat).^2; % distance by chaning the sign of Atem
   cAdiff = sum(Adiff);    % each column summed up
   cAdiffn = sum(Adiffn);  % each column summed up
   cAindx = find(cAdiffn<cAdiff); % index for shorter distance
   A0n = A0u;
   if ~isempty(cAindx)
      A0n(:,cAindx) = -A0u(:,cAindx); % find the shortest or nearest distance
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
elseif IndxNmlr(4)
   a0dpindx = find(diag(A0u)<0);
   A0n = A0u;
   if ~isempty(a0dpindx)
      A0n(:,a0dpindx) = -A0u(:,a0dpindx);   % normalized new draws
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
elseif IndxNmlr(5)
   a0dpindx = find(diag(A0inu)<0);
   A0n = A0u;
   A0inn = A0inu;
   if ~isempty(a0dpindx)
      A0n(:,a0dpindx) = -A0u(:,a0dpindx);
      A0inn(a0dpindx,:) = -A0inu(a0dpindx,:);
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
elseif IndxNmlr(6)         %*** This one has to be MANUALLY handled
   [jnk,nvar]=size(A0u);
   A0dummy=A0u;
   A0dummy(:,1:2)=-A0u(:,1:2);   % keep it to a sign that coincide with Brooking paper
   a0dpindx = find(A0dummy(nvar,:)<0);   % the last row
   A0n = A0u;
   if ~isempty(a0dpindx)
      A0n(:,a0dpindx) = -A0u(:,a0dpindx);   % normalized new draws
      nswitch = nswitch + 1;  %<<>> # of sign switches
   end
end
