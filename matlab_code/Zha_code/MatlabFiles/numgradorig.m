function [g, badg] = numgrad(fcn,x,varargin)
% function [g badg] = numgrad(fcn,xvarargin)
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
delta = 1e-6;  % works for a general case
%delta = 1e-2;  % bigger step -- works for a flat peak
%delta=1e-04;

n=length(x);
tvec=delta*eye(n);
g=zeros(n,1);
%--------------------old way to deal with variable # of P's--------------
%tailstr = ')';
%stailstr = [];
%for i=nargin-2:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end
%f0 = eval([fcn '(x' tailstr]); % Is there a way not to do this?
%---------------------------------------------------------------^yes
f0 = eval([fcn '(x,varargin{:})']);
% disp(' first fcn in numgrad.m ------------------')
%home
% disp('numgrad.m is working. ----') % Jiinil on 9/5/95
% sizex=size(x),sizetvec=size(tvec),x,    % Jinill on 9/6/95
badg=0;
for i=1:n
   scale=1; % originally 1
   % i,tveci=tvec(:,i)% ,plus=x+scale*tvec(:,i) % Jinill Kim on 9/6/95
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end
   g0 = (eval([fcn '(x+scale*tvecv'', varargin{:})']) - f0) ...
         /(scale*delta);
   % disp(' fcn in the i=1:n loop of numgrad.m ------------------')% Jinill 9/6/95
   % disp('          and i is')               % Jinill
   % i                         % Jinill
   % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see below Jinill 9/6/95
% -------------------------- special code to essentially quit here
   % absg0=abs(g0) % Jinill on 9/6/95
   if abs(g0)< 1e15
      g(i)=g0;
      % disp('good gradient') % Jinill Kim
   else
      disp('bad gradient ------------------------') % Jinill Kim
      % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0) %see above
      g(i)=0;
      badg=1;
      % return
      % can return here to save time if the gradient will never be
      % used when badg returns as true.
   end
end
%-------------------------------------------------------------
%     if g0 > 0
%        sided=2;
%        g1 = -(eval([fcn '(x-scale*tvec(:,i)''' tailstr]) - f0) ...
%           /(scale*delta);
%        if g1<0
%           scale = scale/10;
%        else
%           break
%        end
%     else
%        sided=1;
%        break
%     end
%  end
%  if sided==1
%     g(i)=g0;
%  else
%     if (g0<1e20)
%        if (g1>-1e20)
%           g(i)=(g0+g1)/2;
%        else
%           g(i)=0;
%           badg=1;
%           disp( ['Banging against wall, parameter ' int2str(i)] );
%        end
%     else
%        if g1>-1e20
%           if g1<0
%              g(i)=0;
%              badg=1;
%              disp( ['Banging against wall, parameter ' int2str(i)] );
%           else
%              g(i)=g1;
%           end
%        else
%           g(i)=0;
%           badg=1;
%           disp(['Valley around parameter ' int2str(i)])
%        end
%     end
%  end
%end
%save g.dat g x f0
%eval(['save g g x f0 ' stailstr]);
