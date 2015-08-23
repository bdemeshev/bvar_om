      function [ax,h]=subtitle(text)
      %
      %Centers a title over a group of subplots.
      %Returns a handle to the title and the handle to an axis.
      % [ax,h]=subtitle(text)
      %           returns handles to both the axis and the title.
      % ax=subtitle(text)
      %           returns a handle to the axis only.
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


%      ax=axes('Units','Normal','Position',[.05 .05 .9 .9],'Visible','off','FontSize',14,'FontWeight','Bold');
      ax=axes('Units','Normal','Position',[.05 .05 .9 .9],'Visible','off','FontSize',10,'FontWeight','Bold','FontName','Times New Roman');
%      ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
      set(get(ax,'Title'),'Visible','on')
      title(text);
      if (nargout < 2)
        return
      end
      h=get(ax,'Title');

      %%%END CODE%%%
