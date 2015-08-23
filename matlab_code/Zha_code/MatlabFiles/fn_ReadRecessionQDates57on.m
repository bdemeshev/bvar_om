function rec_dates = fn_ReadRecessionQDates57on()
% Recession dates.  For some reasons, Matlab code will only take pairs of rec_dates when using rectangle.
% Note that -1 in the following is for the graph lie in the correct position.
%   For example, 1961+(1-1)/4 means Q1 1961 plotted right at 1961.

rec_dates=[ ...
      1957+(3-1)/4 1958+(2-1)/4
      1960+(2-1)/4 1961+(1-1)/4
      1969+(4-1)/4 1970+(4-1)/4
      1973+(4-1)/4 1975+(1-1)/4
      1980+(1-1)/4 1980+(3-1)/4
      1981+(3-1)/4 1982+(4-1)/4
      1990+(3-1)/4 1991+(1-1)/4
      2001+(1-1)/4 2001+(4-1)/4
      2007+(4-1)/4 2009+(2-1)/4];