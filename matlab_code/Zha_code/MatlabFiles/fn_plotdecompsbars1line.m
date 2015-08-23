function fn_PlotDecompsBars1Line(flags, dates_qm_matlab, nbars, nvars, nrows, ncols, allseries) 
%fn_PlotDecompsBars1Line(flags, dates_qm_matlab, nbars, nvars, nrows, ncols, allseries) 
%
%Ouputs:
%  Plots "nbars" decompositions with "nvars" subplots and with 1 line (or 1 series).
%Inputs:
%  flags.figure: 0: figure; 1: figure('PaperPosition',[0.25 1.0 8 9]).
%  flags.legend: 0: no legend; 1: legend to take allseries.barnames. 
%  flags.nber: 0: no NBER recession bars; 1: NBER recession bars (NOT working yet -- needs someone to help me out).
%  dates_qm_matlab: a vector of Matlab-encoded date numbers.
%  nbars: number of decomposition expressed in bars.
%  nvars: number of variabeles (and thus number of subplots).
%  nrows: number of rows for all subplots.
%  ncols: number of cols for all subplots.
%  allseries.barnames: nbars cells of names for bars.
%  allseries.bars: an fss-by-nbars matrix of decompositions.
%  allseries.lines: an fss-by-nlines matrix of series.
%  allseries.titles: nvars cells of names for the title of each subplot.
%  allseries.ylabels: nvars cells of names for the y-axis label of each subplot.
%  allseries.axis_4by1: a 4-by-1 vector of [yearbeg yearend ymin ymax].
%  allseries.colors: specifiying colors of bars.  Examples are:
%                  colors = [ .1 .1 .75  %blue
%                             .8  0  0   %red
%                              1 .7 .25   %deep yellow
%                              1  1  0    %light yellow
%                              1  0  1    %magenta
%                              0  1  1    %cyan
%                              0  1  0    %green
%                                    ]; %colormap(jet);
%
% Written by T. Zha, October 2011
% Copyright (C) 1997-2013 Tao Zha
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

ndays = 45;  % Number of days used to gauge the width of a bar.  The larger the number, the wider the bar.



if (nvars>(nrows*ncols))
   error('fn_PlotDecompsBars1Line(): nrows*ncols must be greater than or equal to nvars');
end
if ( (size(allseries.bars,2) ~= nbars) || (length(allseries.barnames) ~= nbars) ) 
   error('fn_PlotDecompsBars1Line(): make sure size(allseries.bars,2) = length(allseries.barnames) = nbars');
end
if ( (length(allseries.ylabels) ~= nvars) || (length(allseries.titles) ~= nvars) ) 
   error('fn_PlotDecompsBars1Line(): make sure length(allseries.ylabels) = length(allseries.titles) = nvars');
end
if flags.figure
   figure('PaperPosition',[0.25 1.0 8 9]); %[left bottom witdh height]  Tip: bottom + height = 10.0
else
   figure;
end


fss = length(dates_qm_matlab);
XY = [];
kj = 1;
for (ki=1:nvars)
   subplot(nrows, ncols,kj);
   kj= kj+1;
   
   hold on   
   if (0) (flags.nber)
      %---------------------------
      %-- Plotting with NBER recession bars 
      %---------------------------
      rec_dates = fn_ReadRecessionQDatesMatlabForm57on();   %Reading recession dates.
      %-- Get axis data before assigning it to AxisDat.
      %axis([1975 2011 (min(TightStandardLargeFirms)-0.1*abs(min(TightStandardLargeFirms)))   (max(TightStandardLargeFirms)+0.1*abs(max(TightStandardLargeFirms))) ]);  %1970 will automatically cut off the dates before that in Read_RecessionDates*.m.
      axis(allseries.axis_4by1);  %1970 will automatically cut off the dates before that in Read_RecessionDates*.m.
      AxisDat = axis;
      fn_PlotRecessionShades(rec_dates, AxisDat);
   end
   
   for kk=1:fss
      pos_position = 0; neg_position = 0;
      for jj=1:nbars  %size(sshocks,2)
         st = allseries.bars(kk,jj);
         x = [addtodate(dates_qm_matlab(kk),-ndays,'day') addtodate(dates_qm_matlab(kk),-ndays,'day') addtodate(dates_qm_matlab(kk),ndays,'day') addtodate(dates_qm_matlab(kk),ndays,'day')];
         if st < 0
             yi = st+neg_position;
             y = [neg_position yi yi neg_position];
             neg_position = yi;
         else
             yi = st+pos_position;
             y = [pos_position yi yi pos_position];
             pos_position = yi;
         end
         fill(x,y,allseries.colors(jj,:));
         XY(kk,jj,:) = y;
      end
   end
    
   plot(dates_qm_matlab, allseries.lines(:,1),'k','LineWidth',1.5);
   if (flags.legend)
      if (0)
         if (ki==nvars)
            lh = legend(allseries.barnames,'Location','BestOutside','Orientation','horizontal'); 
         end   
      else   
         legend(allseries.barnames,'Location','Best','Orientation','horizontal'); 
      end   
   end



   hold off

   datetick 'x';
   xlim([min(dates_qm_matlab)-ndays,max(dates_qm_matlab)+ndays])
   grid on
   ylabel(allseries.ylabels{ki})
   title(allseries.titles{ki});
    
end
