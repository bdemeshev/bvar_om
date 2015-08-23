function [yrEnd,qmEnd,forep,forepq,forepy,forelabel] = pararc(q_m)

yrEnd=2000;
qmEnd=12;
%forep = 24;      % 35, allow 1 year to get at least 1 calendar year
if (qmEnd==q_m)     % end of the year
	qmSub=0;         % Sub: substraction or substitute
else
	qmSub=qmEnd;     % Sub: substraction or substitute
end
%forep = 4*q_m-qmSub;
forep = 4*q_m;   %461;   %4*q_m;

%-------------------- below, automatic -------------------
forepq = floor((forep + (3-mod(qmEnd,3))*mod(qmEnd,3)) / 3);
                  % e.g, if qmEnd=5, 2nd quarter belogns to forepq
forepy = floor( (forep + qmEnd*(qmEnd~=12)) / q_m);     % annually
						% e.g., if qmEnd<12, forepy includes this particular year
if forepy == 0
   error('forep must have 12 months to make it one calendar years')
end

ystr=num2str(yrEnd);
forelabel =	[ ystr(3:4) ':' num2str(qmEnd) ' Forecast'];
