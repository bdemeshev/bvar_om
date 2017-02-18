
 function ylag = makelag(y,nlag)

% Creates the nlagth lag of a vector or matrix y
% 
% Usage:  ylag = makelag(y,nlag)
%
% Inputs: y    - vector or matrix of which you want lags
%         nlag - how many period you want to lag
% 
% Output: ylag - nlagth lag of y
%
% By B. Kolb, Oct. 2015

[nr, nc] = size(y);

y1 = y(1:(nr - nlag),:);
ylag = [zeros(nlag,nc); y1];

end