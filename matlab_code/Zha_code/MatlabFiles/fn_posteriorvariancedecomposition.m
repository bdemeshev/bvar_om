function [VarXAll, E_YVarXGivenY, Var_YEXGivenY] = fn_PosteriorVarianceDecomposition(xydraws,probs)
%--- Inputs:
% xydraws: a matrix with rows being x draws and columns being y draws (columns to be conditioned).
% probs: a vector of probabilities for y draws (length(probs) must be the same as size(xydraws,2))
%--- Outputs:
% VarXAll: Var(X) -- overall variance
% E_YVarXGivenY: E_Y Var(X|Y) -- uncertainty of not knowing the correct true value of X for a given Y.
% Var_YEXGivenY: Var_Y E(X|Y) -- varianation in the estimator of X due to uncertainty about the choice of Y.
%
% Mathematically and numerically, it must be that VarXAll = E_YVarXGivenY + Var_YEXGivenY.
%
% Written by T. Zha, October 2009
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

n_ydraws = size(xydraws,2);
if (nargin==1), probs = (1/n_ydraws)*ones(n_ydraws,1); end
                                                                 

%---------------- Overall variance --------------------
VarXAll1st = 0.0; 
VarXAll2nd = 0.0;      
xy1st = mean(mean(xydraws));
xy2nd = mean(mean(xydraws.^2));
for (yi=1:n_ydraws)                    
   xGiveny1st = mean(xydraws(:,yi));
   xGiveny2nd = mean(xydraws(:,yi).^2);
   VarXAll1st = VarXAll1st + probs(yi) * xGiveny1st;
   VarXAll2nd = VarXAll2nd + probs(yi) * xGiveny2nd;
end                                                          
VarXAll = VarXAll2nd - VarXAll1st^2; 


%-------------- The component E_YVarXGivenY --------------------
E_YVarXGivenY = 0.0;
for (yi=1:n_ydraws)                            
   VarXGivenY = mean(xydraws(:,yi).^2) - (mean(xydraws(:,yi)))^2;
   E_YVarXGivenY = E_YVarXGivenY + probs(yi) * VarXGivenY;
end                                            
E_YVarXGivenY = E_YVarXGivenY;

%-------------- The component Var_YEXGivenY --------------------                                               
Var_YEXGivenY1st = 0.0; 
Var_YEXGivenY2nd = 0.0;                        
for (yi=1:n_ydraws)                                  
   EXGivenY = mean(xydraws(:,yi));
   Var_YEXGivenY1st = Var_YEXGivenY1st + probs(yi) * EXGivenY;
   Var_YEXGivenY2nd = Var_YEXGivenY2nd + probs(yi) * EXGivenY^2;
end                                                          
Var_YEXGivenY = Var_YEXGivenY2nd - Var_YEXGivenY1st^2;
