function [x] = bernoullirnd(p)
%-------------------------------------------------------------------
% Generate a scalar x from Bernoulli distribution
%-------------------------------------------------------------------
% x is a 0-1 variable following the distribution:
%                    
%         prob(x) = [(p)^x][(1-p)^(1-x)]
%
%-------------------------------------------------------------------
%  INPUTS:
%   n     - The number of variates you want to generate
%   p     - Associated probability
%
%  OUTPUT:
%   x     - [1x1] Bernoulli variate
%-------------------------------------------------------------------

   
u=rand(1,1);
if u<p   
    x = 0;
else
    x = 1;
end
