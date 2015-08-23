function [a,b,q,z] = qzmoveindex(a,b,q,z,pairs,idx,k,j)
% function [a,b,q,z] = qzmoveindex(a,b,q,z,pairs,idx,k,j)
%
% Takes U.T. matrices a, b, orthonormal matrices q,z, rearranges them
% so that the indices in idx are moved up to the jth position, while 
% preserving U.T. and orthogonal properties and q'az' and q'bz'. 
%
% by Daniel Waggoner based on code by Christopher A. Sims

while k > 0
   [a b q z pairs]=qzslide(a,b,q,z,pairs,idx(k),j);
   if (pairs(j) == -1) 
       j=j-2;
   else
       j=j-1;
   end
   k=k-1;
end