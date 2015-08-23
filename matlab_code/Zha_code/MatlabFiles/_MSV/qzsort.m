function [A,B,Q,Z] = qzsort(A,B,Q,Z)
%function [A,B,Q,Z] = qzsort(stake,A,B,Q,Z)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that abs(B(i,i)/A(i,i)) are increasing, while preserving U.T. and 
% orthonormal properties and Q'AZ' and % Q'BZ'.
%
% by Daniel Waggoner based on code by Christopher A. Sims

n=size(A,1);

for i=2:n
  for j=i:-1:2
    if abs(A(j,j)*B(j-1,j-1)) > abs(A(j-1,j-1)*B(j,j))
      [A B Q Z]=qzswitch(j-1,A,B,Q,Z);
    else
      break;
    end
  end
end

