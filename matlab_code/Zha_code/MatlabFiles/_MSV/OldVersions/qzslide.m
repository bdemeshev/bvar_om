function [a,b,q,z,pairs] = qzslide(a,b,q,z,pairs,i,j)
% function [a,b,q,z,pairs] = qzslide(a,b,q,z,pairs,i,j)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that the ith diagonal element is moved to the jth position, while 
% preserving U.T. and orthogonal properties and Q'AZ' and Q'BZ'.
%
% The array pairs is also rearranged.  Complex conjugate pairs are kept 
% adjacent.  If pairs(i)=1, then position i and i+1 are complex 
% conjugate pairs.  If pairs(i)=-1, then position i and i-1 are complex
% conjugate pairs.  If pairs(i)=0, then position i is real.
% 
%
% by Daniel Waggoner based on code by Christopher A. Sims

n=size(a,1);

if i == j
    return;
end

if i < j
    if pairs(j) == 1
        j=j+1;
    end
    if pairs(i) == -1
        j=j-1;
        i=i-1;
    else
        if (pairs(i) == 1) & (j == n)
            if (i == n-1)
                return;
            end
            j=n-1;
        end
    end
    while i < j
        if (pairs(i) == 1)
            [a b q z]=qzswitch(i+1,a,b,q,z);
            [a b q z]=qzswitch(i,a,b,q,z);
            pairs(i)=pairs(i+2);
            pairs(i+1)=1;
            pairs(i+2)=-1;
        else
            [a b q z]=qzswitch(i,a,b,q,z);
            pairs(i)=pairs(i+1);
            pairs(i+1)=0;
        end
        i=i+1;
    end 
else
    if pairs(j) == -1
        j=j-1;
    end
    if pairs(i) == 1
        j=j+1;
        i=i+1;
    else
        if (pairs(i) == -1) & (j == 1)
            if (i == 2)
                return;
            end
            j=2;
        end
    end
    while i > j
        if (pairs(i) == -1)
            [a b q z]=qzswitch(i-2,a,b,q,z);
            [a b q z]=qzswitch(i-1,a,b,q,z);
            pairs(i)=pairs(i-2);
            pairs(i-2)=1;
            pairs(i-1)=-1;
        else
            [a b q z]=qzswitch(i-1,a,b,q,z);
            pairs(i)=pairs(i-1);
            pairs(i-1)=0;
        end
        i=i-1;
    end
end
