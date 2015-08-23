function [ beta,PROBLEM] = getcoef( mstar,sigma,ixx,MAXTRYS,N,L )
PROBLEM=0;
%
vstar=kron(sigma,ixx);
check=-1;
tryx=1;
   while check<0 && tryx<MAXTRYS
beta=mstar+(randn(1,N*(N*L+1))*chol(vstar))';

CH=stability(beta,N,L);
    if CH==0
    check=10;
    else
      tryx=tryx+1;
    end
   end
   if CH>0
       PROBLEM=1;
   end
   


end

