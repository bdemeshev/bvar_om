r=2;
s=2;
n=6;
explosive=n;

div=1.0;

a=rand(n,n);
b=rand(n,n);
[a b q z]=qz(a,b);
a=q*diag(ones(n,1)+rand(n,1))*z;
d=rand(n,1);
d(1:explosive)=d(1:explosive)+ones(explosive,1);
b=q*diag(d)*z;

psi=rand(n,r);
pi=rand(n,s);
c=zeros(n,1);

[G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(a,b,c,psi,pi,div);
[dwG1,dwimpact,dwgev,dweu]=dw_gensys(a,b,psi,pi,div);

G1
dwG1

impact
dwimpact

eu
dweu