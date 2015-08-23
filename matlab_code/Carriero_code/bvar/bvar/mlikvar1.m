function out=mlikvar1(y,x,yd,xd)
N=size(y,2);
T=size(y,1);
v=rows(yd);
y1=[y;yd];
x1=[x;xd];
%prior moments
xx0=xd'*xd;
invxx0=pinv(xx0); 
b0=invxx0*xd'*yd;
v0=rows(yd);
e0=yd-xd*b0;
sigma0=e0'*e0;
%posterior moments
xx1=x1'*x1;
invxx=pinv(xx1);
b=invxx*x1'*y1;

v1=v0+T; 
e=y1-x1*b; 
sigma1=e'*e;


PP=inv(eye(T)+x*invxx0*x'); 
QQ=sigma0;


lngam_ratio=mgamln(N,v0)-mgamln(N,v1);


py=-(lngam_ratio+(T*N/2)*log(pi))+0.5*N*log(det(PP))+(v0/2)*log(det(QQ))-(v1/2)*log(det(QQ+(y-x*b0)'*PP*(y-x*b0)));
out=py;