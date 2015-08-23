r=2;
s=2;
n=6;
explosive=s;

div=1.0;

cont=1;
while cont > 0
    [q tmp]=qr(rand(n,n));
    [z tmp]=qr(rand(n,n));  
    a=q*diag(ones(n,1)+rand(n,1))*z;
    %a=q*z;
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
    
    if eu(1) == 1
        cont=0;
    end
end

% Simulate
y0=rand(n,1)
n_sim=10;
for i=1:10
    epsilon=rand(r,1);
    y=G1*y0+impact*epsilon;
    dwy=dwG1*y0+dwimpact*epsilon;
    
    y
    dwy
    norm(y-dwy)
    pause
    
    y0=dwy;
end