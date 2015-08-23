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

    [G1,G2,gev,eu]=msv_all(a,b,psi,pi);

    k=size(eu,2);
    
    for i=1:k
        eu{i}
        abs(gev{i}(:,2)./gev{i}(:,1))
        G1{i}
        G2{i}
        i
        k
        pause
    end
end

