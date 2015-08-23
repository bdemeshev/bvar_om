function [ optL,optP,optD,outmlik,tableP,tableL ] = getoptlagprior( data,Lbig,Lamdabig,taubig,deltaP,epsilonP,muP,SigmaP )
N=cols(data);
outmlik=zeros(length(Lbig),length(Lamdabig));
tableL=zeros(length(Lbig),length(Lamdabig));
tableP=tableL;
tableD=tableL;
for i=1:length(Lbig)
    L=Lbig(i);
    
   
    Y=data;
    X=[];
for jj=1:L
X=[X lag0(data,jj) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);

for j=1:length(Lamdabig)
    [yd,xd] = create_dummies(Lamdabig(j),taubig(j),deltaP,epsilonP,L,muP,SigmaP,N);

    %get marginal likelihood
    mlik=mlikvar1(Y,X,yd,xd);
    outmlik(i,j)=mlik;
    tableL(i,j)=L;
    tableP(i,j)=Lamdabig(j);
    tableD(i,j)=taubig(j);
end

end
% outmlik
% tableL
% tableP
id=outmlik==max(max(outmlik));
optL=tableL(id);
optP=tableP(id);
optD=tableD(id);
