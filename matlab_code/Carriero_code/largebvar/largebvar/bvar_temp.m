function fsave=bvar(datain,Lbig,lamdaPbig,tauPbig,epsilonP,REPS,BURN,HORZ,Update,maxtrys,fix,dateid)
data=packr(datain);
N=cols(data);
Y=packr(data);

% Additional priors for VAR coefficients
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:3
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)]
    sigmaP=[sigmaP;stemp]
    e0=[e0 etemp];
end
fsave=0;
