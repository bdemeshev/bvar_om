function fsave=bvar(datain,Lbig,lamdaPbig,tauPbig,epsilonP,REPS,BURN,HORZ,Update,maxtrys,fix,dateid)
data=packr(datain);
N=cols(data);
Y=packr(data);

% Additional priors for VAR coefficients
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
BB=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
    e0=[e0 etemp];
end

if fix==0
[ optL,optP,optD ] = getoptlagprior( Y,Lbig,lamdaPbig,tauPbig,deltaP,epsilonP,muP,sigmaP );
L=optL;
lamdaP=optP;
tauP=optD;
else
L=Lbig;
lamdaP=lamdaPbig;
tauP=tauPbig;
end
%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html


[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

%find optimal prior and lag length
N=cols(Y);
%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);


T=rows(Y);


fsave=zeros(REPS-BURN,HORZ,N); 

% Display empty line (to separate samples in the screen shot)
disp(sprintf(' ')) 


  Y0=[Y;yd];
  X0=[X;xd];
  %conditional mean of the VAR coefficients
  mstar=vec(X0\Y0);  %ols on the appended data
  
  xx=X0'*X0;
  ixx=xx\eye(cols(xx));  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
  sigma=eye(N); %starting value for sigma
  beta0=vec(X0\Y0);
  igibbs=1;
  jgibbs=0;
while jgibbs<REPS-BURN
	
    % Display progress:
    if mod(igibbs,Update)==0 
        disp(sprintf(' Replication %s of %s. Lag %s Prior Tightness %s date %s', ... 
             num2str(igibbs), num2str(REPS),num2str(L),num2str(lamdaP),dateid) );
    end
        
    
     %step 1: Sample VAR coefficients
     [ beta2,PROBLEM] = getcoef( mstar,sigma,ixx,maxtrys,N,L );
     if PROBLEM
         beta2=beta0;
     else
         beta0=beta2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma=iwpq(T+rows(yd),inv(scale));

     if igibbs>BURN && ~PROBLEM
         jgibbs=jgibbs+1;
         yhat = get_pathsVAR(N, HORZ+1, T, L, Y, beta2,sigma);
         
         
         
         
         
         
         fsave(jgibbs,:,:)=[yhat(L+1:end-1,:)];
         
         
     end
     igibbs=igibbs+1;
end



end
  



