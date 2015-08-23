%% Estimates a Large BVAR Model with USA Data using MCMC - Gibbs Sampling Method.

load('usa_data'); %Contains Data for USA
data = usa_data; %Choose all variables

%Inputs
horizon=12;%forecast horizon
reps=10000;%total reps
burn=8000;%burn in
update=1000;%prints every update iter
maxtrys=1000; %max tries for stable draw
L=4;   %max number of lags in the VAR
%lamdap=0.1:0.1:1;% 0.1 %tightness of the prior uses value=max(marginal likelihood)
lamdap= 0.2 %tightness of the prior uses value=max(marginal likelihood)

taup=10*lamdap; %tightness of prior on sum of coefficients
epsilonp=1/1000; %tightness of prior on constant
%fix=0;%1; %if fix is 1 no search over prior and lags is done 
fix=0;%1; %if fix is 1 no search over prior and lags is done 
forecastf=bvar(data,L,lamdap,taup,epsilonp,reps,burn,horizon,update,maxtrys,fix,'1990Q1_usa');