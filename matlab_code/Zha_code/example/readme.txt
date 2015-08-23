                        Guide to Waggoner-Zha Gibbs Sampling (JEDC)


                                     Data Files

xd29.mat:    29 variables assembled in the MAT file.
idenml.mat:  Output from idenml.m.  Needed to run msstart.m and msprob.m
iden6s.prn:  LZ structural identification.
iden6r.prn:  Reduced-form (Choleski).  


                                      Programs

idenml.m:  Finds the ML estimate with a particular identification.
msstart.m: A start file for both data and Bayesian estimation (when IxEstima==0, no Bayesian 
           estimation but only data), which allows you to set (a) data range, (b) sample range, 
           (c) conditions of shocks ‘Cms’, (d) conditions of variables ‘nconstr’, (e) soft 
           conditions ‘nbancon,’  (f)  produce DLS conditional forecast (by setting nconstr and 
           gDLSIx).  
msprob.m:  Generates draws of A0’s (only if impulse responses are commented out) and impulse 
           responses under importance sampling or Gibbs sampling.  Calls msstart.m and exports 
           xdrawsIS50k.mat and xdrawsGibbs50k.mat.


                                      Functions

pararc.m:  Function to set the end period for the sample as well as the forecast horizon.  
           Called by idenml.m and msstart.m.
szasbvar.m:  Function to export SZ asymmetric prior VAR with a prior set at SZ IER.  Called by 
             idenml.m and msstart.m.
ftd_gibbsmp.m:  Function for Gibbs sampling of selected impulse responses.  Check <<>>1 lines 
                for selected variables, shocks, and steps.
ftd_impsmp.m:  Function for importance sampling of selected impulse responses.  Check <<>>1 
               lines for selected variables, shocks, and steps.
 

                           Drawing Impulse Responses
                Drawing A0s with importance-sampling and Gibbs

Change pararc.m (if neccesary)

Run idenml (Only with Rform=0)

Change msstart (must set actup=nSample-lags and Pesudo=0 for get structural shocks)

Run msprob (Choose if 1 for Gibbs and else for importance sampoing) Rename outxdraws.mat to 
           xdrawsGibbs50.mat or xdrawsIS50k.mat.

