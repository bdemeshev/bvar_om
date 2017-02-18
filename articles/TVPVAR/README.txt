
This folder contains some notes and code on a time-varying coefficient VAR à la Primiceri (2005), corrected as in 
Del Negro and Primiceri (2015). 
It is available under http://www.bkolb.eu/codes/TVPVAR.zip. 
By Benedikt Kolb, October 2015.

This folder contains the following files:

==========
run_TVPVAR1.m: Code that replicates the results by Del Negro and Primiceri available at
http://restud.oxfordjournals.org/content/suppl/2015/06/22/rdv024.DC1/Supplementary.zip.

run_TVPVAR2.m: Similar to the above, but with an ordering of the Gibbs sampler as in the note. 
The best for studying with the note, I think.

run_TVPVAR3.m: Similar to the above, but additionally keeping A1(t) and A23(t) separated.

KalmanCarterKohn.m: Routine to run the backward recursions as in Carter and Kohn (1994)on the coefficients.

makelag.m: Simple routine to generate lags of a vector or matrix.

data.mat: The dataset from Del Negro and Primiceri (2015).
==========

For comments, contact me at benedikt@bkolb.eu. Enjoy!