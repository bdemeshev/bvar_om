function out=mgamln(N,v)
constant=(N*(N-1)/4)*log(pi);
term2=sumgamlnx(N,v);
out=constant+term2;