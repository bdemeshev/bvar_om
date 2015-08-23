function out=sumgamlnx(b,a);

out=0;

for i=1:b
out=out+gammaln((a+1-i)/2);
end
