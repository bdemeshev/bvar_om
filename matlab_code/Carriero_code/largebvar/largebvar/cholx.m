function out=cholx(x)
[R,p] = chol(x);
if p==0
    out=R;
else
    out=real(sqrtm(x))';
end