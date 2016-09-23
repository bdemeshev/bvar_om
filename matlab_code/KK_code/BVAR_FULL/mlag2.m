function [Xlag] = mlag2(X,p)
%MLAG2 Create autoregressive lags
%   This function creates lags of the matrix X(t), in the form:
%            Xlag = [X(t-1),...,X(t-p)]
%   Written by Dimitris Korobilis, March 2007
[Traw,N]=size(X);
Xlag=zeros(Traw,N*p);
for ii=1:p
    Xlag(p+1:Traw,(N*(ii-1)+1):N*ii)=X(p+1-ii:Traw-ii,1:N);
end


% %OR:
% [Traw,N]=size(X);
% Xlag=zeros(Traw,N,p);
% for ii=1:p
%     Xlag(p+1:Traw,:,ii)=X(p+1-ii:Traw-ii,:);
% end
% Xlag=Xlag(:,:);