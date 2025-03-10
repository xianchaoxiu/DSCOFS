function [X_0]= intialx(A,B)
A=A-mean(A);
A=A';
maxerr=0;
d=size(A,1);
m=length(unique(B));
for i=1:10
%     [W,~] = qr(randn(d,m),0);
    W=RandOrthMat(d, m);
    AA=A*A';
    err=-trace(W'*AA*W);
    if err<maxerr
        X_0    = W;
        maxerr = err;
    end
end
