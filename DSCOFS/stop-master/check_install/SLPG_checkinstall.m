clear all;

n = 1000;
p = 20;
gamma = 0.02;
m = 200;

A = randn(m,n);
A = A - repmat(mean(A,1),m,1);
scale = norm(A,2);
A =A/scale;
A1 = A;   
AtA = A'*A;

Xinit = randn(n,p);
[Xinit,~] = qr(Xinit,0);


opts.dim = [n,p];
% opts0.maxit = 8000;
% opts0.gtol = 1e-8;
opts.gamma = gamma;
% opts0.X = Xinit;
opts.info = 0;
obj_fun = @fun_PCA_l1_lr;



tic
[out] = stop_slpg(obj_fun,opts,A,opts.gamma);
t_PenC = toc;
feas_penc = norm(out.X'*out.X -eye(p),'fro');

% semilogy(1:out.iter, out.kkts, 1:out.iter, out.feas, 1:out.iter, out.feas_subs)





function [fvalue,grad] = fun_PCA_l1_lr(X,A,gamma)
% min -0.5*(X^TA^TAX) + gamma* ||X||_1 
% low-rank version
% 
[n,p] = size(X);
AX = A*X;
AAX = A'*AX;
grad = -AAX;
fvalue = 0.5*sum(sum(X.*grad)) + gamma* sum(sum(abs(X)))  + 2*norm(X'*X - eye(p),'fro');
end

