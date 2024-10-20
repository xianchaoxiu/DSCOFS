clear all
n = 5000; 
p = 10;

L = gallery('tridiag',n,-1,2,-1); 
alpha = 100; 
[Ll,Lu] = lu(L);
tol_all = 1e-11;
tol_sub = 1e-2;

% intial point
X = randn(n, p);
[U, ~, V] = svd(X, 0); 
X_init = U*V';
tempM2 = alpha*(L\(sum(X_init.^2,2)));
tempM2 = spdiags(tempM2,0,n,n);
tempM = L + tempM2;
[U0, ~, ~] = eigs(tempM, p,'sm'); 
X0 = U0;




optsF_init.gtol = tol_sub;
optsF_init.X = X0;
optsF_init.info = 0;
[outs] = stop_slpg_smooth(@local_fun_KS,optsF_init,L,Ll,Lu,alpha);
Xtemp = outs.X;

[~,gradks] = local_fun_KS(Xtemp,L,Ll,Lu,alpha);
opts.penalparam = 10*norm(gradks,'fro');
opts.Hessf = @Local_Hess_KS;
opts.X = Xtemp;
opts.maxit = 20;
opts.gtol = tol_all;

opts.info = 0;

tic
outN = stop_pencs(@local_fun_KS,opts,L,Ll,Lu,alpha);
t1 = toc;
% fprintf('PenCS|  f: %8.6e, nrmG: %2.1e, cpu: %4.2f, OutIter: %3d, InnerIter: %4d, nfe: %4d,\n',...
%     outN.fval(end),outN.kkts(end), t1, outN.iter, outN.subiter, 000);
CPUs_PenCS = outN.times;

tic;



function [f,g] = local_fun_KS(X,L,Ll,Lu,alpha)
    LX = L*X;
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX); tempa = alpha*tempa;
    f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
    g = LX + bsxfun(@times,tempa,X);
end


function h = Local_Hess_KS(X, U, L,Ll,Lu,alpha)
        rhoX = sum(X.*X,2);
        rhoXdot = 2*sum(X.*U, 2);
        tempa = Lu\(Ll\rhoXdot);
        tempa = alpha*tempa;
        tempb = Lu\(Ll\rhoX);
        tempb = alpha*tempb;
        
        h = L*U + bsxfun( @times,tempa,X) + bsxfun(@times, tempb, U);
     
            
end



