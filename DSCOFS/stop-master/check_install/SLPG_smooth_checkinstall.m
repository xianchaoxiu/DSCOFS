clear all
test_instances = 3;

%% Problem setting
n = 1000;  p = 20; 
A = gallery('tridiag',n,-1,2,-1); 
alpha = 1; 
[Ll,Lu] = lu(A);
[X0,~] = qr(randn(n,p),0); % Set Initial point

switch test_instances
%% Algorithm: all set as default, no initial point
    case 1
    
        opts = [];  % no predefined parameter

        opts.X = [];
        opts.dim = [n,p];

        t1 = tic;
        [Out] = stop_slpg_smooth(@myfun,opts,A, Ll,Lu,alpha);
        t=toc(t1)


%% Algorithm: all set as default, with initial point
    case 2
        opts = [];  % no predefined parameter

        opts.X = X0;

        t1 = tic;
        [Out] = stop_slpg_smooth(@myfun,opts,A, Ll,Lu,alpha);
        t=toc(t1)


%% Algorithm: set all, with initial point
    case 3
    opts = [];  % no predefined parameter
    
    opts.X = X0;
    opts.info_warning = 1;
    opts.stepsize.type = 'ABB';
    opts.stepsize.max_stepsize = 1000;
    opts.stepsize.min_stepsize = 1e-10;
    opts.stepsize.init_stepsize = 1e-2;
    opts.gtol = 1e-8;
    opts.info = 0;
    opts.postprocess = 1;

    t1 = tic;
    [Out] = stop_slpg_smooth(@myfun,opts,A, Ll,Lu,alpha);
    t=toc(t1);
    
    
    case 4 
        L=gallery('tridiag',n,-1,2,-1);
        manifold = stiefelfactory(n,p);
        problem.M = manifold;
        [Ll,Lu] = lu(L);

        problem.cost = @(x, store) myfun_manopt(x, L, Ll,Lu,alpha,store);

        problem.egrad = @(x,store) mygrad_manopt(x, L, Ll,Lu,alpha,store);

        [x, ~] = qr(randn(n,p),0);

        options.maxiter = 1000;
        t1 = tic;
        [outs] =  stop_slpg_smooth_manopt2stop(problem);
        t=toc(t1)
        
        
end
%% Defin objective function
function [funX, F] = myfun(X, A, Ll,Lu,alpha)
    AX = A*X;
    rhoX = sum(X.^2, 2); 
    tempa = Lu\(Ll\rhoX); 
    tempa = alpha*tempa;
    F = A*X + bsxfun(@times,tempa,X);
    funX = 0.5*sum(sum(X.*(AX))) + 1/4*(rhoX'*tempa); 
end

function [fval, store] = myfun_manopt(X, L, Ll,Lu,alpha, store)
    LX = L*X;
    
    if ~isfield(store, 'tempa')
        rhoX = sum(X.^2, 2); % diag(X*X');
        store.tempa = alpha * Lu\(Ll\rhoX); 
    end
    
    tempa = store.tempa; 

    fval = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
end
    
function [g,store] = mygrad_manopt(X, L, Ll,Lu,alpha, store)

        if ~isfield(store, 'tempa')
            rhoX = sum(X.^2, 2); % diag(X*X');
            store.tempa = alpha * Lu\(Ll\rhoX); 
        end
        tempa = store.tempa;
        g = L*X + bsxfun(@times,tempa,X);
end