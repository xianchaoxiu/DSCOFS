clear;
test_instances = 3;

%% Problem setting
n = 1000;  
p = 20;
A = randn(n); 
A = 0.5*(A + A');
B = randn(n,p);
normA = norm(A, 2);
normB = norm(B, 2);
A = A / normA;
B = B / normB / 10;

[X0, ~] = qr(randn(n, p), 0); % Set Initial point

switch test_instances
%% Algorithm: all set as default, no initial point
    case 1
    
        opts = [];  % no predefined parameter
        
        opts.X = [];
        opts.dim = [n, p];

        t1 = tic;
        [Out] = stop_pcal(@myfun, opts, A, B);
        t = toc(t1);


%% Algorithm: all set as default, with initial point
    case 2
        
        opts = [];  % no predefined parameter
        
        opts.X = X0;
        
        t1 = tic;
        [Out] = stop_pcal(@myfun, opts, A, B);        
        t = toc(t1);


%% Algorithm: set all, with initial point
    case 3
        opts = [];
        
        opts.X = X0;
        opts.info_warning = 1;
        opts.stepsize.type = 'ABB';
        opts.solver = 'GP';
        opts.stepsize.max_stepsize = 1000;
        opts.stepsize.min_stepsize = 1e-10;
        opts.stepsize.init_stepsize = 1e-2;
        opts.gtol = 1e-3;
        opts.xtol = 1e-6;
        opts.info = 0;
        opts.linesearch = 0;

        t1 = tic;
        [Out] = stop_pcal(@myfun, opts, A, B);
        t = toc(t1);
    
    
    case 4 
        manifold = stiefelfactory(n, p);
        problem.M = manifold;
        opts.local_info = 0;
        problem.cost = @(x, store) myfun_manopt(x, A, B, store);

        problem.egrad = @(x,store) mygrad_manopt(x, A, B, store);

        [x, ~] = qr(randn(n, p), 0);

        options.maxiter = 1000;
        
        t1 = tic;
        [Out] =  stop_pcal_manopt2stop(problem);
        t = toc(t1);
        
        
end

%% Define objective function
function [fval, grad] = myfun(X, A, B)
    AX = A * X;
    grad = AX + B;
    fval = 0.5 * sum(sum(X .* grad));
end

function [fval, store] = myfun_manopt(X, A, B, store)
    
    if ~isfield(store, 'tempa')
        AX = A * X;
        grad = AX + B;
        store.tempa = grad; 
    end
    
    grad = store.tempa;

    fval = 0.5 * sum(sum(X .* grad));
end
    
function [grad, store] = mygrad_manopt(X, A, B, store)

    if ~isfield(store, 'tempa')
        AX = A*X;
        grad = AX + B;
        store.tempa = grad; 
    end
    
    grad = store.tempa;
end