function [output] = stop_pencf(fun,opts,varargin)
% First-order method for smooth optimization over the Stiefel manifold based on PenC model (PenCF)
% function [output] = stop_pencf(fun,opts)
% function [output] = stop_pencf(fun,opts,varargin)

% The descriptions for stop_pencf solver could also be found at https://stmopt.gitee.io/algorithm_description/PenCF_code.html
% This solver is a part of stmopt package (https://stmopt.gitee.io/)

% The stop_pencf solver starts at opts.X if it is provided (otherwise, generates opts.X whose size is specified in
% opts.dim on the Stiefel manifold). 
% The fun in stop_pencf solver is called by
% [fval, grad] = fun(X) or [fval, grad] = fun(X, varargin)

% The opts in stop_pencf solver is a matlab structure for overwritting the default values. 
% In opts structure, either opts.X or opts.dim should be specified. 
% Any other option have a default value and is hence optional. To assign the value to the opts structure,
% one needs to pass an options structure with a field options.optionname.

% * * opts.X
%     * The initial point for the solver. opts.X should be a matrix with $n$ rows and $p$ columns, where $n\geq p$. We recommend to choose a feasible opts.X in practice. 
%     * When opts.X is empty, the option opts.dim should be specified, and the solver will randomly generate a feasible initial point by [X, ~] = qr(randn(n,p),0); otherwise, the solver will raise an error. 
%   * opts.proj
%     * Projection to a set that contains Stiefel manifold
%     * The default value is set as opts.proj = @(X)  X * min(1.0001*(sqrt(p)/norm(X,'fro')), 1) ;. This is equivalent to project $X$ into a ball with radius $1.001 \sqrt{p}$.
%   * opts.stepsize
%     *  Stepsize arguments for PenCF
%     *  opts.stepsize.type = 0:  
%        *  opts.stepsize.type = 'ABB': Alternating BB stepsize
%        *  opts.stepsize.type = 'BB1': BB1 stepsize
%        *  opts.stepsize.type = 'BB2': BB2 stepsize
%        *  opts.stepsize.type = 'FIX': fixed stepsize
%     *  opts.stepsize.max_stepsize = 1000: Upper-bound for stepsize
%     *  opts.stepsize.min_stepsize = 0: Lower-bound for stepsize
%     *  opts.stepsize.init_stepsize: Initial stepsize
%   * opts.gtol = 1e-8
%     * Tolerance for the kkt condition
%     * PenCF terminates when $\frac{||\nabla f(X_k) - X_k \Lambda(X_k) ||_F}{||\nabla f(X_0) - X_0 \Lambda(X_0) ||_F} \leq gtol$
%   * opts.maxit = 1000
%     * Maximum iteration
%     * PenCF terminates when the iteration exceeds maxit
%   * opts.penalparam
%     * Penalty parameter 
%     * The penalty parameter is important to the performance of PenCF. A small penalty parameter leads to diverge while a large one results in slow convergence rate. 
%     * The default value for opts.penalparam is chosen as $0.001 \cdot ||\nabla f(X_0)||_F$.
%   * opts.safe_proj
%     * Option for safeguards
%     * safe_proj = 1 is a function that pulls $X_k$ close to Stiefel manifold if it is a little far away from Stiefel manifold. It improves the robustness of PenCF, but takes additional computational costs. 
%   * opts.postorth = 1
%     * Option for orthonormalization as postprocess:
%     * The postprocess improves the feasibility but requires additional computational costs. We recommand to set  opts.postorth = 1 when high precision in feasibility is required. 
%   * opts.info = 0
%     * Option for display during the execution
%     * No print or display when opts.info = 0
%   * opts.info_warning = 0
%     * Options for displaying warnings in the solver. 




% All the outputs of stop_pencf is wrapped in the output structure. 

% out: output information

% * out.X
%   * Final solution
% * out.iter
%   * Total iterations taken by PenCF
% * out.fval
%   * Final function value
% * out.kkt
%   * Final substationarity evaluated by $||\nabla f(X_k) - X_k \Lambda(X_k)||_F$
% * out.fea
%   * Final feasibility evaluated by $||X_k^\top X_k - I_p||_F$ 
% * out.fvals
%   * Records of function values in each iteration
% * out.kkts
%   * Record of the substationarities in each iteration
% * out.feas
%   * Record of the feasibility in each iteration
% * out.times
%   * Records of the CPU time in each iteration. 

% The implementation of this solver is based on the following paper:

% N. Xiao, X. Liu and Y. Yuan, A Class of Smooth Exact Penalty Function Methods for Optimization Problems with Orthogonality Constraints, Optimization Methods and Software 
% (DOI:10.1080/10556788.2020.1852236). 
% https://www.tandfonline.com/doi/abs/10.1080/10556788.2020.1852236?journalCode=goms20
% 


%% Warnings 
if ~isfield(opts, 'info_warning')
    info_warning = 0;
else
    info_warning = opts.info_warning;
end

%% Specify the initial point
if isfield(opts, 'X') && (~isempty(opts.X))
    [n,p] = size(opts.X);
    X = opts.X;
elseif isfield(opts, 'dim')
    n = opts.dim(1);
    p = opts.dim(2); % 检查dim的size
    if n < p
        error('Column size should be no less than the size of row. Please check the dimensions of the opts.dim.') % 否定句, n\& p的具体含义
    end
    [X, ~] = qr(randn(n,p),0);
else
    error('No initial point, and the size of initial point is not specified')
end
    
% X = Xinit;

if ~isfield(opts, 'maxit')
    maxlength = 3000;
else
    maxlength = opts.maxit;
end




%% Specify the penalty parameter according to ||\nabla f(X_0)||_F
[fvalue, gradf] = feval(fun,X,varargin{:});
s_upper = norm(gradf,'fro');

if ~isfield(opts, 'penalparam')
    beta = max(min(0.5*s_upper/p, 1e10), 1e-10); % set lower-bound for beta
else
    beta = opts.penalparam;  % warning for too large/small beta
    if info_warning
        if beta > 1e10
            warning('The penalty parameter is too large. Please scale the problem or specify opts.penalparam.');
        elseif beta < 1e-10
            warning('The penalty parameter is too small. Please scale the problem or specify opts.penalparam.');
        end
    end
end



Lambda = X'*gradf;
Lambda = 0.5*(Lambda+Lambda');
XX = X'*X;
Grad = gradf - X*(Lambda  - beta * ( XX- eye(p))  )  ; % approximated gradient
Grad_tmp = gradf - X * Lambda;






%% Specify the options for stepsize
if isfield(opts,'stepsize') 
    if ~isfield(opts.stepsize,'type')
        stepsize_type = 0;
        if info_warning
            warning('No input stepsize type, choose alternating BB stepsize. ')
        end
    else
        if strcmp(opts.stepsize.type, 'ABB')
            stepsize_type = 0;   
        elseif strcmp(opts.stepsize.type, 'BB1')
            stepsize_type = 1;
        elseif strcmp(opts.stepsize.type, 'BB2')
            stepsize_type = 2;
        elseif strcmp(opts.stepsize.type, 'SD') % initial stepsize by difference
            stepsize_type = 3;
        else
            stepsize_type = 0; 
            if info_warning
                warning('Invaild stepsize type, choose ABB instead. ')
            end
        end
        
    end
    
    if ~isfield(opts.stepsize,'max_stepsize')
        max_stepsize = 1000;
    else
        max_stepsize = opts.stepsize.max_stepsize;
    end

    if ~isfield(opts.stepsize,'min_stepsize')
        min_stepsize = 0;
    else
        min_stepsize = opts.stepsize.min_stepsize;
    end
    
    
    if ~isfield(opts.stepsize,'init_stepsize')
        init_stepsize =1/( norm(X,'fro')/2*norm(Grad,'fro') + p * beta);
    else
        init_stepsize = opts.stepsize.init_stepsize;
    end
else
    if info_warning 
        warning('No input stepsize. ')
    end
    stepsize_type = 0;
    max_stepsize = 1000;
    min_stepsize = 0;
    init_stepsize = 1/( norm(X,'fro')/2*norm(Grad,'fro') + p* beta);
end

%% Options to projection to keep feasibility
if ~isfield(opts,'proj')
    proj = @(X)  X * min(1.0001*(sqrt(p)/norm(X,'fro')), 1) ;
else
    proj = opts.proj;
end



%% Options for tolerence

if ~isfield(opts,'gtol') % add gtol
    gtol = 1e-10;
else
    gtol = opts.gtol;
end



if ~isfield(opts,'xtol') % add xtol
    xtol = 1e-15;
else
    xtol = opts.xtol;
end

%% Options for printed infomation
if ~isfield(opts,'info')
    local_info = 2;
else
    local_info = opts.info;
end


%% Options for post-process
if ~isfield(opts,'postorth')
    postprocess = 1;
else
    postprocess = opts.postorth;
end

%% Options for line search

if ~isfield(opts,'linesearch')
    linesearch = 1;
else
    linesearch = opts.linesearch;
end


%% Initialization

feas_init = norm(XX- eye(p),'fro');
kkts_init =  norm(Grad_tmp,'fro');
fval_init = fvalue;

violation_KKT = [kkts_init];
violation_Constraint = [feas_init];
fval = [fval_init];

line_search_record = [0];

time_record = [0];
%t_start = tic;



if local_info
    fprintf('%4s | %12s | %10s | %10s | %8s | %8s\n', 'Iter ', 'F(X) ', 'KKT ', 'Xerr ', 'Fea ', 'stepsize');
    fprintf('%04d \t %1.5e \t %3.2e \t %3.2e \t %3.2e \t %3.2e\n',0, fval_init, kkts_init, 0, feas_init, 0);
end

%% Main loop
for jj = 2:1:maxlength
    X_pre = X;

    % Computing stepsize
    if stepsize_type ~= 3
        if jj < 3 
            stepsize = init_stepsize;
        else
            stepsize_flag = mod(jj,2);
            if stepsize_type == 2
                stepsize_flag = 0;
            elseif stepsize_type == 1
                stepsize_flag = 1;
            end
            if stepsize_flag == 0
                stepsize =  abs(norm(D,'fro')^2 /sum(sum(Grad_Y.*D)));
            else
                stepsize = abs( sum(sum(Grad_Y.*D))/ norm(Grad_Y,'fro')^2);
            end
        end
    else
        if jj < 3
            stepsize = init_stepsize;
        else
            stepsize = 2*abs(fval(end-1)-fval(end))/violation_KKT(end)^2;
        end
    end
    stepsize = max(min_stepsize, min(abs(stepsize),max_stepsize));
    
    
    if linesearch == 0
        X = X - stepsize*Grad; % Update X
        ls_iter = 0;
    else
        local_fval_grad = @(ZZ) feval(fun,ZZ,varargin{:});
        if stepsize_type ~= 3
            [X, ls_iter] = PenCF_lineasarch(X, Grad, stepsize,local_fval_grad,fval, Lambda, XX, beta);
        else
            [X, ls_iter] = PenCF_lineasarch_sd(X, Grad, stepsize,local_fval_grad,fval, Lambda, XX, beta);
        end
    end
    line_search_record = [line_search_record, ls_iter];
    
    
%     X = proj(X);
    X = safe_proj(proj, X, stepsize*Lambda, violation_Constraint(end));
    
    XX = X'*X;
    
    
    D = X-X_pre;
    
    

    [fvalue,gradf] = feval(fun,X,varargin{:});
    Lambda = X'*gradf;
    Lambda = 0.5*(Lambda+Lambda');
    Grad_p = Grad;
    
    Grad_tmp = gradf - X * Lambda;
    
    Grad = gradf - X*  (Lambda  - beta * (XX - eye(p))  );

    
    Grad_Y = Grad - Grad_p;

    violation_KKT =  [violation_KKT,norm(Grad_tmp,'fro')];
    violation_Constraint = [violation_Constraint, norm(XX- eye(p),'fro')];
    %time_record = [time_record, toc(t_start)];
    fval = [fval,fvalue - 0.5* sum(sum((XX- eye(p)).* Lambda)) + 0.25 * beta * norm(XX- eye(p),'fro')^2];
    
    if local_info && mod(jj,20) == 0 
        fprintf('%04d \t %1.5e \t %3.2e \t %3.2e \t %3.2e \t %3.2e\n',jj, fvalue, violation_KKT(jj), norm(D,'fro'), violation_Constraint(jj), stepsize);
    end
    
    
    if violation_KKT(jj) < gtol * kkts_init 
        if local_info
            fprintf('%04d \t %1.5e \t %3.2e \t %3.2e \t %3.2e \t %3.2e\n',jj, fvalue, violation_KKT(jj), norm(D,'fro'), violation_Constraint(jj), stepsize);
            fprintf('The norm of gradient is less than tolerence \n')
        end
        break;
    end
    
    
    
    if norm(D,'fro') < xtol 
        if local_info
            fprintf('%04d \t %1.5e \t %3.2e \t %3.2e \t %3.2e \t %3.2e\n',jj, fvalue, violation_KKT(jj), norm(D,'fro'), violation_Constraint(jj), stepsize);
            fprintf('The norm of X_{k} - X_{k-1} is less than tolerence \n')
        end
        
        break;
    end
    
    

    if violation_KKT(jj) > 1e12
        warning('The algoorithm may diverge \n');
        break;
    end
    
end

if postprocess 
    if local_info 
        fprintf('Postprocess start \n')
    end
    X = post_orth(X);
    [fvalue,gradf] = feval(fun,X,varargin{:});
    Lambda = X'*gradf;
    Lambda = 0.5*(Lambda+Lambda');
    Grad = gradf - X*  (Lambda  - beta * (XX - eye(p))  );


    violation_KKT(jj) =  norm(Grad,'fro');
    violation_Constraint(jj) = norm(X'*X- eye(p),'fro');
    
    if local_info
        fprintf('%04d \t %1.5e \t %3.2e \t %3.2e \t %3.2e \t %3.2e\n',jj, fvalue, violation_KKT(jj), 0, violation_Constraint(jj), 0);
    end
end

%% Output
output.X = X;
output.iter = jj;

output.kkt = violation_KKT(end);
output.fval = fvalue;
output.fea = violation_Constraint(end);

output.kkts = violation_KKT;
output.fvals = fval;
output.feas = violation_Constraint;
%output.times = time_record;

output.line_search_iters = line_search_record;
end



function [out,iter] = PenCF_lineasarch(X,Grad, stepsize,fval_grad,fval_records,Lambda,XX, beta )
if isempty(fval_records)
    out = X - stepsize * Grad;
    iter = 0;
else
    stepsize_local = stepsize;
    [n,p] = size(X);
    for jj = 0:1:2
        Y = X - stepsize_local * Grad;
        [fval,~] = fval_grad(Y);
        fval_pen = fval - 0.5* sum(sum((XX- eye(p)).* Lambda)) + 0.25 * beta * norm(XX- eye(p),'fro')^2;
        if fval_pen <= max(fval_records(max(1,end-50 ):end))
            
            break
        else
            stepsize_local = stepsize_local * 0.5;
        end
        
    end
%     jj
    iter = jj;
    out = Y;
    
end
end




function [out] = post_orth(X)
[n,p] = size(X);

if p > 0.5 * n
    [U, ~, V] = svd(X,0);
    Y = U * V';
else 
    [V, D] = eig(X'*X);
    Y = X*(V*diag(sqrt(1./diag(D) ))*V');
end

[U, ~, V] = svd(X,0);
Z = U * V';



out = Y;
end


function [X] = safe_proj(proj, X, G, tol_feas)
% Sate projection
% When ||X^TX - I_p||_F is large,  project X close to Stiefel manifold


    [n,p] = size(X);
    X_tmp = X;
    
    esti_feas  = (norm(G,'fro')+1) * tol_feas;
    
    
    if esti_feas > 0.5
        XX_tmp = X_tmp'*X_tmp;
        X_tmp = 2*X_tmp / (XX_tmp + eye(p));
    elseif esti_feas > 0.1
        XX_tmp = X_tmp'*X_tmp;
        X_tmp = X_tmp * ( 1.5 * eye(p) - 0.5 *XX_tmp);
    end
    
    X = feval(proj, X_tmp);
    
end


