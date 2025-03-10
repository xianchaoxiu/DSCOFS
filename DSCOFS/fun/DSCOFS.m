  function [out] = DSCOFS(data, options)
%% Get Data Sequence
A = data.X;     % samples data
B = data.Y;     % label

%% Initialization
if isfield(options,'maxiter');   maxiter   = options.maxiter;   else; maxiter   = 1e3;  end
if isfield(options,'mu1');        mu1      = options.mu1;       else; mu1       = 100;  end
if isfield(options,'mu2');        mu2      = options.mu2;       else; mu2       = 100;  end
if isfield(options,'s1');         s1       = options.s1;        else; s1        = 100;  end
if isfield(options,'s2');         s2       = options.s2;        else; s2        = 1;    end
if isfield(options,'X0');        X0        = options.X0;        else; X0        = 0;    end
if isfield(options,'tol');       tol       = options.tol;       else; tol       = 1e-3; end
if isfield(options,'i_x');       i_x       = options.i_x;       else; i_x       = 100; end
[n,~]       = size(A);
ClassLabel  = unique(B);
p           = length(ClassLabel);
Y0          = X0  ;         % Initialize Y
Z0          = X0  ;         % Initialize Z

X         = X0;
Y_current = Y0;
Z_current = Z0;
elements = size(X,1)*size(X,2);
%% Defining Functions
AA=A*A';
f_loss  = @(X)(-1)*trace(X'*AA*X);   
Fnorm   = @(var)norm(var,'fro')^2;

%% Main
for iter = 1:maxiter                                                        
     %% update X 
    opts = [];  % no predefined parameter

    opts.X = X; % specific the initial point
    opts.info_warning = 0; % display costomized warning messages
    opts.stepsize.type = 'ABB'; % set stepsize as alternating Barzilai-Borwein stepsize
    opts.stepsize.max_stepsize = 1000; % specify the maximum stepsize to improve the robustness
    opts.stepsize.min_stepsize = 0; % specify the minimum stepsize to improve the robustness
    opts.stepsize.init_stepsize = 1e-2; % specify the initial stepsize
    opts.gtol = 1e-8; % set the stopping criteria
    opts.local_info = 0; % set the display mode
    opts.postprocess = 1; % turn on the post-process 
    opts.linesearch = 0; % no linsearch
    opts.maxit = i_x;
    opts.info = 0;
    X_current = X;
    [Out] = stop_pencf(@funch1,opts,Y_current,Z_current,mu1,mu2,AA,X_current);
    X = Out.X;        % update X
    %% update Y
    Y_old            = Y_current;
    W                = (X+0.001*Y_old)/(1+0.001);
    [~,Tu]           = maxk(sum(W.^2,2) , s1,'ComparisonMethod','abs'); 
    Y_current        = zeros(n,p);
    Y_current(Tu,:)  = X(Tu,:);  % update Y
   %% update Z
    Z_old            = Z_current;
    V                = (X+0.001*Z_old)/(1+0.001);
    Z1               = reshape(V,1,elements);
    [Zs,~]           = sort(abs(Z1),'descend');
    Z_S              = Zs(floor(s2*elements));
    Z_current        = X;
    Z_current(abs(Z_current)<Z_S)=0; % update Z
%%  check the stop criteria
    f_cost = f_loss(X)+mu1*(Fnorm(X-Y_current))+mu2*(Fnorm(X-Z_current));   
    obj(iter) = f_cost;
    if iter>1
        error_obj(iter) = abs(obj(iter)-obj(iter-1))/(1+abs(obj(iter-1)));      % error_obj
        error_Y(iter)   = Fnorm(Y_current-Y_old)/(1+Fnorm(Y_old));              % error_Y
        error_Z(iter)   = Fnorm(Z_current-Z_old)/(1+Fnorm(Z_old));              % error_Z
        fprintf('%5d\t  %6.2e\t %6.2e\t %6.2e\t %6.2e\n ',iter, obj(iter), error_obj(iter), error_Y(iter), error_Z(iter));
        if (error_obj(iter) <tol); break; end
        if iter == maxiter; fprintf('The number of iterations reaches maxiter.\n'); end
    end
end
out.Tu          = Tu;
  end
function [h_loss, h_grad] = funch1(X,Y,Z,mu1,mu2,AA,X_pre)
%% Calculate the gradient and function values of h                                                                
        h_grad  = -2*AA*X+2*mu1*(X-Y)+2*mu2*(X-Z)+2*0.001*(X-X_pre);
        h_loss =(-1)*trace(X'*AA*X) + mu1*(norm(X-Y,'fro')^2)+mu2*(norm(X-Z,'fro')^2)+0.001*(norm(X-X_pre,'fro')^2);                                 
end