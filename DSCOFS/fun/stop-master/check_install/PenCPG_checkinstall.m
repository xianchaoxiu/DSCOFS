clear all

%% Problem setting
n = 1000;
p = 10;
gamma = 0.045;

% Generate data matrix
m = 50;
A = randn(m,n);
A = A - repmat(mean(A,1),m,1);
scale = [];
for i = 1:n
    scale = [scale norm(A(:,i))];
end
scale = max(scale);
A =A/scale;
A1 = A;   AtA = A'*A;
ddd = eigs(AtA,1);
AtA = AtA/ddd;
A = A/sqrt(ddd);
[X0,~] = qr(randn(n,p),0); % Set Initial point


%% Set options
opts = [];  % no predefined parameter
    
opts.X = X0; % specific the initial point
opts.info_warning = 1; % display costomized warning messages
opts.stepsize.type = 'ABB'; % set stepsize as alternating Barzilai-Borwein stepsize
opts.stepsize.max_stepsize = 1000; % specify the maximum stepsize to improve the robustness
opts.stepsize.min_stepsize = 0; % specify the minimum stepsize to improve the robustness
opts.stepsize.init_stepsize = 1e-2; % specify the initial stepsize
opts.gtol = 1e-4; % set the stopping criteria
opts.info = 0; % set the display mode
opts.postprocess = 1; % turn on the post-process 
opts.linesearch = 0; % no linsearch
opts.gamma = gamma;


%% Algorithm

tic;
[Out] = stop_pencpg(@fun_l21PCA,opts,A,gamma);
t=toc;


%% Defin objective function
function [fvalue,grad] = fun_l21PCA(X,A,gamma)
AX = A'*(A*X);
grad = -AX;
fvalue = 0.5*sum(sum(X.*grad)) + gamma* sum(sqrt(sum(X.*X,2)));
end