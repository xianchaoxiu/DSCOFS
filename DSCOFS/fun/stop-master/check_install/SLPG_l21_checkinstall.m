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


%% Algorithm
opts = [];  % no predefined parameter
opts.X = []; % no initial point
opts.dim = [n,p]; % specify the dimension
opts.gamma = gamma;

opts.info = 0;

t1 = tic;
[Out] = stop_slpg_l21(@fun_l21PCA,opts,A,gamma);
t=toc(t1);


%% Defin objective function
function [fvalue,grad] = fun_l21PCA(X,A,gamma)
AX = A'*(A*X);
grad = -AX;
fvalue = 0.5*sum(sum(X.*grad)) + gamma* sum(sqrt(sum(X.*X,2)));
end