clear

data=load('Isolet'); 
load("W0_Isolet.mat")

s1=100;% the number of selected features
options.maxiter   = 100; % the maxiter of DSCOFS
options.s1        = s1;% global sparsity (L2,0)
options.tol       = 1e-3;% stop criteria
options.i_x       = 1000; % the maxiter of pencf
options.s2=0.1;% local elemental sparsity(L0)
options.mu1=10.^0;% penalty coefficient μ1
options.mu2=10.^6;% penalty coefficient μ2

options.X0=W0_Isolet;% initial solution

% options.X0=intialx(data.X,data.Y);% initial solution

data.X=(data.X-mean(data.X))';
data.Y=data.Y';
[out] =DSCOFS(data, options);% solve
indexs=out.Tu;% selected features

