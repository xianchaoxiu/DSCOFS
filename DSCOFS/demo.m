clear

s1=100;% the number of selected features
options.maxiter   = 100; % the maxiter of DSCOFS
options.s1        = s1;% global sparsity (L2,0)
options.tol       = 1e-3;% stop criteria
options.i_x       = 1000; % the maxiter of pencf

data=load('Isolet'); 
name='W0_Isolet'; 
load("intialw.mat",name)
options.X0=eval(name);% initial solution
T=data.X;
data.X=(data.X-mean(data.X))';
data.Y=data.Y';
options.s2=0.1;% local elemental sparsity(L0)
options.mu1=10.^0;% penalty coefficient μ1
options.mu2=10.^6;% penalty coefficient μ2

indexs=zeros(s1,1);% indexs of selected features
[out] =DSCOFS(data, options);% solve
indexs=out.Tu;

