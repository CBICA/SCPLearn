%> @file  CSolver_lowrank_frob.m
%> @brief Function that solves for the C variable
%======================================================================
%> @brief It takes as input the initial value of C, current value of B,
%> correlation matrices and the solver options. It returns the current
%> solution of C for the given B and data. \n
%> For details see the following paper: \n
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003 \n
%>
%> @param C0 The initial value of B of size D x K
%> @param data The correltaion matrices of size D x D x N
%> @param sample_weights The (optional) subject weights of size N X 1 
%> @param B The current set of coefficients K x N
%> @param E The diagonal estimates of size K x N
%> @param options Structure of solver options containing fields: 
%>               BBMethod = 3 (default)
%>               lambda = sparsity level (set by Block Solver)
%>               K =  number of SCPs (set by Block Solver)
%>               N = number of subjects (set by Block Solver)
%>               D = number of ROIs (set by Block Solver)
%> @retval C SCPs of size D x K
%>
%> @b Author: 
%> Harini Eavani
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
function [C]  = CSolver_lowrank_frob(C0,data,sample_weights,B,E,options)
%
%> Initialize
N = size(C0,2) ;
K = size(C0,1) ;

options.N = N;
options.K = K;
BBMethod = options.BBMethod ;

%> remove diag vals
[~,N] = size(E);
for nn=1:N
    data(:,:,nn) = data(:,:,nn) - diag(E(:,nn));
end

%%
%> Solve using spg
f = @(x)obj_fun(x,data,sample_weights,B,options) ;
g = @(x)obj_gradient(x,data,sample_weights,B,options) ;
projFcn = @(x) x.*(x>0);
x0  = projFcn(C0(:)) ;
sol = spg_05262015(double(x0), f, g, BBMethod, projFcn,[],[],options.verbose) ;
C = reshape(sol.par,K,N) ;

end

function f = obj_fun(C,data,sample_weights,B,options)

K = options.K ;
N = options.N ;
C = reshape(C,K,N) ;


ObjFun = @(X, A) norm(A-X,'fro')^2 ;

f=0;
for i=1:N
    Xi = B*diag(C(:,i))*B' ; 
    f = f + sample_weights(i)*ObjFun(Xi, data(:,:,i));
end


end

function g = obj_gradient(C,data,sample_weights,B,options)

K = options.K ;
N = options.N ;
C = reshape(C,K,N) ;

g = zeros(size(C));

for i=1:N
         g(:,i) =  g(:,i) -  sample_weights(i)*diag(2*B'*(data(:,:,i)-B*diag(C(:,i))*B')*B) ;
 end

g = g(:);

end

