%> @file  BSolver_lowrank_frob.m
%> @brief Function that solves for the B variable
%======================================================================
%> @brief It takes as input the initial value of B, current value of C,
%> correlation matrices and the solver options. It returns the current
%> solution of B for the given C and data. \n
%> For details see the following paper: \n
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003 \n
%>
%> @param B0 The initial value of B of size D x K
%> @param data The correltaion matrices of size D x D x N
%> @param sample_weights The (optional) subject-wise weights
%> @param C The current set of coefficients K x N
%> @param E The diagonal estimates
%> @param options Structure of solver options containing fields: 
%>               BBMethod = 3 (default)
%>               lambda = sparsity level (set by Block Solver)
%>               K =  number of SCPs (set by Block Solver)
%>               N = number of subjects (set by Block Solver)
%>               D = number of ROIs (set by Block Solver)
%> @retval B SCPs of size D x K
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
function [B]  = BSolver_lowrank_frob(B0,data,sample_weights,C,E,options)
%%
%> Initialize
D = size(B0,1) ;
K = size(B0,2) ;

options.D = D;
options.K = K;
lambda = options.lambda ;
BBMethod = options.BBMethod ;

%> remove diag vals
[~,N] = size(C);
for nn=1:N
    data(:,:,nn) = data(:,:,nn) - diag(E(:,nn));
end

%%
%> solving the problem using BBSolver
f = @(x)obj_fun(x,data,sample_weights,C,options) ;
g = @(x)obj_gradient(x,data,sample_weights,C,options) ;
projFcn = @(x)ProjectionBasisVectors(x,lambda,D,K,options.isPositive) ;
x0  = projFcn(B0(:)) ;
sol = spg_05262015(x0, f, g, BBMethod, projFcn,-inf,inf,options.verbose) ;
B = reshape(sol.par,D,K) ;
end% end of function

function f = obj_fun(B, data,sample_weights,C,options) 
K = options.K ;
D = options.D ;
B = reshape(B,D,K) ;
N = size(data,3);

ObjFun = @(X, A) norm(A-X,'fro')^2;

f=0;
for i=1:N
    Xi = B*diag(C(:,i))*B' ; 
    f = f + sample_weights(i)*ObjFun(Xi, data(:,:,i));
end

end% end function 

function g = obj_gradient(B,data,sample_weights,C,options) 

K = options.K ;
D = options.D ;
B = reshape(B,D,K) ;
N = size(data,3);


g = zeros(size(B));

for i=1:N
         g =  g -  4*sample_weights(i)*(data(:,:,i)-B*diag(C(:,i))*B')*B*diag(C(:,i));
 end
g = g(:);

end% end function 

function bp = ProjectionBasisVectors(b,lambda,D,K,isPositive)
    B = reshape(b,D,K) ;
    Bp = zeros(D,K) ;
    
    for ii=1:K
        % this is for individual basis sparsity control
        Bp(:,ii) = ProjectionOnUnitBoxSimplex(abs(B(:,ii)),lambda);
        
        if(isPositive==0)
            % append sign 
            Bp(:,ii) = sign(B(:,ii)).*Bp(:,ii);
        end
    end
    bp = Bp(:) ;  
end% end function 

