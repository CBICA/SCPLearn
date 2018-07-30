%> @file  BlockSolverLowRankFrob.m
%> @brief Function that iteratively calls Bsolver, Csolver and ESolver
%======================================================================
%> @brief It takes as input a structure containing the input data,
%> initilization values, user parameters K and lambda \n
%> (1) a set of SCPs and a set of subject-level coefficients for the data
%> stored in variables B and C \n
%> (2) error value at each iteration, stored in a variables err \n
%> For details see the following paper:
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003
%>
%> @param params Structure containing the mandatory fields:
%>    'data' - Training data of size D x T x N
%>    'initdict' -  Initial dictionary of size D x K 
%>    'lambda' - Sparsity level, a value between 0 and 1
%>
%>     Optional parameters
%>    ------------------------
%>    'iternum' - Number of training iterations. Default: 50
%>    'test' - 0/1; If 1, Calculates C for fixed B. Default: 0 
%> @param verbose integer value if 1 provides verbose messages
%> @retval B SCPs of size D x K
%> @retval C SCP Coefficients of size K x N 
%> @retval err objective function of size iternum X 1 
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
function [B,C,err] = BlockSolverLowRankFrob(params,verbose)
%%
%> read data else error
if (isfield(params,'data'))
    data = params.data;
else
    error('No input data\n')
end% end if

%> read initialized values else error
if (isfield(params,'initdict'))
    B0 = params.initdict;
    K = size(B0,2);
else
    error('SCPs not initialized\n')
end% end if

%> if isPositive not existent make it zero
if (~isfield(params,'isPositive'))
    params.isPositive=0;
end% end if

%> read number of iterations, else default 50
if (isfield(params,'iternum'))
    iternum = params.iternum;
else
    iternum = 51;
end% end if

%%
%> Initialize C variable with random positive values
[~, D, N] = size(data);
C0 = rand(K, N);

%> Set solver options for Csolver
Coptions.BBMethod = 3;
Coptions.K = K;
Coptions.N = N;
Coptions.D = D;
Coptions.verbose=verbose;


%> Set solver options for Bsolver
Boptions.BBMethod = 3;
Boptions.lambda = params.lambda;
Boptions.K = K;
Boptions.N = N;
Boptions.D = D;
Boptions.verbose=verbose;
Boptions.isPositive = params.isPositive;
err = zeros(1,iternum);

B =  B0 ;
C =  C0 ;
E = rand(D,N);
iter = 1;

%%
%> Check if this is only for test (default it is not)
if (isfield(params,'test'))
    if(params.test==1)
        C  = CSolver_lowrank_frob(C0,data,params.sample_weights,B0,E,Coptions) ;
        B = B0;
        err = ComputeObjective(data,params.sample_weights,B,C,E);
        return
    end% end if
end% end if

%%
%> Start Main Loop
while(iter <= iternum)
    
    if mod(iter,3) == 2    
        %> optimization wrt to B
        blockName = 'B' ;
        B0 = B ;
        B  = BSolver_lowrank_frob(B0,data,params.sample_weights,C,E,Boptions) ;
        %> re-normalize the basis (warning:very ad-hoc but common practice)
        [B,C] = renormalize(B,C);
    end
    
    if mod(iter,3) == 1      
        %> optimization wrt to C
        blockName = 'C' ;
        C0 = C ;
        C  = CSolver_lowrank_frob(C0,data,params.sample_weights,B,E,Coptions) ;
    end % end if 
    
    if mod(iter,3) == 0   
        %> optimization wrt to D
        blockName = 'E' ;
        E  = ESolver_lowrank_frob(B,C,data) ;        
    end
   
    curTime = datestr(now) ;
    %> Calculate objective (error)
    err(iter) = ComputeObjective(data,params.sample_weights,B,C,E);
    %> Print update
    if(verbose)
        fprintf(1,'Block %s, (%s)-Iteration %d : Error = %3.2f\n',blockName, curTime,iter,err(iter));
    end

    iter = iter + 1;
end% end for loop
%%
%> sort based on coeff sum
[~,sorted_ind]=sort(sum(C,2),'descend');
B = B(:,sorted_ind);
C = C(sorted_ind,:);

%%
B_norm = bsxfun(@times,abs(B),1./sqrt(diag(B'*B)'));
overlap_coeffs = abs(B_norm'*B_norm);

% starting with the least ranked prune 
scps_to_remove=[];
for kk=K:-1:1
    overlap_coeffs(scps_to_remove,:) = 0;
    overlap_coeffs(:,scps_to_remove) = 0;
    kk_overlap = overlap_coeffs(kk,:);
    high_overlaps = find(abs(kk_overlap)>params.pruningThr);
    num_overlaps = length(high_overlaps)-1;
    if(num_overlaps>=1)
        scps_to_remove = [scps_to_remove, kk];
    end
end

B = B(:,setdiff(1:K,scps_to_remove));
C = C(setdiff(1:K,scps_to_remove),:);

end % end function 

function f = ComputeObjective(data,sample_weights,B,C,E)

N =size(data,3);
ObjFun = @(X, A) norm(A-X,'fro')^2 ;

f=0;
for i=1:N
    Xi = B*diag(C(:,i))*B' + diag(E(:,i)); 
    f = f + sample_weights(i)*ObjFun(Xi, data(:,:,i));
end

end% end function 

function [Bn,Cn] = renormalize(B,C)
    scaling = max(abs(B));
    scaling(scaling==0)=1;%make sure scaling factor is not zero
    Bn = B*diag(1./scaling);
    Cn = diag(scaling.^2)*C;
end% end function 
