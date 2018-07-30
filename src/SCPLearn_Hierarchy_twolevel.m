%> @file  SCPLearn_Hierarchy_twolevel.m
%> @brief Function that runs SCPLearn on time-series data with two hierarchical levels
%======================================================================
%> @brief It takes as input subject level time-series data. It returns
%> (1) a set of SCPs and a set of subject-level coefficients for the data
%> stored in variables B and C
%> (2) Time-series associated with each SCP, stored in B_ts
%> For details see the following paper:
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003
%>
%> @param DataMatFile Input mat file containing correlation matrices shaped
%> as [D,T,N], where D is the number of ROIs, T is the number of time-points,
%> N is the number of subjects
%> @param K Number of SCPs
%> @param lambda Sparsity constraint of SCPs, specified as a positive value
%> @param outprefix prefix for all output files (will be overwritten if exists)
%> @param verbose Integer value if 1 verbose messages will be output
%> @param pruningThr SCPs with inner-product overlap > pruningThr are
%discarded

%> @b Author: 
%> Harini Eavani
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
function [] = SCPLearn_Hierarchy_twolevel(DataMatFile, K, lambda, outprefix,verbose,pruningThr)
%> set up env: set the seed of random number generator to shuffle
rng('shuffle')
set(0,'RecursionLimit',20)

%%
%> load data
load(DataMatFile);
K = str2double(K);
lambda = str2double(lambda);
verbose = str2double(verbose);
pruningThr = str2double(pruningThr);

N = numel(ts);
D = size(ts{1},1);

fprintf('size of data is %d %d\n',D,N)
data = zeros(D,D,N);
for n=1:numel(ts)
    [D1,T] = size(ts{n});
    if(D1 ~= D)
        fprintf('Size of %d matrix is not same as rest\n',n)
    end
    data(:,:,n) =  corrcoef(ts{n}');
end % end of for
data(isnan(data))=0;
if(lambda == -1)
    lambda = D/8;
    fprintf('Sparsity is set to default %1.4f \n',lambda)
end
if(~exist('sample_weights'))
    sample_weights = ones(N,1);
end


%% level 1
%> specify level 1 input parameters
D = size(data,1);
if(K>D)
    multFact = ceil(K/D);
    params.initdict = repmat(eye(D),1,multFact);
    params.initdict = params.initdict(:,1:K);
else
    [IDX, ~] = kmeans(mean(data,3), K);    
    params.initdict = double(bsxfun(@eq,IDX,1:1:K));
end
params.sample_weights = sample_weights;
params.pruningThr = 1;
params.data = data;
params.iternum = 22;
params.lambda = lambda;
%> call block solver that alternately solves for B and C until convergence
%> call block solver that alternately solves for B and C until convergence
[B1,C1,~] = BlockSolverLowRankFrob(params,verbose);
K_pruned = size(B1,2);

%% level 2
%> specify level 2 input parameters
K_second = 50;
lambda_second = 2;
B = zeros(D,K_pruned+K_second^2);
C = zeros(K_pruned+K_second^2,N);
nodes = zeros(1,K_pruned+K_second^2);
B(:,1:K_pruned) = B1;
C(1:K_pruned,:) = C1;
count = K_pruned;
for kk=1:K_pruned
    K_second = 50;
    indices = abs(B(:,kk))>0.05;
    D_second = sum(indices);
    fprintf('Running SCPLearn on SCP %d, size %d\n',kk,D_second);
    data_new=zeros(D_second,D_second,N);
    for n=1:N
        data_new(:,:,n) = data(indices,indices,n);%.*abs(B(indices,kk)*B(indices,kk)');
    end
    data_new(isnan(data_new))=0;
    if(K_second>D_second)
        multFact = ceil(K_second/D_second);
        params.initdict = repmat(eye(D_second),1,multFact);
        params.initdict = params.initdict(:,1:K_second);
    else
        [IDX, ~] = kmeans(mean(data_new,3), K_second);
        params.initdict = double(bsxfun(@eq,IDX,1:1:K_second));
    end
    params.data = data_new;
    params.pruningThr = pruningThr;
    params.lambda = lambda_second;
    params.isPositive=1;
    %> call block solver that alternately solves for B and C until convergence
    [B2,C2,~] = BlockSolverLowRankFrob(params,verbose);
    K_second = size(B2,2);
    B(indices,count+1:count+K_second) = B2;
    C(count+1:count+K_second,:) = C2;
    nodes(count+1:count+K_second) = kk;
    count=count+K_second;
    fprintf('SCP %d split into %d secondary SCPs\n',kk,K_second);
end
%%
B = B(:,1:count);
C=C(1:count,:);
nodes =  nodes(1:count);
%%
%> saving all the results
OutputFile=[outprefix,'_SCPs.mat'];
save(OutputFile,'B','C','nodes');

end % end of function
