%> @file  SCPLearn_SplitSample.m
%> @brief Function that compares results of SCPLearn on split data samples 
%======================================================================
%> @brief It takes as input subject level time-series data. It randomly    
%> splits the data into two subsets. It returns
%> (1) a set of SCPs and a set of subject-level coefficients for both subsets
%> stored in variables B1, C1 and B2, C2.
%> (2) reproducibility (mean and standard deviation ) values
%> stored in variables reprod and reprod_std
%> (3) split-sample error stored in the variable CV_err
%> The above six outputs are saved in the output file specified as
%> the fourth argument of the function. 
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

function [] = SCPLearn_SplitSample(DataMatFile, K, lambda, outprefix,verbose,pruningThr)
%%
%> set up env
rng('shuffle')

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
ts_data=[];
for n=1:numel(ts)
    [D1,T] = size(ts{n});
    if(D1 ~= D)
        fprintf('Size of %d matrix is not same as rest\n',n)
    end
    ts_data = [ts_data,ts{n}];
    data(:,:,n) =  corrcoef(ts{n}');
end % end of for
data(isnan(data))=0;



if(~exist('sample_weights'))
    sample_weights = ones(N,1);
end

%%
%> random permutation of the data
rand_indices = randperm(N);

GroupIndex1 = rand_indices(1:floor(N/2));
GroupIndex2 = rand_indices(floor(N/2)+1:N);

%%
%> Training Run 1
if(K>D)
    [IDX, ~] = kmeans(ts_data, D);
    params.initdict = double(bsxfun(@eq,IDX,1:1:D)); 
    params.initdict = [params.initdict, params.initdict(:,1:K-D)];
else
    [IDX, ~] = kmeans(ts_data, K);
    params.initdict = double(bsxfun(@eq,IDX,1:1:K));
end
params.data = data(:,:,GroupIndex1);
params.iternum = 50;
params.lambda = lambda;
params.pruningThr = pruningThr;
params.sample_weights = sample_weights;
[B1,C1,~] = BlockSolverLowRankFrob(params,verbose);

%%
%> Training Run 2
if(K>D)
    [IDX, ~] = kmeans(ts_data, D);
    params.initdict = double(bsxfun(@eq,IDX,1:1:D)); 
    params.initdict = [params.initdict, params.initdict(:,1:K-D)];
else
    [IDX, ~] = kmeans(ts_data, K);
    params.initdict = double(bsxfun(@eq,IDX,1:1:K));
end
params.data = data(:,:,GroupIndex2);
params.iternum = 50;
params.lambda = lambda;
params.pruningThr = pruningThr;
params.sample_weights = sample_weights;
[B2,C2,~] = BlockSolverLowRankFrob(params,verbose);

%%
%> Testing Run 1
params.data = data(:,:,GroupIndex1);
params.initdict = B2;
params.iternum = 2;
params.lambda = lambda;
params.test = 1;
params.sample_weights = sample_weights;
[~,C1_test,err1_test] = BlockSolverLowRankFrob(params,verbose);

%%
%> Testing Run 2
params.data = data(:,:,GroupIndex2);
params.initdict = B1;
params.iternum = 2;
params.lambda = lambda;
params.test = 1;
params.sample_weights = sample_weights;
[~,C2_test,err2_test] = BlockSolverLowRankFrob(params,verbose);

%%
%> Computing relative error on test set
err1_test = err1_test(err1_test>0);
element_var = (length(GroupIndex1)-1)*var(data(:,:,GroupIndex1),[],3);
mean_sq_var_test1 = sum(element_var(:));
err1_test = err1_test/mean_sq_var_test1 ;

err2_test = err2_test(err2_test>0);
element_var = (length(GroupIndex2)-1)*var(data(:,:,GroupIndex2),[],3);
mean_sq_var_test2 = sum(element_var(:));
err2_test = err2_test/mean_sq_var_test2 ;

CV_err = (1/N)*(length(GroupIndex1)*err1_test(end)+length(GroupIndex2)*err2_test(end));

%%
%> Computing reproducibility between the two sets
[~,reprod,reprod_std] = CompareSCPs(B1, B2);

%%
%> saving all the results
OutputFile=[outprefix,'_SCPs.mat'];
save(OutputFile,'B1','C1','B2','C2','CV_err','reprod','reprod_std');

end % end of function
