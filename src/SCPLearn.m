%> @file  SCPLearn.m
%> @brief Function that runs SCPLearn on time-series data 
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
function [] = SCPLearn(DataMatFile, K, lambda, outprefix,verbose,pruningThr)
%%
%> set up env: set the seed of random number generator to shuffle
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
if(lambda == -1)
    lambda = D/8;
    fprintf('Sparsity is set to default %1.4f \n',lambda)
end
%%
%> specify input parameters
if(K>D)
    multFact = ceil(K/D);
    params.initdict = repmat(eye(D),1,multFact);
    params.initdict = params.initdict(:,1:K);
else
    [IDX, ~] = kmeans(ts_data, K);
    params.initdict = double(bsxfun(@eq,IDX,1:1:K));
end
params.data = data;
params.iternum = 52;
params.lambda = lambda;
params.pruningThr = 1;
if(exist('sample_weights'))
    params.sample_weights = sample_weights;
else
    params.sample_weights = ones(N,1);
end

%> call block solver that alternately solves for B and C until convergence
[B,C,err] = BlockSolverLowRankFrob(params,verbose);
%> update number of SCPs if it is pruned
K = size(B,2);

%%
%> Get SCP time-series
B_ts = {};
for n=1:N
    B_ts{n} = B\ts{n};
end % end of for

%%
%> saving all the results
save(outprefix,'B','C','B_ts','err');

end % end of function
