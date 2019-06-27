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
function [] = SCPLearn(DataMatFile, K, lambda, outprefix,verbose,pruningThr,levels)
%%

%%
%> load data
load(DataMatFile);
K = str2double(K);
lambda = str2double(lambda);
verbose = str2double(verbose);
pruningThr = str2double(pruningThr);
levels = str2double(levels);

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
ts_data(isnan(ts_data))=0;

if(lambda == -1)
    lambda = D/10;
    fprintf('Sparsity is set to default %1.4f \n',lambda)
end
%%
%> specify input parameters
if(K>D)
    multFact = ceil(K/D);
    params.initdict = repmat(eye(D),1,multFact);
    params.initdict = params.initdict(:,1:K);
else
    [IDX, ~] = kmeans(ts_data,K,'Distance','correlation','Replicates',10);  
    params.initdict = double(bsxfun(@eq,IDX,1:1:K));
end
params.data = data;
params.iternum = 51;
params.lambda = lambda;
params.pruningThr = 1;
if(exist('sample_weights','var'))
    params.sample_weights = sample_weights;
else
    params.sample_weights = ones(N,1);
end

%> call block solver that alternately solves for B and C until convergence
[B1,C1,~] = BlockSolverLowRankFrob(params,verbose);
K_pruned = size(B1,2);

%% level 2
if (levels > 0) %% if hierarchical learning is specified
	%> specify level 2 input parameters
	K_second = levels;
	B = zeros(D,K_pruned+K_second*K_pruned);
	C = zeros(K_pruned+K_second*K_pruned,N);
	nodes = zeros(1,K_pruned+K_second*K_pruned);
	B(:,1:K_pruned) = B1;
	C(1:K_pruned,:) = C1;
	count = K_pruned;
	for kk=1:K_pruned
		K_second = levels;
		indices = abs(B(:,kk))>0.3;
		D_second = sum(indices);
		lambda_second = 2;
		fprintf('Running SCPLearn on SCP %d\n',kk);
		data_new=zeros(D_second,D_second,N);
		for n=1:N
		    data_new(:,:,n) = data(indices,indices,n);%.*abs(B(indices,kk)*B(indices,kk)');
		end
		data_new(isnan(data_new))=0;
		ts_data_new = ts_data(indices,:);
		if(K_second>D_second)
		    multFact = ceil(K_second/D_second);
		    params.initdict = repmat(eye(D_second),1,multFact);
		    params.initdict = params.initdict(:,1:K_second);
		else
		    [IDX, ~] = kmeans(ts_data_new,K_second,'Distance','correlation','Replicates',10);  
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
	save(outprefix,'B','C','nodes');
else %% if hierarchical learning not specified
	B = B1;
	C = C1;
	save(outprefix,'B','C');
end

end % end of function
