%> @file  SCPTest.m
%> @brief Function that calculates coefficients for a given basis
%======================================================================
%> @brief It takes as input subject level time-series data. It returns
%> (1) a set of subject-level coefficients for the data
%> stored in variable C
%>
%> @param DataMatFile Input mat file containing correlation matrices shaped
%> as [D,T,N], where D is the number of ROIs, T is the number of time-points,
%> N is the number of subjects
%> @param InputBasisFile Input mat file containing the  basis
%> @param CoeffFile SCP Coeffs corresponding to input data
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
function [] = SCPTest(DataMatFile, InputBasisFile, CoeffFile,verbose)
%%
%> set up env: set the seed of random number generator to shuffle
rng('shuffle')
verbose = str2double(verbose);

%%
%> load data
load(DataMatFile);

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

%> load input basis file
load(InputBasisFile);
if (exist('B','var'))
    fprintf('Reading input basis\n')
else
    fprintf('Input basis not found\n')
    return
end

if (exist('nodes','var'))
    fprintf('Basis is hierarchical\n')
        hierarchy = 1;
else
    fprintf('Basis is not hierarchical\n')
        hierarchy = 0;
end

if(exist('sample_weights','var'))
    params.sample_weights = sample_weights;
else
    params.sample_weights = ones(N,1);
end


%%
%> specify input parameters
if (hierarchy == 0)
        params.initdict = B;
        params.data = data;
        params.iternum = 2;
        params.test = 1;
        [~,C,~] = BlockSolverLowRankFrob(params,verbose);
else
        C = zeros(size(B,2),N);
        count = 0;
        K_primary = sum(nodes==0);
        params.initdict = B(:,nodes==0);
        params.data = data;
        params.iternum = 2;
        params.test = 1;
        [~,C(nodes==0,:),~] = BlockSolverLowRankFrob(params,verbose);
        for kk=1:K_primary
                indices = abs(B(:,kk))>0.3;
                D_second = sum(indices);
   
                fprintf('Computing coeffs for secondary SCP %d\n',kk);
                data_new=zeros(D_second,D_second,N);
                for n=1:N
                    data_new(:,:,n) = data(indices,indices,n);
                end
                data_new(isnan(data_new))=0;

                params.initdict = B(indices,nodes==kk);
                params.data = data_new;
                params.iternum = 2;
                params.test = 1;
                [~,C(nodes==kk,:),~] = BlockSolverLowRankFrob(params,verbose);
        end
end


%%
%> saving the coefficients
save(CoeffFile,'C');

end % end of function
