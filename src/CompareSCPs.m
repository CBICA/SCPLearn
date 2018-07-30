%> @file  CompareSCPs.m
%> @brief Function that compares two sets of SCPs
%======================================================================
%> @brief It takes as input the initial value of C, current value of B,
%> correlation matrices and the solver options. It returns the current
%> solution of C for the given B and data.
%> For details see the following paper:
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003
%>
%> @param B1 First set of SCPs of size D x K
%> @param B2 Second set of SCPs of D x K
%> @retval B_perm Permuted version of B2, that matches B1 the most
%> @retval reprod Average pair-wise inner-product of B1 and B_perm
%> @retval reprod_std Standard deviation of pair-wise inner-product 
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

function [B_perm, reprod,reprod_std] = CompareSCPs(B1, B2)
%%
%> calculating the similarity between the two SCPs
K = size(B1,1);
if(size(B1,1) ~= size(B2,1))
    error('SCPs for comparison do not have same size\n')
end% end if
norm_1 = sqrt(diag(B1'*B1));
norm_2 = sqrt(diag(B2'*B2));
Similarity = abs(B1'*B2 ./ (norm_1*norm_2'));

%%
%> running the hungarian algorithm
[Matching,~] = Hungarian(exp(-1*Similarity));

%%
%> finding the matches in the second set and flipping the sign, if needed
[r,c]=find(Matching==1);
[~,ind]=sort(r,'ascend');

if(isempty(setdiff(ind,1:K)))
    B_perm=B2(:,c(ind));
    norm_1 = sqrt(diag(B1'*B1));
    norm_2 = sqrt(diag(B_perm'*B_perm));
    Similarity = B1'*B_perm ./ (norm_1*norm_2');
    inner_products = diag(Similarity);
    sign_change=ones(1,size(B_perm,2));
    sign_change(inner_products<0)=-1;
    B_perm=bsxfun(@times,B_perm,sign_change);
else
    fprintf('Error with one-to-one Hungarian matching')
end% end if

%%
%> computing similarity between the matched sets
norm_1 = sqrt(diag(B1'*B1));
norm_2 = sqrt(diag(B_perm'*B_perm));
Similarity = B1'*B_perm ./ (norm_1*norm_2');
reprod = mean(diag(Similarity));
reprod_std = sqrt(var(diag(Similarity)));
end% end function 
