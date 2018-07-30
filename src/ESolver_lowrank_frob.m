%> @file  ESolver_lowrank_frob.m
%> @brief Function that solves for the E variable
%======================================================================
%> @brief It takes as input the current value of B, current value of C,
%> correlation matrices. It returns the current
%> solution of D for the given B, C and data.
%> For details see the following paper:
%> http://www.sciencedirect.com/science/article/pii/S1053811914008003
%>
%> @param B The current value of B of size D x K
%> @param data The correltaion matrices of size D x D x N
%> @param C The current set of coefficients K x N
%> @retval E The diagonal estimates of size D x N
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
function [E]  = ESolver_lowrank_frob(B,C,data)
[D,~] = size(B);
[~,N] = size(C);
E = zeros(D,N);
for nn=1:N
    E(:,nn) = diag(data(:,:,nn) - B*diag(C(:,nn))*B');
end
end