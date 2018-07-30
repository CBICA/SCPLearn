%> @file  testSCPLearnCode.m
%> @brief Function that tests all main matlab-based functions
%======================================================================
%> @brief It uses NetSim-generated synthetic data to test run all major
%functions
%> @b Author: 
%> Harini Eavani
%>
%> @b Link: 
%> https://www.cbica.upenn.edu/sbia/software/
%> 
%> @b Contact: 
%> sbia-software@uphs.upenn.edu
%======================================================================
function [] = testSCPLearnCode()
%%
%> clear everything
clear
clc
close all
verbose = '0';

%> input and output test files
inFile = 'netsim_data_05262014_timeseries.mat';
outFile = 'netsim_data_05262014_timeseries_results_SCPs.mat';

%> parameters 
K='20';
lambda='5';
pruningThr='0.7';

%> call main function
fprintf('Starting SCPLearn on simulated data\n ')
SCPLearn(inFile, K, lambda, outFile,verbose,pruningThr);

%> load results
load([outFile]);

%> compare with ground truth SCPs
load(inFile)
[B_perm, reprod,reprod_std] = CompareSCPs(Basis,B);

h=figure('visible','off');
for ii=1:8
subplot(2,8,ii), imagesc(Basis(:,ii)*Basis(:,ii)'), title(sprintf('Simulated%d',ii)),caxis([-1,1]),axis square
subplot(2,8,ii+8), imagesc(B_perm(:,ii)*B_perm(:,ii)'), title(sprintf('Estimated%d',ii)),caxis([-1,1]),axis square
end
filename=['SCPs_simulated.png'];
saveas(h,filename, 'png')

%> call split sample function five times
fprintf('Starting split-sample validation on simulated data\n ')
outFile = ['netsim_data_05262014_SCPs_splitsamples'];
SCPLearn_SplitSample(inFile, K, lambda, outFile,verbose,pruningThr);
load([outFile,'_SCPs.mat']);
fprintf('\n\n\n The accuracy comparing ground truth SCPs and the experimental result is %1.4f pm %1.4f \n\n\n ',reprod, reprod_std)

%> run heirarchical SCP learn 
outFile = ['netsim_data_05262014_heirarchy'];
SCPLearn_Hierarchy_twolevel(inFile, K, lambda, outFile,verbose,pruningThr);
load([outFile,'_SCPs.mat']);

figure('visible','off');
set(gcf,'Renderer','zbuffer');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
treeplot(nodes)
hold on;
[x,y] = treelayout(nodes);
for i=1:length(x)
    text(x(i),y(i),num2str(i),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',8)
end
set(gca,'xtick',[],'ytick',[])
xlabel('')
ylabel('')
title('')
filename=[outFile,'_tree'];
print(filename,'-dpng','-r600')

fprintf('\n\n\n Code run on test data complete\n\n\n')
end


