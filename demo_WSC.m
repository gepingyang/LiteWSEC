%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code for the LiteWWSC algorithm, which is proposed in            %
% the following paper:                                                         %
%LiteWSEC:a Lightweight framework for Web-Scale Spectral Ensemble Clustering   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [label, label_orig]= demo_LiteWSC(1)

dataname = {'./dataset/USPS.mat'};            
s =4000; %the number of data points that are used to generate prototypes and prototype graph.
p = 500; % the number of prototypes
r =5; %the number of nearest prototypes.
k =10; % the number of clusters.

seed.end = 10;
seed.start = 1;
interval = seed.end - seed.start + 1;
ac_sum = 0;
nmi_sum = 0;
time_all = 0;


for i = seed.start : seed.end 
     rand('seed',i);    
     fprintf('Seed No: %d\n',i);
     tic;
     [label, label_orig] = LiteWSC(dataname, s, p, r,  k);
     time_once = toc;
     time_all = time_once + time_all;
     label = bestMap(label_orig,label);
      nmi_result = nmi(label,label_orig);  
      ac_result = length(find(label_orig == label))/length(label);
      fprintf('nmi: %.2f%% +', nmi_result * 100);
      fprintf('ac: %.2f%% + ', ac_result * 100);
      fprintf("runtime: %.2f \n", time_once);
      nmi_sum = nmi_result + nmi_sum;
      ac_sum = ac_result + ac_sum;            
end
    nmi_avg = nmi_sum / interval;
    ac_avg = ac_sum / interval;
    time_avg = time_all / interval;
    fprintf('**************************************************************\n');
    fprintf('avg_nmi: %.2f%% + ', nmi_avg * 100);
    fprintf('avg_ac: %.2f%% + ', ac_avg * 100);
    fprintf('avg_runtime: %.2f s\n\n\n', time_avg);
    fprintf('Algorithm Finished');
             %%%%%%%%%%%%%%%%%%calculate clustering accuaracy%%%%%%%%%%%%%%%%%%%%%%
% end
