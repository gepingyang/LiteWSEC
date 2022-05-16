%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code for the LiteWWSC algorithm, which is proposed in   %
% the following paper:                                                %
%LiteWSSC:a Lightweight framework for Web-Scale Spectral Clustering   %                                             %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ label_orig, label]= LiteWSC(dataname, s, p, r,  k)

% if nargin < 5
%     r = 5; % The number of nearest prototypes.
% end
% if nargin < 4
%     p = 500; % The number of prototypes.
% end
% if nargin < 3
%     s = 0; % tthe number of data points that are used to generate prototypes and prototype graph.
%             % s = 0 denotes that  all data of this btach is used to generate prototypes and
%             % prototype graph.
% end 
% if nargin < 2 
%     disp('you have to input the number of target clustering');
%     return;
% end
% if nargin < 1 | size(dataname) < 1
%     disp('you have to input data');
%     return;
% end
flg = 1;% flg = 1 denotes the first time loading data, we sholud contrcut prototypes graph, otherwise, we just need to assign data points to
        % their nearest prototpyes.
for dataname_batch = dataname
     if flg == 1
         flg = 0;
        dataname_batch = char(dataname_batch);
        load( dataname_batch,'fea','gnd'); %loading batch data from hard disk 
        N = size(fea, 1);
        fea = full(fea); 
        label_orig = gnd; %storage the true label
         %%%%%%%%%%%%%%%%%%%%%%%%%%%samples generation%%%%%%%%%%%%%
        indSmp = randperm(N);
        samples = fea(indSmp(1:s),:);

        %%%%%%%%%%%%Sam2Pro%%%%%%%%%%%%%%%%%%
        [min_Y, Y_p, prototypes]=Sam2Pro(samples, p, r, k);
        %%%%%%%%%%%%%%%%% batch cluster assigment%%%%%%%%%%%%%%%%%%
         if s ~= 0
            fea_remaining = fea(indSmp(s+1:N),:);
            D_remaining = EuDist2(fea_remaining, prototypes, 0);
            [~, min_remaining] = min(D_remaining, [], 2);
         end
         
         label = zeros(N,1);
         label(indSmp(1:s)) = min_Y;
         label(indSmp(s+1:N))= Y_p(min_remaining);

     else
        %%%%%%%%%%%%%%%%% batch cluster assigment%%%%%%%%%%%%%%%%%%
         dataname_batch = char(dataname_batch);
         load( dataname_batch,'fea', 'gnd');
         fea = full(fea);   
         label_orig = [label_orig;gnd];
         D = EuDist2(fea, prototypes, 0);
         [dump, idx2] = min(D, [], 2);
         label_batch = Y_p(idx2);
         label = [label;label_batch];
         clear  dump D fea gnd dump idx2 label_batch
     end
          
end
end


