%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo for the LiteWWSC algorithm, which is proposed in the %
% following paper:                                                  %
%                                                                   %
%LiteWSC:a Lightweight framework for Web-Scale Spectral Clustering %
%                                                                   %
% The code has been tested in Matlab R2018b                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label, label_orig]= LiteWSEC(dataname, s, p, M, r, lowK, upK,  k)
 Ks = randsample(upK-lowK+1,M,true)+lowK-1; %the number of clusters of each ensembles
         flg = 1;
         tic;
         for dataname_batch = dataname
             if flg == 1
                flg = 0;
                dataname_batch = char(dataname_batch);
                load(dataname_batch,'fea','gnd'); %loading batch data from hard disk 
                gnd = double(gnd);
                label_orig = gnd; %storage the true label
                fea = full(fea);
                N = size(fea,1);
                d = size(fea, 2);
                members_s = zeros(s,M); % storage  the clustering labels of samples
                members_batch = zeros(N,M); % storage the clustering labels of this batch(N >= s);
                prototypes =  zeros(p,d,M); % storage the prototypes of m ensembles
                prototypes_label = zeros(p,M); % storage the clustering ;abels of prototypes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%samples generation%%%%%%%%%%%%%
                indSmp = randperm(N);
                samples = fea(indSmp(1:s),:); 
                %%%%%%%% m times Sam2Pro %%%%%%%%%%%%%%%%%%%%%%%
                 for j = 1:M
                     [members_s(:,j), prototypes_label(:,j),  prototypes(:,:,j)]=Sam2Pro(samples,  p, r, Ks(j));
                 end
                 %%%%%%%%%%%%%%%%%%% ensemble Laplacian grpah constrction%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                 maxCls = max(members_s);
                for j = 1:numel(maxCls)-1
                    maxCls(j+1) = maxCls(j+1)+maxCls(j);
                end
                cntCls = maxCls(end);
                members_s(:,2:end) = members_s(:,2:end) + repmat(maxCls(1:end-1),s,1); clear maxCls

               
                Z_e=sparse(repmat([1:s]',1,M),members_s(:),1,s,cntCls); clear members_s % Build the bipartite graph.
                W_e = Z_e'*Z_e; %affinity matrix
                D_e = 1./(sqrt(sum(W_e,2))+10^(-6));%  degree matrix
                L_e = repmat(D_e,1,sum(Ks)).*W_e.*repmat(D_e', sum(Ks), 1); %ensemble Laplacian matrix
                clear A
                 %%%%%%% %%% eigendedecomposition %%%%%%%%%%%%%%%%
               if issparse(L_e)
                    L_e = full(L_e);
               end
               L_e = (L_e+L_e')/2;
              [U, eigvalue] = eig(L_e);
              eigvalue = diag(eigvalue);
              [dump, index] = sort(-eigvalue);
               clear dump
                U = U(:, index);
                U = U(:,2: k+1);
                U = U ./repmat(sqrt(sum(U.^2,2)),1,k);
               %%%%%%%%%%%%%%%% $k$-means%%%%%%%%%%%%%%%%     
               [ensemble_label]=litekmeans(U,k,'MaxIter',100,'Replicates',10);
               %%%%%%%%%%%%%%%%% batch cluster assigment%%%%%%%%%%%%%%%%%%
                label = ensemble_label';
                count = 0;
                for j=1:M
                    label_1 = label(count+1:count+Ks(j));
                    prototypes_label(:,j) = label_1(prototypes_label(:,j));
                    count = count + Ks(j);
                end
                for j=1:M
                    D = EuDist2(fea,prototypes(:,:,j),0);   
                    [dump1, idx1] = min(D, [], 2);
                    current_label = prototypes_label(:,j);
                    members_batch(:,j) = current_label(idx1);

                end
                clear dump1 idx1 fea
                label = mode(members_batch, 2);
         else
                %%%%%%%%%%%%%%%%% batch cluster assigment%%%%%%%%%%%%%%%%%%
                dataname_batch = char(dataname_batch);
                load( dataname_batch,'fea', 'gnd');             
                fea = full(fea);   
                gnd = gnd';
                gnd = double(gnd);
                label_orig = [label_orig;gnd];
                N = size(fea, 1);
                members_batch = zeros(N,M);
                for j=1:M
                    D = EuDist2(fea,prototypes_all(:,:,j),0);   
                    [dump1, idx1] = min(D, [], 2);
                    current_label = prototypes_label(:,j);
                    members_batch(:,j) = current_label(idx1);
                    clear dump1 idx1
                end
                clear  fea
                label_batch =  mode(members_batch, 2);
                label = [label,label_batch];
             end
         end
end

     
