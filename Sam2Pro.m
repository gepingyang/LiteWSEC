%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code for the LiteWWSC algorithm, which is proposed in            %
% the following paper:                                                         %
%LiteWSEC:a Lightweight framework for Web-Scale Spectral Ensemble Clustering   %                                             
%                                                                              %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [label, labels, prototypes] = Sam2Pro(fea,  p, r, k)
        s = size(fea,1);
%%%%%%%%%%%%%%%%%%%%%%%%%prototypes generation%%%%%%%%%%%%%%%%%%%%%%%%%
       [label1, prototypes] = litekmeans(fea, p ,'MaxIter', 3,'Replicates',1);

       clear label1

%%%%%%%prototype Laplacian graph construction%%%%%%%%%%%%%%%%%%%%%%
       D = EuDist2(fea,prototypes,0);
       [dump1, idx1] = min(D, [], 2);
       sigma = mean(mean(D));
       dump = zeros(s,r);
       idx = dump;
       dump(:, 1) = dump1;
       clear dump1
       idx(:,1) = idx1;

       for ii = 2:r
           temp = (idx(:,ii - 1)-1)*s+[1:s]';
           D(temp) = 1e100; 
           [dump(:,ii),idx(:,ii)] = min(D,[],2);
       end
       clear D
       dump = exp(-dump/2/sigma^2);
       Gsdx = dump;
       Gidx = repmat([1:s]',1,r);
       Gjdx = idx;
       clear dump 
       Z=sparse(Gidx(:),Gjdx(:),Gsdx(:),s,p);
       W = Z'*Z; %prototypes affinity matrix
       D_p = 1./(sqrt(sum(W,2))+10^(-6));%  prototypes degree matrix
       L = repmat(D_p,1,p).*W.*repmat(D_p', p, 1); %prototypes Laplacian matrix
       clear A
%%%%%%% %%%eigendedecomposition%%%%%%%%%%%%%%%%
       if issparse(L)
            L = full(L);
       end
        L = (L+L')/2;
        [U, eigvalue] = eig(L);
        eigvalue = diag(eigvalue);
        [dump, index] = sort(-eigvalue);
        clear dump
        U = U(:, index);
        U = U(:,2: k+1) ;
        U = U ./repmat(sqrt(sum(U.^2,2)),1,k);
 %%%%%%%%%%%%%%%% $k$-means%%%%%%%%%%%%%%%%
        [labels]=litekmeans(U,k,'MaxIter',100,'Replicates',10);
        label = labels(idx1);%getting smaples' clustering labels
end


