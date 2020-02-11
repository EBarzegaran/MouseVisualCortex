function [H,Hcc] = HierarchyExtract(M)
% INPUT:
    %   M: connectivity matrix (Nodes x Nodes) with rows as targets and columns as sources

% remove diagonal if any    
M(1:length(M)+1:end)    =   0;
NNode                   =   size(M,1);

%% initialize the connections to FF(=+1) and FB(=-1) and the hierarchyscores
MD = M-M';
% 1- best permutation order
Perms = perms(1:NNode);
MDp = arrayfun(@(x) nansum(nansum(tril(MD(Perms(x,:),Perms(x,:))))),1:size(Perms,1));
[~,BestPerm] = max(MDp);
BestPerm =  1:NNode;%Perms(BestPerm,:); %%%%%%% 

% 2- assign FF and FBs according to order
Mcc                         = ones(size(M));
Mcc(find(triu(Mcc)))        = -1;
Mcc(1:length(Mcc)+1:end)    = 0; 
[~,Ind]                     = sort(BestPerm);
Mcc                         = Mcc(Ind,Ind);

% 3- initialize the node hierarchy score
H(1,:)                     = (nansum(Mcc.*M,2) - nansum((Mcc.*M)',2))/(2*nansum(M(:)));

% 4- initialize the global hierarechy score
Hcc(1)                     =  nansum(nansum((M.*Mcc).*(H(1,:)'-H(1,:))))/nansum(M(:));
%% iteration to optimize for node and global hierarchy
% % err         = 1;
%  for it = 2:5
%   % 2- assign FF and FBs according to order
%     Mcc                         = ones(size(M));
%     Mcc(find(triu(Mcc)))        = -1;
%     Mcc(1:length(Mcc)+1:end)    = 0; 
%     [~,Ind]                     = sort(H(it-1,:));
%     Mcc                         = Mcc(Ind,Ind); 
%     
%     H(it,:)                     = (nansum(Mcc.*M,2) - nansum((Mcc.*M)',2))/(2*nansum(M(:)));
% 
%     Hcc(it)                     =  nansum(nansum((M.*Mcc).*(H(1,:)'-H(1,:))))/nansum(M(:));
     
%     
%     
%     for i = 1:NNode
%         H(it+1,i) = -(nanmean(Mcc(setdiff(1:NNode,i),i)'+H(it,setdiff(1:NNode,i))) - nanmean(Mcc(i,setdiff(1:NNode,i))-H(it,setdiff(1:NNode,i))))/2;
%     end
%     
%     for i = 1:NNode
%         H(it+2,i) = H(it+1,i)-mean(H(it+1,setdiff(1:NNode,i)));
%     end
% end
end