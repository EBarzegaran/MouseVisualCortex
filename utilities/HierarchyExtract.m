function [H,Hcc,BestPerm] = HierarchyExtract(M,Method)
% INPUT:
    %   M: connectivity matrix (Nodes x Nodes) with rows as targets and columns as sources
%%
% References: Bastos et al 2015
% Harris et al, 2019

%%

if ~exist('Method','var')
    Method = 'tracer';
end

%%
% remove diagonal if any    
M(1:length(M)+1:end)    =   0;
NNode                   =   size(M,1);

if strcmpi(Method,'tracer')
    %% initialize the connections to FF(=+1) and FB(=-1) and the hierarchyscores
    MD = (M-M');
    % 1- best permutation order
%      Perms = perms(1:NNode);
%      MDp = arrayfun(@(x) nansum(nansum(tril(MD(Perms(x,:),Perms(x,:))))),1:size(Perms,1));
%      [~,BestPerm] = max(MDp);
    BestPerm =  1:NNode;%
    %BestPerm =  Perms(BestPerm,:); %%%%%%% 

    % 2- assign FF and FBs according to order, Mcc is the mapping function
    Mcc                         = ones(size(M));
    Mcc(find(triu(Mcc)))        = -1;
    Mcc(1:length(Mcc)+1:end)    = 0; 
    [~,Ind]                     = sort(BestPerm);
    Mcc                         = Mcc(Ind,Ind); 

    % 3- initialize the node hierarchy score
    H(1,:)                     = (nansum(Mcc.*M,2) - nansum((Mcc.*M)',2))/(2*nansum(M(:)));

    % 4- initialize the global hierarechy score
    Hcc(1)                     =  nansum(nansum((M.*Mcc).*(H(1,:)'-H(1,:))))/nansum(M(:));
elseif strcmpi(Method,'functional')
    
    DAI = (M'-M)./(M'+M);
    DAI(isnan(DAI))=0;
    %---------------------------- Bastos paper------------------------
    % (1)rescale to -3 to 3
    DAI = (DAI-min(DAI(:)))./(max(DAI(:))-min(DAI(:)));
    DAI = DAI*size(DAI,1)-size(DAI,1)/2;
    %DAI = (DAI)*3;
    % (2) for each target (row) shifte the rescaled mDAI values of all source areas such that the smallest value was one. 
    for roi = 1:size(DAI,1)
        ind = setdiff(1:size(DAI,1),roi);
        %DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
        DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
        DAI(roi,roi) = NaN;
    end
    H = nanmean(DAI);
    %H = nanmean(DAI-nanmean(DAI)');
    [~,BestPerm] = sort(H,'ascend');
    
    Mcc                         = ones(size(M));
    Mcc(find(triu(Mcc)))        = -1;
    Mcc(1:length(Mcc)+1:end)    = 0; 
    [~,Ind]                     = sort(BestPerm);
    Mcc                         = Mcc(Ind,Ind); 
    H2 = nanmean(DAI-nanmean(DAI)');
    Hcc(1)                     =  nansum(nansum((M.*Mcc).*(H2(1,:)'-H2(1,:))))/nansum(M(:));
end
end