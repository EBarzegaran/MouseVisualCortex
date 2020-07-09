function [H,Hcc,BestPerm] = HierarchyExtractTimeFreq(M,Method)
% INPUT:
    %   M: connectivity matrix (Nodes x Nodes X freq x time) with rows as targets and columns as sources
%%
% References: Bastos et al 2015
% Harris et al, 2019

%%

if ~exist('Method','var')
    Method = 'tracer';
end

%%
% remove diagonal if any    
%M(1:length(M)+1:end,:,:)    =   0;
NNode                   =   size(M,1);

if strcmpi(Method,'tracer')
    %% initialize the connections to FF(=+1) and FB(=-1) and the hierarchyscores
    MD = (M-permute(M,[2 1 3 4]));
    % 1- best permutation order
    Perms = perms(1:NNode);
    Filt = repmat(tril(ones(size(MD,1)),-1),[1 1 size(MD,3) size(MD,4)]);
    MDp = arrayfun(@(x) squeeze(nansum(nansum(MD(Perms(x,:),Perms(x,:),:,:).*Filt))),1:size(Perms,1),'uni',false);
    MDp = cat(3,MDp{:});
    [~,BestPerm] = max(MDp,[],3);
    %BestPerm =  1:NNode;%
    BestPerm =  Perms(BestPerm,:); %%%%%%% 

    % 2- assign FF and FBs according to order, Mcc is the mapping function
    Mcc                         = ones(size(M,1));
    Mcc(find(triu(Mcc)))        = -1;
    Mcc(1:length(Mcc)+1:end)    = 0; 
    [~,Ind]                     = arrayfun(@(x) sort(BestPerm(x,:)),1:size(BestPerm),'uni',false);
    Mcc                         = cellfun(@(x) Mcc(x,x),Ind,'uni',false); 

    % 3- initialize the node hierarchy score
    M2                      = reshape(M,[size(M,1) size(M,2) size(M,3)*size(M,4)]);
    H                       = arrayfun(@(x) (nansum(Mcc{x}.*M(:,:,x),2) - nansum((Mcc{x}.*M2(:,:,x)),1)')/(2*nansum(nansum(M2(:,:,x)))),1:numel(Mcc),'uni',false);
   
    % 4- initialize the global hierarechy score
    Hcc                     = arrayfun(@(x) nansum(nansum((M(:,:,x).*Mcc{x}).*(H{x}'-H{x})))/nansum(nansum(M2(:,:,x))),1:numel(Mcc));
    
    H = cat(2,H{:});
    H = reshape(H,size(H,1),size(M,3),size(M,4));
    
    Hcc = reshape(Hcc,size(M,3),size(M,4));
    BestPerm = reshape(BestPerm',size(BestPerm',1),size(M,3),size(M,4));
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
    [~,BestPerm] = sort(H,'ascend');
    
    Mcc                         = ones(size(M));
    Mcc(find(triu(Mcc)))        = -1;
    Mcc(1:length(Mcc)+1:end)    = 0; 
    [~,Ind]                     = sort(BestPerm);
    Mcc                         = Mcc(Ind,Ind); 
    
    Hcc(1)                     =  nansum(nansum((M.*Mcc).*(H(1,:)'-H(1,:))))/nansum(M(:));
end
end