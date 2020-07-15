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

    % 1- best permutation order
    Perms = perms(1:NNode);

    % 2- assign FF and FBs according to order, Mcc is the mapping function
    Mcc                         = ones(size(M));
    Mcc(find(triu(Mcc)))        = -1;
    Mcc(1:length(Mcc)+1:end)    = 0; 
    [~,Ind]                     = arrayfun(@(x) sort(Perms(x,:)),1:size(Perms,1),'uni',false);
    %Mcc                         = Mcc(Ind,Ind); 

    % 3- initialize the node hierarchy score
    H                       = cellfun(@(x) (nansum(Mcc(x,x).*M,2) - nansum((Mcc(x,x).*M)',2))/(2*nansum(M(:))),Ind,'uni',false);

    % 4- initialize the global hierarechy score
    %Hcc                     =  cellfun(@(x) nansum(nansum((M.*Mcc).*(x-x')))/nansum(M(:)),H,'uni',false);
    Hcc                     =  arrayfun(@(x) nansum(nansum((M.*Mcc(Ind{x},Ind{x})).*(H{x}-H{x}')))/nansum(M(:)),1:numel(H),'uni',false);
    Hcc                     = cat(1,Hcc{:});
    [Hcc,Order_init]          = max(Hcc);
    H_init                  = H{Order_init};
    Mcc                     = Mcc(Ind{Order_init},Ind{Order_init});% Correct the FF and FB order
    % 5- iterations: DOES NOT WORK with functional
    H = H_init;
    
%     Mcc(1:length(Mcc)+1:end)=NaN;
%     for i = 1:5
%         Hm = (nanmean(Mcc+H,2) - nanmean((Mcc-H')',2))/(2);
%         for i = 1:numel(Hm)
%             H(i) = Hm(i) - mean(Hm(setdiff(1:numel(Hm),i)));
%         end
%     end
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