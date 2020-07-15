function [H,Hcc,BestPerm] = HierarchyExtract(M,labels,Method)
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
NNode                   =   size(M,2);

if strcmpi(Method,'tracer')
    %% initialize the connections to FF(=+1) and FB(=-1) and the hierarchyscores

    % 2- assign FF and FBs according to order, Mcc is the mapping function
    Mcc                         = ones(size(M)).*squeeze(labels)';

    % 3- initialize the node hierarchy score
    H(1,:)                     = mean((nansum(Mcc.*M,3) - nansum(permute((Mcc.*M),[1 3 2]),3))/(2*nansum(M(:))));

    % 4- initialize the global hierarechy score
    Hcc(1)                     =  0;%nansum(nansum((M.*Mcc).*(H(1,:)'-H(1,:))))/nansum(M(:));
    
   
end
end