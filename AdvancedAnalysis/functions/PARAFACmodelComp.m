function [ord2,sscore,Corrs,model2,CorrPvals,Corrvals_out] = PARAFACmodelComp(model1,model2,modes)

if ~exist('modes','var')
    modes = 1:numel(model1);
end

modenum = numel(model1);
compnum = size(model1{1},2);
if modenum~=numel(model2)
    error('Models are not comparable!')
end

Corrvals = arrayfun(@(x) corr(model1{x},model2{x}),modes,'uni',false);
Corrvals_out = (mean(cat(3,Corrvals{:}),3));

bestords = cellfun(@extractorder,Corrvals);

% then select the most common sequence
uniqord = unique(bestords); % orders
comord  = histc(bestords,uniqord); % hist of the orders
[~,ind] = max(comord);
perm = perms(1:compnum);
ord2 = perm(uniqord(ind),:);

model2 = cellfun(@(x) x(:,ord2),model2,'uni',false);

sscore = mean(cellfun(@(x) mean(diag(x(:,ord2))),Corrvals));
[Corrvals,CorrPvals] = arrayfun(@(x) corr(model1{x},model2{x}),1:modenum,'uni',false);
Corrs = diag(mean(cat(3,Corrvals{:}),3));
CorrPvals = cellfun(@diag,CorrPvals,'uni',false); CorrPvals = cat(2,CorrPvals{:});% comp x mode


end


function [ind,ord2] = extractorder(M)
    perm = perms(1:length(M));
    oval = arrayfun(@(x) sum(diag(M(:,perm(x,:)))),1:size(perm,1));
    [~,ind] = max(oval);
    ord2 = perm(ind,:);
end
