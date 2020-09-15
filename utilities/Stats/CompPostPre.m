function [Stats,FIG] = CompPostPre(PDC,avdims,facdims,plotresults, rdiag, time, freq, figtitle,figpath)

%% initialization

if ~isempty(intersect(avdims,facdims))
    error('avdims and facdims should not contain common dimension')
end
if sum(avdims>2) || sum(facdims>2)
    error('AVerage or factor dimensions should be 1 or 2')
end

avdims = sort(avdims,'descend'); % sort it so it is straightforward for averaging: no problem if not squeeze the matrix

if exist('rdiag','var') && rdiag==1% remove diagonals in case of intra-area PDC
    for c = 1:numel(PDC)
        for i = 1:size(PDC{c},1)
            PDC{c}(i,i,:,:,:) = NaN;
        end
    end
end

if ~exist('figtitle','var')
    figtitle=' ';
end

%% Average data over the requested dimensions
%-------------------(1) average over dimensions-----------------------
if ~isempty(avdims)% average over the requested dimensions, otherwise they are considered as independent samples
    for d=1:numel(avdims)
        PDC = cellfun(@(x) nanmean(x,avdims(d)),PDC,'uni',false);
    end
end
%-------------------(2) pull all the data together---------------------
PDC = cat(5,PDC{:});

%% if no additional factor is needed: it is the simplest case: ttest or ranksum
if isempty(facdims)
    %-------------------(3) there are two ways: ---------------------------
    %either consider the remaining dimensions as
    % indipendent samples (pull them together) or apply the test
    % independently on them!
    TPres = find(time<0 & time>-.2);
    for d1 = 1:size(PDC,1)
        d1
        for d2 = 1:size(PDC,2)
            for d3 = 1:size(PDC,3)
                for d4 = find(round(time,2)==0,1):size(PDC,4)
                    %PDCpres = mean(PDC(d1,d2,d3,TPres,:),4);
                    PDCpres = PDC(d1,d2,d3,TPres,:);
                    [~,P(d1,d2,d3,d4),ci(d1,d2,d3,d4,:),SS] = ttest2(squeeze(PDC(d1,d2,d3,d4,:)),PDCpres(:));
                    Tval(d1,d2,d3,d4) = SS.tstat;
                end
            end
        end
    end
    
    %------------------(4) next step for this simple test is to implement a
    %cluster-based statistics
    
    Stats.P = P;
    Stats.CI = ci;
    Stats.Tval = Tval;


    %--------------------------plot the results---------------------------
    if plotresults
        FIG = PlotStatresults(P,Tval,avdims, time, freq,...
            'figtitle'  ,figtitle,...
            'figpath'   ,figpath,...
            'PThresh'   ,.01/2500,... %benferoni correction!
            'SThresh'   ,20,...
            'Twin'      ,[-.0 1],...
            'ColorM'    ,'jet');
    else
        FIG=[];
    end
end


end