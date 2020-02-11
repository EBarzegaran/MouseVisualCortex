function [Stats,FIG] = CompPDCCells(PDC1,PDC2,avdims,prefactor,plotresults, rdiag, time, freq, figtitle,figpath, DoMean)

% Compares the PDC values in PDC1 and PDC2 using (mixed-effect) anova
% INPUTS:
    % PDC1 (Nx1) and PDC2(Mx1): a cell array containing a signle between/within area
        % or a set of them from the pdc matrix, calculated by
        % ExtractAllROIPDC (it should be applied to all pdc values calculated by
        % LFPF.STOK_analysis function)

        % Attention! Here it is better that the PDC cells are converted to
        % %change  by calculating this PDC_change = (PDC -
        % PDC_prestim)/PDC_prestim if prefactor is false, othrwise use the original PDC
        % values
        
    % avdiums: an empty, single or double element vector, indicating which
        % dimentions of PDC should be averaged over, it can be 1 (input) or 2
        % (output) or both or neither
    
    % prefactor: indicates if pre vs. post time points should be considered
    % as an extra factor in anova analysis, false or true. If true then the
    % stats is a two-way mixed linear model anova
        
        % rdiag is a 2x1 vector, indicating if it should remove the diagonals
        % of PDC1 and PDC2 (1 for remove, 0 for keep)
%% -----------------------------------------------------------------------
% DIAGONALS FOR INTER VS. INTRA

%% initialization
if sum(avdims>2)
    error('AVerage or factor dimensions should be 1 or 2')
end

avdims = sort(avdims,'descend'); % sort it so it is straightforward for averaging: no problem if not squeeze the matrix

if exist('rdiag','var') % remove diagonals in case of intra-area PDC
    dind = find(rdiag);
    for d = 1:numel(dind)
        eval(['PDC = PDC' num2str(dind) ';']);
        for c = 1:numel(PDC)
            for i = 1:size(PDC{c},1)
                PDC{c}(i,i,:,:,:) = NaN;
            end
        end
        eval(['PDC' num2str(dind) ' = PDC;']);
    end
end

if ~exist('figtitle','var')
    figtitle=' ';
end

if ~exist('DoMean','var')
    DoMean=true;
end

%% Average data over the requested dimensions
%-------------------(1) average over dimensions-----------------------
if ~isempty(avdims)% average over the requested dimensions, otherwise they are considered as independent samples
    for d=1:numel(avdims)
        if DoMean
            PDC1 = cellfun(@(x) nanmean(x,avdims(d)),PDC1,'uni',false);
            PDC2 = cellfun(@(x) nanmean(x,avdims(d)),PDC2,'uni',false);
        else
            % bring that dimension to the end and add a one instead of that
            dims = size(PDC1{1});
            %reshape(permute(PDC1{1},[setdiff(1:numel(dims),avdims(d)) avdims(d)]),[dims(1:avdims(d)-1) 1 dims(avdims(d)+1:numel(dims)-1) dims(end)*dims(avdims(d))]);
            PDC1 = cellfun(@(x) reshape(permute(x,[setdiff(1:numel(dims),avdims(d)) avdims(d)]),[dims(1:avdims(d)-1) 1 dims(avdims(d)+1:numel(dims)-1) size(x,5)*dims(avdims(d))]),PDC1,'uni',false);
            PDC2 = cellfun(@(x) reshape(permute(x,[setdiff(1:numel(dims),avdims(d)) avdims(d)]),[dims(1:avdims(d)-1) 1 dims(avdims(d)+1:numel(dims)-1) size(x,5)*dims(avdims(d))]),PDC2,'uni',false);
           
        end
    end
end
%-------------------(2) pull all the data together---------------------
% but first keep the animal ID for anova analysis
PDC1 = cat(5,PDC1{:});
PDC2 = cat(5,PDC2{:});


%% if no additional factor is needed: it is the simplest case: ttest or ranksum   
if ~(prefactor)
    
    %-------------------(3) there are two ways: ---------------------------
    %either consider the remaining dimensions as
    % indipendent samples (pull them together) or apply the test
    % independently on them!
    for d1 = 1:size(PDC1,1)
        d1
        for d2 = 1:size(PDC1,2)
            for d3 = 1:size(PDC1,3)
                for d4 = 1:size(PDC1,4)
                    [~,P(d1,d2,d3,d4),ci(d1,d2,d3,d4,:),SS] = ttest2(PDC1(d1,d2,d3,d4,:),PDC2(d1,d2,d3,d4,:));
                    Tval(d1,d2,d3,d4) = SS.tstat;
                end
            end
        end
    end
    
    %------------------(4) next step for this simple test is to implement a
    %cluster-based statistics --- NOT IMPLEMENTED YET
    Stats.P = P;
    Stats.CI = ci;
    Stats.Tval = Tval;

    %---------------------------plot the results---------------------------
    if plotresults
        FIG = PlotStatresults(P,Tval,avdims, time, freq,...
        'figtitle'  ,figtitle,...
        'figpath'   ,figpath,...
        'PThresh'   ,0.05,...
        'SThresh'   ,8,...
        'Twin'      ,[0 1],...
        'ColorM'    ,'jet');
    end
else
%% %% if additional factor is needed, then I should use anova and linear mixed model, with a within-group factor (pre-post) and between-group factor (PDC1 and PDC2)
    TPres = find(time<0 & time>-.2);
    
    for d1 = 1:size(PDC1,1)
        for d2 = 1:size(PDC1,2)
            for d3 = 1:size(PDC1,3)
                d3
                for d4 = find(round(time,2)==0,1):500%size(PDC1,4)
                    % first make a table with PCD, prepost, interintra, and ID fields
                    % ATTENTION: Careful about categorical factors! double check them!!
%                     PDC1TPre        = PDC1(d1,d2,d3,TPres,:); PDC1TPre = PDC1TPre(:);ID1 = repmat(1:size(PDC1,5),[numel(TPres) 1]);
%                     PDC2TPre        = PDC2(d1,d2,d3,TPres,:); PDC2TPre = PDC2TPre(:);ID2 = repmat((1:size(PDC2,5))+1000,[numel(TPres) 1]);
%                     PDC             = [squeeze(PDC1(d1,d2,d3,d4,:)); PDC1TPre; squeeze(PDC2(d1,d2,d3,d4,:)); PDC2TPre];
%                     prepost         = [ones(size(PDC1,5),1); zeros(numel(PDC1TPre),1); ones(size(PDC2,5),1); zeros(numel(PDC2TPre),1)];
%                     interintra      = [ones(size(PDC1,5),1); ones(numel(PDC1TPre),1); ones(size(PDC2,5),1)+1; ones(numel(PDC2TPre),1)+1];
%                     IDs             = [1:size(PDC1,5) ID1(:)' (1:size(PDC2,5))+1000 ID2(:)']'; % THIS SHOULD BE FIXED
%                     
%                     
                    
                    
                    PDC             = [squeeze(PDC1(d1,d2,d3,d4,:)); squeeze(mean(PDC1(d1,d2,d3,TPres,:),4)); squeeze(PDC2(d1,d2,d3,d4,:)); squeeze(mean(PDC2(d1,d2,d3,TPres,:),4))];
                    prepost         = [ones(size(PDC1,5),1); zeros(size(PDC1,5),1); ones(size(PDC2,5),1); zeros(size(PDC2,5),1)];
                    interintra      = [ones(size(PDC1,5),1); ones(size(PDC1,5),1); ones(size(PDC2,5),1)+1; ones(size(PDC2,5),1)+1];
                    IDs             = [1:size(PDC1,5) 1:size(PDC1,5) (1:size(PDC2,5))+1000 (1:size(PDC2,5))+1000]'; % THIS SHOULD BE FIXED
                    ds              = table(PDC,prepost,interintra,IDs,'VariableNames',{'PDC','prepost','interintra','ID'});
                    ds.prepost      = categorical(ds.prepost);
                    ds.interintra   = categorical(ds.interintra);
                    ds.ID           = categorical(ds.ID);
                    
                    lme             = fitlme(ds,'PDC ~ prepost * interintra + (1|ID)');
                    %lme             = fitlme(ds,'PDC ~ interintra + (1|ID)');
                    % keep the results
                    P(d1,d2,d3,d4,:)    = lme.anova.pValue;
                    Fval(d1,d2,d3,d4,:) = lme.anova.FStat;
                end
            end
        end
    end

    Stats.P = P;
    Stats.Tval = Fval;
    
    %---------------------------plot the results---------------------------
    if plotresults
        FIG = PlotStatresults(permute(P(:,:,:,find(round(time)==0,1)-1:500,:),[1 5 3 4 2]),permute(Fval(:,:,:,find(round(time)==0,1)-1:500,:),[1 5 3 4 2]),...
            1, time(find(round(time)==0,1)-1:500), freq,...
            'figpath'   ,figpath,...
            'figtitle'  ,figtitle,...
            'Pthresh'   ,.01,...
            'SThresh'   ,200,...
            'Twin'      ,[0 1],...
            'Labels'    ,lme.anova.Term,...
            'ColorM'    ,'jet'...
            );
    end

end
end


