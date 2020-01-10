function [Signal_Averaged, FIG] = SignalROIAverage(Session,ROIselect,savefig,Path)

Sessions_ID = fieldnames(Session);

ROIs = cellfun(@(x) Session.(x).ROIs,Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(2,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

for roi = 1:numel(ROIs)
    for S = 1:numel(Sessions_ID)
        if isfield(Session.(Sessions_ID{S}),ROIs{roi})
            Temp{S} =  squeeze(mean(mean(Session.(Sessions_ID{S}).(ROIs{roi}).Y(:,:,:,1:4),1),4));
            Signal_Averaged.(ROIs{roi}).Times = Session.(Sessions_ID{S}).(ROIs{roi}).Times;
        end
    end
    Signal_Averaged.(ROIs{roi}).SM = mean(cat(3,Temp{:}),3);
    Signal_Averaged.(ROIs{roi}).STD = std(cat(3,Temp{:}),[],3);
    Signal_Averaged.(ROIs{roi}).Num = ROI_counts(roi);
    clear Temp;
end
%% plot the results
Colors = brewermap(6,'Spectral');
if ~isempty(ROIselect)
    ROIs = intersect(ROIs,ROIselect);
end
FIG = figure;
for roi = 1:numel(ROIs)
    if numel(ROIs)>5
        subplot(ceil(numel(ROIs)/2),2,roi);
    else
        subplot(numel(ROIs),1,roi)
    end
    hold on;
    Signal = Signal_Averaged.(ROIs{roi}).SM;
    SignalErr = Signal_Averaged.(ROIs{roi}).STD;
    Offset = max(Signal(:))*1;
    
    for l = 1:size(Signal,1)
        X = Signal_Averaged.(ROIs{roi}).Times;
        Y = Signal(l,:)-Offset*l;
        Err = SignalErr(l,:)/sqrt(Signal_Averaged.(ROIs{roi}).Num);
        fill([X X(end:-1:1)],[Y+Err Y(end:-1:1)-Err(end:-1:1)],Colors(l,:),'facealpha',.3,'edgecolor','none');
    end
    
    for l = 1:size(Signal,1)
        plot(Signal_Averaged.(ROIs{roi}).Times,Signal(l,:)-Offset*l,'linewidth',1.5,'color',Colors(l,:));
    end
    xlim([-.1 1])
    title([ROIs{roi} ' - Averaged over ' num2str(Signal_Averaged.(ROIs{roi}).Num) ' sessions'])
    
    if roi == numel(ROIs)
        legend(arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false))
        xlabel('Time(s)')
    end
    vline(0.0,'k--')
    vline(.044,'k--')
end

set(FIG,'color','w','unit','inch','position',[0 0 15 15])
if savefig
    export_fig(FIG,fullfile(Path,['Signal_Averaged_AllROIs']),'-pdf');
end
end