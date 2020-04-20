
clear; clc;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';

load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','Hierarchy');


%% The plan is to load the unaveraged data and bootstrap, average and calculate hierarchy scores
ROIs_Select = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
IDs = fieldnames(StokALL);
ROIsPerID = cellfun(@(x) StokALL.(x).ROIs,IDs,'uni',false);

COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs_Select);
nROIs = numel(ROIs_Select);
load ROInames;

nboots = 150;
bootsize = 5;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs_Select, nboots, bootsize);

%% Average PDCs over bootraps and estimate hierarchy
if exist(['Hierarchy_' FileName '.mat'],'file')
    for b = 1:10%nboots
        % make the sub-sample stok structure
        for s = 1:bootsize
            Stok_sample.(['S' num2str(s)]) = StokALL.(BootIDs{b,s});
        end

        % average over it
        Stok_avg    = STOKAllAverage(Stok_sample,ROIs_Select);
        [Hf(:,:,:,b), Hccf(:,:,b)]  = HierarchyTimeFreq(Stok_avg);
        if mod(b,10)==0
            disp(b)
        end
    end

    save(['Hierarchy_' FileName],'Hf','Hccf','ROIs_Select');
    
else
    load(['Hierarchy_' FileName '.mat']);
end
%% ------------------------------------------------------------------------
% -----------------------Time Domain Hierarchy ----------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Overall hierarchy score
FIG1 = figure;
set(FIG1,'unit','inch','position',[5 0 8 3],'color','w')

Time = StokALL.S766640955.Times(51:end);
HccfN = (Hccf - mean(Hccf(:,Time>-.2 & Time<0,:),2))./mean(Hccf(:,Time>-.2 & Time<0,:),2)*100;
MH = mean(mean(HccfN(:,51:end,:),3),1);
CI(1,:) = quantile(mean(HccfN(:,51:end,:),1),.95,3);
CI(2,:) = quantile(mean(HccfN(:,51:end,:),1),.05,3);
fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],'k','facealpha',.3,'edgecolor','none');
hold on;
plot(Time,MH,'k','linewidth',1.5);
xlim([-.2 2])
hline(0,'k--')
vline(0,'k--')

xlabel('Time (S)')
title('Total Hierarchy Score % Change')
export_fig(FIG1,fullfile(FigPath,['Total_hierarchy_time_' FileName]),'-pdf');close;

%% Hierarchy of ROIs
FIG2 = figure;
set(FIG2,'unit','inch','position',[5 0 8 5],'color','w')
for roi = 1:nROIs
    HfR = squeeze(Hf(:,:,roi,:));
    HccfN = (HfR);
    MH = nanmean(nanmean(HccfN(:,51:end,:),3),1);
    CI(1,:) = quantile(mean(HccfN(:,51:end,:),1),.95,3);
    CI(2,:) = quantile(mean(HccfN(:,51:end,:),1),.05,3);
    fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],[0 0 0],'facealpha',.2,'edgecolor','none');
    hold on;
    LG(roi) = plot(Time,MH,'k','linewidth',1.5,'color',Colors(roi,:));
    xlim([-.2 2])
    %hline(0,'k--')
    vline(0,'k--')
    if roi==nROIs
        xlabel('Time (S)')
        ylabel('Hierarchy Score');
    end
    
    
end

legend(LG,cellfun(@(x) ROI_names.(x),ROIs_Select,'uni',false));
export_fig(FIG2,fullfile(FigPath,['ROIs_hierarchy_time_' FileName]),'-pdf');close;

%% Hierarchy change of ROIs
FIG3 = figure;
set(FIG3,'unit','inch','position',[5 0 8 15],'color','w')
for roi = 1:nROIs
    subplot(6,1,roi)
    HfR = squeeze(Hf(:,:,roi,:));
    HccfN = (HfR - mean(HfR(:,Time>-.2 & Time<0,:),2));
    MH = nanmean(nanmean(HccfN(:,51:end,:),3),1);
    CI(1,:) = quantile(mean(HccfN(:,51:end,:),1),.95,3);
    CI(2,:) = quantile(mean(HccfN(:,51:end,:),1),.05,3);
    fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],[0 0 0],'facealpha',.3,'edgecolor','none');
    hold on;
    plot(Time,MH,'k','linewidth',1.5,'color',Colors(roi,:));
    xlim([-.2 2])
    ylim([-.025 .025])
    %ylim([-30 40])
    hline(0,'k--')
    vline(0,'k--')
    if roi==nROIs
        xlabel('Time (S)')
        ylabel('Hierarchy Score Change');
    end
    title(ROI_names.(ROIs_Select{roi}))
    
end

export_fig(FIG3,fullfile(FigPath,['ROIs_hierarchy_change_time_' FileName]),'-pdf');close;


%% ------------------------------------------------------------------------
% --------------------Time-Frequency Domain Hierarchy ---------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Overall hierarchy score
CIinterv = .05;
FIG1 = figure;
set(FIG1,'unit','inch','position',[5 0 8 3],'color','w')
clear CI
Time = StokALL.S766640955.Times(51:end);
HccfN = (Hccf(:,51:end,:) - mean(Hccf(:,Time>-.2 & Time<0,:),2))./mean(Hccf(:,Time>-.2 & Time<0,:),2)*100;
CI(1,:,:) = quantile((HccfN),1-CIinterv,3);
CI(2,:,:) = quantile((HccfN),CIinterv,3);
Alpha = squeeze((CI(1,:,:)<0) | (CI(2,:,:)>0));
IM = imagesc(Time,[],mean(HccfN,3));
set(IM,'alphadata',double(Alpha)*.8+.2)
caxis([-5 5])
xlim([-.2 2])
vline(0,'k--')
axis xy;
colorbar;

xlabel('Time (S)')
title('Total Hierarchy Score % Change')
colormap('jet')
export_fig(FIG1,fullfile(FigPath,['Total_hierarchy_timefreq_' FileName]),'-pdf');close;

%% Hierarchy of ROIs
FIG2 = figure;
set(FIG2,'unit','inch','position',[5 0 8 15],'color','w')
for roi = 1:nROIs
    subplot(6,1,roi)
    HfR = squeeze(Hf(:,51:end,roi,:));
    HccfN = (HfR);
    
    CI(1,:,:) = quantile((HccfN),.95,3);
    CI(2,:,:) = quantile((HccfN),.05,3);
    Alpha = squeeze((CI(1,:,:)<0) | (CI(2,:,:)>0));
    IM = imagesc(Time,[],mean(HccfN,3));
    set(IM,'alphadata',double(Alpha))
    axis xy
    if roi==nROIs
        xlabel('Time (S)')
        ylabel('Frequency (Hz)');
    end
    xlim([-.2 2])
    vline(0,'w--')
    colorbar;
    title(ROI_names.(ROIs_Select{roi}))
end

export_fig(FIG2,fullfile(FigPath,['ROIs_hierarchy_timefrequency_' FileName]),'-pdf');close;

%% Hierarchy change of ROIs
CIinterv = .025;
FIG3 = figure;
set(FIG3,'unit','inch','position',[5 0 8 15],'color','w')
for roi = 1:nROIs
    subplot(6,1,roi)
    HfR = squeeze(Hf(:,51:end,roi,:));
    HccfN = (HfR - mean(HfR(:,Time>-.2 & Time<0,:),2))*100;
    
    CI(1,:,:) = quantile((HccfN),1-CIinterv,3);
    CI(2,:,:) = quantile((HccfN),CIinterv,3);
    Alpha = squeeze((CI(1,:,:)<0) | (CI(2,:,:)>0));
    IM = imagesc(Time,[],mean(HccfN,3));
    set(IM,'alphadata',double(Alpha)*.8+.2)
    axis xy
    if roi==nROIs
        xlabel('Time (S)')
        ylabel('Frequency (Hz)');
    end
    xlim([-.2 2])
    caxis([-4 4])
    vline(0,'w--')
    colorbar;
    title(ROI_names.(ROIs_Select{roi}))
end
colormap('jet')
export_fig(FIG3,fullfile(FigPath,['ROIs_hierarchy_change_timefrequency_' FileName]),'-pdf');close;
