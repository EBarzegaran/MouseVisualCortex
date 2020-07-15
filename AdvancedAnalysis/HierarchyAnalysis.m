
clear; clc;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';

load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','Hierarchy');
FS=16;

%% The plan is to load the unaveraged data and bootstrap, average and calculate hierarchy scores
ROIs_Select = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
IDs = fieldnames(StokALL);
ROIsPerID = cellfun(@(x) StokALL.(x).ROIs,IDs,'uni',false);

COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs_Select);
nROIs = numel(ROIs_Select);
load ROInames;
ROISN = cellfun(@(x) ROI_names.(x),ROIs_Select,'uni',false);

nboots = 150;
bootsize = 8;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs_Select, nboots, bootsize);

%% Average PDCs over bootraps and estimate hierarchy
% this code is very slow, the first thing to do is to select a shorter time
% window ([-.8 1])
TW = [-.8 1];
Method = 'tracer';
if ~exist(['Hierarchy_' Method '_' FileName '.mat'],'file')
    for b = 1:nboots
        tic
        % make the sub-sample stok structure
        for s = 1:bootsize
            Stok_sample.(['S' num2str(s)]) = StokALL.(BootIDs{b,s});
        end

        % average over it
        
        Stok_avg    = STOKAllAverage(Stok_sample,ROIs_Select);
        
        ind = find((Stok_avg.Time>TW(1)) & (Stok_avg.Time<TW(2)));
        Stok_avg.Time = Stok_avg.Time(ind);
        Time = Stok_avg.Time;
        Stok_avg.PDC = Stok_avg.PDC(:,:,:,ind);
        [Hf(:,:,:,b), Hccf(:,:,b),BestPerm(:,:,:,b)]  = HierarchyTimeFreq(Stok_avg,Method);
        toc
        if mod(b,10)==0
            disp(b)
        end
    end

    %save(['Hierarchy_' FileName '_BastosM'],'Hf','Hccf','ROIs_Select','Time');
    save(['Hierarchy_' Method '_' FileName ],'Hf','Hccf','BestPerm','ROIs_Select','Time','nboots','bootsize','BootIDs');
    
else
    load(['Hierarchy_' Method '_' FileName ]);
    %load(['Hierarchy_' FileName '_BastosM.mat']);
end

% 
%Hf = permute(Hf,[2 3 1 4]);
%% average hierarchy graph plot for gijs
Stok_avg    = STOKAllAverage(StokALL,ROIs_Select);
NROIs = numel(ROIs_Select);

for roi1 = 1:NROIs
   ind1 = (roi1-1)*6+1:roi1*6;
   for roi2 = 1:NROIs
       if roi2 ~=roi1
           ind2 = (roi2-1)*6+1:roi2*6;
           HscoreMR(roi2,roi1,:,:) = nanmean(nanmean(Stok_avg.PDC(ind2,ind1,:,Stok_avg.Time>TW(1) & Stok_avg.Time<TW(2)),1),2);
       end
   end
end
%%
TW_av = [-.3 0; .05 .3; .3 1]; % time windows to average on
Freq = 30:100;

FIG = figure;
set(FIG,'unit','inch','position',[1 1 10 3],'color','w')
Tles = {'[-300 0] msec','[50 300] msec','[300 1000] msec'};
for t = 1:size(TW_av,1)
    inds    = find(Time>TW_av(t,1) & Time<TW_av(t,2));
    CN      = nanmean(nanmean(HscoreMR(:,:,:,inds),4),3);
    HTW     = squeeze(mean(mean(mean(Hf(:,inds,:,:),4),2)));
    
    diags = sum(CN);
    diags=diags-.07;
    CN = (CN-min(CN(:)))./(max(CN)-min(CN(:)))/2;
    CN(1:length(CN)+1:end) = diags;
    subplot(1,size(TW_av,1),t)
    plot_graph(CN,ROISN,Colors,[],.2,round(HTW,2));
    %ylim([-.2 .2]);
    set(gca,'fontsize',14)
    if t==1
        ylabel('Functional Hierarchy Scores')
    end
    title(Tles{t});
end

export_fig(FIG,fullfile(FigPath,['Total_hierarchy_Average_' FileName '_' Method]),'-pdf');close;

%% ------------------------------------------------------------------------
% -----------------------Time Domain Hierarchy ----------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Overall hierarchy score
FIG1 = figure;
set(FIG1,'unit','inch','position',[5 0 8 3],'color','w')
clear CI;
%Time = StokALL.S766640955.Times;
HccfN = (Hccf - mean(Hccf(:,Time>-.2 & Time<0,:),2))./mean(Hccf(:,Time>-.2 & Time<0,:),2)*100;
MH = mean(mean(HccfN(:,:,:),3),1);
CI(1,:) = quantile(mean(HccfN(:,:,:),1),.95,3);
CI(2,:) = quantile(mean(HccfN(:,:,:),1),.05,3);
fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],'k','facealpha',.3,'edgecolor','none');
hold on;
plot(Time,MH,'k','linewidth',1.5);
xlim([-.2 1])
hline(0,'k--')
vline(0,'k--')

xlabel('Time (S)')
title('Total Hierarchy Score % Change')
export_fig(FIG1,fullfile(FigPath,['Total_hierarchy_time_' FileName '_' Method]),'-pdf');close;

%% Hierarchy of ROIs
FIG2 = figure;
set(FIG2,'unit','inch','position',[5 0 8 10],'color','w')
hold on;
clear CI;
for roi = 1:nROIs
    HfR = squeeze(Hf(:,:,roi,:));
    HccfN = (HfR);
    MH = nanmean(nanmean(HccfN(:,:,:),3),1);
    % peastim mean
    MH_prest = mean(MH(Time<0 & Time>-.6));
    line([Time(1) Time(end)],[MH_prest MH_prest],'linestyle','--','color',Colors(roi,:),'linewidth',1.5)
    CI(1,:) = quantile(mean(HccfN(:,:,:),1),.975,3);
    CI(2,:) = quantile(mean(HccfN(:,:,:),1),.025,3);
    fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],[0 0 0],'facealpha',.2,'edgecolor','none');
    LG(roi) = plot(Time,MH,'k','linewidth',1.5,'color',Colors(roi,:));
    xlim([-.2 1])
    %hline(0,'k--')
    vline(0,'k--')
    if roi==nROIs
        xlabel('Time (S)')
        ylabel('Hierarchy Score');
    end
    
    
end
set(gca,'fontsize',FS)
LG = legend(LG,cellfun(@(x) ROI_names.(x),ROIs_Select,'uni',false),'location','northeastoutside');
set(LG,'position',get(LG,'position')+[-.0005 0 0 0])
export_fig(FIG2,fullfile(FigPath,['ROIs_hierarchy_time_' FileName '_' Method]),'-pdf');close;

%% Hierarchy change of ROIs
FIG3 = figure;
set(FIG3,'unit','inch','position',[5 0 6 8],'color','w')
clear CI
for roi = 1:nROIs
    SB = subplot(6,1,roi);
    HfR = squeeze(Hf(:,:,roi,:));
    HccfN = (HfR - mean(HfR(:,Time>-.2 & Time<0,:),2))./abs(mean(HfR(:,Time>-.6 & Time<0,:),2))*100;
    MH = nanmean(nanmean(HccfN(50:100,:,:),3),1);
    CI(1,:) = quantile(mean(HccfN(50:100,:,:),1),.95,3);
    CI(2,:) = quantile(mean(HccfN(50:100,:,:),1),.05,3);
    fill([Time flip(Time)],[CI(1,:) CI(2,end:-1:1)],[0 0 0],'facealpha',.3,'edgecolor','none');
    hold on;
    plot(Time,MH,'k','linewidth',1.5,'color',Colors(roi,:));
    xlim([-.2 1])
    %ylim([-.025 .025])
    %ylim([-60 50])
    %hline(0,'k--')
    line([Time(1) Time(end)],[0 0],'linestyle','--','color',Colors(roi,:),'linewidth',1.5)
    vline(0,'k--')
    if roi==nROIs
        xlabel('Time (msec)')
    end
    if roi==nROIs/2
        ylabel('Hierarchy Score % Change');
    end
    title(ROI_names.(ROIs_Select{roi}))
    set(gca,'fontsize',FS)
    if roi ==nROIs
        set(gca,'xtick',-.2:.2:1,'xticklabel',(-.2:.2:1)*1000)
    else
        set(gca,'xticklabels',[])
    end
    
    set(SB,'position',get(SB,'position')+[0 0 0 .01])
end

export_fig(FIG3,fullfile(FigPath,['ROIs_hierarchy_change_time_' FileName '_' Method]),'-pdf');close;
%% average over prestim
FIG4 = figure;
set(FIG4,'unit','inch','Color','w','position',[1 1 4 3])
Hf_prest = squeeze(mean(mean(Hf(1:end,Time<0 & Time>-.8,:,:),1),2)); 
boxplot(Hf_prest','colors','k','symbol','');%,'Orientation','horizontal')%,'PlotStyle','compact'
box off
h = findobj(gca,'Tag','Box'); 
for j=1:length(h) 
    patch(get(h(j),'XData'),get(h(j),'YData'),Colors(numel(ROIs_Select)-j+1,:),'FaceAlpha',.8);
end 

set(gca,'xtick',1:6,'xticklabel',ROISN,'fontsize',FS-2);

ylabel('Functional Hierarchy Score')
export_fig(FIG3,fullfile(FigPath,['ROIs_hierarchy_Score_Prestimulus_' FileName '_' Method]),'-pdf');close;
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
%caxis([-5 5])
xlim([-.2 2])
vline(0,'k--')
axis xy;
colorbar;

xlabel('Time (S)')
title('Total Hierarchy Score % Change')
colormap('jet')
export_fig(FIG1,fullfile(FigPath,['Total_hierarchy_timefreq_' FileName '_' Method]),'-pdf');close;

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

export_fig(FIG2,fullfile(FigPath,['ROIs_hierarchy_timefrequency_' FileName '_' Method]),'-pdf');close;

%% Hierarchy change of ROIs
CIinterv = .025;
FIG3 = figure;
set(FIG3,'unit','inch','position',[5 0 8 15],'color','w')
for roi = 1:nROIs
    subplot(6,1,roi)
    HfR = squeeze(Hf(:,51:end,roi,:));
    HccfN = (HfR - mean(HfR(:,Time>-.2 & Time<0,:),2))./mean(Hccf(:,Time>-.2 & Time<0,:),2)*100;
    
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
    %caxis([-4 4])
    vline(0,'w--')
    colorbar;
    title(ROI_names.(ROIs_Select{roi}))
end
colormap('jet')
export_fig(FIG3,fullfile(FigPath,['ROIs_hierarchy_change_timefrequency_' FileName '_' Method]),'-pdf');close;
