% This one requires N-way toolbox

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
%addpath(genpath('E:\Elham\Codes\NonGit\nway331'));% add the nway toolbox
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nway331'));% add the nway toolbox
% clear; clc;
FileName = {'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098','drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098'};
%FileName = {'dot_motion__Speed0-01--------0-02--------0-04_iPDC_Mord15_ff098','dot_motion__Speed0-0005-------0-001-------0-005_iPDC_Mord15_ff098'};
DataPath = 'Data_Temp';

%Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
%% Load the iPDCs for all animals
mac = true;
if mac
    Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
else
    Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
end
load(fullfile(Path, ['STOK_ALL_' FileName{1} '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','PARAFAC');

% extract the distance and RF distances
load (fullfile(Path,'Probe_Data_All.mat'));


%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
NROIs = numel(ROIs);
load ROInames;
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
addpath(genpath('../../ARC'));
IDs = fieldnames(StokALL);

% bootstrap
ROIsPerID = cellfun(@(x) StokALL.(x).ROIs,IDs,'uni',false);
nboots = 10;
bootsize = 8;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);

redoanalysis = false;

%mode = '';% just unfolding over inter-area connectivity
mode = 'unfoldedlam';%unfolding over laminar and inter-area connectivity but separately
%mode = 'unfoldedall';
% Time windows for variance explained
TW = [-300 0; 50 150; 150 250; 250 1000];

for cond = 1:2
    load(fullfile(Path, ['STOK_ALL_' FileName{cond} '.mat']));
    StokALL = DistanceEstimate(StokALL,Probe_all);
    StokALL = RFDistanceEstimate(StokALL,Probe_all);
    for NComp = 3:3
        if redoanalysis || ~exist([FileName{cond} ['PARAFAC_covtemp_' mode 'N_'] num2str(NComp) '_ExtraVar.mat'],'file')
            PARAFAC_FC(StokALL,NComp,BootIDs,nboots,ROIs,DataPath,FigPath,FileName{cond},redoanalysis,mode,TW);    
        end
        load(fullfile(DataPath,[FileName{cond} ['PARAFAC_covtemp_' mode num2str(NComp)]]));
        %DistanceBoots = DistanceForBoot(StokALL,BootIDs);
        PARRES{cond} = load(fullfile(DataPath,[FileName{cond} ['PARAFAC_covtemp_' mode 'N_'] num2str(NComp) '_ExtraVar.mat']));
        PARRES{cond}.DistanceBoots = DistanceForBoot(StokALL,BootIDs);
        disp(['NComp = ' num2str(NComp) ', Corcondia = ' num2str(mean(corcondia)) ', Error = ' num2str(mean(err)) ', it = ' num2str(mean(it))]);
        COR(NComp,:) = corcondia;
        ERR(NComp,:) = err;
        IT(NComp,:) = it;
        VAREXP(NComp,:) = 1-err./SSX;
    end
end
%% Number of components selection
FIG = figure;
FS = 12;
set(FIG,'unit','inch','position',[6 6 8 6],'color','w');
X = repmat((1:NComp)',1,nboots)+rand(NComp,nboots)*.2;

subplot(2,2,1)
scatter(X(:),COR(:),'filled');
ylabel('CORCONDIA')
set(gca,'fontsize',FS,'xtick',1:5,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])

subplot(2,2,2)
scatter(X(:),ERR(:),'filled');
ylabel('MSE')
set(gca,'fontsize',FS,'xtick',1:5,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])

subplot(2,2,3)
scatter(X(:),IT(:),'filled');
ylim([0,500])
ylabel('# of Iterations');
xlabel('# of Components')
set(gca,'fontsize',FS,'xtick',1:5,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])

subplot(2,2,4)
scatter(X(:),VAREXP(:),'filled');
ylabel('Variance Explained (%)');
set(gca,'fontsize',FS,'xtick',1:5,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])


export_fig(FIG,fullfile(FigPath,[FileName{1} '_numberofcomoponent']),'-pdf','-r300')

%% variance explained over time
NComp = 3;
redoanalysis = false;
mode = '';
for cond = 1:1
    load(fullfile(Path, ['STOK_ALL_' FileName{cond} '.mat']));    
    load(fullfile(DataPath,[FileName{cond} ['PARAFAC_covtemp_' mode num2str(NComp)]]));
    %DistanceBoots = DistanceForBoot(StokALL,BootIDs);
    PARRES{cond} = load(fullfile(DataPath,[FileName{cond} ['PARAFAC_covtemp_' mode 'N_'] num2str(NComp) '_ExtraVar.mat']));
        
end

addpath(genpath('/Users/elhamb/Documents/Codes/Git/plotting/boxplot2-pkg'))
close;
Freq = 1:100;
con_mode = 1;
%% plotting params
offs = -.06;
FS = 12;
FIG = figure(1);
set(FIG,'unit','inch','position',[0 0 5 6],'color','w')
lstyle = {'-','-'};
lcolor = {'k',[.5 .5 .5]};
connames = {'High Contrast','Low Contrast'};
Compcol = brewermap(10,'Paired');
Compcol = Compcol([8 6 10 2],:);

TLabels = {'[-300  0] msec','[50  150] msec','[150  250] msec','[250  1000] msec'};
NLabels = {'Network1','Network2','Network3'};
% Time windows
TW = [-300 0; 50 150; 150 250; 250 1000; -300 1000];

for cond = 1:2
    subplot(2,1,cond)
    Model_reord =   PARRES{cond}.Model_reord;
    Comp_ord    =   PARRES{cond}.Comp_ord;
    % temporal dynamics: significance
    M_temp = cellfun(@(x) (x{numel(Model_reord{1})}),Model_reord,'uni',false);
    M_temp = cellfun(@(x) ((x)),M_temp,'uni',false);
    M_temp = cat(3,M_temp{:});
    % convert to variance explained
    
    Inds = arrayfun(@(x) find(round(temp_time,2)*1000==x,1),TW);
    Data = arrayfun(@(x) (mean(M_temp(Inds(x,1):Inds(x,2),:,:).^2)./sum(mean(M_temp(Inds(x,1):Inds(x,2),:,:).^2,1),2))*100,1:size(TW,1),'uni',false);
    Total = mean(M_temp(:,:,:).^2)./sum(mean(M_temp(:,:,:).^2,1),2)*100;
    mean(Total,3)
    std(Total,[],3)
    Data = cat(1,Data{:});
    %Data(:,1,:)= Data(:,1,:)-40;
    
    h = boxplot2(permute(log10(Data(1:4,:,:)),[2 1 3]),1:3);
    mean(Data,3)
    %boxplot(reshape(Data,16,500)')
    axis tight
    % correct colors
    for ii = 1:4
        structfun(@(x) set(x(ii,:), 'color', Compcol(ii,:), ...
            'markeredgecolor', Compcol(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-','linewidth',1.5);
    set(h.box, 'linestyle', '-','linewidth',1.5);
    set(h.out, 'marker', 'none');
    
    % axis info
    title(connames{cond})
    set(gca,'ylim',[log10(2) log10(100)]);
    
    if cond ==2
       set(gca,'xtick',1:NComp,'xticklabel',NLabels,'ytick',log10([5 10 20 40 80]),'yticklabel',[5 10 20 40 80]) 
       set(gca,'position',get(gca,'position')+[0 .05 0 0])
       xtickangle(45)
       ylabel('Variance Explained (%)')
       
    else
        set(gca,'xtick',1:NComp,'xticklabel',[],'ytick',log10([5 10 20 40 80]),'yticklabel',[5 10 20 40 80]);
        LG = legend(TLabels);
        legend box off
        set(LG,'position',get(LG,'position')+[.05 0 0 0])
        LG.ItemTokenSize(1)=15;
    end
    %hline(19,'k--')
    set(gca,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])
end
    
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_Variance']),'-pdf','-r200')

