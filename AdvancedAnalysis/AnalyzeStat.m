clear;clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
clear; clc;
FileName = 'dot_motion__Speed0-01--------0-02--------0-04_iPDC_Mord15';%'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%%
load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults');
%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
[PDC,IDs,Time,Freq,ROIs] = ExtractAllRoiPDC(StokALL,ROIs);
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;
clear StokALL;

%% To indicate which part of the evoked PDC is significant
for roi = 1:nROIs
    [~,FIG1] = CompPostPre(PDC(roi,roi),[1],[],1,1,Time,Freq,['Output-intraPDC- ' ROIs{roi} '_' FileName],FigPath);
    [~,FIG2] = CompPostPre(PDC(setdiff(1:nROIs,roi),roi),[1],[],1,1,Time,Freq,['Output-interPDC- ' ROIs{roi} '_' FileName],FigPath);
end

% check separately the evoked FF and FB
ind = 1;
for roi1 = 1:numel(ROIs)
    for roi2 = roi1+2:numel(ROIs)
        %if roi1~=4 && roi2~=4
            PDCFF{ind} = PDC{roi2,roi1};
            PDCFB{ind} = PDC{roi1,roi2};
            ind = ind+1;
        %end
    end
end

[~,FIG2] = CompPostPre(PDCFF,[1],[],1,1,Time,Freq,['Output-interPDC-output -FF_' FileName],FigPath);
[~,FIG2] = CompPostPre(PDCFF,[2],[],1,1,Time,Freq,['Output-interPDC-input -FF_' FileName],FigPath);

[~,FIG2] = CompPostPre(PDCFB,[1],[],1,1,Time,Freq,['Output-interPDC-output -FB_' FileName],FigPath);
[~,FIG2] = CompPostPre(PDCFB,[2],[],1,1,Time,Freq,['Output-interPDC-input -FB_' FileName],FigPath);

%% To compare the evoked PDC in different connections intra vs. inter

% First, the PDCs should be converted to % change
PDC_change = cellfun(@(x) (x - mean(x(:,:,:,Time<0 & Time>-.3,:),4))./(mean(x(:,:,:,Time<0 & Time>-.3,:),4)),PDC,'uni',false);


% layerwise comparison -> intra vs. inter
for roi = 1:nROIs
    CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:nROIs,roi),roi),[1],false,1,[1 0],Time,Freq,['Output-(Intra vs. Inter) - ' ROIs{roi} '_' FileName],FigPath);
end
% Compare all rois included
ind = 1;
ind2 =1;
for roi1 = 1:numel(ROIs)
    PDCintra{ind} = PDC_change{roi1,roi1};
    PDCinter(ind2:ind2+nROIs-2) = PDC_change(setdiff(1:nROIs,roi),roi1);
    ind = ind+1;
    ind2 = ind2 +nROIs-1;
end
% TOTAL intra vs. inter
CompPDCCells(PDCintra,PDCinter,[1],false,1,[1 0],Time,Freq,['Output-Averaged-(Intra vs. Inter)- All_' FileName],FigPath,false);
CompPDCCells(PDCintra,PDCinter,[],false,1,[1 0],Time,Freq,['Outputinput-Averaged-(Intra vs. Inter)- All_' FileName],FigPath,false);


% average comparison -> intra vs. inter
for roi = 1:nROIs
    CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:nROIs,roi),roi),[1 2],false,1,[1 0],Time,Freq,['Output-Averaged-(Intra vs. Inter)- ' ROIs{roi} '_' FileName],FigPath);
end

% all layer comparison -> intra vs. inter
CompPDCCells(PDC_change(1,1),PDC_change(2:end,1),[],false,1,[1 0],Time,Freq,['All-InterIntra - ' ROIs{1} '_' FileName]);

%% compare input vs. output

for roi = 1:nROIs
    CompPDCCells(PDC_change(roi,setdiff(1:nROIs,roi)),PDC_change(setdiff(1:nROIs,roi),roi),[1 2],false,1,[1 0],Time,Freq,['Output-(Input vs. output) - ' ROIs{roi} '_' FileName],FigPath);
end

%% average comparison -> FF vs. FB 
PDC_change2 = PDC_change;
%PDC_change2 = PDC_change2([1:4 6 5 7],[1:4 6 5 7]);
ind = 1;
for roi1 = 1:numel(ROIs)
    for roi2 = roi1+1:numel(ROIs)
        PDCFF{ind} = PDC_change2{roi2,roi1};
        PDCFB{ind} = PDC_change2{roi1,roi2};
        ind = ind+1;
    end
end
% TOTAL FB vs. FF
CompPDCCells(PDCFF,PDCFB,[1 2],false,1,[0 0],Time,Freq,['Output-Averaged-(FF vs. FB)- All_' FileName],FigPath,false);

CompPDCCells(PDCFF,PDCFB,[2],false,1,[0 0],Time,Freq,['Output-Averaged-(FF vs. FB)- layers-input - All_' FileName],FigPath,false);
CompPDCCells(PDCFF,PDCFB,[1],false,1,[0 0],Time,Freq,['Output-Averaged-(FF vs. FB)- layers-output - All_' FileName],FigPath,false);

CompPDCCells(PDCFF,PDCFB,[],false,1,[0 0],Time,Freq,['Output-Averaged-(FF vs. FB)- layers-outputinput - All_' FileName],FigPath,false);


% average comparison -> FF vs. FB
for roi = 2:nROIs-1
    CompPDCCells(PDC_change(roi+1:nROIs,roi),PDC_change(1:roi-1,roi),[1 2],false,1,[0 0],Time,Freq,['Output-Averaged-(FF vs. FB)- ' ROIs{roi} '_' FileName],FigPath);
end

for roi = 2:nROIs-1
    CompPDCCells(PDC_change(roi+1:nROIs,roi),PDC_change(1:roi-1,roi),[1],false,1,[0 0],Time,Freq,['Output-(FF vs. FB)- ' ROIs{roi} '_' FileName],FigPath);
end

%% In this part, I use two-way anova to compare inter-intra , pre-post PDC values

for roi = 2:7
    CompPDCCells(PDC(roi,roi),PDC(setdiff(1:7,roi),roi),[1 2],true,1,[1 0],Time,Freq,['Output-Averaged-(Intra vs. Inter)-ANOVA - ' ROIs{roi}],FigPath);
end


%% ----------------------------------------------------------------------
%------------------------------------------------------------------------
%---------------To correlate with hierarchy scores-----------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%% for each animal calculate the slope of intra-inter separately and then correlate that with hierarchy scores
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_ALL_iPDC.mat');

ids = fieldnames(StokALL);
ROIsizes = contains(ROIs,'VIS')*3+3;
ROIindices = [0 cumsum(ROIsizes)];



for s = 1:numel(ids)
    % (1) for each animal and roi, calculate the intra/intra+inter with prestim
    % removed
    
    % convert to % change
    %PDC_change_roi = (StokALL.(ids{s}).PDC - mean(StokALL.(ids{s}).PDC,4))./mean(StokALL.(ids{s}).PDC,4);
    PDC_change_roi = (StokALL.(ids{s}).PDC);
    
    
    [C,ia] = intersect(ROIs,StokALL.(ids{s}).ROIs);
    [~,I] = sort(ia);
    rois = C(I);
    Labels = ia(I);
    
    indices  =  arrayfun(@(x) ROIindices(x)+1:ROIindices(x+1),Labels,'uni',false);
    indices  = cat(2,indices{:});    
    
    for roi1 = 1:numel(ROIs)
        ind = (roi1-1)*6+1:roi1*6;
        [im]= ismember(indices,ind);
        lind = find(im);
        if isempty(lind)
            % it puts NaN values...
        end
        PDC_out_total{roi1,s} = squeeze(mean(mean(PDC_change_roi(:,lind,:,:))));
        PDC_in_total{roi1,s} = squeeze(mean(mean(PDC_change_roi(lind,:,:,:))));
        PDC_out_intra{roi1,s} = squeeze(mean(mean(PDC_change_roi(lind,lind,:,:))));
        PDC_out_inter{roi1,s} = squeeze(mean(mean(PDC_change_roi(setdiff(1:size(StokALL.(ids{s}).PDC,1),lind),lind,:,:))));
    end
    
    % (2) load hierarchy scores and correlate them
end


load Hierarchyscores;
HS = mean(H(Time>-.3 & Time<0,:));
% PDC_perinter = cellfun(@(x,y) x-y,PDC_out_intra,PDC_out_inter,'uni',false);
% PDC_out_intra = cellfun(@(x) (x-mean(x(:,Time>-.3 & Time<0),2))./mean(x(:,Time>-.3 & Time<0),2),PDC_out_intra,'uni',false);
% PDC_out_inter = cellfun(@(x) (x-mean(x(:,Time>-.3 & Time<0),2))./mean(x(:,Time>-.3 & Time<0),2),PDC_out_inter,'uni',false);


PDC_perinter = cellfun(@(x,y) x./(x+y),PDC_out_intra,PDC_out_inter,'uni',false);
PDC_perinter = cellfun(@(x) (x-mean(x(:,Time>-.3 & Time<0),2))./mean(x(:,Time>-.3 & Time<0),2),PDC_perinter,'uni',false);

% -------------------- also need to fit a line...
TimeWin = Time>.1 & Time<1;
FIG = figure;
for roi = 1:numel(ROIs)
     PDC_temp = PDC_perinter(roi,:);
     PDC_temp = nanmean(cat(3,PDC_temp{:}),3);

    subplot(6,1,roi);

    imagesc(squeeze(PDC_temp));
    %plot(mean(PDC_temp));hold on;
    axis xy;
    %ylim([-.15 .05])
    caxis([-.08 .08])
    xlim([find(round(Time,2)==-.1,1) find(round(Time,2)==1,1)]);
    %hline(0,'k--')
    % estimate coefficients
    PDC_temp = PDC_perinter(roi,:);
    Data = squeeze(nanmean(cat(3,PDC_temp{:})));
    for s = 1:11
        
        coefficients = polyfit(Time(TimeWin)', Data(TimeWin,s), 1);
        Coef(:,s,roi) = coefficients;
    end

end
colormap('jet')
%--------------------------------------------------------------------------
FIG = figure;
set(FIG,'unit','inch','position',[5 5 5 5],'color','w');
ZPoint = squeeze(-Coef(2,:,:)./Coef(1,:,:))*1000;
ZPoint(9:10,5)=NaN;%outlier
Ys = [];
Xs = [];
for roi = 1:nROIs
    scatter(repmat(HS(roi),[1 11]), ZPoint(:,roi),50,Colors(roi,:),'filled');hold on;% color and filled, fit line, correlation
    Xs = [Xs squeeze(repmat(HS(roi),[1 11]))];
    Ys = [Ys squeeze(ZPoint(:,roi))'];
end
%ylim([0 1400])
Xs(isnan(Ys))=[];
Ys(isnan(Ys))=[];
[r,p] = corr(Xs',Ys');
coefficients = polyfit(Xs', Ys', 1);
yFit = polyval(coefficients , -.2:.01:.2);
hold on;
plot(-.2:.01:.2,yFit,'k--','linewidth',1.5);
xlabel('Hierarchy Score');
ylabel('Time to baseline within-between PDC (msec)');
set(gca,'fontsize',18)
title(['r = ' num2str(round(r,2)) ', p = ' num2str(round(p,5))])
%export_fig(fullfile(FigPath,'Time_baseline_interintra'),'-pdf')

%% ------------------------PDC intput and output correlation with HS-----------------------
PDC_out = cellfun(@(x) nanmean(nanmean(x,2),1) , PDC_out_total);
PDC_in  = cellfun(@(x) nanmean(nanmean(x,2),1) , PDC_in_total);
FIG = figure;
set(FIG,'unit','inch','position',[5 5 5 4],'color','w');
Correlate_hierarchy_plot(PDC_out' ,Time ,Colors)
ylabel('Average output PDC') 
export_fig(fullfile(FigPath,'Hierarchy_correlation_AverageOutput'),'-pdf')
 
FIG = figure;
set(FIG,'unit','inch','position',[5 5 5 4],'color','w');
Correlate_hierarchy_plot(PDC_in' ,Time ,Colors);
ylabel('Average input PDC') 
export_fig(fullfile(FigPath,'Hierarchy_correlation_AverageInput'),'-pdf')

%% with all animals pulled together
Colors = COBJ.MatrixColors(ROIs);
Ls = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
% layerwise comparison -> intra vs. inter
tic
for roi = 1:nROIs
    Stats_interintra{roi} = CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:nROIs,roi),roi),[1 2],false,false,[1 0],Time,Freq,[ROI_names.(ROIs{roi})],FigPath);
    %plot(Time,squeeze(mean(Stats{1}.Tval(:,:,30:end,:),3)),colors(roi,:))
end
toc

%% Average fit line plot
TimeWin = Time>.05 & Time<.6;
Xticklabels = -.2:.200:1.000;
Xticks = arrayfun(@(x) find(round(Time,2)==x,1),Xticklabels);

FIG = figure;
for roi = 1:nROIs
    Data = squeeze(mean(Stats_interintra{roi}.Tval(:,:,30:end,:),3));
    subplot(nROIs,1,roi)
    plot(Time*1000,Data,'color',Colors(roi,:),'linewidth',2); hold on;
    
    coefficients = polyfit(Time(TimeWin)', Data(TimeWin), 1);

    yFit = polyval(coefficients , Time(Time>.05 & Time<1));
    hold on;
    plot(Time(Time>.05 & Time<1)*1000, yFit, '--','color','k', 'LineWidth', 2);
    grid on;
    ylim([-8 5])
    xlim([-.1 1]*1000);
    text(600,-7,['Slope = ' num2str(coefficients(1))],'fontsize',14)
    title(ROI_names.(ROIs{roi}))
    %Cof(:,roi) = coefficients;
    set(gca,'fontsize',16);
    if roi==nROIs
        xlabel('Time (msec)');
        ylabel('Tstats');
    end
    set(gca,'xtick',Xticklabels*1000);
end

set(FIG,'unit','inch','position',[1 1 4 16],'color','w')
export_fig(fullfile(FigPath,['Output-(Intra vs. Inter) - fitline - All']),'-pdf','-r200');close;

figure,scatter(HS,-Cof(2,:)./Cof(1,:));
[r,p] = corr(HS',(-Cof(2,:)./Cof(1,:))');

%% aggregate stats and correlate with hierarchy at time and frequency points

for roi  = 1:nROIs
    TVF(:,:,roi) = squeeze(Stats_interintra{roi}.Tval);
end

for f = 1:size(TVF,1)
    for t = 100:size(TVF,2)
        [r(f,t),p(f,t)] = corr(squeeze(TVF(f,t,:)),HS');
    end
end

FIG = figure;
set(FIG,'unit','inch','position',[1 1 5 2],'color','w')
h = imagesc(r);
set(h,'alphadata',(p<0.01)*.8+.2);
Xticklabels = -.2:.200:1.000;
Xticks = arrayfun(@(x) find(round(Time,2)==x,1),Xticklabels);
Yticklabels = 10:20:100;
Yticks = arrayfun(@(x) find(round(Freq,2)==x,1),Yticklabels);
set(gca,'xtick',Xticks,'xticklabel',Xticklabels*1000,'ytick',Yticks,'yticklabel',Yticklabels,'fontsize',16);
xlim([find(round(Time,2)==-.1,1) find(round(Time,2)==1,1)])
xlabel('Time (S)');
ylabel('Freq (Hz)');
axis xy;
caxis([-1 1])
colormap('jet')
colorbar;
export_fig(fullfile(FigPath,['Output-(Intra vs. Inter) - Correlation - All']),'-pdf','-r200');close;


%% --------------------------------------------------------------------------
%---------------------------------layer output - inter intra-------------
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs,'subcolors');
Ls = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);


figure,
TimeWin = Time>-.3 & Time<.0;
S = 8;
for roi = 1:nROIs
    MW = mean(mean(mean(PDC{roi,roi}(:,:,30:end,TimeWin,:),5),4),3);
    MB = cellfun(@(x) mean(mean(mean(x(:,:,30:end,TimeWin,:),5),4),3),PDC(setdiff(1:nROIs,roi),roi),'uni',false);
    MB = nanmean(cat(3,MB{:}),3);
    
    MW(1:length(MW)+1:end)=0;
    MB(1:length(MB)+1:end)=0;
    
    subplot(nROIs,2,(roi-1)*2+1);
    B = bar(sum(MW));
    ylim([0 1]/2);
    ylabel(ROI_names.(ROIs{roi}),'fontweight','bold')
    if roi==1
        title('Within-area PDC')
    end
    B.FaceColor = 'flat';
    for l = 1:6
        B.CData(l,:) = Colors(roi,l,:);
    end
    if roi ==nROIs
        set(gca,'xtick',1:6,'xticklabel',Ls);
    else
        set(gca,'xtick',1:6,'xticklabel',[]);
    end
    
    %----------------------------------------------------------------------
    subplot(nROIs,2,(roi-1)*2+2);
    B = bar(sum(MB));
    ylim([0 1]/2);
    if roi == 1
        title('Between-area PDC')
    end
    B.FaceColor = 'flat';
    for l = 1:6
        B.CData(l,:) = Colors(roi,l,:);
    end
    
    if roi ==nROIs
        set(gca,'xtick',1:6,'xticklabel',Ls);
    else
        set(gca,'xtick',1:6,'xticklabel',[]);
    end
end

FIG = figure;
set(FIG,'unit','inch','position',[5 5 5 3],'color','w')
for roi = 1:nROIs
    MW = mean(mean(mean(PDC{roi,roi}(:,:,:,TimeWin,:),5),4),3);
    SW = std(mean(mean(PDC{roi,roi}(:,:,:,TimeWin,:),3),4),[],5);
    MB = cellfun(@(x) mean(mean(mean(x(:,:,:,TimeWin,:),5),4),3),PDC(setdiff(1:nROIs,roi),roi),'uni',false);
    MB = nanmean(cat(3,MB{:}),3);
    SB = cellfun(@(x) std(mean(mean(x(:,:,:,TimeWin,:),3),4),[],5),PDC(setdiff(1:nROIs,roi),roi),'uni',false);
    SB = nanmean(cat(3,SB{:}),3);
    
    
    MW(1:length(MW)+1:end)=NaN;
    MB(1:length(MB)+1:end)=NaN;
    SW(1:length(MW)+1:end)=NaN;
    SB(1:length(MB)+1:end)=NaN;
    
    MWA(roi) = nanmean(MW(:));
    MBA(roi) = nanmean(MB(:));
    SWA(roi) = nanmean(SW(:));
    SBA(roi) = nanmean(SB(:));
end

Cs = [.2 .2 1; .7 .2 .2];
BR = bar([MBA; MWA]','FaceColor','flat');
for i = 1:numel(BR)
    BR(i).CData = Cs(i,:);%squeeze(Colors(:,(i-1)*5+1,:));
end
hold on;
errorbar((1:6)-.15,MBA,SBA,'.','color','k','linewidth',1.5)
errorbar((1:6)+.15,MWA,SWA,'.','color','k','linewidth',1.5)
%colormap('jet');
legend('Between-area','Within-area');
ylabel('Output PDC')
title('Averaged pre-stimuli PDCs');
set(gca,'xticklabel',arrayfun(@(x) ROI_names.(ROIs{x}),1:nROIs,'uni',false));
set(gca,'fontsize',14)

export_fig(fullfile(FigPath,['Output-(Intra vs. Inter) - Average']),'-pdf','-r200');close;



%% I was trying to look at the intra-laminar connectivity patterns in different time points of FF abd FB -> what about decomopsition?
Colors = COBJ.MatrixColors(ROIs,'subcolors');

figure,
for roi = 1:6
    subplot(2,3,roi)
    PDC_temp = squeeze(mean(mean(mean(PDC{roi,roi}(:,:,:,Time>-.3 & Time<0,:),3),4),5));
    PDC_temp(1:7:end)=0;
    PDC_temp = (PDC_temp-min(PDC_temp(:)))./(max(PDC_temp(:))-min(PDC_temp(:)));
    PDC_temp(1:7:end)=sum(PDC_temp,1);
    plot_graph(PDC_temp,{'L1','L2','L3','L4','L5','L6'},squeeze(Colors(roi,:,:)),'down',0.7)
    
    PDC_temp_boot = bootstrap_PDC(PDC_change{roi,roi},100);
    PDC_boot_base = mean(PDC_temp_boot(:,:,:,Time>0.05 & Time<=0.1,:),4);
    for f = 1:100
        f
        for t = 250:numel(Time)
            for i = 1:6
                for j = 1:6
                    [h,p(i,j,f,t),~,Stats] = ttest(PDC_boot_base(i,j,f,1,:),PDC_temp_boot(i,j,f,t,:));
                    Tval(i,j,f,t) = Stats.tstat;
                end
            end
        end
    end
   
end


    
dynet_connplot(-1*Tval,Time,Freq)
 
M           = mean(mean(Tval(:,:,:,Time>0.05 & Time<.1),3),4)*-1/20;
M(1:7:end)  = 0;
M(M<0)      = 0;
plot_graph(M,{'L1','L2','L3','L4','L5','L6'},squeeze(Colors(roi,:,:)),'down',0.7)



M           = mean(mean(Tval(:,:,1:20,Time>0.1 & Time<.3),3),4)*-1/20;
M(1:7:end)  = 0;
M(M<0)      = 0;
figure,plot_graph(M,{'L1','L2','L3','L4','L5','L6'},squeeze(Colors(roi,:,:)),'down',0.7)


