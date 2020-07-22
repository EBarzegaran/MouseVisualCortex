
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098';%'_dot_motion__Speed0-01--------0-02--------0-04_iPDC_Mord15';%'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord10';%
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%%
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/dynet_toolbox-master'));
fileparts(mfilename('fullpath'))

load([Path 'STOK_Average_' FileName '.mat']);
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
SavePath = Path;
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs);
load ROInames;
NROIs = numel(STOK_avg.ROIs);
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
%% plot each ROI separately

MODE = 2;
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);
if MODE==1
    TW = [0 2];
else
    TW = [-.2 1];
end

for roi = 1:NROIs
    ind = [(roi-1)*6+1:roi*6];
    FIG = dynet_connplot(PDC(ind,ind,:,Time>TW(1) & Time<TW(2)),Time(Time>TW(1) & Time<TW(2)),Freq,cellfun(@(x) strrep(x,'_','-'),STOK_avg.labels(ind),'uni',false),[],[],[],1,STOK_avg.ROIs{roi});
    set(FIG,'unit','inch','position',[0,0,20,12],'color','w')
    if MODE==1
        export_fig(FIG,fullfile(SavePath,['STOK_individual_sustained_' STOK_avg.ROIs{roi} FileName]),'-pdf'); %close; 
    else
        export_fig(FIG,fullfile(SavePath,['STOK_individual_transient_' STOK_avg.ROIs{roi} FileName]),'-pdf'); %close; 
    end
end
%% Plot averaged oevr layers of all ROIs connectivity
FS = 12;
MODE = 2;
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = (STOK_avg.PDC);
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4)*100;
if MODE==1
    TW = [0 2];
else
    TW = [-.2 1];
end
TimeInd = (Time>TW(1) & Time<TW(2));
clear PDC_avg;
for roi = 1:NROIs
    for roi2 = 1:NROIs
        Ind = ((roi-1)*6+1):(roi*6);
        Ind2 = ((roi2-1)*6+1):(roi2*6);
        if roi~=roi2
            PDC_avg(roi,roi2,:,:) = squeeze(sum(sum(PDC(Ind,Ind2,:,TimeInd),1),2)./(numel(Ind)*numel(Ind2)));
        else
            Temp=zeros(size(PDC,3),sum(TimeInd));
            for i=1:numel(Ind)
                Temp = Temp+squeeze(PDC(Ind(i),Ind(i),:,TimeInd));
            end
            Temp = Temp./numel(Ind);
            PDC_avg(roi,roi2,:,:) = Temp;
        end
    end
end

for l1 = 1:6
    for l2 = 1:6
            Ind = (0:(NROIs-1))*6+l1;
            Ind2 = (0:(NROIs-1))*6+l2;
        if l1~=l2
            PDC_layer_avg(l1,l2,:,:) = squeeze(sum(sum(PDC(Ind,Ind2,:,TimeInd),1),2)./(numel(Ind)*numel(Ind2)));
        else
            Temp=zeros(size(PDC,3),sum(TimeInd));
            for i=1:numel(Ind)
                Temp = Temp+squeeze(PDC(Ind(i),Ind(i),:,TimeInd));
            end
            Temp = Temp./numel(Ind);
            PDC_layer_avg(l1,l2,:,:) = Temp;
        end
    end
end

FIG1 = dynet_connplot(PDC_avg,Time(TimeInd),Freq,cellfun(@(x) ROI_names.(x), STOK_avg.ROIs,'uni',false),[],[],[],0);
colormap('jet')
set(FIG1,'unit','inch','position',[1 1 6 3],'color','w')
export_fig(FIG1,fullfile(SavePath,['STOK_AllROIs' FileName]),'-pdf');

FIG1 = figure;
imagesc(Time(TimeInd),Freq,squeeze(mean(mean(PDC_avg,1),2)));
axis xy
caxis([-10 50])
colormap('jet')
set(FIG1,'unit','inch','position',[1 1 5 2.8],'color','w')
set(gca,'fontsize',FS)
xlabel('Time(S)')
ylabel('Frequency(Hz)')
vline(0,'w-')
CL = colorbar;
set(get(CL,'title'),'string','%Change')
export_fig(FIG1,fullfile(SavePath,['STOK_AverageAll' FileName]),'-pdf');


FIG1 = dynet_connplot(PDC_layer_avg,Time(TimeInd),Freq,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false),[],[],[],1);
colormap('jet')
set(FIG1,'unit','inch','position',[1 1 6 3],'color','w')
export_fig(FIG1,fullfile(SavePath,['STOK_AllLayers' FileName]),'-pdf');


%PDC_avg_diff = PDC_avg - permute(PDC_avg,[2 1 3 4]);
%% plot PSDs for ROIs and layers
FIG = figure;
FS=12;
set(FIG,'unit','centimeters','position',[1 1 6 8]*2,'color','w')
M = max(PDC_avg(:))/200;

m = min(PDC_avg(:))/1;
for roi = 1:NROIs
    subplot(NROIs,2,(roi-1)*2+1)
    imagesc(Time(TimeInd),Freq,squeeze((PDC_avg(roi,roi,:,:))))
    axis xy;
    caxis([m M])
 
    set(gca,'position',get(gca,'position')+[-0.03 -.02 0 .03])
    set(gca,'fontsize',FS)
    if roi~=NROIs
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','lineWidth',1.5,'yticklabel',[])

    end
    text(-.46,50,ROISN{roi},'color',Colors(roi,:),'fontsize',FS+2,'fontweight','bold')
    vline(0,'--w')
    
end

for l = 1:6
    subplot(6,2,(l-1)*2+2)

    %yyaxis right
    set(gca,'YColor','k');
    imagesc(Time(TimeInd),Freq,squeeze((PDC_layer_avg(l,l,:,:))))
    %yyaxis left
    axis xy;
    caxis([m M])
    if l==1
        Pos_temp = get(gca,'position'); 
        colorbar
        set(gca,'position',Pos_temp)
    end
    set(gca,'position',get(gca,'position')+[-0.03 -.02 0 .03])
    set(gca,'fontsize',FS)
    if l~=6
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','yticklabel',[],'lineWidth',1.5)
        yyaxis right
        set(gca,'YColor','k');
        set(gca,'ytick',.2:.2:1,'yticklabel',(.2:.2:1)*100)
        
        YL=ylabel('Frequency(Hz)');
        set(YL,'position',get(YL,'position')+[-0.07 0 0])
        xlabel('Time(sec)');
        yyaxis left
    end
    %
    if l==1
       
        TL = title('PSD');
        set(TL,'position',get(TL,'position')+[-.8 10 0])
    end
    
    text(-.35,50,['l' num2str(l)],'color','k','fontsize',FS+2,'fontweight','bold')
    vline(0,'--w')
    
end

colormap('jet')
export_fig(FIG,fullfile(SavePath,['STOK_PSD_AllLayers_AllROIs' FileName]),'-pdf');

%% plot output iPDC
FIG = figure;
FS=12;
set(FIG,'unit','centimeters','position',[1 1 6 8]*2,'color','w')
M = max(PDC_avg(:))/300;
for roi = 1:NROIs
    subplot(NROIs,2,(roi-1)*2+1)
    imagesc(Time(TimeInd),Freq,squeeze(mean(PDC_avg(:,roi,:,:))))
    axis xy;
    caxis([0 M])
 
    set(gca,'position',get(gca,'position')+[-0.03 -.02 0 .03])
    set(gca,'fontsize',FS)
    if roi~=NROIs
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','lineWidth',1.5,'yticklabel',[])

    end
    text(-.46,50,ROISN{roi},'color',Colors(roi,:),'fontsize',FS+2,'fontweight','bold')
    vline(0,'--w')
    
end

for l = 1:6
    subplot(6,2,(l-1)*2+2)

    %yyaxis right
    set(gca,'YColor','k');
    imagesc(Time(TimeInd),Freq,squeeze(mean(PDC_layer_avg(:,l,:,:))))
    %yyaxis left
    axis xy;
    caxis([0 M])
    if l==1
        Pos_temp = get(gca,'position'); 
        colorbar
        set(gca,'position',Pos_temp)
    end
    set(gca,'position',get(gca,'position')+[-0.03 -.02 0 .03])
    set(gca,'fontsize',FS)
    if l~=6
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','yticklabel',[],'lineWidth',1.5)
        yyaxis right
        set(gca,'YColor','k');
        set(gca,'ytick',.2:.2:1,'yticklabel',(.2:.2:1)*100)
        
        YL=ylabel('Frequency(Hz)');
        set(YL,'position',get(YL,'position')+[-0.07 0 0])
        xlabel('Time(sec)');
        yyaxis left
    end
    %
    if l==1
       
        TL = title('iPDC output');
        set(TL,'position',get(TL,'position')+[-.8 10 0])
    end
    
    text(-.35,50,['l' num2str(l)],'color','k','fontsize',FS+2,'fontweight','bold')
    vline(0,'--w')
    
end

colormap('jet')
export_fig(FIG,fullfile(SavePath,['STOK_iPDCout_AllLayers_AllROIs' FileName]),'-pdf');

%% plot laminar iPDCs for an examples
FIG = figure;
FS=12;
set(FIG,'unit','centimeters','position',[1 1 6 10]*2,'color','w')
M = max(PDC_avg(:))/300;
ind1 = 1:6;
ind2 = 7:12;
for l = 1:6
    subplot(NROIs,2,(l-1)*2+1)
    imagesc(Time(TimeInd),Freq,squeeze(mean(PDC(ind2,ind1(l),:,Time>TW(1) & Time<TW(2)))))
    axis xy;
    caxis([0 M])
 
    set(gca,'position',get(gca,'position')+[-0.03 -.03 0 .02])
    set(gca,'fontsize',FS)
    if l~=6
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','lineWidth',1.5,'yticklabel',[])

    end
    if l==1
        title('source(V1)')
    end
    text(-.45,50,['L' num2str(l)],'color','k','fontsize',FS+2,'fontweight','normal')
    vline(0,'--w')
    
end

for l = 1:6
    subplot(6,2,(l-1)*2+2)

    %yyaxis right
    set(gca,'YColor','k');
    imagesc(Time(TimeInd),Freq,squeeze(mean(PDC(ind2(l),ind1,:,Time>TW(1) & Time<TW(2)),2)))
    %yyaxis left
    axis xy;
    caxis([0 M])
    if l==1
        Pos_temp = get(gca,'position'); 
        colorbar
        set(gca,'position',Pos_temp)
    end
    set(gca,'position',get(gca,'position')+[-0.03 -.03 0 .02])
    set(gca,'fontsize',FS)
    if l~=6
        set(gca,'TickDir','out','xticklabel',[],'yticklabel',[],'lineWidth',1.5)
    else
        set(gca,'TickDir','out','yticklabel',[],'lineWidth',1.5)
        yyaxis right
        set(gca,'YColor','k');
        set(gca,'ytick',.2:.2:1,'yticklabel',(.2:.2:1)*100)
        
        YL=ylabel('Frequency(Hz)');
        set(YL,'position',get(YL,'position')+[-0.07 0 0])
        xlabel('Time(sec)');
        yyaxis left
    end
    %
    if l==1
        title('Target(LM)');
        text(-1.45,155,'Evoked Laminar iPDCs (V1->LM)','fontsize',FS,'fontweight','bold');
        set(TL,'position',get(TL,'position')+[-.8 10 0])
    end
    
    text(-.45,50,['L' num2str(l)],'color','k','fontsize',FS+2,'fontweight','normal')
    vline(0,'--w')
    
end

colormap('jet')
export_fig(FIG,fullfile(SavePath,['STOK_iPDCLaminar_V1_LM' FileName]),'-pdf');
%%
for i = 2:6
ind2 = (i-1)*6+1:i*6;
FIG1 = dynet_connplot(PDC(ind2,ind1,:,Time>TW(1) & Time<TW(2)),Time(TimeInd),Freq,arrayfun(@(x) ['L' num2str(x)], 1:6,'uni',false),[],[],[],1);
colormap('jet')
set(FIG1,'unit','inch','position',[1 1 12 8],'color','w')
export_fig(FIG1,fullfile(SavePath,['STOK_V1' ROISN{i} FileName]),'-pdf');
end

%%
FIG1 = dynet_connplot(PDC_avg_diff,Time(TimeInd),Freq,cellfun(@(x) ROI_names.(x), STOK_avg.ROIs,'uni',false),[],[],[],1);
set(FIG1,'unit','inch','position',[0 0 25 15],'color','w');
if MODE==1
    export_fig(FIG1,fullfile(SavePath,['STOK_AllROIs_diff_sustained' FileName]),'-pdf'); close; 
else
    export_fig(FIG1,fullfile(SavePath,['STOK_AllROIs_diff_transient' FileName]),'-pdf'); close; 
end

FIG2 = dynet_connplot(PDC_avg,Time(TimeInd),Freq,cellfun(@(x) ROI_names.(x), STOK_avg.ROIs,'uni',false),[],[],[],1);
set(FIG2,'unit','inch','position',[0 0 25 15],'color','w');
if MODE==1
    export_fig(FIG2,fullfile(SavePath,['STOK_AllROIs_sustained' FileName]),'-pdf'); close; 
else
    export_fig(FIG2,fullfile(SavePath,['STOK_AllROIs_transient' FileName]),'-pdf'); close; 
end
%% Inter-Intra connectivity

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
% Convert to percentage
%PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);
% Remove diagonal
for i = 1:size(PDC,1)
    PDC(i,i,:,:)=0;
end

%PDC= PDC(1:42,1:42,:,:);

Cnd_var = {'PDCInter','PDCIntra','PDCInterIntra','PDCInterIntraPercent','PDCInterIntraPercent_prestimNorm'};

M = max(PDC(:))/5;
for Cnd =5:5
    FIG1 = figure(1);
    for roi = 1:NROIs
        figure(1);
        ind = [(roi-1)*6+1:roi*6];
        IND = setdiff(1:size(PDC,1),ind);
        switch Cnd
            case 1
                PDCInterIntra = squeeze(nanmean(PDC(IND,ind,:,:),1));%
            case 2
                PDCInterIntra = squeeze(nanmean(PDC(ind,ind,:,:),1));
            case 3
                PDCInterIntra = squeeze(nanmean(PDC(ind,ind,:,:),1) - nanmean(PDC(IND,ind,:,:),1));
            case 4
                PDCInterIntra = squeeze(nansum(PDC(ind,ind,:,:),1))./squeeze(nansum(PDC(ind,ind,:,:),1) + nansum(PDC(IND,ind,:,:),1));
                
            case 5
                PDCInterIntra = squeeze(nansum(PDC(ind,ind,:,:),1))./squeeze(nansum(PDC(ind,ind,:,:),1) + nansum(PDC(IND,ind,:,:),1));
                PDCInterIntra = PDCInterIntra - mean(PDCInterIntra(:,:,Time<0 & Time>-.3),3);
            case 6
                PDCInterIntra = squeeze(nanmean(PDC(ind,ind,:,:),1)-nanmean(PDC(IND,ind,:,:),1))./squeeze(nanmean(PDC(ind,ind,:,:),1) + nanmean(PDC(IND,ind,:,:),1));
                %PDCInterIntra = PDCInterIntra - mean(PDCInterIntra(:,:,Time<0 & Time>-.3),3);
                
        end
        subplot(NROIs,1,roi),
        imagesc(squeeze(mean(PDCInterIntra,1)));
        if Cnd==4
            %caxis([0 1])
        else
            caxis([-M M])
        end
        axis xy;
        
        title(STOK_avg.ROIs{roi});
        set(gca,'xtick',50:100:numel(Time),'xticklabel',round(Time(50:100:numel(Time)),2),...
            'ytick',5:20:100,'yticklabel',10:20:100);
        xlim([200 501]);
        vline(find(round(Time,2)==0,1),'k--');
        if roi==7
            xlabel('Time (msec)');
            ylabel('Frequency (Hz)');
            SUBPOS = get(gca,'position');
            CB = colorbar;
            set(gca,'position',SUBPOS);
            set(CB,'position',get(CB,'position')-[.04 0 0 0])
        end
        
        % This is the time and frequency distribution
    FIG2 = figure(2);
    plot(Time,squeeze(mean(mean(PDCInterIntra,1),2))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;

    FIG3 = figure(3);
    subplot(2,1,1);
    plot(Freq,squeeze(nanmean(nanmean(PDCInterIntra(:,:,Time>.045 & Time<.1),1),3))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;
    
    subplot(2,1,2);
    plot(Freq,squeeze(mean(mean(PDCInterIntra(:,:,Time>.3 & Time<1),1),3))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;
 
    %-----------------------------------------------------------------
    if Cnd ==4
        PDC_temp = PDCInterIntra(:,:,Time>-.3 & Time<0);
        disp([STOK_avg.ROIs{roi} ' - intra% ' num2str(nanmean(PDC_temp(:))*100)]);
    end
    %-----------------------------------------------------------------
    end
    
    
    figure(1);
    colormap('jet')
    set(FIG1,'unit','inch','position',[0,0,6,15],'color','w');
    export_fig(FIG1,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent' FileName]),'-pdf'); close; 
    
    figure(2);
    legend(STOK_avg.ROIs,'location','northeast')
    xlim([-.2 1])
    ylim([-8 4])
    vline(0,'k--');
    hline(0,'k--')
    xlabel('Time (S)');
    ylabel('%change');
    set(FIG2,'unit','inch','position',[0,0,6,4],'color','w');
    export_fig(FIG2,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_time' FileName]),'-pdf'); 
    
    
    figure(3);
    legend(STOK_avg.ROIs,'location','northeast')
    subplot(2,1,1);hline(0,'k--'); title('45-100 MS');xlim([5 100]);ylim([-8 2]);
    subplot(2,1,2);hline(0,'k--'); title('300-1000 MS');xlim([5 100]);ylim([-8 2]);
    ylabel('%change');
    xlabel('Frequency (Hz)');
    set(FIG3,'unit','inch','position',[0,0,5,10],'color','w');
    export_fig(FIG3,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_freq' FileName]),'-pdf'); 
    
    close all
end

    
%% Now based on those hierarchy orders, estimate FF and FB

% should it be percent change or not
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC_sub     = STOK_avg.PDC;

% Convert to percentage
PDC_sub     = (PDC_sub - mean(PDC_sub(:,:,:,Time>-.3 & Time<0),4))./mean(PDC_sub(:,:,:,Time>-.3 & Time<0),4);

M = max(abs(PDC_sub(:)));


Cnd_var = {'PDCFF','PDCFB','PDCFFFB','PDCFFFBpercent','PDCFFFBpercent_PrestimNorm'};

for Cnd = 1:2
    
    for roi = 1:NROIs
        FIG = figure(1);
        indO        = (roi-1)*6+1:roi*6; % indices for original ROI
        indFF       = setdiff(roi*6+1:size(PDC_sub,1),indO); % indices for FF ROI
        indFB       = setdiff(1:(roi-1)*6,indO); % indices for FF ROI
        PDCFF       = squeeze((nanmean(PDC_sub(indFF,indO,:,:),1)));
        PDCFB       = squeeze((nanmean(PDC_sub(indFB,indO,:,:),1)));
        
        PDCFF(isnan(PDCFF))=0;
        PDCFB(isnan(PDCFB))=0;
        switch Cnd
            case 1
                PDC_temp    = PDCFF;
            case 2
                PDC_temp    = PDCFB;
            case 3
                PDC_temp    = PDCFF-PDCFB;
            case 4
                PDCFF       = squeeze((nansum(PDC_sub(indFF,indO,:,:),1)));
                PDCFB       = squeeze((nansum(PDC_sub(indFB,indO,:,:),1)));
                PDCFF(isnan(PDCFF))=0;
                PDCFB(isnan(PDCFB))=0;
                PDC_temp    = PDCFF./(PDCFF+PDCFB);
                
            case 5
                PDCFF       = squeeze((nansum(PDC_sub(indFF,indO,:,:),1)));
                PDCFB       = squeeze((nansum(PDC_sub(indFB,indO,:,:),1)));
                PDCFF(isnan(PDCFF))=0;
                PDCFB(isnan(PDCFB))=0;
                PDC_temp    = PDCFF./(PDCFF+PDCFB);
                PDC_temp    = PDC_temp - mean(PDC_temp(:,:,Time<0 & Time>-.3),3);
                
        end
        
        subplot(7,1,roi)
        imagesc(squeeze(nanmean(PDC_temp(:,:,:),1)));
        hold on;
        
        if Cnd>3
            caxis([-1 1]/20)
        else
            plot(squeeze(nanmean(nanmean(PDC_temp(:,:,:),1),2))*30,'color','k','linewidth',1.5)
            %caxis([-M M])
            caxis([-2 2]/2)
        end
        axis xy
        title(STOK_avg.ROIs{roi});
        set(gca,'xtick',1:50:numel(Time),'xticklabel',round(Time(1:50:numel(Time)),2)*1000,...
            'ytick',5:20:100,'yticklabel',10:20:100);
        xlim([200 501]);
        vline(find(round(Time,2)==0,1),'k--');
        if roi==7
            xlabel('Time (msec)');
            ylabel('Frequency (Hz)');
            SUBPOS = get(gca,'position');
            CB = colorbar;
            set(gca,'position',SUBPOS);
            set(CB,'position',get(CB,'position')-[.04 0 0 0])
        end
        
            % This is the time and frequency distribution
    FIG2 = figure(2);
    plot(Time,squeeze(mean(mean(PDC_temp,1),2))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;

    FIG3 = figure(3);
    subplot(2,1,1);
    plot(Freq,squeeze(nanmean(nanmean(PDC_temp(:,:,Time>.045 & Time<.1),1),3))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;
    
    subplot(2,1,2);
    plot(Freq,squeeze(mean(mean(PDC_temp(:,:,Time>.3 & Time<1),1),3))*100,'color',Colors(roi,:),'linewidth',1.5);
    hold on;
 
    end
    figure(1);
    colormap('jet');
    set(FIG,'unit','inch','position',[0,0,6,15],'color','w')
    export_fig(FIG,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers' FileName]),'-pdf','-r200'); close; 
    
    figure(2);
    legend(STOK_avg.ROIs,'location','northeast')
    xlim([-.2 1])
    ylim([-20 110])
    vline(0,'k--');
    hline(0,'k--')
    xlabel('Time (S)');
    ylabel('%change');
    set(FIG2,'unit','inch','position',[0,0,6,4],'color','w');
    export_fig(FIG2,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_time' FileName]),'-pdf'); 
    
    
    figure(3);
    legend(STOK_avg.ROIs,'location','northeast')
    subplot(2,1,1);hline(0,'k--'); title('45-100 MS');xlim([5 100]);ylim([-20 120]);
    subplot(2,1,2);hline(0,'k--'); title('300-1000 MS');xlim([5 100]);ylim([-20 120]);
    ylabel('%change');
    xlabel('Frequency (Hz)');
    set(FIG3,'unit','inch','position',[0,0,5,10],'color','w');
    export_fig(FIG3,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_freq' FileName]),'-pdf'); 
    
    close all

end


%% Now make a graph plot for prestim, transient and steady-state response


Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;

PDC_perchange    = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);


TW =  [-.3 0; .05 .08; .3 1;.3 1];
FW = [1 100; 1 100; 30 55; 1 30];
subtitles = {'Pre-stimulus PDC','%change (50-100 msec)','%change (300-1000 msec,30-100 Hz)','%change (300-1000 msec,1-30 Hz)'};
clear PDCM;

for roi1 = 1:NROIs
    for roi2 = 1:NROIs
        ind1    = (roi1-1)*6+1:roi1*6; % indices for source ROI
        ind2    = (roi2-1)*6+1:roi2*6; % indices for target ROI
        
        for TC = 1:size(TW,1)
            if TC==1
                PDCM{TC}(roi2,roi1,:) = mean(mean(mean(PDC (ind2,ind1,FW(TC,1):FW(TC,2),Time>TW(TC,1) & Time<TW(TC,2)),1),2),4);
            else
                PDCM{TC}(roi2,roi1,:) = mean(mean(mean(PDC_perchange (ind2,ind1,FW(TC,1):FW(TC,2),Time>TW(TC,1) & Time<TW(TC,2)),1),2),4);

            end
        end
        
    end
end

%
close all;
FIG1 = figure(1);
%FIG2 = figure(2);
set(FIG1,'unit','inch','position',[2 10 12 5*size(TW,1)],'color','w');
%set(FIG2,'unit','inch','position',[2 0 15 5],'color','w');
load Hierarchyscores.mat;
HS = mean(H(Time<0 & Time>-.3,:));
LG = (cellfun(@(x) ROI_names.(x),STOK_avg.ROIs,'uni',false));

Mcc                         = ones((NROIs));
Mcc(find(triu(Mcc)))        = -1;
Mcc(1:length(Mcc)+1:end)    = 0; 

for Cnd = 1:size(TW,1)
    figure(1);subplot(size(TW,1),2,(Cnd-1)*2+1);
    Matrix = squeeze(mean(PDCM{Cnd},3));
    m1 = sum(sum(triu(Matrix.*Mcc*100)))/(NROIs*(NROIs-1)/2);
    m2 = sum(sum(tril(Matrix.*Mcc*100)))/(NROIs*(NROIs-1)/2);
    imagesc(Matrix.*Mcc*100);
    set(gca,'xtick',1:NROIs,'xticklabels',LG,'ytick',1:NROIs,'yticklabels',LG,'fontsize',14)
    %title(subtitles{Cnd})
    caxis([-1 1]*100);
    colormap(jmaColors('coolhot2'))
    CB = colorbar;
    CB.TickLabels = {'100','50','0','50','100'};
    xlabel('Source');ylabel('Target')
    title(subtitles{Cnd})
    
    figure(1);subplot(size(TW,1),2,(Cnd-1)*2+2);
    if Cnd ==1
        %Matrix = (Matrix-min(Matrix(:)))./(max(Matrix(:))-min(Matrix(:)));
        Matrix(1:length(Matrix)+1:end)=NaN;
        Matrix = Matrix - nanmean(Matrix(:));
        Matrix(Matrix<0)=0;
        Matrix = (Matrix-min(Matrix(:)))./(max(Matrix(:))-min(Matrix(:)));
        Matrix(1:length(Matrix)+1:end)=1;
       
    else
        %Matrix = Matrix./(max(Matrix(:)));
        Matrix(1:length(Matrix)+1:end)=NaN;
        Matrix = Matrix - nanmean(Matrix(:));
        Matrix(Matrix<0)=0;
        Matrix = (Matrix-min(Matrix(:)))./(max(Matrix(:))-min(Matrix(:)));
        Matrix(1:length(Matrix)+1:end)=1;
        Matrix = Matrix;
    end
    plot_graph(Matrix,LG,Colors,[],0.0,HS);
    %
    box off
    axis off;
    set(gca,'fontsize',14)
end

export_fig(FIG1,fullfile(SavePath,'PremResults',['Matrix_graph_PDCs' FileName]),'-pdf','-r300');  
%export_fig(FIG2,fullfile(SavePath,'PremResults',['Graph_PDCs']),'-pdf','-r200'); 
close all;

%% Hierarchical organization based on siegle paper

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
LG = (cellfun(@(x) ROI_names.(x),STOK_avg.ROIs,'uni',false));
% hierarchy scores based on input and output
%PDCdiff = (permute(PDC,[2 1 3 4])-PDC)./((permute(PDC,[2 1 3 4])+PDC)/2);
%PDC   = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

% remove inter connections
%PDC= PDC(1:42,1:42,:,:);

% average the scores

for roi1 = 1:NROIs
   ind1 = (roi1-1)*6+1:roi1*6;
   for roi2 = 1:NROIs
       if roi2 ~=roi1
           ind2 = (roi2-1)*6+1:roi2*6;
           HscoreMR(roi2,roi1,:,:) = nanmean(nanmean(PDC(ind2,ind1,:,:),1),2);
       end
   end
end
%
%HscoreMR = PDC;

for t = 50:numel(Time)
    M = (nanmean(HscoreMR(:,:,:,t),3));
    %M = M./nansum(M(:));
    [H(t,:),Hcc(t)] = HierarchyExtract(M);
end 

for f = 1:100
    f
    for t = 50:numel(Time)
        M = squeeze((HscoreMR(:,:,f,t)));
        %M = M./nansum(M(:));
        [Hf(f,t,:),Hccf(f,t)] = HierarchyExtract(M);
    end 
end


FIG = figure;
set(FIG,'unit','inch','position',[0,0,5,4.5],'color','w')
CMM = mean(H)'-mean(H);%squeeze(mean(CM(Time>-.3 & Time<0,:,:)));
imagesc(CMM)
colormap(jmaColors('coolhot2'))
set(gca,'xtick',1:NROIs,'xticklabel',LG,'ytick',1:NROIs,'yticklabel',LG,'fontsize',18,'fontweight','bold');
xlabel('Source');
ylabel('Target');
CB = colorbar;
CP = get(gca,'position');
caxis([-max(abs(CMM(:))) max(abs(CMM(:)))]);
set(CB,'position',get(CB,'position')+[0 .4 0.03 -.4]);
set(gca,'position',CP);
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_Prestim' FileName]),'-pdf'); close; 

%----------------------------Graph style-----------------------------------
HS = mean(H(Time<0 & Time>-.3,:));
Edges = mean(mean(HscoreMR(:,:,:,Time<0 & Time>-.3),4),3);
Edges = Edges ./max(Edges(:));
FIG = figure;
set(FIG,'unit','inch','position',[1,1,4,3.5],'color','w')
plot_graph(Edges.^2,LG,Colors,[],.0,HS);
axis on;
box off;
set(gca,'xtick',[],'ytick',round(HS,2));
ylabel('Functional Hierachy score');
set(gca,'fontsize',14)
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_Prestim_graph' FileName]),'-pdf'); close; 

%--------------------------------------------------------------------------
% OVER TIME
load ROInames;
Xticklabels = -.2:.200:1.000;
Xticks = arrayfun(@(x) find(round(Time,2)==x,1),Xticklabels);


HM = (H - mean(H(Time<0 & Time>-.3,:)))./abs(mean(H(Time<0 & Time>-.3,:)))*100;
FIG = figure;
set(FIG,'unit','inch','position',[0,0,6.5,6.5],'color','w');
subplot(2,1,1),
for roi = 1:NROIs
    plot(Time(100:end)*1000,squeeze(HM(100:end,roi)),'Color',Colors(roi,:),'linewidth',1.5);
    hold on
end
lg = legend(LG);
set(lg,'position',get(lg,'position')+[.05 0 0 0])
xlim([-.1 1]*1000)
title('ROI Hierarchy Scores % change');

vline(0,'k--')
set(gca,'fontsize',14)

subplot(2,1,2),%plot(Time(100:end)*1000,Hcc(100:end),'color','k','linewidth',1.5)
HccfPC = (Hccf - mean(Hccf(:,Time<0 & Time>-.3),2))./mean(Hccf(:,Time<0 & Time>-.3),2);
imagesc(HccfPC(:,:)*100); axis xy;
xlim([find(round(Time,2)==-.1,1) find(round(Time,2)==1,1)])
ylabel('Frequency (Hz)')
xlabel('Time (msec)')
title('Total Hierarchy Score % change')
%caxis([.14 .16])
caxis([-.05 .05]*100)
set(gca,'xtick',Xticks,'xticklabel',Xticklabels*1000)
SPP = get(gca,'position');
CB = colorbar;
set(gca,'position',SPP);
set(CB,'position',get(CB,'position')+[-.03 0 0 0])
%set(get(CB,'title'),'string','% change')
colormap('jet')
hold on
vline(find(round(Time,2)==0,1),'k--')
set(gca,'fontsize',14)


export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_PDC' FileName]),'-pdf'); close; 
ROIs = STOK_avg.ROIs;
save('Hierarchyscores','H','ROIs');


%--------------------------------------------------------------------------
% OVER TIME and FREQUENCY
load ROInames;
Xticklabels = -.2:.200:1.000;
Xticks = arrayfun(@(x) find(round(Time,2)==x,1),Xticklabels);

FIG = figure;
set(FIG,'unit','inch','position',[0,0,6.5,10.5],'color','w');

HMf = (Hf - mean(Hf(:,Time<0 & Time>-.3,:),2))./abs(mean(Hf(:,Time<0 & Time>-.3,:),2))*100;

for roi = 1:NROIs
    subplot(NROIs,1,roi),
    imagesc(squeeze(HMf(:,:,roi)));
    hold on
    caxis([-30 30])
    xlim([find(round(Time,2)==-.1,1) find(round(Time,2)==1,1)])
    axis xy;
    title(STOK_avg.ROIs{roi})
    set(gca,'xtick',Xticks,'xticklabel',Xticklabels*1000,'fontsize',14)
    if roi ==NROIs
        ylabel('Frequency (Hz)')
        xlabel('Time (msec)')
    end
end
colormap('jet')


export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_PDC_time_freq' FileName]),'-pdf'); close; 


