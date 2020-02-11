
clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
%load('/Users/elhamb/switchdrive/EB/AllenBrainData/STOK_Average_iPDC_ff.98_MOrd10_Thalamus.mat');
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_Average_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs);
load ROInames;
NROIs = numel(STOK_avg.ROIs);
%% plot each ROI separately
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

for roi = 1:NROIs
    ind = [(roi-1)*6+1:roi*6];
    FIG = dynet_connplot(PDC(ind,ind,:,Time>-.2 & Time<1),Time(Time>-.2 & Time<1),Freq,cellfun(@(x) strrep(x,'_','-'),STOK_avg.labels(ind),'uni',false),[],[],[],[],STOK_avg.ROIs{roi});
    set(FIG,'unit','inch','position',[0,0,20,12],'color','w')
    export_fig(FIG,fullfile(SavePath,['STOK_individual_sustained_' STOK_avg.ROIs{roi} '_iPDC']),'-pdf'); close; 
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
    export_fig(FIG1,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent']),'-pdf'); close; 
    
    figure(2);
    legend(STOK_avg.ROIs(1:7),'location','northeast')
    xlim([-.2 1])
    ylim([-20 110])
    vline(0,'k--');
    hline(0,'k--')
    xlabel('Time (S)');
    ylabel('%change');
    set(FIG2,'unit','inch','position',[0,0,6,4],'color','w');
    export_fig(FIG2,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_time']),'-pdf'); 
    
    
    figure(3);
    legend(STOK_avg.ROIs(1:7),'location','northeast')
    subplot(2,1,1);hline(0,'k--'); title('45-100 MS');xlim([5 100]);ylim([-20 110]);
    subplot(2,1,2);hline(0,'k--'); title('300-1000 MS');xlim([5 100]);ylim([-20 110]);
    ylabel('%change');
    xlabel('Frequency (Hz)');
    set(FIG3,'unit','inch','position',[0,0,5,10],'color','w');
    export_fig(FIG3,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_freq']),'-pdf'); 
    
    close all
end

    
%% Now based on those hierarchy orders, estimate FF and FB

% should it be percent change or not
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;

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
    export_fig(FIG,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers']),'-pdf','-r200'); close; 
    
    figure(2);
    legend(STOK_avg.ROIs(1:7),'location','northeast')
    xlim([-.2 1])
    ylim([-20 110])
    vline(0,'k--');
    hline(0,'k--')
    xlabel('Time (S)');
    ylabel('%change');
    set(FIG2,'unit','inch','position',[0,0,6,4],'color','w');
    export_fig(FIG2,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_time']),'-pdf'); 
    
    
    figure(3);
    legend(STOK_avg.ROIs(1:7),'location','northeast')
    subplot(2,1,1);hline(0,'k--'); title('45-100 MS');xlim([5 100]);ylim([-20 120]);
    subplot(2,1,2);hline(0,'k--'); title('300-1000 MS');xlim([5 100]);ylim([-20 120]);
    ylabel('%change');
    xlabel('Frequency (Hz)');
    set(FIG3,'unit','inch','position',[0,0,5,10],'color','w');
    export_fig(FIG3,fullfile(SavePath,'PremResults',[Cnd_var{Cnd} '_AveragedLayers_percent_freq']),'-pdf'); 
    
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

export_fig(FIG1,fullfile(SavePath,'PremResults',['Matrix_graph_PDCs']),'-pdf','-r300');  
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

for t = 100:numel(Time)
    M = (nanmean(HscoreMR(:,:,:,t),3));
    %M = M./nansum(M(:));
    [H(t,:),Hcc(t)] = HierarchyExtract(M);
end 

for f = 1:100
    f
    for t = 100:numel(Time)
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
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_Prestim']),'-pdf'); close; 

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
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_Prestim_graph']),'-pdf'); close; 

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


export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_PDC']),'-pdf'); close; 
ROIs = STOK_avg.ROIs;
save('Hierarchyscores','H','ROIs');


