
clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
load('/Users/elhamb/switchdrive/EB/AllenBrainData/STOK_Average_iPDC_ff.98_MOrd10_Thalamus.mat');
%load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_Average_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs([1:7]));

%% plot each ROI separately
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

for roi = 1:8
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

PDC= PDC(1:42,1:42,:,:);

Cnd_var = {'PDCInter','PDCIntra','PDCInterIntra','PDCInterIntraPercent','PDCInterIntraPercent_prestimNorm'};

M = max(PDC(:))/5;
for Cnd =5:5
    FIG1 = figure(1);
    for roi = 1:7
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
        subplot(7,1,roi),
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

PDC_sub = PDC(1:42,1:42,:,:);
% Convert to percentage
PDC_sub     = (PDC_sub - mean(PDC_sub(:,:,:,Time>-.3 & Time<0),4))./mean(PDC_sub(:,:,:,Time>-.3 & Time<0),4);

M = max(abs(PDC_sub(:)));


Cnd_var = {'PDCFF','PDCFB','PDCFFFB','PDCFFFBpercent','PDCFFFBpercent_PrestimNorm'};

for Cnd = 1:2
    
    for roi = 1:7
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

PDC         = PDC(1:42,1:42,:,:);
PDC_perchange    = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);


TW = [-.3 0; .05 .09; .1 1];
TName = {'Prestim','Transient','SteadyState'};

for roi1 = 1:7
    for roi2 = 1:7
        ind1    = (roi1-1)*6+1:roi1*6; % indices for source ROI
        ind2    = (roi2-1)*6+1:roi2*6; % indices for target ROI
        
        for TC = 1:3
            if TC==1
                PDCM(roi2,roi1,:,TC) = mean(mean(mean(PDC (ind2,ind1,:,Time>TW(TC,1) & Time<TW(TC,2)),1),2),4);
            else
                PDCM(roi2,roi1,:,TC) = mean(mean(mean(PDC_perchange (ind2,ind1,:,Time>TW(TC,1) & Time<TW(TC,2)),1),2),4);

            end
        end
        
    end
end

%
close all;
FIG1 = figure(1);
FIG2 = figure(2);
set(FIG1,'unit','inch','position',[2 10 20 5],'color','w');
set(FIG2,'unit','inch','position',[2 0 15 5],'color','w');
subtitles = {'Prestim PDC','%change PDC (transient)','%change PDC (Steady)'};

for Cnd = 1:3
    figure(1);subplot(1,3,Cnd);
    Matrix = squeeze(mean(PDCM(:,:,:,Cnd),3));
    imagesc(Matrix);
    set(gca,'xtick',1:7,'xticklabels',STOK_avg.ROIs(1:7),'ytick',1:7,'yticklabels',STOK_avg.ROIs(1:7))
    title(subtitles{Cnd})
    colorbar
    figure(2);subplot(1,3,Cnd);
    if Cnd ==1
        Matrix = (Matrix-min(Matrix(:)))./(max(Matrix(:))-min(Matrix(:)));
    end
    plot_graph(Matrix,STOK_avg.ROIs(1:7),Colors,'up',0.6);
    title(subtitles{Cnd})
end

export_fig(FIG1,fullfile(SavePath,'PremResults',['Matrix_PDCs']),'-pdf','-r200');  
export_fig(FIG2,fullfile(SavePath,'PremResults',['Graph_PDCs']),'-pdf','-r200'); 
close all;

%% Hierarchical organization based on input and output

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
% hierarchy scores based on input and output
%PDCdiff = (permute(PDC,[2 1 3 4])-PDC)./((permute(PDC,[2 1 3 4])+PDC)/2);
%PDC   = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

% remove inter connections
PDC= PDC(1:42,1:42,:,:);
% average the scores

for roi1 = 1:7
   ind1 = (roi1-1)*6+1:roi1*6;
   for roi2 = 1:7
       if roi2 ~=roi1
           ind2 = (roi2-1)*6+1:roi2*6;
           HscoreMR(roi2,roi1,:,:) = nanmean(nanmean(PDC(ind2,ind1,:,:),1),2);
       end
   end
end
%
for t = 100:numel(Time)
    M = (nanmean(HscoreMR(:,:,:,t),3));
    %M = M./nansum(M(:));
    [H(t,:),Hcc(t)] = HierarchyExtract(M);
end 

%

FIG = figure;
subplot(2,1,1),
for roi = 1:7
    plot(Time(100:end),squeeze(H(100:end,roi)),'Color',Colors(roi,:),'linewidth',1.5);
    hold on
end
legend(STOK_avg.ROIs(1:7))
xlim([-.3 1])
ylabel('Hierarchy score')

subplot(2,1,2),plot(Time(100:end),Hcc(100:end),'color','k','linewidth',1.5)
xlim([-.3 1])
ylabel('Total Hierarchy')
xlabel('Time (msec)')

set(FIG,'unit','inch','position',[0,0,8,10],'color','w');
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_percentChange']),'-pdf'); close; 


%% Hierarchical organization based on input and output

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
% Convert to percentage
%PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);

% hierarchy scores based on input and output
PDCdiff = (permute(PDC,[2 1 3 4])-PDC)./((permute(PDC,[2 1 3 4])+PDC)/2);
% remove inter connections

PDCdiff= PDCdiff(1:42,1:42,:,:);
% average the scores

for roi = 1:7
    if roi<9
        ind1 = (roi-1)*6+1:roi*6;
    else
        switch roi
            case 9
                ind1 = 49:51;
            case 10
                ind1 = 52:54;
        end
    end
   ind2 = setdiff(1:size(PDCdiff,1),(roi-1)*6+1:roi*6);
   HscoreM(roi,:,:) = nanmean(nanmean(PDCdiff(ind2,ind1,:,:),1),2);
   for roi2 = 1:7
       ind2 = (roi2-1)*6+1:roi2*6;
       HscoreMR(roi2,roi,:,:) = nanmean(nanmean(PDCdiff(ind2,ind1,:,:),1),2);

   end
end


Htotal = arrayfun(@(x) nanmean(nanmean(tril(squeeze(nanmean(PDCdiff(:,:,:,x),3))),1),2),1:size(PDCdiff,4));
FIG = figure;
plot(Time,Htotal,'color','k','linewidth',1.5);
xlim([-.3 1]);
set(FIG,'unit','inch','position',[0,0,8,5],'color','w');
export_fig(FIG,fullfile(SavePath,'PremResults',['Total_hierarchy_Score']),'-pdf'); close; 


%--------------------------------------------------------------------------
% ROIs relative hierarchy averaged
ROIName = STOK_avg.ROIs(1:7);
TW = [-.3 0; .05 .09; .1 1];
TName = {'Prestim','Transient','SteadyState'};
for TC = 1:3
    FIG = figure;
    HS = nanmean(nanmean(HscoreMR(:,:,:,Time>TW(TC,1) & Time<TW(TC,2)),3),4)';
    HierarchyExtract(HS);
    HSroi = nanmean(HS,2);
    HS(isnan(HS))=0;
    imagesc(HS)
    colormap(jmaColors('coolhot'))
    set(gca,'xtick',1:7,'xticklabel',STOK_avg.ROIs([1:7]),'ytick',1:7,'yticklabel',STOK_avg.ROIs([1:7]));
    colorbar
    set(FIG,'unit','inch','position',[0,0,5,4],'color','w')
    export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_' TName{TC}]),'-pdf'); close; 

    % check the best hierarchy
    Perm = perms(1:length(HS));
    HSp = arrayfun(@(x) sum(sum(tril(HS(Perm(x,:),Perm(x,:))))),1:size(Perm,1));
    [~,BestPerm] = max(HSp);
    BestPerm = Perm(BestPerm,:);
    disp ('Optimum order of ROIs for TName{TC}:' );
    disp(ROIName(BestPerm));
end


%--------------------------------------------------------------------------
% Hierarchy over time
FIG = figure;
COBJ = LFPF.RColors();
%Colors = COBJ.MatrixColors(STOK_avg.ROIs([1:7]));
for roi = 1:7
    plot(Time*1000,squeeze(mean(HscoreM(roi,:,:),2))','linewidth',1.5,'color',Colors(roi,:));
    hold on;
end
xlim([-500 1000])
legend(STOK_avg.ROIs([1:7]))
xlabel('Time (msec)')
ylabel('Hierarchy score')
set(FIG,'unit','inch','position',[0,0,10,6],'color','w')
export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_percentChange']),'-pdf'); close; 


%% Output/input of layers

Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
% Convert to percentage
PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);
% Remove diagonal
for i = 1:size(PDC,1)
    PDC(i,i,:,:)=0;
end
%PDC = permute(PDC,[2 1 3 4]);
%--------------------------------------------------------------------------
%for each area separetly
M = max(PDC(:));
for roi1 = 1:7
    FIG = figure;
    for roi2 = 1:7
        ind1 = [(roi1-1)*6+1:roi1*6];
        ind2 = [(roi2-1)*6+1:roi2*6];
        
        PDC_temp = squeeze(nanmean(PDC(ind2,ind1,:,:),1));%,
        
        for l = 1:6
            subplot(6,8,(l-1)*8+roi2),imagesc(squeeze(PDC_temp(l,:,:)));
            set(gca,'xtick',50:100:numel(Time),'xticklabel',round(Time(50:100:numel(Time)),2),...
            'ytick',5:20:100,'yticklabel',10:20:100);
           % xtickangle(90);
            xlim([200 501]);
            vline(find(round(Time,2)==0,1),'k--');
            axis xy
            caxis([-M M]/10)
            if l ==1
                title(STOK_avg.ROIs{roi2})
            end
        end
        if l==6 && roi2==1
            xlabel('Time (msec)');
            ylabel('Frequency (Hz)');
        end
    end
    
    colormap('jet')
    axes('position',[.5 .98 .1 .05]);
    text(0,0,STOK_avg.ROIs{roi1},'fontsize',12);axis off
    set(FIG,'unit','inch','position',[0,0,13,8],'color','w')
    export_fig(FIG,fullfile(SavePath,'PremResults',['STOK_individual_' STOK_avg.ROIs{roi1} '_iPDC_OutPut']),'-pdf','-r200'); close; 
end

% --------------------------------------------------------------------------
% Averaged over layers
Time    = STOK_avg.Time;
Freq    = STOK_avg.Freq;
PDC     = STOK_avg.PDC;
% Convert to percentage
%PDC     = (PDC - mean(PDC(:,:,:,Time>-.3 & Time<0),4))./mean(PDC(:,:,:,Time>-.3 & Time<0),4);
% Remove diagonal
for i = 1:size(PDC,1)
    PDC(i,i,:,:)=0;
end
M = max(PDC(:));

FIG = figure;
for roi1 = 1:7
    for roi2 = 1:7
        ind1 = (roi1-1)*6+1:roi1*6;
        ind2 = (roi2-1)*6+1:roi2*6;
        
        PDC_temp = squeeze(nanmean(nanmean(PDC(ind2,ind1,:,:),1),2));%
        PDCt(roi2,roi1,:,:)= PDC_temp;
        subplot(7,7,(roi2-1)*7+roi1),imagesc(PDC_temp);
               
        set(gca,'xtick',50:100:numel(Time),'xticklabel',round(Time(50:100:numel(Time)),2),...
        'ytick',5:20:100,'yticklabel',10:20:100);
        xlim([200 501]);
        vline(find(round(Time,2)==0,1),'k--');
        axis xy
        caxis([-M M]/2)
        if roi1 ==1
            ylabel(STOK_avg.ROIs{roi2},'fontweight','bold')
        end
        
        if roi2 ==1
            title(STOK_avg.ROIs{roi1})
        end
       
    end
end
colormap('jet');
set(FIG,'unit','inch','position',[0,0,20,10],'color','w')
export_fig(FIG,fullfile(SavePath,'PremResults',['STOK_individual_AllROIs_iPDC_nonnorm']),'-pdf','-r200'); close; 
