clear; clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
%load('/Users/elhamb/switchdrive/EB/AllenBrainData/STOK_Average_iPDC_ff.98_MOrd10_Thalamus.mat');
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_Average_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(STOK_avg.ROIs);
load ROInames;
NROIs = numel(STOK_avg.ROIs);

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
doperm = false;
for TC = 1:1
    FIG = figure;
    set(FIG,'unit','inch','position',[0,0,5,4.5],'color','w')
    HS = nanmean(nanmean(HscoreMR(:,:,:,Time>TW(TC,1) & Time<TW(TC,2)),3),4)';
    HierarchyExtract(HS);
    HSroi = nanmean(HS,2);
    HS(isnan(HS))=0;
    imagesc(HS)
    colormap(jmaColors('coolhot2'))
    set(gca,'xtick',1:7,'xticklabel',cellfun(@(x) ROI_names.(x),STOK_avg.ROIs([1:7]),'uni',false),'ytick',1:7,'yticklabel',cellfun(@(x) ROI_names.(x),STOK_avg.ROIs([1:7]),'uni',false),'fontsize',14,'fontweight','bold');
    xlabel('Source');
    ylabel('Target');
    CB = colorbar;
    CP = get(gca,'position');
    set(CB,'position',get(CB,'position')+[0 .4 0.03 -.4]);
    set(gca,'position',CP);
    %export_fig(FIG,fullfile(SavePath,'PremResults',['Hierarchy_Score_ROIs_' TName{TC}]),'-pdf'); close; 

    % check the best hierarchy
    if doperm
        Perm = perms(1:length(HS));
        HSp = arrayfun(@(x) sum(sum(tril(HS(Perm(x,:),Perm(x,:))))),1:size(Perm,1));
        [~,BestPerm] = max(HSp);
        BestPerm = Perm(BestPerm,:);
        disp ('Optimum order of ROIs for TName{TC}:' );
        disp(ROIName(BestPerm));
    end
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

