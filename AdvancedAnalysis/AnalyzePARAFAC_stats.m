addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%%
load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','PARAFAC');
%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
addpath(genpath('../../ARC'));
IDs = fieldnames(StokALL);

% bootstrap
ROIsPerID = cellfun(@(x) StokALL.(x).ROIs,IDs,'uni',false);
nboots = 100;
bootsize = 7;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);

    

for b = 1:nboots
    clear Stok_sample;
    % make the sub-sample stok structure
    for s = 1:bootsize
        Stok_sample.(['S' num2str(s)]) = StokALL.(BootIDs{b,s});
    end

    disp(num2str(b))
    [PDC,~,Time,Freq,ROIs] = ExtractAllRoiPDC(Stok_sample,ROIs);
    %
    clear PDCFF PDCFB indFF indFB
    ind = 1;
    indl = reshape(1:36,[6 6]);
    for roi1 = 1:numel(ROIs)
        for roi2 = roi1+1:numel(ROIs)
            PDCFF{ind} = PDC{roi2,roi1};
            indFF(ind) = indl(roi2,roi1);
            PDCFB{ind} = PDC{roi1,roi2};
            indFB(ind) = indl(roi1,roi2);
            ind = ind+1;
        end
    end

    PDCFFtemp = cellfun(@(x) mean(x,5),PDCFF,'uni',false);
    PDCFBtemp = cellfun(@(x) mean(x,5),PDCFB,'uni',false);
    PDCPulled = cat(2,PDCFFtemp,PDCFBtemp);
    PDCPulled = cat(5,PDCPulled{:});
    indTotal = [indFF indFB];

    NComp = 5;
    TW = Time>-0.5 & Time<1;
    [model{b}]=parafac(permute(PDCPulled(:,:,:,TW,:),[5 1 2 3 4]),NComp,[0 0 0 0 0],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    temp_time = Time(TW);
end
save([FileName 'PARAFAC_temp'],'model','NComp','temp_time','BootIDs','ROIs','IDs','Freq','indTotal')
load([FileName 'PARAFAC_temp.mat'])
%% cluster the components  500 x 500 similarity
% for example here I calculate correlation for the first mode : this does
% not work

% We have different modes of data ...
for i = 1:nboots
    for j = 1:nboots
        [maps(i,j,:),sscore(i,j),Corrs(i,j,:),~,CorrPvals(i,j,:,:)] = PARAFACmodelComp(model{i},model{j});
    end
end

% sscore(100,100)=0;
% sscore = max(sscore,sscore');
[~,RefInd] = max(mean(sscore)); % find the reference bootstrap
Corrsref = squeeze(Corrs(RefInd,:,:));

% reordering the models according to the reference component
for i = 1:nboots
    [~,sscore_reord(i),Corrs_reord(i,:),Model_reord{i}] = PARAFACmodelComp(model{RefInd},model{i});
end

%% if you want to check the components' stability, you should use the non-ordered ones: corrs and maps

% (1) Everything according to the reference bootstrap
for i = 1:nboots
    for j = 1:nboots
        map1 = squeeze(maps(RefInd,i,:));
        map2 = squeeze(maps(i,j,:));
        Corrs_reord2(i,j,:) = Corrs(i,j,map1);% or squeeze(maps(RefInd,j,:))
        CorrPvals_reord2(i,j,:,:) = CorrPvals(i,j,map1,:);% or squeeze(maps(RefInd,j,:))
    end
end

Comp_PerSig = squeeze(sum(reshape(CorrPvals_reord2,nboots*nboots,5,5)<.05)/(nboots.^2));

Comp_ord = [1 4 3 2 5];

FIG = figure;
set(FIG,'unit','inch','position',[0 0 2*NComp 4],'color','w')
bar(1:5,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,5)),.5);hold on;
errorbar(1:NComp,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,5)),var(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,5)),'.','color','k')
xlim([.5 5.5])
ylim([0.0 1.2])
xlabel('Component');
ylabel('Average inter boostrap Correlations')
for c = 1:NComp
    text(c,1,['Consistency=' num2str(round(mean(Comp_PerSig(Comp_ord(c),:)),2))],'HorizontalAlignment','center')
end
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'Consistency']),'-pdf','-r200')

%% mean variance of the models
Freq = 1:100;
model_reord =   cat(1,Model_reord{:});
%
Comp_ord = [1 4 3 2 5];
for id = RefInd%numel(IDs)
    FIG = figure;
    set(FIG,'unit','inch','position',[0 0 3.5*NComp 10],'color','w')

    model_temp_M = arrayfun(@(x) mean(cat(3,model_reord{:,x}),3),1:size(model_reord,2),'uni',false);
    model_temp_S = arrayfun(@(x) std(cat(3,model_reord{:,x}),[],3),1:size(model_reord,2),'uni',false);
    for c = 1:NComp
        cc = Comp_ord(c);
        %%PLOT SOURCE    
        MS = subplot(8,NComp,(c-1)+1); hold on;
        bar(model_temp_M{3}(:,cc)); colormap(MS,jmaColors('coolhot'))
        errorbar(1:6,model_temp_M{3}(:,cc),model_temp_S{3}(:,cc),'.','color','k')
        set(MS,'position',get(MS,'position')+[-.11+(.017*c) .05 .0 -.02],'xtick',1:6);
        ylim([0 1])
        
        xlabel('Source layers');
        title(['Consistency=' num2str(round(Comp_PerSig(cc,3),2))])
        

        % PLOT TARGET
        MT = subplot(8,NComp,(c-1)+1+NComp);hold on;
        bar(model_temp_M{2}(:,cc)); colormap(MS,jmaColors('coolhot'))
        errorbar(1:6,model_temp_M{2}(:,cc),model_temp_S{2}(:,cc),'.','color','k')
        set(MT,'position',get(MT,'position')+[-.11+(.017*c) .01 .0 -.02],'xtick',1:6);
        ylim([0 1])
        
        xlabel('Target layers')
        title(['Consistency=' num2str(round(Comp_PerSig(cc,2),2))])
        

        % PLOT TIME MODE
        M = model_temp_M{5}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
        TP = subplot(4,NComp,(c-1)+NComp+1); hold on;
        h = fill([temp_time flip(temp_time)],[M+model_temp_S{5}(:,cc); flip(M-model_temp_S{5}(:,cc))]','b','edgecolor','none');
        set(h,'facealpha',.5)
        plot(temp_time,M,'linewidth',1.5,'color','k')
        ylim([0 .1])
        xlim([-.2 1])
        hline(0,'k--')
        if c==1
            xlabel('Time(Sec)')
            ylabel('% change')
        end
        set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
        title(['Consistency=' num2str(round(Comp_PerSig(cc,5),2))])

        % PLOT FREQUENCY DISTRIBUTION
        TP = subplot(4,NComp,(c-1)+NComp*2+1); hold on;
        h = fill([Freq flip(Freq)],[model_temp_M{4}(:,cc)+model_temp_S{4}(:,cc); flip(model_temp_M{4}(:,cc)-model_temp_S{4}(:,cc))],'b','edgecolor','none');
        set(h,'facealpha',.5)
        plot(model_temp_M{4}(:,cc),'linewidth',1.5,'color','k')
        if c==1
            xlabel('Frequency(Hz)')
            %ylabel('% change')
        end
        set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
        title(['Consistency=' num2str(round(Comp_PerSig(cc,4),2))])
         ylim([0 .25])
         
        % PLOT CONNECTIVITY PATTERN
        CN = zeros(36,1);
        CN(indTotal) = model_temp_M{1}(:,cc);
        CN = reshape(CN,6,6);
        CNP = subplot(4,NComp,(c-1)+NComp*3+1);
        imagesc(CN);
        colormap(CNP,jmaColors('hotcortex'))
        set(CNP,'position',get(CNP,'position')+[-.1+(.017*c) -.02 0 0]);
        set(gca,'xtick',1:6,'xticklabel',ROISN,'ytick',1:6,'yticklabel',ROISN)
        colorbar
        CNT = CN-CN';
        % 1- best permutation order
        Perms = perms(1:length(CN));
        MDp = arrayfun(@(x) nansum(nansum(tril(CNT(Perms(x,:),Perms(x,:))))),1:size(Perms,1));
        [~,BestPerm] = max(MDp);
        xlabel(['Order = ' join(ROISN(Perms(BestPerm,:)),'-')],'fontweight','bold')
        title(['Consistency=' num2str(round(Comp_PerSig(cc,1),2))])
        
    end

end

export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_Bootstrap']),'-pdf','-r200')
%%

FIG = figure;
set(FIG,'unit','inch','position',[0 0 3.5*NComp 10],'color','w')


for c = 1:NComp
    cc = Comp_ord(c);
    %%PLOT SOURCE    
    MS = subplot(4,NComp*2,(c-1)*2+1);
    imagesc(model_temp_M{3}(:,cc),[-1 1]); colormap(MS,jmaColors('coolhot'))
    set(MS,'position',get(MS,'position')+[-.098+(.017*c) 0 -.015 0],'xticklabel',[]);
    title('Source')

    % PLOT TARGET
    MT = subplot(4,NComp*2,(c-1)*2+2);
    imagesc(model_temp_M{2}(:,cc),[-1 1]); colormap(MT,jmaColors('coolhot'))
    set(MT,'position',get(MT,'position')+[-.115+(.017*c) 0 -.015 0],'xticklabel',[],'yticklabel',[]);
    title('Target')

    % PLOT TIME MODE
    M = model_temp_M{5}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
    TP = subplot(4,NComp,(c-1)+NComp+1);
    plot(temp_time,M,'linewidth',1.5,'color','k')
    xlim([-.2 1])
    hline(0,'k--')
    if c==1
        xlabel('Time(Sec)')
        ylabel('% change')
    end
    set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 0]);
    ylim([0 .1])

    % PLOT FREQUENCY DISTRIBUTION
    TP = subplot(4,NComp,(c-1)+NComp*2+1);
    plot(model_temp_M{4}(:,cc),'linewidth',1.5,'color','k')
    if c==1
        xlabel('Frequency(Hz)')
        %ylabel('% change')
    end
    set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 0]);
    ylim([0 .24])

    % PLOT CONNECTIVITY PATTERN
    CN = zeros(36,1);
    CN(indTotal) = model_temp_M{1}(:,cc);
    CN = reshape(CN,6,6);
    CNP = subplot(4,NComp,(c-1)+NComp*3+1);
    imagesc(CN);
    colormap(CNP,jmaColors('hotcortex'))
    set(CNP,'position',get(CNP,'position')+[-.1+(.017*c) 0 0 0]);
    set(gca,'xtick',1:6,'xticklabel',ROISN,'ytick',1:6,'yticklabel',ROISN)
    colorbar
    CNT = CN-CN';
    % 1- best permutation order
    Perms = perms(1:length(CN));
    MDp = arrayfun(@(x) nansum(nansum(tril(CNT(Perms(x,:),Perms(x,:))))),1:size(Perms,1));
    [~,BestPerm] = max(MDp);
    xlabel(['Order = ' join(ROISN(Perms(BestPerm,:)),'-')],'fontweight','bold')
end

export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_Bootstrap_mean']),'-pdf','-r200')
