addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nway331'));
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/dynet_toolbox-master'));
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
bootsize = 8;
BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);

NComp = 3;

for roi = 2:numel(ROISN)

for b = 1:nboots
    clear Stok_sample;
    % make the sub-sample stok structure
    for s = 1:bootsize
        Stok_sample.(['S' num2str(s)]) = StokALL.(BootIDs{b,s});
    end

    disp(num2str(b))
    [PDC,~,Time,Freq,ROIs] = ExtractAllRoiPDC(Stok_sample,ROIs);
    %
    clear PDCFF PDCFB PDCLam indFF indFB 
    ind = 1;
    indl = reshape(1:36,[6 6]);
    for roi1 = 1:numel(ROIs)
        PDCLam{roi1} = PDC{roi1,roi1};
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
    
    PDCLam = cellfun(@(x) mean(x,5),PDCLam,'uni',false);
    PDCLam = cat(5,PDCLam{:});
    for x = 1:6
        PDCLam(x,x,:,:,:)=0;
    end
    indTotal = [indFF indFB];
    
    TW = Time>-0.5 & Time<1;
    %[ssX,Corco,It] = pftest(4,reshape(permute(PDCLam(:,:,:,TW,:),[4 1 2 5 3]),[375 36*6 100]),3,[2 2 2 0 0]);
    %[model_laminar{b},it(b),err(b),corcondia(b)]=parafac(squeeze(permute(PDCLam(:,:,:,TW,:),[4 1 2 3 5])),NComp,[0 1 0 0 NaN],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    %[model_laminar{b},it(b),err(b),corcondia(b)]=parafac(reshape(permute(PDCLam(:,:,:,TW,:),[4 1 2 5 3]),[375 36*6 100]),NComp,[0 1 0 0 NaN],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    [model_laminar{b},it(b),err(b),corcondia(b)]=parafac(reshape(permute(PDCLam(:,:,:,TW,roi),[4 1 2 5 3]),[375 36 100]),NComp,[0 1 0 0 NaN],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    
    temp_time = Time(TW);
end

save([FileName 'PARAFAC_covtemp_' num2str(NComp) '_laminar2_' ROISN{roi}],'model_laminar','NComp','temp_time','BootIDs','nboots','ROIs','IDs','Freq','indTotal','it','err','corcondia')
%load([FileName 'PARAFAC_covtemp_' num2str(NComp) '_laminar2_' ROISN{roi}])

%% cluster the components  500 x 500 similarity
% for example here I calculate correlation for the first mode : this does
% not work

% We have different modes of data ...
for i = 1:nboots
    for j = 1:nboots
        [maps(i,j,:),sscore(i,j),Corrs(i,j,:),~,CorrPvals(i,j,:,:),Corrvals(i,j,:,:)] = PARAFACmodelComp(model_laminar{i},model_laminar{j});
    end
end

[~,RefInd] = max(mean(sscore)); % find the reference bootstrap
Corrsref = squeeze(Corrs(RefInd,:,:));

% reordering the models according to the reference component
for i = 1:nboots
    [~,sscore_reord(i),Corrs_reord(i,:),Model_reord{i}] = PARAFACmodelComp(model_laminar{RefInd},model_laminar{i});
end

%% how much of variance are explained by each component
%[1 4 3 2 5];

clear Model_vars;
for y = 1:NComp
    Model_vars{y} = cellfun(@(x) norm(x{1}(:,y)),Model_reord);
end



Model_vars=cat(1,Model_vars{:})';
Model_var_percent = Model_vars./sum(Model_vars,2);

[~,Comp_ord] = sort(mean(Model_var_percent),'descend');%[1 3 5 2 4];

FIG = figure;
set(FIG,'unit','inch','position',[0 0 2*NComp 4],'color','w')
bar(1:NComp,mean(Model_var_percent(:,Comp_ord)),.5);hold on;
errorbar(1:NComp,mean(Model_var_percent(:,Comp_ord)),var(Model_var_percent(:,Comp_ord)),'.','color','k')
xlabel('Component');
ylabel('Variance explained')
xlim([.5 NComp+.5])
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_VarianceExplained_laminar_' ROISN{roi}]),'-pdf','-r200')

%% if you want to check the components' Consistency, you should use the non-ordered ones: corrs and maps

% (1) Everything according to the reference bootstrap
for i = 1:nboots
    for j = 1:nboots
        map1 = squeeze(maps(RefInd,i,:));
        map2 = squeeze(maps(i,j,:));
        Corrs_reord2(i,j,:) = Corrs(i,j,map1);% or squeeze(maps(RefInd,j,:))
        CorrPvals_reord2(i,j,:,:) = CorrPvals(i,j,map1,:);% or squeeze(maps(RefInd,j,:))
    end
end

Comp_PerSig = squeeze(sum(reshape(CorrPvals_reord2,nboots*nboots,NComp,3)<.05)/(nboots.^2));


FIG = figure;
set(FIG,'unit','inch','position',[0 0 2*NComp 4],'color','w')
bar(1:NComp,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),.5);hold on;
errorbar(1:NComp,mean(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),var(reshape(Corrs_reord2(:,:,Comp_ord),nboots*nboots,NComp)),'.','color','k')
xlim([.5 NComp+.5])
ylim([0.0 1.2])
xlabel('Component');
ylabel('Average inter boostrap Correlations')
for c = 1:NComp
    text(c,1,['Consistency=' num2str(round(mean(Comp_PerSig(Comp_ord(c),:)),2))],'HorizontalAlignment','center')
end
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'Consistency_laminar_' ROISN{roi}]),'-pdf','-r200')


%% mean variance of the models
Freq = 1:100;
model_reord =   cat(1,Model_reord{:});
%
for id = RefInd%numel(IDs)
    FIG = figure;
    set(FIG,'unit','inch','position',[0 0 3.5*NComp 10],'color','w')

    model_temp_M = arrayfun(@(x) mean(cat(3,model_reord{:,x}),3),1:size(model_reord,2),'uni',false);
    %model_temp_M= Breakload2cells(model_temp_M,2, [6 6 6]);
    model_temp_S = arrayfun(@(x) std(cat(3,model_reord{:,x}),[],3),1:size(model_reord,2),'uni',false);
    %model_temp_S= Breakload2cells(model_temp_S,2, [6 6 6]);
    for c = 1:NComp
        cc = Comp_ord(c);
        %%PLOT SOURCE    
        source_dim = 2;
        MS = subplot(4,NComp,(c-1)+1); 
        imagesc(reshape(model_temp_M{source_dim}(:,cc),6,6))
        set(MS,'position',get(MS,'position')+[-.11+(.017*c) .0 .0 .04],'xtick',1:6,'ytick',1:6);
        %ylim([0 1]/2)
        
        xlabel('Source layers');
        ylabel('Target Layer')
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,3),2))])
        

        % PLOT TIME MODE
        temp_dim = 1;
        M = model_temp_M{temp_dim}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
        TP = subplot(4,NComp,(c-1)+NComp+1); hold on;
        h = fill([temp_time flip(temp_time)],[M+model_temp_S{temp_dim}(:,cc); flip(M-model_temp_S{temp_dim}(:,cc))]','b','edgecolor','none');
        set(h,'facealpha',.5)
        plot(temp_time,M,'linewidth',1.5,'color','k')
        %ylim([0 .1])
        xlim([-.2 1])
        hline(0,'k--')
        if c==1
            xlabel('Time(Sec)')
            ylabel('% change')
        end
        set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,5),2))])

        % PLOT FREQUENCY DISTRIBUTION
        freq_dim = 3;
        TP = subplot(4,NComp,(c-1)+NComp*2+1); hold on;
        h = fill([Freq flip(Freq)],[model_temp_M{freq_dim}(:,cc)+model_temp_S{freq_dim}(:,cc); flip(model_temp_M{freq_dim}(:,cc)-model_temp_S{freq_dim}(:,cc))],'b','edgecolor','none');
        set(h,'facealpha',.5)
        plot(model_temp_M{freq_dim}(:,cc),'linewidth',1.5,'color','k')
        if c==1
            xlabel('Frequency(Hz)')
            %ylabel('% change')
        end
        set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,4),2))])
         ylim([0 .25])
         
        
    end

end

export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_Bootstrap_laminar_' ROISN{roi}]),'-pdf','-r200')


end

%% mean variance of the models
% Freq = 1:100;
% model_reord =   cat(1,Model_reord{:});
% %
% for id = RefInd%numel(IDs)
%     FIG = figure;
%     set(FIG,'unit','inch','position',[0 0 3.5*NComp 10],'color','w')
% 
%     model_temp_M = arrayfun(@(x) mean(cat(3,model_reord{:,x}),3),1:size(model_reord,2),'uni',false);
%     %model_temp_M= Breakload2cells(model_temp_M,2, [6 6 6]);
%     model_temp_S = arrayfun(@(x) std(cat(3,model_reord{:,x}),[],3),1:size(model_reord,2),'uni',false);
%     %model_temp_S= Breakload2cells(model_temp_S,2, [6 6 6]);
%     for c = 1:NComp
%         cc = Comp_ord(c);
%         %%PLOT SOURCE    
%         source_dim = 3;
%         MS = subplot(8,NComp,(c-1)+1); hold on;
%         bar(model_temp_M{source_dim}(:,cc)); colormap(MS,jmaColors('coolhot'))
%         errorbar(1:6,model_temp_M{source_dim}(:,cc),model_temp_S{source_dim}(:,cc),'.','color','k')
%         set(MS,'position',get(MS,'position')+[-.11+(.017*c) .05 .0 -.02],'xtick',1:6);
%         ylim([0 1]/2)
%         
%         xlabel('Source layers');
%         %title(['Consistency=' num2str(round(Comp_PerSig(cc,3),2))])
%         
% 
%         % PLOT TARGET
%         target_dim =2;
%         MT = subplot(8,NComp,(c-1)+1+NComp);hold on;
%         bar(model_temp_M{target_dim}(:,cc)); colormap(MS,jmaColors('coolhot'))
%         errorbar(1:6,model_temp_M{target_dim}(:,cc),model_temp_S{target_dim}(:,cc),'.','color','k')
%         set(MT,'position',get(MT,'position')+[-.11+(.017*c) .01 .0 -.02],'xtick',1:6);
%         ylim([0 1]/2)
%         
%         xlabel('Target layers')
%         %title(['Consistency=' num2str(round(Comp_PerSig(cc,2),2))])
%         
% 
%         % PLOT TIME MODE
%         temp_dim = 1;
%         M = model_temp_M{temp_dim}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
%         TP = subplot(4,NComp,(c-1)+NComp+1); hold on;
%         h = fill([temp_time flip(temp_time)],[M+model_temp_S{temp_dim}(:,cc); flip(M-model_temp_S{temp_dim}(:,cc))]','b','edgecolor','none');
%         set(h,'facealpha',.5)
%         plot(temp_time,M,'linewidth',1.5,'color','k')
%         %ylim([0 .1])
%         xlim([-.2 1])
%         hline(0,'k--')
%         if c==1
%             xlabel('Time(Sec)')
%             ylabel('% change')
%         end
%         set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
%         %title(['Consistency=' num2str(round(Comp_PerSig(cc,5),2))])
% 
%         % PLOT FREQUENCY DISTRIBUTION
%         freq_dim = 4;
%         TP = subplot(4,NComp,(c-1)+NComp*2+1); hold on;
%         h = fill([Freq flip(Freq)],[model_temp_M{freq_dim}(:,cc)+model_temp_S{freq_dim}(:,cc); flip(model_temp_M{freq_dim}(:,cc)-model_temp_S{freq_dim}(:,cc))],'b','edgecolor','none');
%         set(h,'facealpha',.5)
%         plot(model_temp_M{freq_dim}(:,cc),'linewidth',1.5,'color','k')
%         if c==1
%             xlabel('Frequency(Hz)')
%             %ylabel('% change')
%         end
%         set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 -.04]);
%         %title(['Consistency=' num2str(round(Comp_PerSig(cc,4),2))])
%          ylim([0 .25])
%          
%          % PLOT ROIs
% %         ROI_dim =5;
% %         MT = subplot(4,NComp,(c-1)+NComp*3+1);hold on;
% %         bar(model_temp_M{ROI_dim}(:,cc)); colormap(MS,jmaColors('coolhot'))
% %         errorbar(1:6,model_temp_M{ROI_dim}(:,cc),model_temp_S{ROI_dim}(:,cc),'.','color','k')
% %         set(MT,'position',get(MT,'position')+[-.11+(.017*c) .01 .0 -.02],'xtick',1:6,'xticklabel',cellfun(@(x) ROI_names.(x),ROIs,'uni',false));
% %         ylim([0 .5])
% %         
% %         xlabel('ROIs')
%         
%         %title(['Consistency=' num2str(round(Comp_PerSig(cc,2),2))])
%         
%     end
% 
% end
% 
% export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_Bootstrap_laminar']),'-pdf','-r200')
% 
% %% Taking the data back to data space to compare temporal dynamics in four frequency bands
% Compcolors = [1 0 0; .5 0 0; .7 .4 .3;.3 .4 .7; 0 0 1];
% 
% Fbins =[0 5;5 10; 10 20; 20 100];
% FIG = figure;
% set(FIG,'unit','inch','position',[0 0 10 5],'color','w');
% hold on;
% 
% cc= Comp_ord;%[1 3 2 5 4];
% 
%     
% MT = cellfun(@(x) (x{1}),Model_reord,'uni',false);
% MT = permute(cat(3,MT{:}),[3 2 1]); 
% 
% clear M S;
% 
% for l = NComp:-1:1
%     CIu = squeeze(quantile(MT(:,cc(l),:),.95,1));
%     CId = squeeze(quantile(MT(:,cc(l),:),.05,1));
%     MTT = (squeeze(MT(:,cc(l),:))>repmat(CId',[nboots,1])) & (squeeze(MT(:,cc(l),:))<repmat(CIu',[nboots,1]));
%     MTS = squeeze(MT(:,cc(l),:));
%     for t = 1:size(MTS,2)
%         M(t) = squeeze(mean(MTS(MTT(:,t),t)))';
%     end
%     subp(l) = plot(temp_time  ,M,'linewidth',1.5,'color',Compcolors(l,:));
%     
%     for t = 1:size(MTS,2)
%         S(t) = squeeze(std(MTS(MTT(:,t),t)))'./sqrt(sum(MTT(:,t)));
%     end
%     %S = squeeze(std(MT(:,cc(l),:)))'./sqrt(nboots);
%     h = fill([temp_time flip(temp_time)],[M+S flip(M-S)],Compcolors(l,:),'edgecolor','none');
%     %h = fill([temp_time flip(temp_time)],[CIu; flip(CId)],Compcolors(l,:),'edgecolor','none');
%     set(h,'facealpha',.3)
% end
% legend(subp,{'C1','C2','C3','C4','C5'})
% 
% xlim([-.3 1])
% vline(0,'k--')
% xlabel('Time (S)')
% ylabel('Component Amplitude');
% set(gca,'fontsize',14)
% export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_TemporalDynamics_laminar']),'-pdf','-r200')
% 
% %% Now indicate the significant evoked components
% % compare clusters with average prestim data...
% 
% FIG = figure;
% set(FIG,'unit','inch','position',[0 0 10 15],'color','w');
% 
% pre_win = find(temp_time<0.0 & temp_time>-0.2);
% 
% for c = 1:NComp
%     subplot(NComp,1,c); hold on;
%     comp = Comp_ord(c);
%     Pres = MT(:,comp,pre_win);   
%     
%     for t = find(temp_time>-0.3)
%         %[~,PVal(t,c),~,stat_temp] = ttest(MT(:,comp,t),mean(MT(:,comp,pre_win)));
%         [~,PVal(t,c),~,stat_temp] = ttest2(MT(:,comp,t),Pres(:));
%         Tstat(t,c) = stat_temp.tstat;
%     end
%     
%      M = squeeze(mean(MT(:,comp,:)))';
%      S = squeeze(std(MT(:,comp,:)))'./sqrt(nboots);
%      minY = min(M-S);
%      maxY = max(M+S);
%     % plot the significant results as shaded
%     SigRes = find(PVal(:,c)<.01);
%     plot(temp_time(SigRes),repmat(maxY*1.1,[1 size(SigRes,1)]),'.','color','k','linewidth',2)
%     
%     % The average temporal loading
%     subp(c) = plot(temp_time  ,squeeze(mean(MT(:,comp,:)))','linewidth',1.5,'color',Compcolors(c,:));
%     h = fill([temp_time flip(temp_time)],[M+S flip(M-S)],Compcolors(c,:),'edgecolor','none');
%     set(h,'facealpha',.3)
%     
%     xlim([-0.2 1])
%     %ylim([0 11])
%     if c==NComp
%         xlabel('time (s)')
%         ylabel('Amplitude')
%     end
%     title(['C' num2str(c)])
% end
% 
% % Next step will be to include cluster level analysis
% export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_TemporalDynamics_Stat_laminar']),'-pdf','-r200')
% 
% %% Now indicate the significant evoked components
% % compare clusters with average prestim data...
% 
% FIG = figure;
% set(FIG,'unit','inch','position',[0 0 10 15],'color','w');
% 
% pre_win = find(temp_time<0.0 & temp_time>-0.2);
% 
% for c = 1:NComp
%     subplot(NComp,1,c); hold on;
%     comp = Comp_ord(c);
%     
% 
%     MTN = (MT - mean(MT(:,:,temp_time<0 & temp_time>-.3),3));%./mean(MT(:,:,temp_time<0 & temp_time>-.3),3);
%     M = squeeze(mean(MTN(:,comp,:)))';
%     X = squeeze(MTN(:,comp,:));
%     Y = repmat(temp_time,[nboots,1]);
%     scatter(Y(:),X(:),5,'filled');
% 
%     xlim([-0.2 1])
%     %ylim([0 11])
%     if c==NComp
%         xlabel('time (s)')
%         ylabel('Amplitude')
%     end
%     title(['C' num2str(c)])
% end
% 
% export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_TemporalDynamics_All_laminar']),'-pdf','-r200')
% 
% 
% 
