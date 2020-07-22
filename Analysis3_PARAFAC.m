% This one requires N-way toolbox

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
%addpath(genpath('E:\Elham\Codes\NonGit\nway331'));% add the nway toolbox
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nway331'));% add the nway toolbox
clear; clc;
FileName = {'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098','drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098'};
DataPath = 'Data_Temp';

%Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
%%
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
load(fullfile(Path, ['STOK_ALL_' FileName{1} '.mat']));
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
nboots = 200;
bootsize = 8;

redoanalysis = false;
NComp = 4;

for cond = 1:2
    if redoanalysis || ~exist([FileName{cond} 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar.mat'],'file')
        BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);
        PARAFAC_FC(StokALL,NComp,BootIDs,nboots,ROIs,DataPath,FigPath,FileName{cond});    
    end
    load(fullfile(DataPath,[FileName{cond} 'PARAFAC_covtemp_' num2str(NComp)]));
    PARRES{cond} = load(fullfile(DataPath,[FileName{cond} 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar.mat']));
end
%% New figure
Freq = 1:100;
con_mode = 1;
% plotting params
offs = -.06;
FS = 14;
FIG = figure(1);
set(FIG,'unit','inch','position',[0 0 2.5*NComp 9],'color','w')
lstyle = {'-','-'};
lcolor = {'k',[.5 .5 .5]};
connames = {'C = .8','C = .1'};
        
for cond = 1:2
    
    Model_reord =   PARRES{cond}.Model_reord;
    Comp_ord    =   PARRES{cond}.Comp_ord;
    % temporal dynamics: significance
    pre_win = find(temp_time<-0.0 & temp_time>-0.5);
    M_temp = cellfun(@(x) (x{5}),Model_reord,'uni',false);
    M_temp = permute(cat(3,M_temp{:}),[3 2 1]); 
    % Bastos based hierarchy
    clear Hscore;
    model_reord =   cat(1,Model_reord{:});


    model_temp_M = arrayfun(@(x) mean(abs(cat(3,model_reord{:,x})),3),1:size(model_reord,2),'uni',false);
    model_temp_S = arrayfun(@(x) std(abs(cat(3,model_reord{:,x})),[],3),1:size(model_reord,2),'uni',false);
    % different for temporal dimension: first %change and then mean and std
    mod_temp =5; % dimension where temporal mode is located
    temporal_temp = cat(3,model_reord{:,mod_temp});
    temporal_temp = (temporal_temp-mean(temporal_temp(temp_time<0,:,:)))./mean(temporal_temp(temp_time<0,:,:))*100;
    model_temp_M{mod_temp} = mean(temporal_temp,3);
    model_temp_S{mod_temp} = std(temporal_temp,[],3);

    for c = 1:NComp
        figure(1);
        cc = Comp_ord(c);
        %-----------------------------PLOT SOURCE--------------------------   
        if cond==1
            MS = subplot(5,NComp*2,(2*c-1)); hold on;
            bar(model_temp_M{3}(:,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 .9]); %colormap(MS,jmaColors('coolhot'))
            errorbar(1:6,model_temp_M{3}(:,cc),model_temp_S{3}(:,cc),'.','color','k')
            view([-90 -90])
            set(MS,'position',get(MS,'position')+[offs+(.017*c) -.02 .0 .0],'xtick',[],'fontsize',FS);
            ylim([0 1.2])
            title(['Source']);
            box on;
            if cc==1
                XL = xlabel('Laminar Layers');
                set(XL,'position',get(XL,'position')+[0 .67 0])
            end        
            text(-3,.75,['SubNetwork' num2str(c)],'fontsize',FS)
            %--------------------------PLOT TARGET-----------------------------
            MT = subplot(5,NComp*2,(2*c));hold on;
            bar(model_temp_M{2}(:,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 .9]); colormap(MS,jmaColors('coolhot'))
            errorbar(1:6,model_temp_M{2}(:,cc),model_temp_S{2}(:,cc),'.','color','k')
            view([90 90])
            set(MT,'position',get(MT,'position')+[offs+(.017*c) -.02 .0 .0],'xtick',1:6,'fontsize',FS);
            ylim([0 1.2])
            title('Target')
            if cc ==1
               % XL = ylabel('Strength');
                %set(XL,'position',get(XL,'position')-[0.12 .5 0])
            end
            box on;

            annotation('arrow',[(cc-1)*.217+.163 (cc-1)*.217+.163+.03], [0.917 0.917],'color','b','linewidth',1.5 )
        end
        %------------------------PLOT TIME MODE----------------------------
        temp_time_ms = temp_time;%*1000;
        M = model_temp_M{5}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:))*100;
        TP = subplot(5,NComp,(c-1)+NComp*2+1+((cond-1)*NComp)); hold on;
        h = fill([temp_time_ms flip(temp_time_ms)],[M+model_temp_S{5}(:,cc); flip(M-model_temp_S{5}(:,cc))]',lcolor{cond},'edgecolor','none');
        set(h,'facealpha',.5)
        PLT = plot(temp_time_ms,M,'linewidth',1.5,'color',lcolor{cond});
        
        ylim([-40 115])
        xlim([-.1 .85]);
        vline(0,'k--');
        hline(0,'c--');
        set(TP,'position',get(TP,'position')+[offs+(.012*c) 0 .02 -.02],'fontsize',FS)
        xtickangle(30)
        if c==1 && cond==2
            xlabel('Time (sec)');
            YL = ylabel('Evoked Change (%)');
            set(YL,'position',get(YL,'position')+[0 70 0])
        end
            
        if cc~=1
            set(TP,'ytick',[]);
        end
        if cond==1
            set(TP,'xticklabel',[]);
        else
            set(TP,'position',get(TP,'position')+[0 .068 .0 0],'fontsize',FS)
        end
        box on
        
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,5),2))])
        %---------------- STATS ------------------
        maxY = 90;
        Pres = mean(M_temp(:,cc,pre_win),3);   

        for t = find(temp_time>-0.5)
            [~,PVal(t,cc),~,stat_temp] = ttest2(M_temp(:,cc,t),Pres(:));
            Tstat(t,cc) = stat_temp.tstat;
        end
        PVal(temp_time<.03)=1;
        SigRes1 = PVal(:,cc)<.01;
        CC = bwconncomp(SigRes1);
        for i = 1:numel(CC.PixelIdxList)
            if numel(CC.PixelIdxList{i})>10
                SigRes = CC.PixelIdxList{i};
                SIGF = fill([temp_time(SigRes) flip(temp_time(SigRes))],[ones(1,numel(SigRes))*-40 ones(1,numel(SigRes))*110],'r','linestyle','none');
                set(SIGF,'facealpha',.1)
                %plot(temp_time(SigRes)*1000,repmat(maxY*1.1,[1 size(SigRes,1)]),'.','color','r','linewidth',2)
            end
        end
        
        if cc==NComp
            lgnd = legend(PLT,connames{cond});
            set(lgnd,'color','none','box','off');
        end
        %-------------------PLOT FREQUENCY DISTRIBUTION--------------------
%         if cond==1
            TP = subplot(5,NComp,(c-1)+NComp*1+1); hold on;
%         else
%             subplot(5,NComp,(c-1)+NComp*1+1); hold on;
%         end
        h = fill([Freq flip(Freq)],[model_temp_M{4}(:,cc)+model_temp_S{4}(:,cc); flip(model_temp_M{4}(:,cc)-model_temp_S{4}(:,cc))],lcolor{cond},'edgecolor','none');
        set(h,'facealpha',.75-.25*cond)
        PL(cond) = plot(model_temp_M{4}(:,cc),'linewidth',1.5,'color',lcolor{cond},'linestyle',lstyle{cond});
        TPP = get(TP,'position');
        if c==1
            
            xlabel('Freq (Hz)')
            YL = ylabel('Frequency');
            set(YL,'position',get(YL,'position')+[-2.1 -.03 0])
        end
        if cc~=1
            set(TP,'ytick',[]);
        end
        box on
        if cond==2
            set(TP,'position',TPP+[offs+(.012*c) 0 .02 -.02],'fontsize',FS,'ytick',0:.1:.2);
            if c==NComp
                lgnd = legend(PL,connames);
                set(lgnd,'color','none','box','off');
            end
        end
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,4),2))])
        % ylim([0 .25])

        % ----------------PLOT Hierarchy scores----------------------------

        Con_temp = abs(cat(5,model_reord{:,con_mode}));
        for perm = 1:size(Con_temp)
            CN = zeros(36,1);
            CN(indTotal) = Con_temp(:,cc,perm);
            CN = reshape(CN,6,6);

            if true
                DAI = (CN'-CN)./(CN+CN');
                DAI(isnan(DAI))=0;
                %---------------------------- Bastos paper------------------------
                % (1)rescale to -3 to 3
                 DAI = (DAI-min(DAI(:)))./(max(DAI(:))-min(DAI(:)));
                 DAI = DAI*6-3;
                %DAI = (DAI)*3;
                % (2) for each target (row) shifte the rescaled mDAI values of all source areas such that the smallest value was one. 
                for roi = 1:size(DAI,1)
                    ind = setdiff(1:size(DAI,1),roi);
                    %DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
                    DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
                    DAI(roi,roi) = NaN;
                end
                Hscore(:,:,perm,cond) = DAI;

            else

                Hscore(:,perm,cond)   = (HierarchyExtract(CN));

            end
        end

        %-----------------plot the hierarchies-------------------------
        if cond==2
            CNP = subplot(5,NComp,(c-1)+NComp*4+1);
            if true
                boxplot(squeeze(nanmean(nanmean(Hscore,1),4))','colors','k','symbol','')%,'PlotStyle','compact'
            else
                boxplot(squeeze(Hscore)','colors','k','symbol','')%,'PlotStyle','compact'
            end
            box on
            if cc~=1
                set(gca,'ytick',[]);
            end
            box on
            xtickangle(45)
            set(gca,'xtick',1:6,'xticklabel',ROISN);

            h = findobj(gca,'Tag','Box'); 
            for j=1:length(h) 
                patch(get(h(j),'XData'),get(h(j),'YData'),Colors(numel(ROIs)-j+1,:),'FaceAlpha',.8);
            end 
            if cc==1
                set(CNP,'position',get(CNP,'position')+[offs+(.01*c) -0.0 .03 .025],'fontsize',FS);
            else
                set(CNP,'position',get(CNP,'position')+[offs+(.01*c)+.015 -0.0 0.015 .025],'fontsize',FS);
            end
            if cc==1
                YL = ylabel('Hierarchy Score');
                set(YL,'position',get(YL,'position')+[-.5 0 0])
            end
        end
    end
end
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_HarrisM']),'-pdf','-r200')

