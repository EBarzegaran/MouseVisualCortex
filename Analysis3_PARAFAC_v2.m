% This one requires N-way toolbox

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
%addpath(genpath('E:\Elham\Codes\NonGit\nway331'));% add the nway toolbox
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nway331'));% add the nway toolbox
clear; clc;
FileName = {'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098','drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098'};
DataPath = 'Data_Temp';

%Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
%% Load the iPDCs for all animals
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
load(fullfile(Path, ['STOK_ALL_' FileName{1} '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults','PARAFAC');

% extract the distance and RF distances
load (fullfile(Path,'Probe_Data_All.mat'));
StokALL = DistanceEstimate(StokALL,Probe_all);
StokALL = RFDistanceEstimate(StokALL,Probe_all);

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
    %DistanceBoots = DistanceForBoot(StokALL,BootIDs);
    PARRES{cond} = load(fullfile(DataPath,[FileName{cond} 'PARAFAC_covtemp_' num2str(NComp) '_ExtraVar.mat']));
    PARRES{cond}.DistanceBoots = DistanceForBoot(StokALL,BootIDs);
end

%% Main figure of all components
close all;
Freq = 1:100;
con_mode = 1;
% plotting params
offs = -.06;
FS = 11;
FIG = figure(1);
set(FIG,'unit','inch','position',[0 0 2.5*NComp 9],'color','w')
lstyle = {'-','-'};
lcolor = {'k',[.5 .5 .5]};
connames = {'High Contrast','Low Contrast'};
        
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
    % all target and source layers to indicate the significant layers
    S_all = cat(3,model_reord{:,3});
    T_all = cat(3,model_reord{:,2});
    
    for c = 1:NComp
        figure(1);
        cc = Comp_ord(c);
        %---------------------significant layers of target and source------
        S_all_temp = squeeze(S_all(:,c,:));
        T_all_temp = squeeze(T_all(:,c,:));
        for l = 1:6
            TH= .025;
            %[H(l), P(l)] = ttest(S_all_temp(l,:)-sqrt(1./6));
            temp = sort(S_all_temp(l,:)); temp(temp==0)=[];
            CI = [temp(round(numel(temp)*TH)) temp(round(numel(temp)*(1-TH)))];
            H_S(l,c,cond) = ~(CI(1)<sqrt(1./6) & CI(2)>sqrt(1./6)) & CI(1)>sqrt(1./6);
            
            temp = sort(T_all_temp(l,:)); temp(temp==0)=[];
            CI = [temp(round(numel(temp)*TH)) temp(round(numel(temp)*(1-TH)))];
            H_T(l,c,cond) = ~(CI(1)<sqrt(1./6) & CI(2)>sqrt(1./6)) & CI(1)>sqrt(1./6);
        end
        
        %-----------------------Cond1 Source and Target-------------------- 
        if cond==1        
            MS1(cc,:) = model_temp_M{3}(end:-1:1,cc);
            SS1(cc,:)= model_temp_S{3}(end:-1:1,cc);
            MT1(cc,:) = model_temp_M{2}(:,cc);
            ST1(cc,:) = model_temp_S{2}(:,cc);
        end
        
        if cond==2
           
            axes('position',[.17+(c-1)*.23 .98 .18 .13])
            text(0, 0,['Network' num2str(c)],'fontsize',FS+1,'HorizontalAlignment' ,'center');
            axis off
                 
            MS = axes('position',[.07+(c-1)*.23 .81 .08 .13]);%axes('position',[.07+(c-1)*.23 .6 .18 .13]);
            hold on;
           
            MS2(cc,:) = model_temp_M{3}(end:-1:1,cc);
            SS2(cc,:) = model_temp_S{3}(end:-1:1,cc);
            
            B = bar([MS2(cc,:); MS1(cc,:)]',1,'EdgeColor',[0 0 0]); %colormap(MS,jmaColors('coolhot'))
            B(2).FaceColor = [.5 .5 .5];
            B(1).FaceColor = [.9 .9 .9];
            for C_temp = 1:2
                xData = B(C_temp).XData+(1.5-C_temp)*-.3;%B(C_temp).XOffset;
                MS_temp=eval(['MS' num2str(3-C_temp)]);
                SS_temp=eval(['SS' num2str(3-C_temp)]);
                errorbar(xData,MS_temp(cc,:),SS_temp(cc,:),'.','color','k')
            end
            camroll(90)
            %set(MS,'yDir','reverse')
            set(MS,'xtick',[],'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            ylim([0 1.2])
            title('Source','fontweight','normal');
            box off;      
            
             %--------------------plot noise level--------------------------
            xlim([.2 6.8])
            noisel = sqrt(1./6);
            line([.2 6.8],[noisel noisel],'color','b','linestyle','--','linewidth',1);
            %---------------------larger than noise------------------------
            for l = 1:6
                if sum(H_S(l,c,:))>0
                    %text(7-l-.3, 1.2,'*','fontsize',20,'color','b');
                    text(7-l-.3, MS_temp(c,7-l)+SS_temp(c,7-l)+.3,'*','fontsize',20,'color','b');
                end
            end
            %--------------------------------------------------------------
            
            %--------------------------PLOT TARGET-----------------------------
            MT = axes('position',[.18+(c-1)*.23 .81 .08 .13]);hold on;
            MT2(cc,:) = model_temp_M{3}(:,cc);
            ST2(cc,:) = model_temp_S{3}(:,cc);
            
            B = bar([MT2(cc,:); MT1(cc,:)]',1,'EdgeColor',[0 0 0]); %colormap(MS,jmaColors('coolhot'))
            B(1).FaceColor = [.5 .5 .5];
            B(2).FaceColor = [.9 .9 .9];
            for C_temp = 1:2
                xData = B(C_temp).XData+(1.5-C_temp)*-.3;
                MT_temp=eval(['MT' num2str(3-C_temp)]);
                ST_temp=eval(['ST' num2str(3-C_temp)]);
                errorbar((xData),MT_temp(cc,:),ST_temp(cc,:),'.','color','k')
            end
            
            view(90,90)
            set(MT,'xtick',1:6,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            ylim([0 1.2])
            title('Target','fontweight','normal')
            if cc ==1
                XL = ylabel('Strength');
                set(XL,'position',get(XL,'position')-[.05 .8 0])
            end
            box off;
            
            if cc==NComp
%                 lgnd = legend(B,connames);
%                  set(lgnd,'color','none','box','off','fontsize',FS-2,'position',get(lgnd,'position')+[-.035 -.04 0 0]);
%                  lgnd.ItemTokenSize=[8 6]
            end
            %annotation('arrow',[(cc-1)*.217+.163 (cc-1)*.217+.163+.03], [0.917 0.917],'color','b','linewidth',1.5 )
            
             %--------------------plot noise level--------------------------
            xlim([.2 6.8])
            noisel = sqrt(1./6);
            L = line([.2 6.8],[noisel noisel],'color','b','linestyle','--','linewidth',1);
            %---------------------larger than noise------------------------
            for l = 1:6
                if sum(H_T(l,c,:))>0
                    %text(l+.3, 1.1,'*','fontsize',20,'color','b');
                    text(l+.3, MT_temp(c,l)+ST_temp(c,l)+.2,'*','fontsize',20,'color','b');
                end
            end
            %--------------------------------------------------------------
            if c==1
                lgnd = legend(L,'Uniform loading distribution');
                set(lgnd,'color','none','box','off','fontsize',FS,'position',get(lgnd,'position')+[-0.12 -0.90 0 0]);
            end
            
        end
        %------------------------PLOT TIME MODE----------------------------
        temp_time_ms = temp_time;%*1000;
        M = model_temp_M{5}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:))*100;
        if cond ==1
            TP = axes('position',[.07+(c-1)*.23 .44 .18 .13]);
        else
            TP = axes('position',[.07+(c-1)*.23 .23 .18 .13]);
        end
       % TP = subplot(5,NComp,(c-1)+NComp*2+1+((cond-1)*NComp)); 
        hold on;
        h = fill([temp_time_ms flip(temp_time_ms)],[M+model_temp_S{5}(:,cc); flip(M-model_temp_S{5}(:,cc))]',lcolor{cond},'edgecolor','none');
        set(h,'facealpha',.5)
        PLT = plot(temp_time_ms,M,'linewidth',1.5,'color',lcolor{cond});
        
        ylim([-40 130])
        xlim([-.1 .85]);
        vline(0,'k--');
        hline(0,'k--');
        set(TP,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])
        ytickangle(90);
        %xtickangle(30)
        if c==1 && cond==2
            XL = xlabel('Time (s)');
            set(XL,'position',get(XL,'position')+[0 8 0])
            
        end
            
        if cc~=1
            %set(TP,'yticklabel',[]);
        end
        if cond==1
            set(TP,'xticklabel',[]);
        else
            set(TP,'position',get(TP,'position')+[0 .068 .0 0],'fontsize',FS)
        end
        box off
        
        
        %title(['Consistency=' num2str(round(Comp_PerSig(cc,5),2))])
        %---------------- STATS ------------------
        maxY = 90;
        Pres = mean(M_temp(:,cc,pre_win),3);   

        for t = find(temp_time>-0.0)
            [~,PVal(t,cc),~,stat_temp] = ttest2(M_temp(:,cc,t),Pres(:));
            Tstat(t,cc) = stat_temp.tstat;
        end
        PVal(temp_time<.03,:)=1;
        SigRes1 = PVal(:,cc)<.001;
        CC = bwconncomp(SigRes1);
        for i = 1:numel(CC.PixelIdxList)
            if numel(CC.PixelIdxList{i})>10
                SigRes = CC.PixelIdxList{i};
                x = temp_time(SigRes);
                if cond ==1
                    y = ones(numel(x),1)'*120;
                else
                    y = ones(numel(x),1)'*120;
                end
                z = zeros(size(x));
                col = abs(Tstat(SigRes,cc))';  % This is the color, vary with x in this case.

                F = fill([x flip(x)],[zeros(size(y))-40 y],'b','edgecolor','none');
                set(F,'facealpha',.13)    
            end
        end
        
        if c==NComp
%             lgnd = legend(PLT,connames{cond});
%             
%              set(lgnd,'color','none','box','off','fontsize',FS-2,'position',get(lgnd,'position')+[.03 0 0 0]);
%              lgnd.ItemTokenSize=[0 18]   
            text(.7,100,connames{cond},'fontsize',FS-2,'HorizontalAlignment' ,'center')
        elseif c==3
            if cond==1
                lgnd = legend(F,'Significant evoked response');
                set(lgnd,'color','none','box','off','fontsize',FS,'position',get(lgnd,'position')+[.2 -0.525 0 0]);
            end
        end
        %-------------------PLOT FREQUENCY DISTRIBUTION--------------------
        if cond==1
            FP(c) = axes('position',[.07+(c-1)*.23 .63 .18 .13]);%subplot(5,NComp,(c-1)+NComp*1+1); hold on;
        else
            axes(FP(c));
        end
        hold on;
        h = fill([Freq flip(Freq)],[model_temp_M{4}(:,cc)+model_temp_S{4}(:,cc); flip(model_temp_M{4}(:,cc)-model_temp_S{4}(:,cc))],lcolor{cond},'edgecolor','none');
        set(h,'facealpha',.75-.25*cond)
        PL(cond) = plot(model_temp_M{4}(:,cc),'linewidth',1.2,'color',lcolor{cond},'linestyle',lstyle{cond});
        
        %--------------------plot noise level--------------------------
        %noisel = sqrt(1./100);
        %L = line([0 100],[noisel noisel],'color','b','linestyle','--','linewidth',1);
        %--------------------------------------------------------------
         
        if cc~=1
            %set(FP(c),'yticklabel',[]);
        end

        if cond==2
            %ylim([0.05 .25])
            YLIM = ylim;

            TICKS = 0:round(YLIM(2)/3,2):YLIM(2)*1.5;
            set(FP(c),'fontsize',FS,'ytick', TICKS(1:end-1));
            ylim([0 YLIM(2)-eps])
            
            if c==NComp
                lgnd = legend(PL,connames);
                set(lgnd,'color','none','box','off','fontsize',FS-2,'position',get(lgnd,'position')+[.03 0 0 0]);
                lgnd.ItemTokenSize=[20 18];
                
            elseif c==1
%                 lgnd = legend(L,'Uniform loading distribution');
%                 set(lgnd,'color','none','box','off','fontsize',FS,'position',get(lgnd,'position')+[0.1 -0.71 0 0]);
            end    
        end
        set(gca,'TickDir','out','linewidth',1.5,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],...
            'xtick',0:20:100);
        ytickangle(90)
        
        if cc==1
            
            XL = xlabel('Frequency (Hz)');
            set(XL,'position',get(XL,'position')+[0 .005 0])
            
        end
         
        box off;
           
        
        % ----------------PLOT Hierarchy scores----------------------------

        Con_temp = squeeze(abs(cat(5,model_reord{:,con_mode})));
        for perm = 1:size(Con_temp,3)
            CN = zeros(36,1);
            CN(indTotal) = Con_temp(:,cc,perm);
            CN = reshape(CN,6,6);

            if true
%                 CN = CN-1./sqrt(30);
%                 CN(CN<0)=0;
                DAI = (CN'-CN)./(CN+CN');                
                DAI(isnan(DAI))=0;
                DAII = DAI;
                %-------------surrogate data-------------------------------
                CN_sur = rand(size(DAI));
                DAI_sur = (CN_sur'-CN_sur)./(CN_sur+CN_sur');                
                DAI_sur(isnan(DAI_sur))=0;
                %---------------------------- Bastos paper------------------------
                % (1)rescale to -3 to 3
                 DAI = (DAI-min(DAI(:)))./(max(DAI(:))-min(DAI(:)));
                 DAI = DAI*6-3;
                 
                 DAI_sur = (DAI_sur-min(DAI_sur(:)))./(max(DAI_sur(:))-min(DAI_sur(:)));
                 DAI_sur = DAI_sur*6-3;
                 
                 
                %DAI = (DAI)*3;
                % (2) for each target (row) shifte the rescaled mDAI values of all source areas such that the smallest value was one. 
                for roi = 1:size(DAI,1)
                    ind = setdiff(1:size(DAI,1),roi);
                    %DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
                    DAI(roi,ind) = DAI(roi,ind)-min(DAI(roi,ind))+1;
                    DAI(roi,roi) = NaN;
                    
                    DAI_sur(roi,ind) = DAI_sur(roi,ind)-min(DAI_sur(roi,ind))+1;
                    DAI_sur(roi,roi) = NaN;
                end
                Hscore(:,:,perm,cond) = DAI;
                Hscore_sur(:,:,perm,cond) = DAI_sur;

            else

                Hscore(:,perm,cond)   = (HierarchyExtract(CN));

            end
        end

        %-----------------plot the hierarchies-------------------------
        if cond==1
            CNP = axes('position',[.07+(c-1)*.23 .075 .18 .15]);%subplot(5,NComp,(c-1)+NComp*4+1);
            if true
                hold on;
                %line(1:6,squeeze(nanmean(nanmean(Hscore_sur,1),3)),'color','b','linestyle','--','linewidth',1);
                boxplot(squeeze(nanmean(nanmean(Hscore,1),4))','colors','k','symbol','')%,'PlotStyle','compact'
                
            else
                boxplot(squeeze(Hscore)','colors','k','symbol','')%,'PlotStyle','compact'
            end
             
            %-----------------Alternative----------------------
%             imagesc(DAII*-1);caxis([-1 1])
%             colormap(jmaColors2('coolhotlight'))
%             set(gca,'xtick',1:6,'xticklabel',ROISN,'ytick',1:6,'yticklabel',ROISN);
            %-----------------Alternative----------------------
            box off
            xtickangle(45)
            ytickangle(90)
            set(gca,'xtick',1:6,'xticklabel',ROISN);
            ylim([.2 6])

            h = findobj(gca,'Tag','Box'); 
            for j=1:length(h) 
                patch(get(h(j),'XData'),get(h(j),'YData'),Colors(numel(ROIs)-j+1,:),'FaceAlpha',.8);
            end 

            if cc==1
                set(CNP,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            else
                set(CNP,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            end
            if cc==1
                
                axes('position',[.02 .81 .02 .15]);
                text(0, .5,'Laminar Layers','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                axes('position',[.02 .63 .02 .15]);
                text(0, .5,'Frequency Loading','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                axes('position',[.02 .35 .02 .18]);
                text(0, .5,'Evoked change (%)','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                
                axes('position',[.02 .08 .02 .18]);
                text(0, .5,'Hierarchy Score','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                %text(0, .5,'DAI','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                axes('position',[.4 .028 .02 .18]);
                text(-1, -0.03,'*','fontsize',20,'color','b');
                text(0, 0,'Larger than uniform distribution','fontsize',FS);
                %text(0, .5,'DAI','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
            end
            
        end
    end
end
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_HarrisM']),'-pdf','-r200')
%% variance explained over time
close;
Freq = 1:100;
con_mode = 1;
% plotting params
offs = -.06;
FS = 12;
FIG = figure(1);
set(FIG,'unit','inch','position',[0 0 9 4],'color','w')
lstyle = {'-','-'};
lcolor = {'k',[.5 .5 .5]};
connames = {'Contrast = 0.8','Contrast = 0.1'};
Compcol = brewermap(10,'Paired');
Compcol = Compcol([8 6 10 2],:);

pre_win = find(temp_time<-0.0 & temp_time>-0.5);
for cond = 1:2
    
    Model_reord =   PARRES{cond}.Model_reord;
    Comp_ord    =   PARRES{cond}.Comp_ord;
    % temporal dynamics: significance
    M_temp = cellfun(@(x) (x{5}),Model_reord,'uni',false);
    M_temp = cellfun(@(x) ((x)),M_temp,'uni',false);
    M_temp = cat(3,M_temp{:});
    %M_temp = M_temp-mean(M_temp(pre_win,:,:));

    M_vars = mean(M_temp(:,:,:),3);
    S_vars = var(M_temp(:,:,:),[],3);
    subplot(1,2,cond),hold on;
    
    for c = 1:NComp
         
        if c==1
            MData = ((M_vars(:,c)./sum(M_vars,2)))*100-20;
            
        else
           MData = ((M_vars(:,c)./sum(M_vars,2)))*100;
        end
        h = fill([temp_time flip(temp_time)],...
             [MData+S_vars(:,c); flip(MData-S_vars(:,c))],Compcol(c,:),'edgecolor','none');
             set(h,'facealpha',.5)
        P(c) = plot(temp_time,MData,'linewidth',1.5,'color',Compcol(c,:));
    end
    
    ylim([7 47])
    vline(0,'k--')
    hline(25,'k-')
    set(gca,'ytick',10:5:45,'yticklabel',[10:5:20 45:5:65])
    if cond ==1
        LG = legend(P,'SubNet 1','SubNet 2','SubNet 3','SubNet 4');
        set(LG,'box','off')
        xlabel('Time (s)')
        ylabel('Variance Explained %')
        %set(gca,'position',get(gca,'position')+[0 .04 0 0])
    end
    
    set(gca,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])
    %set(gca,'position',get(gca,'position')+[-0.04 .0 0 0])
    box off
    axis tight;
    title(connames{cond})
    xlim([-.2 1])
end
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_Variance']),'-pdf','-r200')

%% define different distribution functions to fit to frequency spectrum
clear model_fun;
close all;
model_fun.powerlaw = @(p,x) ( p(2)*(x).^-p(1)) ;
%model_fun.expo = @(p,x)( p(2) *exp(-p(1)*x) );
model_fun.lognormal = @(p,x) (p(3)+(p(4)./(x*p(2)*sqrt(2*pi))).* exp(-(log(x)-p(1)).^2./(2*(p(2).^2))));
model_fun.normal = @(p,x) (p(3)+(p(4)./(p(2)*sqrt(2*pi)))*exp(-.5*((x-p(1))./p(2)).^2)) ;

FIG = figure;
set(FIG,'unit','inch','position',[5 5 12 6],'color','w')
FS = 12;

Model_names = fieldnames(model_fun);
nb=500;
Conds = {'Contrast = 0.8','Contrast = 0.1'};
for comp = 1:NComp
    for C = 1:2
        clear P_test MSE this_p good;
        for b = 1:nb
            
            P_test(b,:) = PARRES{C}.Model_reord{b}{4}(:,comp);
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'bisquare';
            if comp>1 & comp<5
                MSE(b,1) = Inf;
                this_p{1}(:,b)=[0 0 0];
            else
                try
                    [this_p{1}(:,b),~,~,~,MSE(b,1)] =nlinfit(Freq(1:end),P_test(b,1:end),model_fun.powerlaw,[1.,1,0],opts);
                catch
                    MSE(b,1) = nan;
                    this_p{1}(:,b)=[nan nan nan];
                end
            end
            %[this_p{2}(:,b),~,~,~,MSE(b,2)] =nlinfit(Freq(1:end),P_test(b,1:end),model_fun.expo,[1.,1],opts);
            Pdf = P_test(b,:);%./sum(P_test(b,:));
            [this_p{2}(:,b),~,~,~,MSE(b,2)] =nlinfit( Freq(1:end),Pdf,model_fun.lognormal,[4.,5,.1,1],opts);
            [this_p{3}(:,b),~,~,~,MSE(b,3)] =nlinfit(Freq(1:end),Pdf,model_fun.normal,[50,10,.1,1],opts);

            

        end

         % Plot the results
         subplot(2,4,comp+(C-1)*NComp)
         %
         M = mean((P_test));
         SD = std((P_test));
         PL(1)=plot(Freq,mean(P_test),'k','linewidth',1.5);hold on;

         h = fill(([Freq flip(Freq)]),[M+SD flip(M-SD)],'k','edgecolor','none');
         set(h,'facealpha',.3)
         [~,ind_fun] = min(nanmean(MSE));
         %-----------------Remove the bad fits?----------------------------
         
         for i = 1:numel(Model_names)
            %[~,G] = rmoutliers(MSE(:,i));
            
            RM = this_p{i}>100 | this_p{i}<-100;
            good(:,i) = sum(RM)';
         end
         %-----------------------------------------------------------------
         
         Is = 1:nb;
         y = arrayfun(@(x) model_fun.(Model_names{ind_fun})(this_p{ind_fun}(:,x),Freq),Is(good(:,ind_fun)==0),'uni',false);
         y = cat(1,y{:});
         if comp<4
             PL(2) = plot(Freq,nanmean(y),'--r','linewidth',1.5);hold on;
         else
             PL(2) = plot(Freq(5:end),mean(y(:,5:end)),'--r','linewidth',1.5);hold on;
             ylim([0 .3])
         end
         if comp>2
             HFreq = .1:.1:100;
             y = arrayfun(@(x) model_fun.(Model_names{ind_fun})(this_p{ind_fun}(:,x),HFreq),Is(good(:,ind_fun)==0),'uni',false);
             y = nanmean(cat(1,y{:}));
             HM = mean([mean(y),max(y)]);
             FWHM = HFreq(find(diff(sign(y-HM))))
         end
        %title(Conds{C})
        Params = nanmean(this_p{ind_fun}(:,good(:,ind_fun)==0),2);
        if ind_fun==1
            title([Model_names{ind_fun} ', \beta=' num2str(round(Params(1),1))],'fontweight','normal')
        else
            title([Model_names{ind_fun} ', \mu=' num2str(round(Params(1),1)) ', \sigma=' num2str(round(Params(2),1))],'fontweight','normal')
        end
        
        box off
         set(gca,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],...
             'xtick',0:20:100)
        if C==2 && comp==1
            
            xlabel('Frequency (Hz)')
            ylabel('Normalized Loading')
            %set(gca,'position',get(gca,'position')+[0 .04 0 0])
            
            %title(['P1=' num2str(Params(1)) ' P2=' num2str(Params(2))])
        end
        
        if C==2 && comp==NComp
            L = legend(PL,{'Average frequency loading','Fitted distribution'},'fontsize',FS-3,'box','off');

            L.ItemTokenSize=[8 6];
        end
    end
end

axes('position',[.07 .7 .1 .1])
text(0,0,'Contrast=0.8','rotation',90,'fontsize',FS,'fontweight','bold')
axis off;
axes('position',[.07 .2 .1 .1])
text(0,0,'Contrast=0.1','rotation',90,'fontsize',FS,'fontweight','bold')
axis off

for comp=1:4
    axes('position',[-.02+.2*comp .99 .1 .1])
    text(0,0,['Subnetwork' num2str(comp)],'fontsize',FS,'fontweight','bold')
    axis off;
end

export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_FrequencyFit']),'-pdf','-r200')
%% FFT on the temporal loadings
FIG = figure;
set(FIG,'unit','inch','color','w','position',[5 5 10 4])

for comp = 1:NComp
    for C = 1:2
        clear P_test MSE this_p good;
        for b = 1:nb
            
            P = PARRES{C}.Model_reord{b}{5}(:,comp);
            P_test(b,:) = P./norm(P);
        end
        Pfft = fft(P_test,[],2);
        L = 375;
        Fs = 250;
        f = Fs*(0:(L/2))/L;
        subplot(2,4,comp+(C-1)*4);hold on;
        M = mean(abs(Pfft(:,3:100)).^2);
        plot(f(3:100),M,'linewidth',1.5,'color','k');
        %bar(f(3:100),mean(abs(Pfft(:,3:100)).^2));
        S = std(abs(Pfft(:,3:100)).^2);
        h = fill([f(3:100) flip(f(3:100))],[M+S flip(M-S)],'k','edgecolor','none');
        set(h,'facealpha',.5)
        
        YL = ylim;
        ylim([0 .4])
        if C==1
            title(['Network' num2str(comp)])
        end
        set(gca,'fontsize',FS)
        if comp==1
            ylabel(connames{C})
            if C==2
                xlabel('Frequency(Hz)')
            end
        end
        %
        xlim([2 20])
    end
end

axes('position',[0.05 ,.5, .1,.1]);
text(.0,.5,'FFT of Temporal Loading','rotation',90,'HorizontalAlignment','Center','fontsize',FS)
axis off
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_TemporalFFT']),'-pdf','-r200')

%% Reconstruct the data

for cond = 1:2
    
    Model_reord =   PARRES{cond}.Model_reord;
    for m = 1:numel(Model_reord)
        m
        Data = Loading2Data(Model_reord{m});
        Data = cellfun(@(x) (mean(mean(mean(mean(x(:,:,:,:,temp_time>0),2),3),4),5)),Data,'uni',false);
        Data_con{cond}(:,:,m) = cat(2,Data{:});
    end
end
%% RF distance analysis, carefull: keep the connections from the same layer
 % first extract the loadings of between area connections
 close all;
 Freq = 1:100;
con_mode = 1;
Compcol = brewermap(10,'Paired');
Compcol = Compcol([8 6 10 2],:);
%close ;
 FIG1 = figure(1);
 set(FIG1,'unit','inch','position',[0 0 10 4],'color','w')
 
 FIG2 = figure(2);
 set(FIG2,'unit','inch','position',[0 0 10 10],'color','w')
 for cond = 1:2
    
    Model_reord =   PARRES{cond}.Model_reord;
    model_reord =   cat(1,Model_reord{:}); 
    Con_temp    =   squeeze(abs(cat(5,model_reord{:,con_mode})));
    for c=1:NComp
        for b   =   1:size(Con_temp,3)
            CN  =   zeros(36,1);          
              CN(indTotal)    = Data_con{cond}(:,c,b); 
              CN              = reshape(CN,6,6);
              CNL(:,:,b)      = CN;
              
%               DAI =   ((CN-CN')./(CN+CN'));
%               DAIL(:,:,b)     = DAI;
%               CNL(:,:,b)      =DAIL (:,:,b);
        end
        dists   = PARRES{cond}.DistanceBoots.RFDists;
        RFNan   = PARRES{cond}.DistanceBoots.RFDists_nans(:,:,4,:);
        d       = squeeze((dists(:,:,4,:)));
        XY      =[];
        XYM     =[];
        for roi1 = 1:NROIs
%             if c<3
                r2 = roi1+1:NROIs;
%             else

%             r2 = 1:roi1;%NROIs
%                 r2 = 1:NROIs
%             end
            for roi2=r2 
                if roi1~=roi2
                    x   = squeeze(d(roi1,roi2,~isnan(d(roi1,roi2,:)) & CNL(roi1,roi2,:)~=0 & RFNan(roi1,roi2,:)<4))*10;
                    y   = squeeze(CNL(roi1,roi2,~isnan(d(roi1,roi2,:)) & CNL(roi1,roi2,:)~=0 & RFNan(roi1,roi2,:)<4));
                    XY  = [XY; [x(:) y(:)]];
                    XYM = [XYM; [nanmean(x) nanmean(y)]];
                    % individual connections
                     [CR(roi1,roi2,c,cond),PR(roi1,roi2,c,cond)] = corr(x,y);
                     p = polyfit(x,y,1);
                     x1 = linspace(0,60);
                     y1 = polyval(p,x1);
%                      if cond==1
%                          figure(2)
%                          subplot(NROIs,NROIs,roi1+(roi2-1)*NROIs);hold on;
%     %                     scatter(x,y,15,'filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.2);
%                          plot(x1,y1,'--','linewidth',1.5,'color',Compcol(c,:))
%     %                     xlim([0 5])
%                      end
                end
            end
        end
        figure(1);
        subplot(2,NComp,c+(cond-1)*NComp);hold on;
        %-----------------figure 1-all point-------------------------------
        hist3(XY,'CdataMode','auto','nbins',[20 20])
%         xlabel('Distance')
%         ylabel('iPDC')
        colorbar
        view(2)
        axis xy;
        
        %------------------------------------------------------------------
        
%         XYM((isnan(XYM(:,2))),:)=[];
%         scatter(XYM(:,1),XYM(:,2),10,'filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.2);
       [C,P] = corr(XY(:,1),XY(:,2));
       if isnan(C)
           disp(C);
       end
        title(['Corr=' num2str(round(C,2)) ]);
        p = polyfit(XY(:,1),XY(:,2),1);
        x1 = linspace(0,60);
        y1 = polyval(p,x1);
        plot(x1,y1,'w--','linewidth',1.8)
        plot(y1,x1,'w--','linewidth',1.8)
        if c==1
            ylim([0 0.025])
            if cond==2
                xlabel('Distance')
                ylabel('iPDC')
            end
        else    
           
           ylim([0 0.01])
        end
        xlim([0 60])
%         ylim([-1 1]) 
%         hline(0,'k--')
        %------------------------------------------------------------------
    end
 end
%squeeze(mean(mean(squeeze(CR).*squeeze(PR<.05),1),2))
 colormap('jet')
export_fig(FIG1,fullfile(FigPath,['PARAFAC_RFDistance_Corr_Overall_poststim_iPDC_evoked']),'-pdf','-r200');
%export_fig(FIG1,fullfile(FigPath,['PARAFAC_Distance_Corr_Overall_poststim_iPDC']),'-pdf','-r200');