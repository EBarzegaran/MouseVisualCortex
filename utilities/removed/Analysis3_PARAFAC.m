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
connames = {'Contrast=.8','Contrast=.1'};
        
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
            axes('position',[.17+(c-1)*.23 .98 .18 .13])
            text(0, 0,['SubNetwork' num2str(c)],'fontsize',FS+1,'HorizontalAlignment' ,'center');
            axis off
                
            %MS = subplot(5,NComp*2,(2*c-1)); 
            MS = axes('position',[.07+(c-1)*.23 .81 .08 .13]);%axes('position',[.07+(c-1)*.23 .6 .18 .13]);
            hold on;
            bar(model_temp_M{3}(end:-1:1,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]); %colormap(MS,jmaColors('coolhot'))
            errorbar(1:6,model_temp_M{3}(end:-1:1,cc),model_temp_S{3}(end:-1:1,cc),'.','color','k')
            %view([-90 -90])
            
            camroll(90)
            %set(MS,'yDir','reverse')
            set(MS,'xtick',[],'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            ylim([0 1.2])
            title('Source','fontweight','normal');
            box off;      
            
            %--------------------------PLOT TARGET-----------------------------
            MT = axes('position',[.18+(c-1)*.23 .81 .08 .13]);hold on;
            bar(model_temp_M{2}(:,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]); colormap(MS,jmaColors('coolhot'))
            errorbar(1:6,model_temp_M{2}(:,cc),model_temp_S{2}(:,cc),'.','color','k')
            %camroll(90)
            view(90,90)
            %set(MT,'yDir','reverse')
            set(MT,'xtick',1:6,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2]);
            ylim([0 1.2])
            title('Target','fontweight','normal')
            if cc ==1
                XL = ylabel('Strength');
                set(XL,'position',get(XL,'position')-[.5 .8 0])
            end
            box off;

            %annotation('arrow',[(cc-1)*.217+.163 (cc-1)*.217+.163+.03], [0.917 0.917],'color','b','linewidth',1.5 )
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
        hline(0,'c--');
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

        for t = find(temp_time>-0.5)
            [~,PVal(t,cc),~,stat_temp] = ttest2(M_temp(:,cc,t),Pres(:));
            Tstat(t,cc) = stat_temp.tstat;
        end
        PVal(temp_time<.03)=1;
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
%                 surface([x;x],[y;y],[z;z],[col;col],...
%                         'facecol','no',...
%                         'edgecol','interp',...
%                         'linew',8); 
%                 colormap(jmaColors('hotcortex'));
%                 caxis([-max(col)*.2 max(col)])
                F = fill([x flip(x)],[zeros(size(y))-40 y],'k','edgecolor','none');
                set(F,'facealpha',.13)    
            end
        end
        
        if cc==NComp
            lgnd = legend(PLT,connames{cond});
            %set(lgnd,'color','none','box','off');
             set(lgnd,'color','none','box','off','fontsize',FS-2,'position',get(lgnd,'position')+[.03 0 0 0]);
             lgnd.ItemTokenSize=[20 18]   
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
        
        
        if cc~=1
            %set(FP(c),'yticklabel',[]);
        end

        if cond==2
            YLIM = ylim;

            TICKS = 0:round(YLIM(2)/3,2):YLIM(2)*1.5;
            set(FP(c),'fontsize',FS,'ytick', TICKS(1:end-1));
            ylim([0 YLIM(2)-eps])
            
            if c==NComp
                lgnd = legend(PL,connames);
                set(lgnd,'color','none','box','off','fontsize',FS-2,'position',get(lgnd,'position')+[.03 0 0 0]);
                lgnd.ItemTokenSize=[20 18]
            end
        end
        set(gca,'TickDir','out','linewidth',1.5,'TickLength',[0.03 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],...
            'xtick',0:20:100);
        ytickangle(90)
        %ylim([0.05 .25])
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
                DAI = (CN'-CN)./(CN+CN');                
                DAI(isnan(DAI))=0;
                DAII = DAI;
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
        if cond==1
            CNP = axes('position',[.07+(c-1)*.23 .075 .18 .15]);%subplot(5,NComp,(c-1)+NComp*4+1);
            if true
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
                text(0, .5,'Frequency','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                axes('position',[.02 .35 .02 .18]);
                text(0, .5,'Evoked change (%)','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
                axis off;
                
                
                axes('position',[.02 .08 .02 .18]);
                text(0, .5,'Hierarchy Score','fontsize',FS,'HorizontalAlignment' ,'center','rotation',90);
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
set(FIG,'unit','inch','position',[0 0 10 8],'color','w')
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
    subplot(2,2,cond*2-1),hold on;
    
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
    if cond ==2
        LG = legend(P,'SubNet 1','SubNet 2','SubNet 3','SubNet 4');
        set(LG,'box','off')
        xlabel('Time (s)')
        ylabel('Variance Explained %')
        set(gca,'position',get(gca,'position')+[0 .04 0 0])
    end
    
    set(gca,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])
    set(gca,'position',get(gca,'position')+[-0.04 .0 0 0])
    box off
    axis tight;
    title(connames{cond})
    xlim([-.2 1])
end
%export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_Variance']),'-pdf','-r200')


% Estimate beta exponent for the component #1
% FIG = figure;
% set(FIG,'unit','inch','position',[1 1 4 6],'color','w')
%FS = 14;

% define the power law function to fit
%model_fun.power_law = @(p,x) ( (p(2)*x+1).^-p(1)) ;
model_fun.power_law = @(p,x) ( p(2)*(x).^-p(1)) ;
model_fun.expo = @(p,x)( p(1) *exp(-p(1)*x) );
model_fun.log_norm = @(p,x) ((1./(x*p(2)*sqrt(2*pi))).* exp(-(log(x)-p(1)).^2./(2*(p(2).^2))));

Conds = {'Contrast = 0.8','Contrast = 0.1'};
for comp = 1
    for C = 1:2
        for b = 1:420
            P_test(b,:) = PARRES{C}.Model_reord{b}{4}(:,comp);
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'bisquare';
            [this_p(:,b),R,J,CovB,MSE(b)] =nlinfit(Freq(1:end),P_test(b,1:end),model_fun.power_law,[1.,1],opts);
            %[this_p(:,b),R,J,CovB,MSE(b)] =nlinfit(Freq(1:end),P_test(b,1:end),model_fun.log_norm,[4.,3],opts);

        end

        % Plot the results
        subplot(2,2,C*2)
        M = mean(log10(P_test));
        SD = std(log10(P_test));
        plot(log10(Freq),mean(log10(P_test)),'k','linewidth',1.5);hold on;
        %plot((Freq),mean((P_test)),'k','linewidth',1.5);hold on;
        h = fill(log10([Freq flip(Freq)]),[M+SD flip(M-SD)],'k','edgecolor','none');
        set(h,'facealpha',.3)

        %y = (mean(this_p(2,:),2)*Freq+1).^-(mean(this_p(1,:),2)) ;
        y = arrayfun(@(x) model_fun.power_law(this_p(:,x),Freq),1:420,'uni',false);
        y = cat(1,y{:});

        PL = plot(log10(Freq),log10(mean(y)),'--r','linewidth',1.5);hold on;
       % PL = plot((Freq),(mean(y)),'--r','linewidth',1.5);hold on;
        F = [1 2 5 10 20 40 70 100];
        set(gca,'xtick',log10(F),'xticklabel',F,'ytick',-1.5:.5:-.5,'yticklabel',round(10.^(-1.5:.2:-.5),2));
        LG = legend(PL,['\beta = ' num2str(round(mean(this_p(1,:)),2))]);
        set(LG,'box','off')
        
        
        title(Conds{C})
        axis tight
        box off
        set(gca,'fontsize',FS,'TickDir','out','linewidth',1.2,'TickLength',[0.02 0.035],'XColor',[.2 .2 .2],'YColor',[.2 .2 .2])
        if C==2
            xlabel('Frequency (Hz)')
            ylabel('Normalized Loading')
            set(gca,'position',get(gca,'position')+[0 .04 0 0])
        end
   
    end
end
% subplot(2,1,2)
% plot((Freq),mean((P_test)));hold on;
% plot((Freq),y,'r');hold on;
export_fig(FIG,fullfile(FigPath,[FileName{1} '_PARAFAC_N' num2str(NComp) '_Bootstrap_BetaExponentVariance']),'-pdf','-r200')
% Now let's fit an exponential


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