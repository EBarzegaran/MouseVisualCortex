% This one requires N-way toolbox

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
%addpath(genpath('E:\Elham\Codes\NonGit\nway331'));% add the nway toolbox
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nway331'));% add the nway toolbox
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
%Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
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
nboots = 500;
bootsize = 8;

redoanalysis = false;
NComp = 4;

if redoanalysis || ~exist([FileName 'PARAFAC_covtemp_' num2str(NComp) '.mat'],'file')
    BootIDs = BootAllROIs(IDs, ROIsPerID, ROIs, nboots, bootsize);
    for NComp = 2:6

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

            TW = Time>-0.5 & Time<1;
            %[ssX,Corco,It] = pftest(5,permute(PDCPulled(:,:,:,TW,:),[4 1 2 3 5]),4,[2 2 2 0 0]);
            [model_temp,it(b),err(b),corcondia(b)]=parafac(permute(PDCPulled(:,:,:,TW,:),[4 1 2 3 5]),NComp,[1e-7 10 0 0 NaN],repmat(2,1,NComp));% dimensions: connections x in x out x freq x time
            model{b} = model_temp([5 2 3 4 1]);
            temp_time = Time(TW);
        end
        save([FileName 'PARAFAC_covtemp_' num2str(NComp)],'model','NComp','temp_time','BootIDs','nboots','ROIs','IDs','Freq','indTotal','it','err','corcondia')
    end
else
    load([FileName 'PARAFAC_covtemp_' num2str(NComp)])
end
%% cluster the components  500 x 500 similarity
% for example here I calculate correlation for the first mode : this does
% not work

% We have different modes of data ...
for i = 1:nboots
    i
    for j = 1:nboots
        [maps(i,j,:),sscore(i,j),Corrs(i,j,:),~,CorrPvals(i,j,:,:),Corrvals(i,j,:,:)] = PARAFACmodelComp(model{i},model{j});
    end
end

[~,RefInd] = max(mean(sscore)); % find the reference bootstrap
Corrsref = squeeze(Corrs(RefInd,:,:));

% reordering the models according to the reference component
for i = 1:nboots
    [~,sscore_reord(i),Corrs_reord(i,:),Model_reord{i}] = PARAFACmodelComp(model{RefInd},model{i});
end

%% how much of variance are explained by each component
%[1 4 3 2 5];

clear Model_vars;
for y = 1:NComp
    Model_vars{y} = cellfun(@(x) norm(x{5}(:,y)),Model_reord);
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
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_VarianceExplained']),'-pdf','-r200')

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

Comp_PerSig = squeeze(sum(reshape(CorrPvals_reord2,nboots*nboots,NComp,5)<.05)/(nboots.^2));


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
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'Consistency']),'-pdf','-r200')


%% New figure
% temporal dynamics: significance
pre_win = find(temp_time<-0.0 & temp_time>-0.5);
M_temp = cellfun(@(x) (x{5}),Model_reord,'uni',false);
M_temp = permute(cat(3,M_temp{:}),[3 2 1]); 
% Bastos based hierarchy
clear Hscore;
Freq = 1:100;
model_reord =   cat(1,Model_reord{:});
con_mode = 1;
% plotting params
offs = -.08;
FS = 14;
for id = RefInd%numel(IDs)
    FIG = figure(1);
    set(FIG,'unit','inch','position',[0 0 2.3*NComp 10],'color','w')
    
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
        MS = subplot(5,NComp*2,(2*c-1)); hold on;
        bar(model_temp_M{3}(:,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 .9]); %colormap(MS,jmaColors('coolhot'))
        errorbar(1:6,model_temp_M{3}(:,cc),model_temp_S{3}(:,cc),'.','color','k')
        view([-90 -90])
        set(MS,'position',get(MS,'position')+[offs+(.017*c) .0 .0 .0],'xtick',[],'fontsize',FS);
        ylim([0 1])
        title(['Source']);
        box on;
        if cc==1, xlabel('Layers');end        

        %--------------------------PLOT TARGET-----------------------------
        MT = subplot(5,NComp*2,(2*c));hold on;
        bar(model_temp_M{2}(:,cc),.6,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 .9]); colormap(MS,jmaColors('coolhot'))
        errorbar(1:6,model_temp_M{2}(:,cc),model_temp_S{2}(:,cc),'.','color','k')
        view([90 90])
        set(MT,'position',get(MT,'position')+[offs+(.017*c) .0 .0 .0],'xtick',1:6,'fontsize',FS);
        ylim([0 1])
        title('Target')
        if cc ==1
            XL = ylabel('Strength');
            set(XL,'position',get(XL,'position')-[0.12 .5 0])
        end
        box on;
        
        annotation('arrow',[(cc-1)*.217+.14 (cc-1)*.217+.14+.03], [0.937 0.937],'color','b','linewidth',1.5 )

        %------------------------PLOT TIME MODE----------------------------
        temp_time_ms = temp_time*1000;
        M = model_temp_M{5}(:,cc); %M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:))*100;
        TP = subplot(5,NComp,(c-1)+NComp*2+1); hold on;
        h = fill([temp_time_ms flip(temp_time_ms)],[M+model_temp_S{5}(:,cc); flip(M-model_temp_S{5}(:,cc))]','k','edgecolor','none');
        set(h,'facealpha',.5)
        plot(temp_time_ms,M,'linewidth',1.5,'color','k');
        ylim([-40 107])
        xlim([-.2 1]*1000);
        vline(0,'k--');
        hline(0,'c--');
        if c==1
            xlabel('Time (msec)');
            %ylabel('% change');
        end
        set(TP,'position',get(TP,'position')+[offs+(.017*c) 0 .0 -.02],'fontsize',FS);
        if cc~=1
            set(TP,'ytick',[]);
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
        SigRes1 = PVal(:,cc)<.01;
        CC = bwconncomp(SigRes1);
        for i = 1:numel(CC.PixelIdxList)
            if numel(CC.PixelIdxList{i})>10
                SigRes = CC.PixelIdxList{i};
                SIGF = fill([temp_time(SigRes) flip(temp_time(SigRes))]*1000,[ones(1,numel(SigRes))*-40 ones(1,numel(SigRes))*100],'r','linestyle','none');
                set(SIGF,'facealpha',.1)
                %plot(temp_time(SigRes)*1000,repmat(maxY*1.1,[1 size(SigRes,1)]),'.','color','r','linewidth',2)
            end
        end
        %-------------------PLOT FREQUENCY DISTRIBUTION--------------------
        TP = subplot(5,NComp,(c-1)+NComp*1+1); hold on;
        h = fill([Freq flip(Freq)],[model_temp_M{4}(:,cc)+model_temp_S{4}(:,cc); flip(model_temp_M{4}(:,cc)-model_temp_S{4}(:,cc))],'k','edgecolor','none');
        set(h,'facealpha',.5)
        plot(model_temp_M{4}(:,cc),'linewidth',1.5,'color','k')
        if c==1
            xlabel('Frequency (Hz)')
            %ylabel('% change')
        end
        set(TP,'position',get(TP,'position')+[offs+(.017*c) 0 .0 -.02],'fontsize',FS,'ytick',0:.1:.2);
        if cc~=1
            set(TP,'ytick',[]);
        end
        box on
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
                Hscore(:,:,perm) = DAI;
                
            else
                
                Hscore(:,perm)   = (HierarchyExtract(CN));
            
            end
        end
        CNP = subplot(5,NComp,(c-1)+NComp*3+1);
        if true
            boxplot(squeeze((nanmean(Hscore,1)))','colors','k','symbol','')%,'PlotStyle','compact'
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
            %CNP1 = get(CNP,'position');
            set(CNP,'position',get(CNP,'position')+[offs+(.017*c) 0 .0 -.02],'fontsize',FS);
        else
            set(CNP,'position',get(CNP,'position')+[offs+(.017*c)+.02 0 -0.02 -.02],'fontsize',FS);
        end
         %ylim([-.25 .24])
         if cc==1
            
            %ylabel('Hierarchy Score')
         end
        
         % ----------------PLOT CONNECTIVITY PATTERN-----------------------
         
         CNPp = subplot(5,NComp,(c-1)+NComp*4+1);
         CN = zeros(36,1);
         CN(indTotal) = model_temp_M{1}(:,cc);
         CN = reshape(CN,6,6);
         DAI = (CN-CN')./(CN+CN');
         DAI(isnan(DAI))=0;
         DAI(DAI<0)=0;
         for roi = 1:length(DAI)
            DAI(roi,roi)=sum(DAI(:,roi))*2;
            if DAI(roi,roi)==0
                DAI(roi,roi)=eps;
            end
         end
        
        plot_graph(DAI,ROISN,Colors,[],0);
        set(CNPp,'position',get(CNPp,'position')+[offs+(.017*c) -.03 .0 0.02],'fontsize',FS);
        if cc==1
        colororder({'k','k'})
        yyaxis left
        ylabel('FeedForward');
        yyaxis right;
        set(gca,'ytick',[])
        ylabel('Feedback');
        end
        
    end

end

export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_Bootstrap_HarrisM']),'-pdf','-r200')

%% convert back the components into data space

for b = 1:2%nboots
    Model_temp  = Model_reord{nboots};
    modes       = [1 4 5];% remove laminar modes
    Model_temp  = Model_temp(modes);
    Sizes       = cellfun(@(x) size(x,1),Model_temp);
    ord = [1 3 4; 3 1 4; 3 4 1];
    for m = 1:numel(Model_temp)
        Model_temp{m} = permute(Model_temp{m},[2 ord(m,:)]);
        if m==1
            MT = Model_temp{m};
        else
            MT = Model_temp{m}.*MT;
        end
    end
    
    % NOW we need to do hierarchy analysis B)
    FF_labels = [1 1 -1 -1];
    for f = 1:size(MT,3)
        for t = 1:size(MT,4)
            M = zeros(NComp,numel(ROIs)^2);
            M(:,indTotal) = (MT(:,:,f,t));%
            M = reshape(M,NComp,numel(ROIs),numel(ROIs));
            for c = 1:NComp
                [H(f,t,:,c,b)] = HierarchyExtract(squeeze(M(c,:,:)),'functional');
            end
        end
    end
    b
end
%%
figure,
for roi = 1:6
    plot(squeeze(mean(mean(H(30:end,:,roi,4,:)),4)),'color',Colors(roi,:));hold on
end
