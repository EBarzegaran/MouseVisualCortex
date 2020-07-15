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
[PDC,IDs,Time,Freq,ROIs] = ExtractAllRoiPDC(StokALL,ROIs);
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;
clear StokALL;
ROISN = cellfun(@(x) ROI_names.(x),ROIs,'uni',false);
%%
addpath(genpath('../../ARC'));

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
%%
% [~,Corco,~] = pftest(1,PDCPulled(:,:,:,Time>-.5 & Time<1,:),5,[0 0 0 0 0],[0 2 2 0 2]);% COMP = 4
%PDCPulled2 = PDCPulled - mean(PDCPulled(:,:,:,Time>-0.5 & Time<0,:),4);%kõ
NComp = 5;
if ~exist(fullfile(Path,[FileName '_PARAFAC_' num2str(NComp) '.mat']),'file')
    
    TW = Time>-0.5 & Time<1;
    [model]=parafac(permute(PDCPulled(:,:,:,TW,:),[5 1 2 3 4]),NComp,[0 0 0 0 0],[2 2 2 2 2]);% dimensions: connections x in x out x freq x time
    temp_time = Time(TW);
    save(fullfile(Path,[FileName '_PARAFAC_' num2str(NComp)]),'model','temp_time','NComp');
else
    load(fullfile(Path,[FileName '_PARAFAC_' num2str(NComp)]));
end

% for contrast .1
Comp_ord = [1 2 5 3 4];
% for contrast .8
Comp_ord = [1 5 3 2 4];

%%
FIG = figure;
set(FIG,'unit','inch','position',[0 0 3.5*NComp 10],'color','w')


for c = 1:NComp
    cc = Comp_ord(c);
    %%PLOT SOURCE    
    MS = subplot(4,NComp*2,(c-1)*2+1);
    imagesc(model{3}(:,cc),[-1 1]); colormap(MS,jmaColors('coolhot'))
    set(MS,'position',get(MS,'position')+[-.098+(.017*c) 0 -.015 0],'xticklabel',[]);
    title('Source')

    % PLOT TARGET
    MT = subplot(4,NComp*2,(c-1)*2+2);
    imagesc(model{2}(:,cc),[-1 1]); colormap(MT,jmaColors('coolhot'))
    set(MT,'position',get(MT,'position')+[-.115+(.017*c) 0 -.015 0],'xticklabel',[],'yticklabel',[]);
    title('Target')

    % PLOT TIME MODE
    M = model{5}(:,cc); M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
    TP = subplot(4,NComp,(c-1)+NComp+1);
    plot(temp_time,M,'linewidth',1.5,'color','k')
    xlim([-.2 1])
    hline(0,'k--')
    if c==1
        xlabel('Time(Sec)')
        ylabel('% change')
    end
    set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 0]);

    % PLOT FREQUENCY DISTRIBUTION
    TP = subplot(4,NComp,(c-1)+NComp*2+1);
    plot(model{4}(:,cc),'linewidth',1.5,'color','k')
    if c==1
        xlabel('Frequency(Hz)')
        %ylabel('% change')
    end
    set(TP,'position',get(TP,'position')+[-.11+(.017*c) 0 .02 0]);

    % PLOT CONNECTIVITY PATTERN
    CN = zeros(36,1);
    CN(indTotal) = model{1}(:,cc);
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


export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) '_ordered']),'-pdf','-r200')
%print(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp)]),'-dtiff')
%% plot the timings and frequencies overlapped

colors = [.4 .1 0; 1 .3 0; .8 0 0; 0 .0 .5; 0 .4 1];
FIG = figure;
set(FIG,'unit','inch','position',[0 0 8 8],'color','w')
subplot(2,1,1); hold on
for i = 1:NComp
     M = model{5}(:,Comp_ord(i)); M = (M-mean(M(temp_time<0,:)))./mean(M(temp_time<0,:));
    plot(temp_time,M,'color',colors(i,:),'linewidth',1.5);
end

xlim([-.2 1])
%ylim([0 1])
set(gca,'xtick',-.2:.1:1)
vline(0,'--k')
grid on;
xlabel('time (S)')
ylabel('%change')

subplot(2,1,2); hold on
for i = 1:NComp
     M = model{4}(:,Comp_ord(i));
     M = (M -min(M))./(max(M) -min(M));
    plot(M,'color',colors(i,:),'linewidth',1.5);
end

set(gca,'xtick',0:20:100,'ytick',[])
grid on;
xlabel('Frequency (Hz)')

legend({'C1','C2','C3','C4','C5'})
export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'TimeFreq_ordered']),'-pdf','-r200')
%% plot graph !
FIG = figure;
set(FIG,'unit','inch','position',[0 0 5*NComp 3.5],'color','w')


for c = 1:NComp
    cc = Comp_ord(c);
    
    % PLOT CONNECTIVITY PATTERN
    CN = zeros(36,1);
    CN(indTotal) = model{1}(:,cc);
    CN = reshape(CN,6,6);
    CNP = subplot(1,NComp,(c));
    CNT = (CN-CN')./(CN+CN');
    CNT(CNT<0)=0;
    %CN = CN./15; %max(CN(:));
    plot_graph(CNT,ROISN,Colors,[],0)
end

export_fig(FIG,fullfile(FigPath,[FileName '_PARAFAC_N' num2str(NComp) 'Graph_ordered']),'-pdf','-r200')

%% compare .1 and .8
NComp = 5;
FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
load(fullfile(Path,[FileName '_PARAFAC_' num2str(NComp)]));
model08 = model;
FileName = 'drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15_ff098';
load(fullfile(Path,[FileName '_PARAFAC_' num2str(NComp)]));
model01 = model;

% for contrast .8
Comp_ord08 = [1 5 3 2 4];

% for contrast .1
Comp_ord01 = [1 2 5 3 4];

Recon_PDC08 = arrayfun(@(c) model08{1}(:,c).*model08{2}(:,c)'.*permute(model08{3}(:,c),[3 2 1]).*permute(model08{4}(:,c),[4 3 2 1]).*permute(model08{5}(:,c),[5 4 3 2 1]),Comp_ord08,'uni',false);
Recon_PDC01 = arrayfun(@(c) model01{1}(:,c).*model01{2}(:,c)'.*permute(model01{3}(:,c),[3 2 1]).*permute(model01{4}(:,c),[4 3 2 1]).*permute(model01{5}(:,c),[5 4 3 2 1]),Comp_ord01,'uni',false);

for c = 1:NComp
    temp_PDC08  = permute(Recon_PDC08{c},[2 3 4 5 1]);
    temp_PDC01  = permute(Recon_PDC01{c},[2 3 4 5 1]);
    CompPDCCells({temp_PDC08},{temp_PDC01},[1 2],false,1,[0 0],temp_time,Freq,['PARAFAC-CompNum-' num2str(c) ],FigPath,false);
end


%% plot reconstructed data
Recon_PDC = arrayfun(@(c) model{1}(:,c).*model{2}(:,c)'.*permute(model{3}(:,c),[3 2 1]).*permute(model{4}(:,c),[4 3 2 1]).*permute(model{5}(:,c),[5 4 3 2 1]),1:NComp,'uni',false);


for c = 1:NComp
    Dims = size(Recon_PDC{c});
    temp = zeros([36,Dims(2:end)]);
    temp(indTotal,:,:,:,:) = Recon_PDC{c};
    temp = squeeze(mean(mean(temp,2),3));
    temp = reshape(temp,6,6,size(temp,2),size(temp,3));
    Temp{c} = temp - mean(temp(:,:,:,temp_time<0),4);
    dynet_connplot(Temp{c},temp_time,1:100,ROISN);
end


%% to confirm that the sum of 

PDC_avg = cellfun(@(x) mean(x(:,:,:,Time>-0.5 & Time<1,:),5),PDC,'uni',false);
for i = 1:6
    for j = 1:6
        PDCP1(i,j,:,:) = squeeze(mean(mean(PDC_avg{i,j},1),2));
    end
end

PDC_recon = cell(6);
for c = 1:NComp
    Dims = size(Recon_PDC{c});
    temp = zeros([36,Dims(2:end)]);
    temp(indTotal,:,:,:,:) = Recon_PDC{c};
    temp = reshape(temp,[6 6 Dims(2:end)]);
    for i = 1:6
        for j = 1:6
            if ~isempty(PDC_recon{i,j})
                PDC_recon{i,j} = PDC_recon{i,j}+squeeze(temp(i,j,:,:,:,:));
            else
                PDC_recon{i,j} = squeeze(temp(i,j,:,:,:,:));
            end
        end
    end
end

for i = 1:6
    for j = 1:6
        PDCP2(i,j,:,:) = squeeze(mean(mean(PDC_recon{i,j},1),2));
    end
end

dynet_connplot(PDCP1,temp_time,1:100,ROISN)

dynet_connplot(PDCP2,temp_time,1:100,ROISN)
