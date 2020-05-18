clear;clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';

%% load stok data
load(fullfile(Path, ['STOK_ALL_' FileName '.mat']));
SavePath = Path;
FigPath = fullfile(Path,'StatResults');

%% load probe data
load (fullfile(Path,'Probe_Data_All.mat'));

%% plot the centroid and area distributions of the RFs
close all
RFDataplot(StokALL,Probe_all,.05,fullfile(Path,'DistancePlot'));

%% call the function to estimate the distance

StokALL = DistanceEstimate(StokALL,Probe_all);
StokALL = RFDistanceEstimate(StokALL,Probe_all);

%% reshape all the PDCs and distances and then correlate them

DistancePDCcorr(StokALL,1,fullfile(Path,'DistancePlot'),FileName);

RFPDCcorr(StokALL,.05,fullfile(Path,'DistancePlot'),FileName);

%% Part 2: correlate individual connections
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
[PDC,IDs,Time,Freq,ROIs,RFDist,RFPval] = ExtractAllRoiPDC(StokALL,ROIs,true);
%%
for roi1 = 1:numel(ROIs)
    for roi2 = 1:numel(ROIs)
        if roi1==roi2
            diag=true;
        else
            diag = false;
        end
        [R(roi1,roi2,:,:), P(roi1,roi2,:,:)] = ConnectPDCCorr(PDC{roi1,roi2},RFDist{roi1,roi2},RFPval{roi1,roi2},.05,diag,true,Time);
    end
end
%%
% plot the results
FIG = figure;
set(FIG,'unit','inch','position',[1 1 20,15],'color','w')
load ROInames.mat

R_p = R(:,:,:,Time>0);
P_p = P(:,:,:,Time>0);
M = max(R_p(:))*.5;
for roi1 = 1:numel(ROIs)
    for roi2 = 1:numel(ROIs)
        subplot(numel(ROIs),numel(ROIs),(roi1-1)*numel(ROIs)+roi2)
        IM = imagesc(Time(Time>0),[],squeeze(R_p(roi1,roi2,:,:)));
        set(IM,'alphadata',(squeeze(P_p(roi1,roi2,:,:))<0.01)*.8+.2)
        axis xy;
        caxis([-M M])
        colormap('jet')
        if roi2==1
            ylabel(ROI_names.(ROIs{roi1}),'fontweight','bold');
        end
        
        if roi1==1
            title(ROI_names.(ROIs{roi2}));
        end
        xlim([0 1])
        if roi1 == numel(ROIs) && roi2 == numel(ROIs)
            CS = get(gca,'position');
            colorbar
            set(gca,'position',CS)
        end
    end
end


export_fig(FIG,fullfile(Path,'DistancePlot',[FileName '_All_Connections']),'-pdf','-r200')

%%

Coords_all = [Probe_all.AP_CCF Probe_all.ML_CCF Probe_all.DV_CCF];
[RF_PVal  Cnd] = min(Probe_all.RF_PValue,[],2);
RF_cent    = arrayfun(@(x) Probe_all.RF_Centroid(x,:,Cnd(x)),1:330,'uni',false);
RF_cent     = cat(1,RF_cent{:});

RF_Dist = squareform(pdist(RF_cent));
Co_Dist = squareform(pdist(Coords_all));

scatter(Co_Dist(:),RF_Dist(:))
xlim([10000 17000])