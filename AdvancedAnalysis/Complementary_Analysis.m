clear;clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_ALL_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
FigPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/ComplementaryResults';

%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
[PDC,IDs,Time,Freq,ROIs] = ExtractAllRoiPDC(StokALL,ROIs);
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;
clear StokALL;
NumSub = 11;
%% first prepare the relative hierarchy scores (RHC)
load Hierarchyscores; % Attention to use the correct order of ROIs for hierarchy scores
HS = mean(H(Time>-.3 & Time<0,:));
RHC = HS'-HS;

%% SNR scores and its relashiionship with output iPDCs
% ideas: first to regress frequencies of output PDC with signal SNR, then
% we can correlate PDC with multiplication of SNR of target and source
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/Signal_Averaged.mat');
Signal_Averaged.('VISrl').SNR(11)=0;
SNR = arrayfun(@(x) Signal_Averaged.(ROIs{x}).SNR,1:numel(ROIs),'uni',false);
SNR = cat(1,SNR{:});
SNR(SNR==0) = nan;

% Output PDC ------> with and without intra PDCS...
PDC_prestim_all = nan(nROIs,nROIs,numel(Freq),NumSub);
relSNR          = nan(nROIs,nROIs,NumSub);
relHC           = nan(nROIs,nROIs,NumSub);
for roi1 = 1:nROIs
    for roi2 = 1:nROIs
        if roi1~=roi2
            SIDs    = ~strcmpi(IDs{roi1,roi2},'NaN');
            relSNR(roi1,roi2,:)  = SNR(roi1,:); % Relative SNR
            relHC(roi1,roi2,SIDs)   = (repmat(RHC(roi1,roi2),1,sum(SIDs)));
            PDC_temp    = squeeze(mean(mean(PDC{roi1,roi2},1),2));
            PDC_prestim = squeeze(mean(PDC_temp(:,Time<0 & Time<-.3,:),2));
            PDC_prestim_all(roi1,roi2,:,SIDs) = PDC_prestim;
            
            PDC_evoked  = (PDC_temp - (mean(PDC_temp(:,Time<0 & Time<-.3,:),2)))./(mean(PDC_temp(:,Time<0 & Time<-.3,:),2));
            PDC_evoked_all(roi1,roi2,:,SIDs) = squeeze(mean(PDC_evoked(:,Time>.05 & Time<1,:),2));
        end
    end
end

PDC_prestim_all2 = reshape(permute(PDC_evoked_all,[1 2 4 3]),nROIs*nROIs*NumSub,numel(Freq));
Regresor = relSNR(:);
KeepInd = ~isnan(Regresor);

[Regresor, Ind] = sort(Regresor(KeepInd));
PDC_prestim_all2 = PDC_prestim_all2(KeepInd,:);

subplot(2,1,1);
plot(mean(PDC_prestim_all2((Ind),30:100),2));xlim([1 numel(Regresor)])
%imagesc(PDC_prestim_all2((Ind),:)');axis xy;
subplot(2,1,2); plot(Regresor,'linewidth',2,'color','k');xlim([1 numel(Regresor)])


