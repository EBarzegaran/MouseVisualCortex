clear; clc;

load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/Probe_Data_All.mat')
addpath(genpath('/Users/elhamb/Documents/Codes/Git/AllenMouseAtlas'));

%%
ROIs    = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};

RCs     = LFPF.RColors();
Colors  = RCs.MatrixColors(ROIs);

plotBrainGrid_adjusted(); axis on

hold on;
S = arrayfun(@(i) strcmp(Probe_all.Session_ID(i,:),'766640955'),1:size(Probe_all,1));
%Probe_all = Probe_all(S,:);
for roi = 1:numel(ROIs)
    Ind = cellfun(@(x) ~isempty(x),strfind(Probe_all.Labels,[ROIs{roi} '_']));
    scatter3(Probe_all.AP_CCF(Ind)/10,Probe_all.ML_CCF(Ind)/10,Probe_all.DV_CCF(Ind)/10,[],Colors(roi,:),'filled');
end

% for layer = 1:6
%     Ind = cellfun(@(x) ~isempty(x),strfind(Probe_all.Labels,['_L' num2str(layer)]));
%     scatter3(Probe_all.AP_CCF(Ind)/10,Probe_all.ML_CCF(Ind)/10,Probe_all.DV_CCF(Ind)/10,Probe_all.RF_Area(Ind)*2,Colors(layer,:),'filled');
% end

