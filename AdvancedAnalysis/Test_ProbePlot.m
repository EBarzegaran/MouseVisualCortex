clear; clc;

load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/Probe_Data_All.mat')
addpath(genpath('/Users/elhamb/Documents/Codes/Git/AllenMouseAtlas'));

%%
ROIs    = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};

RCs     = LFPF.RColors();
Colors  = RCs.MatrixColors(ROIs);

%plotBrainGrid_adjusted(); axis on
FIG = figure;
set(FIG,'unit','inch','color','w','position',[5 5 4 4]);
hold on;
S = arrayfun(@(i) strcmp(Probe_all.Session_ID(i,:),'766640955'),1:size(Probe_all,1));
%Probe_all = Probe_all(S,:);
for roi = 1:numel(ROIs)
    figure(1);
    Ind = cellfun(@(x) ~isempty(x),strfind(Probe_all.Labels,[ROIs{roi} '_']));
    scatter3(Probe_all.AP_CCF(Ind)/10,Probe_all.ML_CCF(Ind)/10,Probe_all.DV_CCF(Ind)/10,10*2,Colors(roi,:),'filled');
%     figure(2);hold on;
%     scatter(Probe_all.AP_CCF(Ind)/1000,Probe_all.ML_CCF(Ind)/1000,20,Colors(roi,:),'filled');
    ylim([6.5 10.5]*100);
    xlim([6 11]*100)
    view(-90,-90)
    set(gca,'ytick',(6.5:1:10.5)*100,'yticklabel',(6.5:1:10.5),'xtick',(6:1:11)*100,'xticklabel',(6:1:11))
    grid on;
    xlabel('CCFv3 M\L coordinate (mm)');
    ylabel('CCFv3 A\P coordinate (mm)');
    set(gca,'fontsize',12)
end

export_fig(FIG,'Probes_plotss','-pdf');
print(FIG,'Probes_plotss','-dtiff')