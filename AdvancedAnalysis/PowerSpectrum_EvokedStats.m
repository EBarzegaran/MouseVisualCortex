%% check the power spectrum maybe?
%Question why nlinfit does not work?
FileNames2 = {'Signal_PSD_drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15.mat',...
    'Signal_PSD_drifting_gratings_75_repeats__contrast0-1_iPDC_Mord15.mat'};
Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';

%%
close all
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))))
%% Organize PDC values for statistics
ROIs = {'VISp','VISl','VISrl','VISal','VISpm','VISam'};
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
nROIs = numel(ROIs);
load ROInames;

%%
TW = [-.3 0; .05 .5];
Colors = {'r','b'};
FIG = figure;
set(FIG,'unit','inch','color','w','position',[3 0 8 20 ]);
for i = 1:2
    PSD = load(fullfile(Path,FileNames2{i}));
    stokpsd = PSD.STOKPSD;
    for roi = 1:numel(ROIs)
        PSD_temp = stokpsd.(ROIs{roi});
        PSD_temp = PSD_temp(~cellfun(@isempty,PSD_temp));
        PSD_temp = cellfun(@(x) permute(x,[4 5 2 3 1]),PSD_temp,'uni',false);
        [Stats,~] = CompPostPre(PSD_temp,[],[],0,0,PSD.t,PSD.Freqs,[],[]);
        subplot(6,2,(roi-1)*2+(i))
        if i==2 && roi==1
            pcolorb=true;% plot color bar
        else
            pcolorb=false;% plot color bar   
        end
        
        if i==1 && roi==6
            xylabels = true;% x and y labels
        else
            xylabels = false;
        end
        
        FIG = PlotStatresults(Stats.P,Stats.Tval,[1 2], PSD.t, PSD.Freqs,...
                'figtitle'  ,ROI_names.(ROIs{roi}),...
                'newfig'    ,false,...
                'PThresh'   ,.01,... %benferoni correction!
                'SThresh'   ,15,...
                'Twin'      ,[-.0 1],...
                'ColorM'    ,'jet',...
                'colorb'    ,pcolorb,...
                'xylabels'  ,xylabels);
    end
end

axes('position',[.3 .96 .1 .05]);
axis off
text(0,0,'High Contrast','horizontalalignment','center','fontsize',12,'fontweight','bold')

axes('position',[.75 .96 .1 .05]);
axis off
text(0,0,'Low Contrast','horizontalalignment','center','fontsize',12,'fontweight','bold')

export_fig(fullfile(Path,['PSD_stats_Gratings']),'-pdf');