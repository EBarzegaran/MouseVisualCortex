clear;clc;

addpath(genpath(fileparts(mfilename('fullpath'))));
load('/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/STOK_ALL_iPDC.mat');
SavePath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel';
FigPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/StatResults';
%% Organize PDC values for statistics
[PDC,Time,Freq,ROIs] = ExtractAllRoiPDC(StokALL);
COBJ = LFPF.RColors();
Colors = COBJ.MatrixColors(ROIs);
clear StokALL;
%%
CellNums = cellfun(@(x) size(x,5),PDC);
CellNums(CellNums==1)=0;

    
%% To indicate which part of the evoked PDC is significant
for roi = 2:7
    [~,FIG1] = CompPostPre(PDC(roi,roi),[1],[],1,1,Time,Freq,['Output-intraPDC- ' ROIs{roi}],FigPath);
    [~,FIG2] = CompPostPre(PDC(setdiff(1:7,roi),roi),[1],[],1,1,Time,Freq,['Output-interPDC- ' ROIs{roi}],FigPath);
end

%% To compare the evoked PDC in different connections

% First, the PDCs should be converted to % change
PDC_change = cellfun(@(x) (x - mean(x(:,:,:,Time<0 & Time>-.3,:),4))./(mean(x(:,:,:,Time<0 & Time>-.3,:),4)),PDC,'uni',false);
% layerwise comparison -> intra vs. inter
for roi = 1:7
    CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:7,roi),1),[1],false,1,[1 0],Time,Freq,['Output-(Intra vs. Inter) - ' ROIs{roi}],FigPath);
end

% average comparison -> intra vs. inter
for roi = 1:7
    CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:7,roi),1),[1 2],false,1,[1 0],Time,Freq,['Output-Averaged-(Intra vs. Inter)- ' ROIs{roi}],FigPath);
end

% all layer comparison -> intra vs. inter
CompPDCCells(PDC_change(1,1),PDC_change(2:end,1),[],false,1,[1 0],Time,Freq,['All-InterIntra - ' ROIs{1}]);

% average comparison -> FF vs. FB
PDC_change2 = PDC_change;
PDC_change2 = PDC_change2([1:4 6 5 7],[1:4 6 5 7]);
ind = 1;
for roi1 = 1:numel(ROIs)
    for roi2 = roi1+1:numel(ROIs)
        %if roi1~=4 && roi2~=4
            PDCFF{ind} = PDC_change2{roi2,roi1};
            PDCFB{ind} = PDC_change2{roi1,roi2};
            ind = ind+1;
        %end
    end
end
CompPDCCells(PDCFF,PDCFB,[],[],1);

%% In this part, I use two-way anova to compare inter-intra , pre-post PDC values

for roi = 1:7
    CompPDCCells(PDC(roi,roi),PDC(setdiff(1:7,roi),1),[1 2],true,1,[1 0],Time,Freq,['Output-Averaged-(Intra vs. Inter)-ANOVA - ' ROIs{roi}],FigPath);
end

% roi = 1;
% CompPDCCells(PDC_change(roi,roi),PDC_change(setdiff(1:7,roi),1),[1 2],true,1,[1 0],Time,Freq,['Output-(Intra vs. Inter) - ' ROIs{roi}],FigPath);
% 

