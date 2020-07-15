
% This one requires N-way toolbox

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath('E:\Elham\Codes\NonGit\nway331'));% add the nway toolbox
clear; clc;
FileName = 'drifting_gratings_75_repeats__contrast0-8_iPDC_Mord15_ff098';
%Path = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results/Averaged/Fullmodel/';
Path = 'E:\Elham\Data\AllenBrain\preliminary_results\Averaged\Fullmodel\';
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
[PDC,~,Time,Freq,ROIs] = ExtractAllRoiPDC(StokALL,ROIs);

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
%[model_temp,it(b),err(b),corcondia(b)]=parafac(permute(PDCPulled(:,:,:,TW,:),[4 1 2 3 5]),NComp,[0 10 0 0 NaN],repmat(2,1,NComp));% dimensions: connections x in x out x freq x time

XvalResult = ncrossdecomp('parafac',permute(PDCPulled(:,:,:,TW,:),[4 1 2 3 5]),2,3,2,1,1);
 
 
 