function [STOK_avg] = STOKAllAverage(STOKAll,ROIs_Select,Path,SaveName)

% INPUT:
        % STOKALL: is a structure with animal ID as field names and each
            % STOKALL.Session_ID is a structure itself with the following fields: PDC, f,
            % Times, ROIs

        % ROI_Select: A cell array, with name of ROIs (put them in correct order to have the averaged STOk in correct order)
        % savefig: true or false
        % Path: Path to save the file
        % SaveName: name for the file to save the results
        
%% Organize ROIs and their order
Sessions_ID = fieldnames(STOKAll);

ROIs = cellfun(@(x) (STOKAll.(x).ROIs),Sessions_ID,'uni',false);
[ROIs,~,ic] = unique(cat(2,ROIs{:})); % ROI names
ROI_counts = accumarray(ic,1); % ROI counts

[C,ia,ib] = intersect(ROIs_Select,ROIs);
[~,I] = sort(ia);
ROIs = C(I);
ROI_counts = ROI_counts(ib(I));

%% Make the grand average matrix with all ROIS/Full model
Dims = size(STOKAll.(Sessions_ID{1}).PDC);
%
ROIsizes = contains(ROIs_Select,'VIS')*3+3;
ROIindices = [0 cumsum(ROIsizes)];
%
STOK_AV = zeros(sum(ROIsizes),sum(ROIsizes),Dims(3),Dims(4));
ROI_num = zeros(sum(ROIsizes),sum(ROIsizes));
for S = 1:numel(Sessions_ID)
   PDC = STOKAll.(Sessions_ID{S}).PDC;
   rois = STOKAll.(Sessions_ID{S}).ROIs;
   [C,ia] = intersect(ROIs_Select,rois);
   [~,I] = sort(ia);
   rois = C(I);
   Labels = ia(I);
   
   %indices = arrayfun(@(x) [(x-1)*6+1:x*6],Labels,'uni',false);
   %indices = cat(2,indices{:});
   indices  =  arrayfun(@(x) ROIindices(x)+1:ROIindices(x+1),Labels,'uni',false);
   indices = cat(2,indices{:});
   
   STOK_AV(indices,indices,:,:) = STOK_AV(indices,indices,:,:) + PDC;
   ROI_num(indices,indices,:,:) = ROI_num(indices,indices,:,:) + 1;
   
end

STOK_AV = squeeze(STOK_AV)./ROI_num;
%--------------------------------------------------------------------------
saveresults = false;

STOK_avg.PDC    = STOK_AV;
STOK_avg.Time   = STOKAll.(Sessions_ID{1}).Times;
STOK_avg.Freq   = STOKAll.(Sessions_ID{1}).f;
STOK_avg.ROIs   = ROIs_Select;
STOK_avg.ROIs_count = ROI_counts;
labels = arrayfun(@(y) arrayfun(@(x) [ROIs_Select{y} '_L' num2str(x)],1:ROIsizes(y),'uni',false),1:numel(ROIs_Select),'uni',false);
labels = cat(2,labels{:});
STOK_avg.labels  = labels;

if saveresults
    save(fullfile(Path,'Fullmodel',['STOK_Average_' SaveName]),'STOK_avg');
    clear STOK_avg;
end

end


