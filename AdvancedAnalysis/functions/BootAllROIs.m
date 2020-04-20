function [BootIDs, ROIs_Select] = BootAllROIs(IDs, ROIsPerID, ROIs_Select, nboots, bootsize)

%% Organize ROIs 
ROIind = zeros(numel(IDs),numel(ROIs_Select));
for id = 1:numel(IDs)
    [~,ia] = intersect(ROIs_Select,ROIsPerID{id});
    ROIind(id,ia)=1;
end

%% select bootstrap from the animals so that it contains all ROIs

for b = 1:nboots
    sample_temp = randi(numel(IDs),1,bootsize);
    while sum(sum(ROIind(sample_temp,:),1)>0)==numel(IDs)
        sample_temp = randi(numel(IDs),1,bootsize);
    end
    boot_idx(b,:) = sample_temp;
end

%% convert to IDs
BootIDs = IDs(boot_idx);
end