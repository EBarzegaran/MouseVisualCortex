function VisualizeProbes(ProjectPath,Sessions_subset)

%%
Sessions_ID     = subfolders(ProjectPath,0); % find all the available sessions

if ~isempty(Sessions_subset)
    Sessions_ID     = intersect(Sessions_ID,Sessions_subset); % you can select a subset of the sessions
end

Savepath = fullfile(ProjectPath,'Averaged');
if ~exist(Savepath)
    mkdir(Savepath);
end

%%
tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
st = loadStructureTree('structure_tree_safe_2017.csv'); % a table of what all the labels mean


figure;
plotBrainGrid_adjusted();hold on;

%% for each session and stimuli load the probes data
for S = 1:numel(Sessions_ID)
    
    disp(['Loading Probe locations from #' Sessions_ID{S}])
    % read probe information and save them in session information
    Probes_info= subfiles(fullfile(ProjectPath,Sessions_ID{S},'MatlabData','*_ProbeInfo.mat'),0);
    Probes_ID   = cellfun(@(x) x(1:subsref(strfind(x,'_'),struct('type','()','subs',{{1}}))-1),Probes_info,'uni',false);
    Probes_ID   = unique(Probes_ID);
    session.(['S' Sessions_ID{S}]).Probes_ID = Probes_ID;
    for p = 1:numel(Probes_ID)
        load(fullfile(ProjectPath,Sessions_ID{S},'MatlabData',[Probes_ID{p} '_ProbeInfo.mat']));
        CO = [Coords.AP_CCF; Coords.DV_CCF; Coords.ML_CCF];

        [IDS, ind] = sort(Coords.id);

        CO = CO(:,ind);
        scatter3(CO(1,:)/10,CO(3,:)/10,CO(2,:)/10,10,'filled');

%         AN = arrayfun(@(x) av(round(CO(1,x)/10),round(CO(2,x)/10),round(CO(3,x)/10)),1:size(CO,2));
% 
%         for x = 1:numel(AN)
%             name{x} = st.safe_name(st.sphinx_id==AN(x));
%             ac{x}   = st.acronym(st.sphinx_id==AN(x));
%         end
%         Summary = [name' ac'];

        
        
        
    end
end



        
end