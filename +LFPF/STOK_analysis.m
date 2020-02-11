function STOK_analysis(ProjectPath,Stimuli,varargin)


opt = ParseArgs(varargin,   ... % parse the inputs and indicate default variables
        'MOrd'              ,15,...
        'ff'                ,.98, ...
        'PDCMethod'         ,'iPDC',...
        'ReReadData'        ,false,...
        'ParamEstimate'     ,false,...
        'TimeWin'           ,[],...
        'Sessions_subset'   ,[],...
        'Freqs'             ,1:80,...
        'ROIs'              ,{'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam','VISmma','LP','LGd'}...
        );

%%
st = [];av=[];

Sessions_ID     = subfolders(ProjectPath,0); % find all the available sessions

if ~isempty(opt.Sessions_subset)
    Sessions_ID     = intersect(Sessions_ID,opt.Sessions_subset); % you can select a subset of the sessions
end

Savepath = fullfile(ProjectPath,'Averaged');
if ~exist(Savepath)
    mkdir(Savepath);
end

%% for each session and stimuli load the probes data
if ~exist(fullfile(ProjectPath,'Averaged','Fullmodel',['STOK_ALL_' opt.PDCMethod '.mat']))
    for S = 1:numel(Sessions_ID)

        disp(['Analyzing data from session #' Sessions_ID{S}])
        if exist(fullfile(ProjectPath,Sessions_ID{S},'MatlabData','session_Prep.mat')) && ~opt.ReReadData
            load(fullfile(ProjectPath,Sessions_ID{S},'MatlabData','session_Prep.mat'));
            session.(['S' Sessions_ID{S}]) = session_cur;
        else
            % read probe information and save them in session information
            Probes_Stim = subfiles(fullfile(ProjectPath,Sessions_ID{S},'MatlabData',['*' Stimuli '250.mat']),0);
            Probes_ID   = cellfun(@(x) x(1:subsref(strfind(x,'_'),struct('type','()','subs',{{1}}))-1),Probes_Stim,'uni',false);
            Probes_ID   = unique(Probes_ID);
            session.(['S' Sessions_ID{S}]).Probes_ID = Probes_ID;
            session.(['S' Sessions_ID{S}]).Stimuli = Stimuli;

            % read data for each probe and estimate optimum model order
            try % load Laminar info first
                    LayerInfo = readtable(fullfile(ProjectPath,Sessions_ID{S},'MatlabData','Laminar_Layer3.xlsx'));
                catch
                    error('Indicate the Laminar layers first...')
            end
            ROIsNames = {};
            for p =1:numel(Probes_ID)
                disp(['optimizing parameters for probe #' Probes_ID{p}]);
                disp('This may take a few minutes...')
                %-------------------------Read Data----------------------------
                % load the LFP data
                Data = load(fullfile(ProjectPath,Sessions_ID{S},'MatlabData',[Probes_ID{p} '_' Stimuli '250.mat']));
                if ~isempty(opt.TimeWin)% select a time window
                    T1          = opt.TimeWin(1);
                    T2          = opt.TimeWin(2);
                    Data.Times   = Data.Times(Data.Times>T1 & Data.Times<T2);
                    Data.Y      = Data.Y(:,:,Data.Times>T1 & Data.Times<T2,:);
                end


                % load the probe coordinate and label data
                probeinfo = load(fullfile(ProjectPath,Sessions_ID{S},'MatlabData',[Probes_ID{p} '_Probeinfo']));
                % use the allen atlas to label the layers and ROIs
                [Summary,av,st] = GetAllenLabels(probeinfo.Coords,av,st);
                % segment the LFP data according to the ROIs and sub(laminar) labels
                session.(['S' Sessions_ID{S}]) = extractROIs(Data,probeinfo,LayerInfo.Layers,LayerInfo.(['P' Probes_ID{p}]),Summary,opt.ROIs,session.(['S' Sessions_ID{S}]),fullfile(Savepath,Sessions_ID{S}));

                %--------------------plot bipolar maps---------------------
                %plot_bipolars(Data,LayerInfo.(['P' Probes_ID{p}]),Probes_ID{p},fullfile(ProjectPath,Sessions_ID{S}));
                %-----------------------estimate AR params---------------------
                %if opt.ParamEstimate
                % [Error(p) session.(['S' Sessions_ID{S}]).(Data.ROI).Y]= LFPF.ParamEstimation(Data,LayerInfo.(['P' Probes_ID{p}]),opt.Freqs);                
                %else
                %    session.(['S' Sessions_ID{S}]).(Data.ROI).Y = -1*Data.Y(:,LayerInfo.(['P' Probes_ID{p}]),:,Data.cnd_info.contrast==.8); 
                %end

            end

            % session.(['S' Sessions_ID{S}]).ROIs = ROIsNames;
            if opt.ParamEstimate
                % check the MSE for all probes and choose a general model order
                Res = squeeze(cat(3,Error.R));
                [~,MI] = max(mean(Res,2));
                session.(['S' Sessions_ID{S}]).Mord = Error(1).pp(MI); % optimum model order according to MSE of spectrum
            end

            % save the result of analysis, so not evey time you need to estimate
            % model orders
            session_cur = session.(['S' Sessions_ID{S}]);
            save(fullfile(ProjectPath,Sessions_ID{S},'MatlabData','session_Prep.mat'),'session_cur')
        end

        %% Run STOK analysis with indicated model order, first on single ROIs
    %     ROIs = session.(['S' Sessions_ID{S}]).ROIs;
    %     for roi = 1:numel(ROIs)
    %         [Temp.PDC,Temp.f,Temp.Times] = LFPF.STOKEstimate(session.(['S' Sessions_ID{S}]).(ROIs{roi}),opt.MOrd, opt.ff,opt.PDCMethod,opt.Freqs);
    %         Temp.ROI = ROIs{roi};
    %         Stok.(['S' Sessions_ID{S}]).(Temp.ROI) = Temp;
    %     end

        %% Estimate STOK over all ROIs
         tic
              [Temp.PDC,Temp.f,Temp.Times,Temp.ROIs] = LFPF.STOKEstimate_All(session.(['S' Sessions_ID{S}]),opt.MOrd, opt.ff,opt.PDCMethod,opt.ROIs,opt.Freqs);
              StokALL.(['S' Sessions_ID{S}]) = Temp;
         toc

    end

    save(fullfile(ProjectPath,'Averaged','Fullmodel',['STOK_ALL_' opt.PDCMethod]),'StokALL','-v7.3');
else
    load(fullfile(ProjectPath,'Averaged','Fullmodel',['STOK_ALL_' opt.PDCMethod]));
end

savefig = true;

%% plot the average

% STOK_Averaged = STOKROIAverage(Stok,[],savefig,Savepath,opt.PDCMethod);
SpecEstim =false;
SignalROIAverage(session,opt.ROIs,savefig,Savepath,SpecEstim,opt.Freqs);    

%% average the whole visual cortex networks

STOKAllAverage(StokALL,opt.ROIs,savefig,Savepath,opt.PDCMethod);
STOKAllLayerStats(StokALL,opt.ROIs,savefig,Savepath,opt.PDCMethod);
STOKAllDirectStats(StokALL,opt.ROIs,savefig,Savepath,opt.PDCMethod);

end
    
