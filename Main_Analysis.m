clear;
clc;

%% Add required repos and set the project path, where to read the data
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/'))
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'))
addpath(genpath(pwd));

ProjectPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results';

% Alan atlas
addpath(genpath('/Users/elhamb/Documents/Codes/Git/AllenMouseAtlas'));
%% Load the data
Animal_List = {'766640955','767871931','768515987','771160300','771990200','774875821','778240327','778998620','779839471','781842082','821695405'};
    
LFPF.STOK_analysis(ProjectPath,'drifting_gratings_75_repeats',...
    'Sessions_subset'   ,Animal_List,...
    'ReReadData'        ,true,...% not necessary to set to true everytime you run the code
    'ParamEstimate'     ,false,...% not really necessary -> not efficient
    'MOrd'              ,15,...
    'ff'                ,.98,...
    'Freqs'             ,1:100,...
    'ROIs'              ,{'VISp','VISl','VISrl','VISli','VISal','VISpm','VISam'},...%,{'LGd','VISp','VISl','VISrl','VISli','LP','VISal','VISpm','VISam','VISmma'},...
    'PDCMethod'         ,'iPDC');

%% 



