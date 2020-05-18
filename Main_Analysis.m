% clear;
% clc;

%% Add required repos and set the project path, where to read the data
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/'))
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'))
addpath(genpath(pwd));

ProjectPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results';

% Alan atlas
addpath(genpath('/Users/elhamb/Documents/Codes/Git/AllenMouseAtlas'));
%% Load the data
Animal_List = {'766640955','767871931','768515987','771160300','771990200','774875821','778240327','778998620','779839471','781842082','821695405'};

% in case of grating condition
Condition = 'drifting_gratings_75_repeats';
StimParams.contrast = .1;
%StimParams.orientation = [90 45];

%Condition = 'dot_motion';
%StimParams.Speed = [0.0100 0.0200 0.0400];% 0.0005 0.0010 0.0050 0.0100 0.0200 0.0400

LFPF.STOK_analysis(ProjectPath,Condition,...
    'Sessions_subset'   ,Animal_List,...
    'StimParams'         ,StimParams,...
    'ReReadData'        ,true,...% not necessary to set to true everytime you run the code
    'ParamEstimate'     ,false,...% not really necessary -> not efficient
    'MOrd'              ,15,...
    'ff'                ,.98,...
    'Freqs'             ,1:100,...
    'ROIs'              ,{'VISp','VISl','VISrl','VISal','VISpm','VISam'},...%{'VISp','VISl','VISrl','VISal','VISpm','VISam'},...%
    'PDCMethod'         ,'iPDC');

%% 



