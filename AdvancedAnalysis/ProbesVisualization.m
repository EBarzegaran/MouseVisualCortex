clear; clc;

clear;
clc;

%% Add required repos and set the project path, where to read the data
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/'))
addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'))
addpath(genpath(pwd));
addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/AllenMouseAtlas'));
ProjectPath = '/Users/elhamb/Documents/Data/AllenBrainObserver/preliminary_results';

%%
Sessions_subset = {'766640955','767871931','768515987','771160300','771990200','774875821','778240327','778998620','779839471','781842082','821695405'};

VisualizeProbes(ProjectPath,Sessions_subset)