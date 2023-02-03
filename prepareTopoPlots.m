function [] = prepareTopoPlots()

% Load montage and related stuff 

% Add EEG Lab to path 
addpath(genpath('E:\MSU PhD\MATLAB\Toolboxes\eeglab13_6_5b'))

% Electrode names for montage 
elecn=strvcat('FP1','FPz','FP2','AF3','AF4','F7','F5','F3','F1','Fz','F2','F4','F6','F8',...
    'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8',...
    'T7','C5','C3','C1','Cz','C2','C4','C6','T8',...
    'TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8',...
    'P7','P5','P3','P1','Pz','P2','P4','P6','P8',...
    'PO7','PO3','POz','PO4','PO8',...
    'O1','Oz','O2');

% Extract montage 
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',cellstr(elecn));
[G,H] = GetGH(M);

electrodes = cellstr(elecn);
% Convert locations 
ConvertLocations ( '10-5-System_Mastoids_EGI129.csd', '10-5-System_Mastoids_EGI129.locs', electrodes);
