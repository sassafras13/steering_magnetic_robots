% MAIN_02_ExpDataProcessing.m
% Purpose: Process the experimental data and save the average values and
% standard deviations to a csv file for future use.
% Author: Emma Benjaminson

%% Motion Primitive A
fignum = 1 ; 
directory = fullfile(localDir, '/experimental_data') ; 
filename = 'mpA' ; 

[indices, origins] = extractIndOr(directory,filename) ; 

filename = 'tp12_mpA_iter' ; 
outputname = 'mpA-CompositeData' ; 

[meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
    directory, filename, outputname, fignum) ; 


%% Motion Primitive B
% fignum = 2 ; 
% directory = fullfile(localDir, '/experimental_data') ; 
% filename = 'mpB' ; 
% 
% [indices, origins] = extractIndOr(directory,filename) ; 
% 
% filename = 'tp12_mpB_iter' ; 
% outputname = 'mpB-CompositeData' ; 
% 
% [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
%     directory, filename, outputname, fignum) ; 


%% Case 1
% fignum = 3 ; 
% directory = fullfile(localDir, '/experimental_data') ; 
% filename = 'c1' ; 
% 
% [indices, origins] = extractIndOr(directory,filename) ; 
% 
% filename = 'tp12_c1_v' ; 
% outputname = 'Case1-CompositeData' ; 
% 
% [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
%     directory, filename, outputname, fignum) ; 

