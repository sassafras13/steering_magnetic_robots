% MAIN_ExpDataProcessing.m
% Purpose: Process the experimental data and save the average values and
% standard deviations to a csv file for future use.
% Author: Emma Benjaminson
% References: 

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% MPA

directory = '/home/emma/Documents/Research/TP12/Motion Primitives/For Analysis' ; 
filename = 'mpA' ; 

[indices, origins] = extractIndOr(directory,filename) ; 

filename = 'tp12_mpA_iter' ; 
outputname = 'mpA-CompositeData' ; 

[meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
    directory, filename, outputname) ; 


%% MPB

directory = '/home/emma/Documents/Research/TP12/Motion Primitives/For Analysis' ; 
filename = 'mpB' ; 

[indices, origins] = extractIndOr(directory,filename) ; 

filename = 'tp12_mpB_iter' ; 
outputname = 'mpB-CompositeData' ; 

[meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
    directory, filename, outputname) ; 


%% Case 1

directory = '/home/emma/Documents/Research/TP12/Composites/Case 1/Mod' ; 
filename = 'c1' ; 

[indices, origins] = extractIndOr(directory,filename) ; 

filename = 'tp12_c1_v' ; 
outputname = 'Case1-CompositeData' ; 

[meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins, ...
    directory, filename, outputname) ; 

%% Primitives 1

% directory = '/home/emma/Documents/Research/TP11/Sequence 1/5 sec_frame Primitives/P1' ; 
% 
% [indices, origins] = extractIndOr(directory) ; 
%         
% filename = 's1-p1' ; 
% 
% [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins,...
%     directory, filename, ntraj, wp, 's1-P1Data') ; 

%% Primitives 2
% 
% directory = '/home/emma/Documents/Research/TP11/Sequence 1/5 sec_frame Primitives/P2' ; 
% 
% [indices, origins] = extractIndOr(directory) ; 
%        
% filename = 's1-p2' ; 
% 
% [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins,...
%     directory, filename, ntraj, wp, 's1-P2Data') ; 

%% Primitives 3
% 
% directory = '/home/emma/Documents/Research/TP11/Sequence 1/5 sec_frame Primitives/P3' ; 
% 
% [indices, origins] = extractIndOr(directory) ; 
%        
% filename = 's1-p3' ; 
% 
% [meanposx,meanposy,stdposx,stdposy] = processExpData(indices, origins,...
%     directory, filename, ntraj, wp, 's1-P3Data') ; 
