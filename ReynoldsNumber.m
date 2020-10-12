% MAIN_ReynoldsNumber.m
% Purpose: Compute the Reynolds number for this work and compare to
% literature.
% Author: Emma Benjaminson

%% 

clc ; clear all ; close all ; 

%% Directories

localDir = fileparts(mfilename('fullpath')) ;
restoredefaultpath ;
addpath(fullfile(localDir, 'active_functions')) ;

%% Current Work

rho = 1260 ; % [kg/m3] density of glycerin
v = 2.9E-04 ; % [m/s] velocity of swimmer
D = 0.010 ; % [m] length characteristic of swimmer
mu = 0.950 ; % [N s / m2] dynamic viscosity of glycerin

Re = ReNumber(rho, v, D, mu)

%% Steager et al

rho_s = 1000 ; % [kg/m3] density of water
v_s = 10E-06 ; % [m/s] velocity of robot
D_s = 30E-06 ; % [m] length characteristic of robot
mu_s = 0.00089 ; % [N s / m2] dynamic viscosity of water

Re_steager = ReNumber(rho_s, v_s, D_s, mu_s)

%% Functions

function Re = ReNumber(rho, v, D, mu)
    Re = (rho * v * D) / mu ;
end