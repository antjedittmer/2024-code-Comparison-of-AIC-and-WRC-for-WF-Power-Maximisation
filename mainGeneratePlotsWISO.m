clc; clear; close all;

% Add necessary paths
addpath('WFSimCode');
addpath('KoopmanIODMD');

plotAll = 0;

%% Generate plots for the analytical calculations, for FLORIS and FASTFarm
generatePlots_AICWRC_WISO(plotAll);

%% Generate plots for WFSim
generatePlots_WFSim(plotAll);