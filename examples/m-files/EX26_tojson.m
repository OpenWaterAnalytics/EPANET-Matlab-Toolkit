%% Exports values to json format.
%% This example contains: 
%% 
% * Load a network. 
% * Run complete analysis. 
% * Create json text variables from single or many variables. 
% * Create json files containing specific variables. 
% * Create json files containing all variables from the analysis. 
% * Unload library.
%% Clear - Start Toolkit

clear; close('all'); clc;
start_toolkit;
%% Load a network

d = epanet('Net1.inp');
%% Run complete analysis.

values = d.getComputedTimeSeries;
%% Create json text variables from single or many variables.

demsJsonTxt = d.toJson(values.Demand)
flowJsonTxt = d.toJson(values.Flow)
allJsonTxt = d.toJson(values)
%% Create json files containing all variables from the analysis.

d.toJsonFile(values.Flow, 'Flows'); % input: value, filename (.json)
d.toJsonFile(values, 'AllValues');
%% Unload library.

d.unload;