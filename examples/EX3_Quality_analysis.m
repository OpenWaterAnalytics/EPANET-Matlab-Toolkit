%% Quality analysis example
% This example contains: 
% Load a network
% Compute Quality without MSX
% Load EPANET-MSX files
% Compute Quality with MSX (specify type)
% Get quality of specific nodes
% Get quality of specific links
% Get species names
% Get quality for specific species type (nodes and links)
% Plot concentration for specific node and all species 
% Plot concentration for specific link and all species 

%%
clear; close('all'); clc;
start_toolkit;

%% Run Water Quality analysis of a network

% Load a network
d = epanet('net2-cl2.inp');

% Compute Quality without MSX
% (This function contains events)
qual_res = d.getComputedQualityTimeSeries %Value x Node, Value x Link

% Compute Quality step by step
% d.solveCompleteHydraulics #needed
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T=[];Q=[];  
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    P=[P; d.getNodePressure];
    Q=[Q; d.getNodeActualQuality];
    T=[T; t];
    tleft = d.stepQualityAnalysisTimeLeft;
end
d.closeQualityAnalysis;

% Load EPANET-MSX files
d.loadMSXFile('net2-cl2.msx')

% Compute Quality with MSX (specify type)
qual_res_MSX = d.getMSXComputedQualitySpecie('CL2')

% Get quality of specific nodes
sensor_index = [2, 3, 5];
sensors_names = d.getNodeNameID(sensor_index)
QN = d.getMSXComputedQualityNode(sensor_index)

% Get quality of specific links
QL = d.getMSXComputedQualityLink(sensor_index)

% Get species names
type = d.getMSXSpeciesNameID;

% Get quality for specific species type (nodes and links)
for species = type
    MSX_comp = d.getMSXComputedQualitySpecie(species{1})
end

% Plot concentration for specific node and all species 
d.plotMSXSpeciesNodeConcentration(1,1:d.MSXSpeciesCount)

% Plot concentration for specific link and all species 
d.plotMSXSpeciesLinkConcentration(1,1:d.MSXSpeciesCount)

% Unload libraries
d.unloadMSX;
d.unload;