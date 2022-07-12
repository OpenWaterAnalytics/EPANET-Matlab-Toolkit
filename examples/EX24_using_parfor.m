% Using PARFOR
clc; clear; close all; clear class;
start_toolkit;

tic

d = epanet('Net1.inp'); 
iterations = 100;

parfor i = 1:iterations
% Uncomment section for MATLAB R2020 and previous versions.  
    d.loadlibrary;
    d.loadEPANETFile(d.TempInpFile); 

    % set parameters
    elevations = d.getNodeElevations-d.getNodeElevations*rand(1)*.5;
    d.setNodeElevations(elevations*10);
    
    % Computed Hydraulics
    H{i} = d.getComputedHydraulicTimeSeries;
    Q{i} = d.getComputedQualityTimeSeries;
end
d.unload;
toc

