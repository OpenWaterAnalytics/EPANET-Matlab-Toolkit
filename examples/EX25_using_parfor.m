% Using PARFOR
% More fast function with parfor is coming..

tic
start_toolkit;

d = epanet('Net1.inp'); 

clear H;clc;

number_scenarios = 100;

parfor i = 1:number_scenarios
%   d.loadEPANETFile(d.TempInpFile);
    % set parameters
    elevations = d.getNodeElevations-d.getNodeElevations*rand(1)*.5;
    d.setNodeElevations(elevations*10);
    % 
    
    % Computed Hydraulics
    H{i} = d.getComputedHydraulicTimeSeries;
end
d.unload;
toc


