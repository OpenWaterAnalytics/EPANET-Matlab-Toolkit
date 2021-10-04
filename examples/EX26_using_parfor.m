% Using PARFOR
% More fast function with parfor is coming..

tic
start_toolkit;
try
    unloadlibrary('epanet2')
catch
end

d = epanet('Net1.inp');%, 'loadfile'); 

clear H;clc;

number_scenarios = 100;

parfor i = 1:number_scenarios
    loadlibrary('epanet2', 'epanet2.h')
    d.loadEPANETFile(d.TempInpFile);
    
    % set parameters
    elevations = d.getNodeElevations-d.getNodeElevations*rand(1)*.5;
    d.setNodeElevations(elevations);
    % 
    
    % Computed Hydraulics
    H{i} = d.getComputedHydraulicTimeSeries;
    d.closeNetwork;
end
d.unload;
toc


