%% Closing pipes during simulation
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Link index for change the status
link_index=2;

Status = [0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 1 1]';

% Run step by step hydraulic analysis
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
i=1;tstep=1; F=[];
while (tstep>0)
    t=d.runHydraulicAnalysis;
    d.setLinkStatus(link_index, Status(i)); i=i+1;
    F = [F; d.getLinkFlows];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis;

% Get flows for the specific link index
Flows = F(:, link_index);
T = table(Flows, Status);

fprintf('\nFlows and status for node index 2:\n')
T

% Unload library
d.unload
