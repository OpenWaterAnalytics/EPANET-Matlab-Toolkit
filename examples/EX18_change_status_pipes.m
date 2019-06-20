%% Change status randomly at pipes during simulation
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Get pipe count
pipe_count = d.getLinkPipeCount;
% Get pipe indices
pipe_indices = d.getLinkPipeIndex;

% Run step by step hydraulic analysis
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
i=1;tstep=1; F=[];
while (tstep>0)
    % Set status random 0/1 for pipes
    Status = round(rand(1,pipe_count))';

    t=d.runHydraulicAnalysis;
    d.setLinkStatus(pipe_indices, Status); i=i+1;
    F = [F; d.getLinkFlows];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis;

% Plot flows
plot(F)


% Unload library
d.unload
