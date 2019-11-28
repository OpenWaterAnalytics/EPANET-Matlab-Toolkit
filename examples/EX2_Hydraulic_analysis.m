%% Runs the hydraulic analysis of a network.
% This example contains:
%   Load a network.
%   Set simulation time duration.
%   Hydraulic analysis using ENepanet binary file.
%   Hydraulic analysis using epanet2d.exe binary file.
%   Hydraulic analysis.
%   Hydraulic analysis step-by-step.
%   Unload library.

%%
% Clear
clear; close('all'); clc;
start_toolkit;

%% Run hydraulic analysis of a network

% Load a network.
d = epanet('Net1.inp');

% Set simulation time duration.
hrs = 100;
d.setTimeSimulationDuration(hrs*3600);

% Hydraulic analysis using epanet2d.exe binary file.
% (This function ignores events)
tic
hyd_res_1 = d.getBinComputedAllParameters 
time_1 = toc;
tic;
hyd_res_2 = d.getComputedTimeSeries
time_2 = toc;

% Hydraulic analysis using ENepanet binary file (fastest).
% (This function ignores events)
tic;
hyd_res_3 = d.getComputedTimeSeries_ENepanet
time_3 = toc;

% Hydraulic analysis using the functions ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH.
% (This function contains events)
tic;
hyd_res_4 = d.getComputedHydraulicTimeSeries 
time_4 = toc;

% Hydraulic analysis step-by-step using the functions ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH.
% (This function contains events)
tic;
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1;P=[];T_H=[];D=[];H=[];F=[];
while (tstep>0)
    t=d.runHydraulicAnalysis;
    P=[P; d.getNodePressure];
    D=[D; d.getNodeActualDemand];
    H=[H; d.getNodeHydaulicHead];
    F=[F; d.getLinkFlows];
    T_H=[T_H; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis;
time_5 = toc;

% Unload library.
d.unload

disp(['Elapsed time for the function `getComputedTimeSeries` is: ' num2str(time_1)])
disp(['Elapsed time for the function `getBinComputedAllParameters` is: ' num2str(time_2)])
disp(['Elapsed time for the function `getComputedTimeSeries_ENepanet` is: ' num2str(time_3)])
disp(['Elapsed time for the function `getComputedHydraulicTimeSeries` is: ' num2str(time_4)])
disp(['Elapsed time for the function `step-by-step` is: ' num2str(time_5)])
