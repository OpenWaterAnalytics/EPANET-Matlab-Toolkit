%% Hydraulic analysis example
% This function contains:
% Load a network
% Set simulation time duration
% Hydraulic analysis using ENepanet binary file
% Hydraulic analysis using epanet2d.exe binary file
% Hydraulic analysis 
% Hydraulic analysis step-by-step

%% 
clear; close('all'); clc;
start_toolkit;

%% Run hydraulic analysis of a network

% Load a network
d = epanet('Net1.inp');

% Set simulation time duration
hrs = 100;
d.setTimeSimulationDuration(hrs*3600);

% Hydraulic analysis using ENepanet binary file (fastest)
% (This function ignores events)
hyd_res_1 = d.getComputedTimeSeries 

% Hydraulic analysis using epanet2d.exe binary file
% (This function ignores events)
hyd_res_2 = d.getBinComputedAllParameters 

% Hydraulic analysis using the functions ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH
% (This function contains events)
hyd_res_3 = d.getComputedHydraulicTimeSeries 

% Hydraulic analysis step-by-step using the functions ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH
% (This function contains events)
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

% Unload library
d.unload