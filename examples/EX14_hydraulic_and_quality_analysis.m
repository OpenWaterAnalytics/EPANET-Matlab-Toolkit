%% Hydraulic and Quality analysis
% This function contains:
% Load a network
% Hydraulic and Quality analysis STEP-BY-STEP
%% 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Set time hydraulic and quality steps
tstep = 3600;
d.setTimeHydraulicStep(tstep);
d.setTimeQualityStep(tstep);

% Hydraulic and Quality analysis STEP-BY-STEP
d.openHydraulicAnalysis;
d.openQualityAnalysis;
d.initializeHydraulicAnalysis(0);
d.initializeQualityAnalysis(d.ToolkitConstants.EN_NOSAVE);

T = []; P = []; F = []; QN = [];
while (tstep>0)
    t = d.runHydraulicAnalysis;
    qt = d.runQualityAnalysis;

    P = [P; d.getNodePressure];
    F = [F; d.getLinkFlows];
    
    QN = [QN; d.getNodeActualQuality];
    T = [T; t];
    
    qtstep = d.nextQualityAnalysisStep;
    tstep = d.nextHydraulicAnalysisStep;
end
d.closeQualityAnalysis;
d.closeHydraulicAnalysis;

P
F
QN

% Unload library 
d.unload;
