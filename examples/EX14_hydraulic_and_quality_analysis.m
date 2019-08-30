%% Runs the Hydraulic and Quality analysis of a network.
% This example contains:
%   Load a network.
%   Hydraulic and Quality analysis STEP-BY-STEP.
%   Display nodes pressures, links flows, nodes actual qualities, links actual qualities.
%   Unload library.

%% 
clear; close('all'); clc;
start_toolkit;

% Load a network.
d = epanet('Net1.inp');

% Set time hydraulic and quality steps
% etstep = 300;
% d.setTimeReportingStep(etstep);
% d.setTimeHydraulicStep(etstep);
% d.setTimeQualityStep(etstep);
% Hstep = min(Pstep,Hstep)
% Hstep = min(Rstep,Hstep)
% Hstep = min(Qstep,Hstep)

% Hydraulic and Quality analysis STEP-BY-STEP.
d.openHydraulicAnalysis;
d.openQualityAnalysis;
d.initializeHydraulicAnalysis(0);
d.initializeQualityAnalysis(d.ToolkitConstants.EN_NOSAVE);

tstep = 1;
T = []; P = []; F = []; QN = []; QL = [];
while (tstep>0)
    t = d.runHydraulicAnalysis;
    qt = d.runQualityAnalysis;

    P = [P; d.getNodePressure];
    F = [F; d.getLinkFlows];
    
    QN = [QN; d.getNodeActualQuality];
    QL = [QL; d.getLinkActualQuality];
    T = [T; t];

    tstep = d.nextHydraulicAnalysisStep;
    qtstep = d.nextQualityAnalysisStep;
end
d.closeQualityAnalysis;
d.closeHydraulicAnalysis;

% Display nodes pressures, links flows, nodes actual qualities, links actual qualities.
P
F
QN
QL

% Unload library.
d.unload;
