%% Set pump curves
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Computed Hydraulic & Quality Time Series via ENepanet binary file
Results = d.getComputedTimeSeries;

nodeID = '10';
indexNode = d.getNodeIndex(nodeID);

% Plot pressures
figure;
subplot(2,1,1);
plot(Results.Time, Results.Pressure(:,indexNode));
title('Before set curve');
xlabel('Time (hrs)'); 
ylabel(['Pressure (', d.NodePressureUnits,')'])
    
% Get head curve
headCurve = d.getLinkPumpHeadCurveIndex;

% set new head curve values
d.setCurve(headCurve,[2000 250]); 

% Computed hydraulics
Results = d.getComputedTimeSeries;  
subplot(2,1,2);
plot(Results.Time, Results.Pressure(:,indexNode));
title('After set curve');
xlabel('Time (hrs)'); 
ylabel(['Pressure (', d.NodePressureUnits,')'])

% Unload library
d.unload

