%% Plots node pressure for a specific time.
% This function contains:
%   Load a network.
%   Hydraulic analysis using EN functions.
%   Plot pressure after d.plot.
%   Plot pressure only nodes.
%   Unload library.

%% Clear - Start Toolkit
clear; close('all'); clc;
start_toolkit;

%% Load a network.
inpname = 'net2.inp';
d = epanet(inpname);

%% Hydraulic analysis using EN functions.
values = d.getComputedHydraulicTimeSeries;
press = values.Pressure;
flow = values.Flow;
coords = d.getNodeCoordinates;

%% Plot pressure after d.plot
d.plot('legendposition', 'best')
hold on
for i=1:length(press(1,:))
    scatter(coords{1}(i),coords{2}(i),35, press(1,i), 'filled')
end
title(['Pressure at hour 1 ,', inpname])
c = colorbar('southoutside');
c.Label.String = ['Pressure (', d.NodePressureUnits, ')'];

%% Plot pressure only nodes
d.plot('point', 'no', 'legend', 'hide');
hold on
for i=1:length(press(1,:))
    scatter(coords{1}(i),coords{2}(i),35, press(1,i), 'filled')
end
title(['Pressure at hour 1, ', inpname])
c = colorbar('southoutside');
c.Label.String = ['Pressure (', d.NodePressureUnits, ')'];

%% Unload library.
d.unload