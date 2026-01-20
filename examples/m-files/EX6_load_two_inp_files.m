%% Loads 2 different networks.
%% This example contains: 
%% 
% * Load 2 Input files. 
% * Load networks. 
% * Disp elevations for the two networks. 
% * Set new elevations. 
% * Disp new elevations. 
% * Disp computed values for the two networks.
% * Plot networks. 
% * Unload libraries.
%% Clear - Start Toolkit 

clear; close('all'); clc;
start_toolkit;
%% Load networks.

d1 = epanet('Net1.inp', 'ph'); 
d2 = epanet('Net2.inp', 'ph');
%% Disp elevations for the two networks.

disp('Net1 - Elevations:')
disp('--------------------')
d1_Elevs = d1.getNodeElevations;
disp(d1_Elevs)
disp('Net2 - Elevations:')
disp('------------------')
d2_Elevs = d2.getNodeElevations;
disp(d2_Elevs)
%% Set new elevations.

disp('Net1 - New Elevations:'); 
disp('----------------------')
d1.setNodeElevations(d1_Elevs + 200);
disp(d1.getNodeElevations)
disp('Net2 - New Elevations:')
disp('----------------------')
d2.setNodeElevations(d2_Elevs + 200);
disp(d2.getNodeElevations)
%% Disp computed values for the two networks.

disp('Net1 - Computed Time Series Values:')
disp('-----------------------------------')
disp(d1.getComputedTimeSeries)
disp('Net2 - Computed Time Series Values:')
disp('-----------------------------------')
disp(d2.getComputedTimeSeries)
%% Plot networks.

d1.plot;
d2.plot;
%% Unload libraries.

d1.unload;
d2.unload;