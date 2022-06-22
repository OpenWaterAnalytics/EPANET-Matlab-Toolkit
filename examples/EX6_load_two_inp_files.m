%% Loads 2 different networks.
% This example contains:
%   Load 2 Input files.
%   Load networks.
%   Disp elevations for the two networks.
%   Plot networks.
%   Unload libraries.

%%
% Clear 
clear; close('all'); clc;
start_toolkit;

% Load networks.
d1 = epanet('Net1.inp');
d2 = epanet('Net2.inp');


% Disp elevations for the two networks.
disp('Net1 - Elevations:')
disp('------------------')
disp(d1.getNodeElevations)
disp('Net2 - Elevations:')
disp('------------------')
disp(d2.getNodeElevations)

% Plot networks.
d1.plot;
d2.plot;

% Unload libraries.
d1.unload;
d2.unload;
