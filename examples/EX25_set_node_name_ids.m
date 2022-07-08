%% Set node name IDs
% This example contains:
%   Load a network.
%   Get node name ids.
%   Set new node name ids.
%   Get new node name ids.
%   Unload library.

%%
% Clear - Start Toolkit
clear; close('all'); clc;

start_toolkit;

% Load network
d = epanet('Net1.inp');

disp('Node name ids:')
d.getNodeNameID'

% Set your prefix 
junction_prefix = 'J';
reservoir_prefix = 'R';
tank_prefix = 'T';

% Update node names 
for i=d.getNodeJunctionIndex
    d.setNodeNameID(i, [junction_prefix, '-', num2str(i)]);
end

for i=d.getNodeReservoirIndex
    d.setNodeNameID(i, [reservoir_prefix, '-', num2str(i)]);
end

for i=d.getNodeTankIndex
    d.setNodeNameID(i, [tank_prefix, '-', num2str(i)]);
end

disp('Node name ids after:')
d.getNodeNameID'

d.unload;
