%% ADD a CV Pipe in a network
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

disp('Add a CV pipe')

% Plot network components
d.plot;

% Get from and to nodes for add cv pipe
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
node_index = d.addLinkPipeCV('CVPipe',fromNode,toNode);

% Plot new network with changes
d.plot('links','yes');

% Unload library
d.unload