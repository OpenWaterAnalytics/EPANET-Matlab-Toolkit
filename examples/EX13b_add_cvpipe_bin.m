%% Add CV pipe in a network with Bin Functions
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp', 'bin');

% Get from and to nodes for add cv pipe
fromNode = d.getBinNodeNameID.BinNodeNameID{2};
toNode = d.getBinNodeNameID.BinNodeNameID{6};

% Plot network
d.Binplot();
errcode = d.addBinCVPipe('CV-P1',fromNode,toNode,1000,10,100);

% Plot network with changes
d.Binplot('links','yes');

% Unload library
d.unload