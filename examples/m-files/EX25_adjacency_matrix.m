%% Compute the adjacency matrix (connectivity graph) considering the flows,at different time steps or the mean flow. 
% **** For convenience, use EPANET class function getAdjacencyMatrix **** 
%% This example contains: 
%% 
% * Load network. 
% * Compute the new adjacency matrix based on the mean flow in the network. 
% * Unload library.
%% Clear - Start Toolkit

clear; close('all'); clc;
start_toolkit;
%% Load network

G = epanet('Net1.inp');
Fmean = mean(G.getComputedHydraulicTimeSeries('flow').Flow,1);
%% Compute the new adjacency matrix based on the mean flow in the network

Fsign = sign(Fmean);
Nidx = G.getLinkNodesIndex;
for i = 1:size(Nidx,1)
    if Fsign(i) == 1
        A(Nidx(i,1),Nidx(i,2)) = 1;
    else
        A(Nidx(i,2),Nidx(i,1)) = 1;
    end
end
A
%% Unload library

G.unload;