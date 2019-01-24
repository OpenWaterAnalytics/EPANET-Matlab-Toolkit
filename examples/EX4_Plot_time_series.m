%% Visualise/plot time series example
% This example contains:
% Load a network
% Hydraulic analysis using ENepanet binary file
% Change time-stamps from seconds to hours
% Plot node pressures for specific nodes 
% Plot water velocity for specific links
% Plot water flow for specific links

%%
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Hydraulic analysis using ENepanet binary file
% (This function ignore events)
hyd_res = d.getComputedTimeSeries;

% Change time-stamps from seconds to hours
hrs_time = hyd_res.Time/3600;

% Plot node pressures for specific nodes 
node_indices = [1, 3, 5];
node_names = d.getNodeNameID(node_indices)
for i=node_indices
    figure;
    plot(hrs_time, hyd_res.Pressure(:,i));
    title(['Pressure for the node id "', d.getNodeNameID{i},'"']);
    xlabel('Time (hrs)'); 
    ylabel(['Pressure (', d.NodePressureUnits,')'])
end

% Plot water velocity for specific links
link_indices = [4, 8, 10];
link_names = d.getNodeNameID(link_indices)
for i=link_indices
    figure;
    plot(hrs_time, hyd_res.Velocity(:,i));
    title(['Velocity for the link id "', d.getLinkNameID{i},'"']);
    xlabel('Time (hrs)'); 
    ylabel(['Velocity (', d.LinkVelocityUnits,')'])
end

% Plot water flow for specific links
link_indices = [2, 3, 9];
for i=link_indices
    figure;
    plot(hrs_time, hyd_res.Flow(:,i));
    title(['Flow for the link id "', d.getLinkNameID{i},'"']);
    xlabel('Time (hrs)'); 
    ylabel(['Flow (', d.LinkFlowUnits,')'])
end

% Unload library
d.unload
