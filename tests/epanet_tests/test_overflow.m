%% Test overflow
% https://github.com/OpenWaterAnalytics/EPANET/blob/dev/tests/test_overflow.cpp

d = epanet('Net1.inp');

% Get index of the tank and its inlet/outlet pipe
Nindex = d.getNodeIndex('2');
Lindex = d.getLinkIndex('110');

% Set initial & maximum level to 130
d.setNodeTankInitialLevel(Nindex, 130);
d.setNodeTankMaximumWaterLevel(Nindex, 130)

% Set duration to 1 hr
d.setTimeSimulationDuration(3600);

% Solve hydraulics with default of no tank spillage allowed
d.solveCompleteHydraulics

% Check that tank remains full
level = d.getNodeTankInitialLevel(Nindex);
if (abs(level - 130.0) < 0.0001)
    fprintf('\nTrue');
else
    fprintf('\nFalse');
end

% Check that there is no spillage
spillage = d.getNodeActualDemand(Nindex);
if (abs(spillage) < 0.0001)
    fprintf('\nTrue');
else
    fprintf('\nFalse');
end

% Check that inflow link is closed
inflow = d.getLinkFlows(Lindex);
if abs(inflow) < 0.0001
    fprintf('\nTrue');
else
    fprintf('\nFalse');
end

% Turn tank overflow option on
d.setNodeTankCanOverFlow(Nindex, 1);

% Solve hydraulics again
d.solveCompleteHydraulics

% Check that tank remains full
level = d.getNodeTankInitialLevel(Nindex);
if (abs(level - 130.0) < 0.0001)
    fprintf('\nTrue');
else
    fprintf('\nFalse');
end

% Check that there is spillage equal to tank inflow
% (inflow has neg. sign since tank is start node of inflow pipe)
spillage = d.getNodeActualDemand(Nindex);
if (abs(spillage) > 0.0001)
    fprintf('\nTrue');
else
    fprintf('\nFalse');
end

inflow = d.getLinkFlows(Lindex);
if (abs(-inflow - spillage) < 0.0001)
    fprintf('\nTrue\n\n');
else
    fprintf('\nFalse\n\n');
end

% % Save project to file and then close it
d.saveInputFile('net1_overflow.inp');
d.unload;

% % Re-open saved file & run it
d = epanet('net1_overflow.inp');
d.solveCompleteHydraulics;
% 
% % Check that tank spillage has same value as before
spillage2 = d.getNodeActualDemand(Nindex);
if (abs(spillage - spillage2) < 0.0001)
    fprintf('\nTrue\n\n');
else
    fprintf('\nFalse\n\n');
end
 
% % Clean up
d.unload;