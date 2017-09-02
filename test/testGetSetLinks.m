%% EPANET-Matlab Toolkit Test Part 3
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc;
clear;
close all;

% Create EPANET object using the INP file
inpname='Net1.inp';  
% Net1 Net2 Net3 BWSN_Network_1 
d=epanet(inpname);
% d=epanet(inpname, 'epanet2');

%% *Get Links Data (EXAMPLES)*
all_diameters = d.getLinkDiameter;
disp(all_diameters);

diameterForIndex10 = d.getLinkDiameter(10);
disp(diameterForIndex10);

all_diameters2 = d.getLinkDiameter(1:d.getLinkCount);
disp(all_diameters2);

diametersSp = d.getLinkDiameter([1 5 10]);
disp(diametersSp);

% similar..
d.getLinkLength
d.getLinkLength(5)

d.getLinkRoughnessCoeff
d.getLinkRoughnessCoeff(3)

d.getLinkMinorLossCoeff
d.getLinkMinorLossCoeff(2)

d.getLinkInitialStatus
d.getLinkInitialStatus(2)

d.getLinkInitialSetting
d.getLinkInitialSetting(2)

d.getLinkBulkReactionCoeff
d.getLinkBulkReactionCoeff(2)

d.getLinkWallReactionCoeff
d.getLinkWallReactionCoeff(2)

d.getLinkType
d.getLinkType(13)

d.getLinkTypeIndex
d.getLinkTypeIndex(5:d.getLinkCount)

% Runs hydraulics Step-by-step
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1; T=[]; V=[]; H=[];F=[];
index=1;
while (tstep>0)
    t=d.runHydraulicAnalysis;
    V=[V; d.getLinkVelocity(index)];
    H=[H; d.getLinkHeadloss(index)];
    F=[F; d.getLinkFlows(index)];
    T=[T; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis
disp(V);
disp(H);
disp(F);
d.getLinkStatus
d.getLinkStatus(2)

d.getLinkSettings
d.getLinkSettings(2)

% d.getLinkQuality    % bug in epanet
% d.getLinkQuality(2)

%% Set links info
d.getLinkDiameter
d.setLinkDiameter(2*d.getLinkDiameter);
d.getLinkDiameter
d.getLinkDiameter(2)
d.setLinkDiameter(2,200); %index, value
d.getLinkDiameter(2)

d.getLinkLength
d.setLinkLength(2*d.getLinkLength)
d.getLinkLength
d.getLinkLength(2)
d.setLinkLength(2,500)%index, value
d.getLinkLength(2)

d.getLinkRoughnessCoeff
d.setLinkRoughnessCoeff(2*d.getLinkRoughnessCoeff)
d.getLinkRoughnessCoeff
d.getLinkRoughnessCoeff(2)
d.setLinkRoughnessCoeff(2,150)%index, value
d.getLinkRoughnessCoeff(2)

d.getLinkMinorLossCoeff
d.setLinkMinorLossCoeff(d.getLinkMinorLossCoeff+1.1)
d.getLinkMinorLossCoeff
d.getLinkMinorLossCoeff(2)
d.setLinkMinorLossCoeff(2,1.01)%index, value
d.getLinkMinorLossCoeff(2)

d.getLinkInitialStatus
d.setLinkInitialStatus(0*d.getLinkInitialStatus)
d.getLinkInitialStatus
d.getLinkInitialStatus(2)
d.setLinkInitialStatus(2,1)
d.getLinkInitialStatus(2)

d.getLinkBulkReactionCoeff
d.setLinkBulkReactionCoeff(d.getLinkBulkReactionCoeff-0.055)
d.getLinkBulkReactionCoeff
d.getLinkBulkReactionCoeff(1)
d.setLinkBulkReactionCoeff(1,0.2)
d.getLinkBulkReactionCoeff(1)

d.getLinkWallReactionCoeff
d.setLinkWallReactionCoeff(-1.1*d.getLinkWallReactionCoeff)
d.getLinkWallReactionCoeff
d.getLinkWallReactionCoeff(2)
d.setLinkWallReactionCoeff(2,-2)
d.getLinkWallReactionCoeff(2)

linkset=d.getLinkInitialSetting
if d.LinkValveCount
    linkset(d.LinkValveIndex)=0;
end
d.setLinkInitialSetting(linkset*10)
d.getLinkInitialSetting
d.getLinkInitialSetting(2)
d.setLinkInitialSetting(2,10)
d.getLinkInitialSetting(2)

d.getLinkStatus %dynamic
d.setLinkStatus(0*d.getLinkStatus)
d.getLinkStatus 
d.getLinkStatus(2)
d.setLinkStatus(2,1)
d.getLinkStatus(2) 

d.getLinkSettings %dynamic
d.setLinkSettings(d.getLinkSettings+10)
d.getLinkSettings
d.getLinkSettings(2)
d.setLinkSettings(2,121)
d.getLinkSettings(2)

d.unload
fprintf('Test finished.\n')