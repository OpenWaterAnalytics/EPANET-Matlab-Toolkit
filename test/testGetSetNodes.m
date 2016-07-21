%% EPANET-Matlab Class Test Net1 **01/07/2016**
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc;
clear;
close all;

% Create EPANET object using the INP file
inpname='networks/Net1_Rossman2000.inp'; % Net1_Rossman2000
% Net2_Rossman2000 Net3_Rossman2000 BWSN1_Ostfeld2008 
d=epanet(inpname);

%% *Get Nodes Data (EXAMPLES)*
all_elevations = d.getNodeElevations;
disp(all_elevations);

elevationForIndex10 = d.getNodeElevations(10);
disp(elevationForIndex10);

all_elevations2 = d.getNodeElevations(1:d.getNodeCount);
disp(all_elevations2);

elevationsSp = d.getNodeElevations([1 5 10]);
disp(elevationsSp);


d.getNodeDemandCategoriesNumber
d.getNodeDemandCategoriesNumber(2)

numCategories=1;nodeindex=2;
d.getNodeDemandPatternIndex{numCategories}
d.getNodeDemandPatternIndex{numCategories}(nodeindex)
d.getNodeDemandPatternNameID{numCategories}
d.getNodeDemandPatternNameID{numCategories}(nodeindex)

d.getNodeBaseDemands{numCategories}
d.getNodeBaseDemands{numCategories}(nodeindex)

d.getNodePatternIndex
d.getNodePatternIndex(2)

d.getNodeEmitterCoeff
d.getNodeEmitterCoeff(3)

d.getNodeInitialQuality
d.getNodeInitialQuality(1)

d.getNodeSourceQuality
d.getNodeSourceQuality(1)

d.getNodeSourcePatternIndex
d.getNodeSourcePatternIndex(2)

d.getNodeSourceTypeIndex
d.getNodeSourceTypeIndex(2)

d.getNodeSourceType
d.getNodeSourceType(3)

%% Tanks
d.getNodeTankInitialLevel
d.getNodeTankInitialLevel(11)

d.getNodeTankInitialWaterVolume
d.getNodeTankInitialWaterVolume(11)

d.getNodeTankMixZoneVolume
d.getNodeTankMixZoneVolume(11)

d.getNodeTankDiameter
d.getNodeTankDiameter(11)

d.getNodeTankMinimumWaterVolume
d.getNodeTankMinimumWaterVolume(11)

d.getNodeTankVolumeCurveIndex
d.getNodeTankVolumeCurveIndex(11)

d.getNodeTankMinimumWaterLevel
d.getNodeTankMinimumWaterLevel(11)

d.getNodeTankMaximumWaterLevel
d.getNodeTankMaximumWaterLevel(11)

d.getNodeTankMinimumFraction
d.getNodeTankMinimumFraction(11)

d.getNodeTankBulkReactionCoeff
d.getNodeTankBulkReactionCoeff(11)

d.getNodeTankVolume
d.getNodeTankVolume(11)

d.getNodeTankMaxVolume
d.getNodeTankMaxVolume(11)

d.getNodeType
d.getNodeType(11)

d.getNodeNameID
d.getNodeNameID(11)

d.getNodeCoordinates
d.getNodeCoordinates{1}
d.getNodeCoordinates{2}
d.getNodeCoordinates(2)
d.getNodeCoordinates(1:d.NodeCount)
% Runs hydraulics Step-by-step
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1; T=[];P=[];H=[];D=[];Q=[];
index=2;
while (tstep>0)
    t=d.runHydraulicAnalysis;
    D=[D; d.getNodeActualDemand(index)];
    H=[H; d.getNodeHydaulicHead(index)];
    P=[P; d.getNodePressure(index)];
    Q=[Q; d.getNodeActualQuality(index)];
    T=[T; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis
disp(D);
disp(Q);
disp(H);
disp(P);

%% Set nodes info
d.getNodeElevations
d.setNodeElevations(2*d.getNodeElevations);
d.getNodeElevations
d.getNodeElevations(2)
d.setNodeElevations(2,200); %index, value
d.getNodeElevations(2)

d.getNodeEmitterCoeff
d.setNodeEmitterCoeff(2*ones(1,d.NodeCount));
d.getNodeEmitterCoeff
d.getNodeEmitterCoeff(2)
d.setNodeEmitterCoeff(2,1.5); %index, value
d.getNodeEmitterCoeff(2)

d.getNodeInitialQuality
d.setNodeInitialQuality(2*d.getNodeInitialQuality);
d.getNodeInitialQuality
d.getNodeInitialQuality(2)
d.setNodeInitialQuality(2,1.5); %index, value
d.getNodeInitialQuality(2)

d.getNodeCoordinates(2)
d.setNodeCoordinates(2,[10 10]);
d.getNodeCoordinates(2)

d.getNodeBaseDemands{1}
d.setNodeBaseDemands(3,20);
d.getNodeBaseDemands{1}

d.getNodeDemandPatternIndex{1}
d.setNodeDemandPatternIndex(3,0); %remove pattern..
d.getNodeDemandPatternIndex{1}

d.getNodeSourceType
d.setNodeSourceType(1,'MASS')
d.setNodeSourceType(2,'CONCEN')
d.setNodeSourceType(3,'SETPOINT')
d.setNodeSourceType(4,'FLOWPACED')
d.getNodeSourceType

d.getNodeSourcePatternIndex
d.setNodeSourcePatternIndex(1,1);
d.getNodeSourcePatternIndex

d.getNodeSourceQuality
d.setNodeSourceQuality(2*d.getNodeSourceQuality);
d.getNodeSourceQuality
d.setNodeSourceQuality(3,20);
d.getNodeSourceQuality

%% Set tanks info
indTank = d.getNodeTankIndex

d.getNodeTankInitialLevel
d.setNodeTankInitialLevel(d.getNodeTankInitialLevel+20);
v=d.getNodeTankInitialLevel
d.setNodeTankInitialLevel(indTank(1),v(indTank(1))+10);
d.getNodeTankInitialLevel

d.getNodeTankDiameter
d.setNodeTankDiameter(d.getNodeTankDiameter+20);
d.getNodeTankDiameter
d.setNodeTankDiameter(indTank(1),100);
d.getNodeTankDiameter

d.getNodeTankBulkReactionCoeff
d.setNodeTankBulkReactionCoeff(d.getNodeTankBulkReactionCoeff+1);
d.getNodeTankBulkReactionCoeff
d.setNodeTankBulkReactionCoeff(indTank(1),-1);
d.getNodeTankBulkReactionCoeff

d.getNodeTankMaximumWaterLevel
d.setNodeTankMaximumWaterLevel(d.getNodeTankMaximumWaterLevel+21);
d.getNodeTankMaximumWaterLevel
d.setNodeTankMaximumWaterLevel(indTank(1),200);
d.getNodeTankMaximumWaterLevel

d.getNodeTankMinimumWaterLevel
d.setNodeTankMinimumWaterLevel(d.getNodeTankMinimumWaterLevel-21);
n=d.getNodeTankMinimumWaterLevel
d.setNodeTankMinimumWaterLevel(indTank(1),n(indTank(1))+20);
d.getNodeTankMinimumWaterLevel

d.getNodeTankMinimumFraction
d.setNodeTankMinimumFraction(d.getNodeTankMinimumFraction+0.1);
d.getNodeTankMinimumFraction
d.setNodeTankMinimumFraction(indTank(1),0.2);
d.getNodeTankMinimumFraction

d.getNodeTankMinimumWaterVolume
d.setNodeTankMinimumWaterVolume(d.getNodeTankMinimumWaterVolume+10000);
d.getNodeTankMinimumWaterVolume
d.setNodeTankMinimumWaterVolume(indTank(1),20000);
d.getNodeTankMinimumWaterVolume

values = d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
values{end}='MIX2';
d.setNodeTankMixingModelType(values);
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
d.setNodeTankMixingModelType(indTank(1),'FIFO');
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode

d.unload

fprintf('Test finished.\n')