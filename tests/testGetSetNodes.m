%% EPANET-Matlab Toolkit Test Part 2
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

%% *Get Nodes Data (EXAMPLES)*
all_elevations = d.getNodeElevations;
disp(all_elevations);

elevationForIndex10 = d.getNodeElevations(10);
disp(elevationForIndex10);

all_elevations2 = d.getNodeElevations(1:d.getNodeCount);
disp(all_elevations2);

elevationsSp = d.getNodeElevations([1 5 10]);
disp(elevationsSp);


newFunctionsDev2_1 = {'ENgetpumptype' , 'ENgetheadcurveindex', 'ENsetcurvevalue',...
            'ENsetcurve', 'ENaddcurve', 'ENgetcurvevalue', 'ENgetcurve',...
            'ENgetcurvelen', 'ENgetcurveid', 'ENgetcurveindex', 'ENsetcoord',...
            'ENgetcoord', 'ENgetstatistic', 'ENgetnumdemands', 'ENgetbasedemand',...	
            'ENgetdemandpattern', 'ENsetbasedemand', 'ENgetaveragepatternvalue'};
        
%% New Functions 2.1
nF=0; % old dll
for i=1:length(newFunctionsDev2_1)
    if sum(strcmp(d.libFunctions,newFunctionsDev2_1(i)))
        nF=1; % new dll
        break;
    end
end
if nF==1
    d.getNodeDemandCategoriesNumber
    d.getNodeDemandCategoriesNumber(2)

    numCategories=1;nodeindex=2;
    d.getNodeDemandPatternIndex{numCategories}
    d.getNodeDemandPatternIndex{numCategories}(nodeindex)
    d.getNodeDemandPatternNameID{numCategories}
    d.getNodeDemandPatternNameID{numCategories}(nodeindex)
    
    d.getNodeBaseDemands{numCategories}
    d.getNodeBaseDemands{numCategories}(nodeindex)
end

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
tankInd = d.NodeTankIndex(1);

d.getNodeTankInitialLevel
d.getNodeTankInitialLevel(tankInd)

d.getNodeTankInitialWaterVolume
d.getNodeTankInitialWaterVolume(tankInd)

d.getNodeTankMixZoneVolume
d.getNodeTankMixZoneVolume(tankInd)

d.getNodeTankDiameter
d.getNodeTankDiameter(tankInd)

d.getNodeTankMinimumWaterVolume
d.getNodeTankMinimumWaterVolume(tankInd)

d.getNodeTankVolumeCurveIndex
d.getNodeTankVolumeCurveIndex(tankInd)

d.getNodeTankMinimumWaterLevel
d.getNodeTankMinimumWaterLevel(tankInd)

d.getNodeTankMaximumWaterLevel
d.getNodeTankMaximumWaterLevel(tankInd)

d.getNodeTankMinimumFraction
d.getNodeTankMinimumFraction(tankInd)

d.getNodeTankBulkReactionCoeff
d.getNodeTankBulkReactionCoeff(tankInd)

d.getNodeTankVolume
d.getNodeTankVolume(tankInd)

if nF==1
    d.getNodeTankMaximumWaterVolume
    d.getNodeTankMaximumWaterVolume(tankInd)
end

d.getNodeType
d.getNodeType(tankInd)

d.getNodeNameID
d.getNodeNameID(tankInd)

d.getNodeCoordinates
d.getNodeCoordinates{1}
d.getNodeCoordinates{2}
if nF==1
    d.getNodeCoordinates(2)
    % d.getNodeCoordinates(1:d.NodeCount)
end

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

if nF==1
    d.getNodeCoordinates(2)
    d.setNodeCoordinates(2,[10 10]);
    d.getNodeCoordinates(2)

    d.getNodeBaseDemands{1}
    d.setNodeBaseDemands(3,20);
    d.getNodeBaseDemands{1}

    d.getNodeDemandPatternIndex{1}
    d.setNodeDemandPatternIndex(3,0); %remove pattern..
    d.getNodeDemandPatternIndex{1}
end

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
d.getNodeTankInitialLevel
d.setNodeTankInitialLevel(d.getNodeTankInitialLevel+20);
v=d.getNodeTankInitialLevel
d.setNodeTankInitialLevel(tankInd,v(tankInd)+10);
d.getNodeTankInitialLevel

d.getNodeTankDiameter
d.setNodeTankDiameter(d.getNodeTankDiameter+20);
d.getNodeTankDiameter
d.setNodeTankDiameter(tankInd,100);
d.getNodeTankDiameter

d.getNodeTankBulkReactionCoeff
d.setNodeTankBulkReactionCoeff(d.getNodeTankBulkReactionCoeff+1);
d.getNodeTankBulkReactionCoeff
d.setNodeTankBulkReactionCoeff(tankInd,-1);
d.getNodeTankBulkReactionCoeff

d.getNodeTankMaximumWaterLevel
d.setNodeTankMaximumWaterLevel(d.getNodeTankMaximumWaterLevel+21);
d.getNodeTankMaximumWaterLevel
d.setNodeTankMaximumWaterLevel(tankInd,200);
d.getNodeTankMaximumWaterLevel

d.getNodeTankMinimumWaterLevel
d.setNodeTankMinimumWaterLevel(d.getNodeTankMinimumWaterLevel-21);
n=d.getNodeTankMinimumWaterLevel
d.setNodeTankMinimumWaterLevel(tankInd,n(tankInd)+20);
d.getNodeTankMinimumWaterLevel

d.getNodeTankMinimumFraction
d.setNodeTankMinimumFraction(d.getNodeTankMinimumFraction+0.1);
d.getNodeTankMinimumFraction
d.setNodeTankMinimumFraction(tankInd,0.2);
d.getNodeTankMinimumFraction

d.getNodeTankMinimumWaterVolume
d.setNodeTankMinimumWaterVolume(d.getNodeTankMinimumWaterVolume+10000);
d.getNodeTankMinimumWaterVolume
d.setNodeTankMinimumWaterVolume(tankInd,20000);
d.getNodeTankMinimumWaterVolume

values = d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
values{end}='MIX2';
d.setNodeTankMixingModelType(values);
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
d.setNodeTankMixingModelType(tankInd,'FIFO');
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode

d.unload
fprintf('Test finished.\n')