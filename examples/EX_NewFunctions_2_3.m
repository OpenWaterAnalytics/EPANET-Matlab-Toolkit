%% 
%% Clear and start
clear; close all; clc;
start_toolkit;

%% Load Network
inpfile = 'Net1.inp';
%inpfile = 'Richmond_standard.inp';
d = epanet(inpfile);

%% --- Leakage / Expansion / Demand ---
d.solveCompleteHydraulics();

d.setLinkLeakArea(2, 10.5)
disp('Link Leak Area:');
d.getLinkLeakArea(2)

d.setLinkExpansionProperties(5, 2)
disp('Link Expansion Properties:');
d.getLinkExpansionProperties(5)

disp('Link Leakage Rate:');
d.getLinkLeakageRate()

disp('Consumer Demand Requested:');
d.getConsumerDemandRequested(5)

disp('Consumer Demand Delivered:');
d.getConsumerDemandDelivered(5)

disp('Node Emitter Flow:');
d.getNodeEmitterFlow(1)

%% --- Curve Types ---
curveIndex = 1;

disp('--- Curve Types ---');
disp('Original curve type:');
d.getCurveType(curveIndex)

disp('setCurveTypePump');
d.setCurveTypePump(curveIndex)
d.getCurveType(curveIndex)

disp('setCurveTypeVolume');
d.setCurveTypeVolume(curveIndex)
d.getCurveType(curveIndex)

disp('setCurveTypeGeneral');
d.setCurveTypeGeneral(curveIndex)
d.getCurveType(curveIndex)

disp('setCurveTypeHeadloss');
d.setCurveTypeHeadloss(curveIndex)
d.getCurveType(curveIndex)

disp('setCurveTypeEfficiency');
d.setCurveTypeEfficiency(curveIndex)
d.getCurveType(curveIndex)

disp('setCurveTypeValveCurve');
d.setCurveTypeValveCurve(curveIndex)
d.getCurveType(curveIndex)

%% --- Set Vertices ---
linkID = '10';  
xCoords = [22, 24, 28];
yCoords = [30, 68, 69];
d.setLinkVertices(linkID, xCoords, yCoords);
disp('Vertices after setLinkVertices:');
d.getLinkVertices(linkID)

% Example using setVertex
d.setVertex(1, 1, 1, 1);
disp('Vertices after setVertex:');
d.getLinkVertices('10')

%% --- Link/Valve Curves ---
linkIndex = 1;
condition = 1;

d.setLinkTypeValveGPV(linkIndex, condition)
d.setLinkValveCurveGPV(linkIndex, 1)
disp('GPV Curve:');
d.getLinkValveCurveGPV(linkIndex)

d.setLinkTypeValvePCV(linkIndex, condition)
d.setLinkValveCurvePCV(linkIndex, 1)
disp('PCV Curve:');
d.getLinkValveCurvePCV(linkIndex)

%% --- Demand Pattern ---
disp('Demand Pattern:');
d.getOptionsDemandPattern()
d.setOptionsDemandPattern(0)
d.getOptionsDemandPattern()

%% --- Link / Valve Types ---
disp('Link Type:');
d.getLinkType(1)
linkID = 1;
d.setLinkTypeValvePCV(linkID)
d.getLinkType(linkID)

% Add new PCV valve
valveID = 'newValvePCV';
fromNode = '10';
toNode = '21';
d.addLinkValvePCV(valveID, fromNode, toNode);

%% --- Emitter Backflow Options ---
disp('Emitter Backflow:');
d.getOptionsEmitterBackFlow()
d.setOptionsEmitterBackFlowDisallowed()
d.getOptionsEmitterBackFlow()
d.setOptionsEmitterBackFlowAllowed()
d.getOptionsEmitterBackFlow()

%% --- Controls ---
disp('Link/Node Controls:');
d.getLinkInControl()
d.getLinkInControl(10)
d.getNodeInControl()
d.getNodeInControl(11)

%% --- Flow Units & Pressure Units ---
d.setFlowUnitsCMS()
disp('Flow Units:');
d.getFlowUnits()

d.setOptionsPressureUnitsPSI()
disp('Pressure Units:');
d.getOptionsPressureUnits()

d.setFlowUnitsIMGD()
d.setOptionsPressureUnitsPSI()
disp('Pressure Units after IMGD + PSI:');
d.getOptionsPressureUnits()

d.setOptionsPressureUnitsKPA()
disp('Pressure Units after KPA:');
d.getOptionsPressureUnits()

% Reset to defaults
d.setFlowUnitsGPM()
d.setOptionsPressureUnitsPSI()

%% --- Node Emitters ---
d.setNodeEmitterCoeff([2,3], [0.5,0.6])
disp('Node Emitter Coefficients:');
d.getNodeEmitterCoeff()

%% --- Status Report ---
d.setOptionsStatusReportNo()
disp('Status Report:');
d.getOptionsStatusReport()
d.setOptionsStatusReportNormal()
d.getOptionsStatusReport()
d.setOptionsStatusReportFull()
d.getOptionsStatusReport()

%% --- Unload ---
d.unload()
