clear; close('all'); clc;
disp('==================== Starting EPANET Toolkit Script ====================');
start_toolkit;

%% --- Load Network ---
disp('==================== Loading Network: BWSN_Network_1.inp ====================');
d = epanet('BWSN_Network_1.inp');

%% --- Links / Valves ---
disp('==================== Valve Creation and Type ====================');
d.addLinkValvePCV('V1', '1', '2');
disp('Converting existing link to PCV');
d.setLinkTypeValvePCV(1);  
disp('Valve setup complete.');
disp(' ');

%% --- Curve Types ---
disp('==================== Curve Types ====================');
curveIndex = 1;
disp('Original curve type:');
disp(d.getCurveType(curveIndex));

disp('Changing curve types');
d.setCurveType(curveIndex, d.ToolkitConstants.EN_VOLUME_CURVE);
disp('After setCurveType(EN_VOLUME_CURVE):'); disp(d.getCurveType(curveIndex));

d.setCurveTypeVolume(curveIndex);
disp('After setCurveTypeVolume:'); disp(d.getCurveType(curveIndex));

d.setCurveTypePump(curveIndex);
disp('After setCurveTypePump:'); disp(d.getCurveType(curveIndex));

d.setCurveTypeEfficiency(curveIndex);
disp('After setCurveTypeEfficiency:'); disp(d.getCurveType(curveIndex));

d.setCurveTypeHeadloss(curveIndex);
disp('After setCurveTypeHeadloss:'); disp(d.getCurveType(curveIndex));

d.setCurveTypeGeneral(curveIndex);
disp('After setCurveTypeGeneral:'); disp(d.getCurveType(curveIndex));

d.setCurveTypeValveCurve(curveIndex);
disp('After setCurveTypeValveCurve:'); disp(d.getCurveType(curveIndex));
disp(' ');

%% --- Link / Valve Curve Accessors ---
disp('==================== Link / Valve Curve Accessors ====================');
linkID = d.getLinkPipeNameID{1};  
condition = 1;                     % example condition

% --- GPV Valve ---
disp('GPV Curve BEFORE:'); disp(d.getLinkValveCurveGPV());
indexGPV = d.setLinkTypeValveGPV(linkID, condition);  % assign GPV type
d.setLinkValveCurveGPV(indexGPV, 1);                  % set curve value
disp('GPV Curve AFTER:'); disp(d.getLinkValveCurveGPV());

% --- PCV Valve ---
disp('PCV Curve BEFORE:'); disp(d.getLinkValveCurvePCV());
indexPCV = d.setLinkTypeValvePCV(linkID, condition);  % assign PCV type
d.setLinkValveCurvePCV(indexPCV, 1);                  % set curve value
disp('PCV Curve AFTER:'); disp(d.getLinkValveCurvePCV());
disp(' ');

%% --- Controls & Rules ---
disp('==================== Controls and Rules ====================');
disp('Control state all:'); disp(d.getControlState());
disp('Control state 1:'); disp(d.getControlState(1));

disp('Rule enabled all:'); disp(d.getRuleEnabled());
disp('Rule enabled 1:'); disp(d.getRuleEnabled(1));

disp('Disabling rule 1');
d.setRuleEnabled(1, 0);
disp('Rule enabled 1 after disable:'); disp(d.getRuleEnabled(1));

disp('Links under control:'); disp(d.getLinkInControl());
disp('Nodes under control:'); disp(d.getNodeInControl());

disp('Time to next event:'); disp(d.getTimeToNextEvent());
d.unload;

d = epanet('Net1.inp');
%% --- Set Intermediate Vertex ---
disp('==================== Set Vertex ====================');
xCoords = [100.0, 150.0, 200.0];
yCoords = [200.0, 250.0, 300.0];
numVertices = numel(xCoords);

[~, countBefore] = d.apiENgetvertexcount(1, d.LibEPANET, d.ph);
disp('Vertex count before setting:'); disp(countBefore);

d.apiENsetvertices(1, xCoords, yCoords, numVertices, d.LibEPANET, d.ph);

[~, countAfter] = d.apiENgetvertexcount(1, d.LibEPANET, d.ph);
disp('Vertex count after setting:'); disp(countAfter);

[~, xOld, yOld] = d.apiENgetvertex(1, 1, d.LibEPANET, d.ph);
disp('Vertex 1 before manual set:'); disp([xOld, yOld]);

disp('Setting vertex 1 to (300, 800)');
d.setVertex(1, 1, 300.0, 800.0);

[~, xNew, yNew] = d.apiENgetvertex(1, 1, d.LibEPANET, d.ph);
disp('Vertex 1 after manual set:'); disp([xNew, yNew]);
disp(' ');

%% --- Patterns ---
disp('==================== Patterns ====================');
disp('Loading pattern file abc.pat as pattern 4');
d.loadPatternFile('abc.pat','4');
disp('Pattern average default value:'); disp(d.getPatternAverageDefaultValue());
disp(' ');

%% --- Options & Units ---
disp('==================== Options & Units ====================');
disp('Current pressure units:'); disp(d.getOptionsPressureUnits());
d.setOptionsPressureUnits(d.ToolkitConstants.EN_METERS);
d.setOptionsPressureUnitsPSI();
d.setOptionsPressureUnitsKPA();

disp('Status report:'); disp(d.getOptionsStatusReport());
d.setOptionsStatusReport(d.ToolkitConstants.EN_NO_REPORT);
d.setOptionsStatusReportNormal();
disp('Status report after normal:'); disp(d.getOptionsStatusReport());
d.setOptionsStatusReportFull();
disp('Status report after full:'); disp(d.getOptionsStatusReport());

d.setOptionsDemandPattern(2);
disp('Demand pattern set to:'); disp(d.getOptionsDemandPattern());

disp('Emitter backflow before:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowAllowed();
disp('Emitter backflow allowed:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowDisallowed();
disp('Emitter backflow disallowed:'); disp(d.getOptionsEmitterBackFlow());

disp('Flow units before change:'); disp(d.getFlowUnits());
disp('Setting flow units to CMS');
d.setFlowUnitsCMS();
disp('Flow units after change:'); disp(d.getFlowUnits());
disp(' ');

%% --- Leakage Expansion Demand ---
% See EX28_Leakage 

%% --- Project I/O & Bulk Retrieval ---
disp('==================== Project I/O & Bulk Retrieval ====================');
d.openX('Net1.inp', 'Net1.rpt', 'Net1.out');
disp('Link info (flow):'); disp(d.getLinkInfo(d.ToolkitConstants.EN_FLOW));
disp(' ');

%% --- Statistics ---
disp('==================== Statistics ====================');
d.runsCompleteSimulation();
disp('Iterations:'); disp(d.getStatisticIterations());
disp('Relative error:'); disp(d.getStatisticRelativeError());
disp('Deficient nodes:'); disp(d.getStatisticDeficientNodes());
disp('Demand reduction:'); disp(d.getStatisticDemandReduction());
disp('Total leakage loss:'); disp(d.getStatisticTotalLeakageLoss());

d.unload();
