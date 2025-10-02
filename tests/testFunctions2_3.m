clear; close('all'); clc;
disp('==================== Starting EPANET Toolkit Script ====================');
start_toolkit;

%% --- Load Network ---
disp('==================== Loading Network: BWSN_Network_1.inp ====================');
d = epanet('BWSN_Network_1.inp');

disp('pattern count before loading abc.pat: '); disp(d.getPatternCount);
disp('Loading pattern file abc.pat');
d.loadPatternFile('abc.pat','5');
disp('pattern count after loading abc.pat'); disp(d.getPatternCount);

%% --- Curve Types ---
disp('==================== Curve Types ====================');
curveIndex = 1;
disp('Original curve type:');
disp(d.getCurveType(curveIndex));

disp('Changing curve types');
d.setCurveType(curveIndex, d.ToolkitConstants.EN_VOLUME_CURVE);
disp('After setCurveType(EN_VOLUME_CURVE):'); disp(d.getCurveType(curveIndex)');

d.setCurveTypePump(curveIndex);
disp('After setCurveTypePump:'); disp(d.getCurveType(curveIndex));

d.setCurveTypeVolume(curveIndex);
disp('After setCurveTypeVolume:'); disp(d.getCurveType(curveIndex)');

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
disp('GPV Curve BEFORE:'); disp(d.getLinkValveCurveGPV()');
indexGPV = d.setLinkTypeValveGPV(linkID, condition);  % assign GPV type
d.setLinkValveCurveGPV(indexGPV, 1);                  % set curve value
disp('GPV Curve AFTER:'); disp(d.getLinkValveCurveGPV()');

% --- PCV Valve ---
disp('PCV Curve BEFORE:'); disp(d.getLinkValveCurvePCV()');
indexPCV = d.setLinkTypeValvePCV(linkID, condition);  % assign PCV type
d.setLinkValveCurvePCV(indexPCV, 1);                  % set curve value
disp('PCV Curve AFTER:'); disp(d.getLinkValveCurvePCV()');
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

disp('Links under control:'); disp(d.getLinkInControl()');
disp('Nodes under control:'); disp(d.getNodeInControl()');

disp('Time to next event:'); disp(d.getTimeToNextEvent());
d.unload;

d = epanet('Net2.inp');
%% --- Links / Valves ---
disp('==================== Valve Creation and Type ====================');
d.getLinkType()
disp('Adding new pcv valve')
d.addLinkValvePCV('testvalve', '1', '2');
idx = d.getLinkIndex('testvalve');
d.getLinkType()
%% 
disp('Converting new link at index 42 from GPV to PCV');
d.addLinkValveGPV('testvalvegpv', '1', '3');
d.getLinkType(42)
disp('setLinkTypeValvePCV');
d.setLinkTypeValvePCV(42);
d.getLinkType(42)
disp('Valve setup complete.');

%% --- Set Intermediate Vertex ---
disp('==================== Set Vertex ====================');
xCoords = [100.0, 150.0, 200.0];
yCoords = [200.0, 250.0, 300.0];
numVertices = numel(xCoords);

[countBefore] = d.getLinkVerticesCount(1);
disp('Vertex count before setting:'); disp(countBefore);

d.apiENsetvertices(1, xCoords, yCoords, numVertices, d.LibEPANET, d.ph);

[countAfter] = d.getLinkVerticesCount(1);
disp('Vertex count after setting:'); disp(countAfter);

[old] = d.getLinkVertices(1, 1);
x1 = old{1}.x(1);             % x of first vertex
y1 = old{1}.y(1);             % y of first vertex
disp('X and Y of Vertex 1 before manual set:'); disp([x1 ,y1]);

disp('Setting vertex 1 to (300, 800)');
d.setVertex(1, 1, 300.0, 800.0);

[new] = d.getLinkVertices(1, 1);
x2 = new{1}.x(1);             % x of first vertex
y2 = new{1}.y(1);             % y of first vertex
disp('X and Y of Vertex 1 after manual set:'); disp([x2 ,y2]);

%% --- Options & Units ---
disp('==================== Options & Units ====================');
disp('Current pressure units:'); disp(d.getOptionsPressureUnits());
disp('Changing to KPA to show it is not allowed')
d.setOptionsPressureUnitsKPA();
disp('Current pressure units:'); disp(d.getOptionsPressureUnits());
disp('Changing flow units to LPS so we can change to pressure units to KPA');d.setFlowUnitsLPS();
disp('Setting pressure units to KPA ')
d.setOptionsPressureUnitsKPA();
disp('Current pressure units:'); disp(d.getOptionsPressureUnits());
d.setOptionsPressureUnitsMeters(); disp(d.getOptionsPressureUnits());

disp('Status report:'); disp(d.getOptionsStatusReport());
d.setOptionsStatusReport(d.ToolkitConstants.EN_NO_REPORT);
d.setOptionsStatusReportNormal();
disp('Status report after normal:'); disp(d.getOptionsStatusReport());
d.setOptionsStatusReportFull();
disp('Status report after full:'); disp(d.getOptionsStatusReport());

disp('Flow units before change:'); disp(d.getFlowUnits());
disp('Setting flow units to CMS');
d.setFlowUnitsCMS();
disp('Flow units after change:'); disp(d.getFlowUnits());
disp(' ');

d.setOptionsDemandPattern(2);
disp('Demand pattern set to:'); disp(d.getOptionsDemandPattern());

disp('Emitter backflow before:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowDisallowed();
disp('Emitter backflow disallowed:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowAllowed();
disp('Emitter backflow allowed:'); disp(d.getOptionsEmitterBackFlow());

d.unload();

%% Load Net1 with openX 
d = epanet();
disp('==================== Project I/O & Bulk Retrieval ====================');
d.openX('Net1.inp', 'dummy.rpt', 'dummy.out');

disp('Link info FLOW:') 
d.openHydraulicAnalysis();
d.initializeHydraulicAnalysis();

Time_H = [];       % Time steps
LinkFlows = [];    % Link flows
tStep = 1;
stepCount = 0;
maxSteps = 2;

% Run hydraulic analysis
while tStep > 0 && stepCount < maxSteps
    % Run next time step
    t = d.runHydraulicAnalysis();
    
    % Get link flow 
    flows = d.getLinkInfo(d.ToolkitConstants.EN_FLOW);
    
    % Save current time
    Time_H = [Time_H; t];
    
    % Display flows 
    fprintf('Time: %.2f sec\n', t);
    disp('Link flows:'); disp(flows');
    
    tStep = d.nextHydraulicAnalysisStep();
    stepCount = stepCount + 1;
end

d.closeHydraulicAnalysis();

%% --- Statistics ---
disp('==================== Statistics ====================');
disp('Iterations:'); disp(d.getStatisticIterations());
disp('Relative error:'); disp(d.getStatisticRelativeError());
disp('Deficient nodes:'); disp(d.getStatisticDeficientNodes());
disp('Demand reduction:'); disp(d.getStatisticDemandReduction());
disp('Total leakage loss:'); disp(d.getStatisticTotalLeakageLoss());

d.unload;

%% --- Leakage Expansion Demand ---
% See EX28_Leakage 