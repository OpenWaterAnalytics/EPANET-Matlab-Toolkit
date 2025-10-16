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
disp('Rule enabled 1 after disable:'); disp(d.getRuleEnabled());

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
index = d.addLinkValveGPV('testvalvegpv', '1', '3');
d.getLinkType(index)
disp('setLinkTypeValvePCV');
d.setLinkTypeValvePCV(index);
d.getLinkType(index)
disp('Valve setup complete.');

%% --- Set Intermediate Vertex ---
disp('==================== Set Vertex ====================');
xCoords = [100.0, 150.0, 200.0];
yCoords = [200.0, 250.0, 300.0];
numVertices = numel(xCoords);

index_link = 1;
[countBefore] = d.getLinkVerticesCount(index_link);
disp('Vertex count before setting:'); disp(countBefore);
linkID = d.getLinkNameID(index_link);
d.setLinkVertices(linkID, xCoords, yCoords)

[countAfter] = d.getLinkVerticesCount(index_link);
disp('Vertex count after setting:'); disp(countAfter);

[old] = d.getLinkVertices{1};
x1 = old.x(1);             % x of first vertex
y1 = old.y(1);             % y of first vertex
disp('X and Y of Vertex 1 before manual set:'); disp([x1 ,y1]);

disp('Setting vertex 1 to (300, 800)');
d.setVertex(1, 1, 300.0, 800.0);

[new] = d.getLinkVertices{1};
x2 = new.x(1);             % x of first vertex
y2 = new.y(1);             % y of first vertex
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

disp('Demand pattern:'); disp(d.getOptionsDemandPattern());
d.setOptionsDemandPattern(3);
disp('Demand pattern set to:'); disp(d.getOptionsDemandPattern());

disp('Emitter backflow before:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowDisallowed();
disp('Emitter backflow disallowed:'); disp(d.getOptionsEmitterBackFlow());
d.setOptionsEmitterBackFlowAllowed();
disp('Emitter backflow allowed:'); disp(d.getOptionsEmitterBackFlow());

d.unload();

%% Load Net1
d = epanet('Net1.inp');
disp('==================== Project I/O & Bulk Retrieval ====================');

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
    flows = d.getLinkValues(d.ToolkitConstants.EN_FLOW);
    
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

%% open broken file using ENopenX in ENopen 
d = epanet('Net1_broken.inp');
d.unload


%% --- Leakage Expansion Demand ---
% See EX28_Leakage

%% Test functions where getLinkValue -> getLinkValues
d = epanet('Net1.inp');

d.getLinksInfo

d.runsCompleteSimulation;
disp('getLinkFlows(1):');              disp(d.getLinkFlows(1));
disp('getLinkVelocity(1):');           disp(d.getLinkVelocity(1));
disp('getLinkHeadloss(1):');           disp(d.getLinkHeadloss(1));
disp('getLinkStatus(1):');             disp(d.getLinkStatus(1));
disp('getLinkPumpState(1):');          disp(d.getLinkPumpState(1));
disp('getLinkSettings(1):');           disp(d.getLinkSettings(1));
disp('getLinkEnergy(1):');             disp(d.getLinkEnergy(1));
disp('getLinkActualQuality(1):');      disp(d.getLinkActualQuality(1));
disp('getLinkPumpEfficiency(1):');     disp(d.getLinkPumpEfficiency(1));
disp('getLinkDiameter(1):');           disp(d.getLinkDiameter(1));
disp('getLinkLength(1):');             disp(d.getLinkLength(1));
disp('getLinkQuality(1):');            disp(d.getLinkQuality(1));
disp('getLinkRoughnessCoeff(1):');     disp(d.getLinkRoughnessCoeff(1));
disp('getLinkMinorLossCoeff(1):');     disp(d.getLinkMinorLossCoeff(1));
disp('getLinkInitialStatus(1):');      disp(d.getLinkInitialStatus(1));
disp('getLinkInitialSetting(1):');     disp(d.getLinkInitialSetting(1));
disp('getLinkBulkReactionCoeff(1):');  disp(d.getLinkBulkReactionCoeff(1));
disp('getLinkWallReactionCoeff(1):');  disp(d.getLinkWallReactionCoeff(1));
disp('getLinkLeakArea(1):');           disp(d.getLinkLeakArea(1));
disp('getLinkExpansionProperties(1):'); disp(d.getLinkExpansionProperties(1));
disp('getLinkLeakageRate(1):');        disp(d.getLinkLeakageRate(1));

%% --- all links ---
disp('--- all links ---');
disp('getLinkFlows:');                  disp(d.getLinkFlows');
disp('getLinkVelocity:');               disp(d.getLinkVelocity');
disp('getLinkHeadloss:');               disp(d.getLinkHeadloss');
disp('getLinkStatus:');                 disp(d.getLinkStatus');
disp('getLinkPumpState:');              disp(d.getLinkPumpState');
disp('getLinkSettings:');               disp(d.getLinkSettings');
disp('getLinkEnergy:');                 disp(d.getLinkEnergy');
disp('getLinkActualQuality:');          disp(d.getLinkActualQuality');
disp('getLinkPumpEfficiency:');         disp(d.getLinkPumpEfficiency');
disp('getLinkDiameter:');               disp(d.getLinkDiameter');
disp('getLinkLength:');                 disp(d.getLinkLength');
disp('getLinkQuality:');                disp(d.getLinkQuality');
disp('getLinkRoughnessCoeff:');         disp(d.getLinkRoughnessCoeff');
disp('getLinkMinorLossCoeff:');         disp(d.getLinkMinorLossCoeff');
disp('getLinkInitialStatus:');          disp(d.getLinkInitialStatus');
disp('getLinkInitialSetting:');         disp(d.getLinkInitialSetting');
disp('getLinkBulkReactionCoeff:');      disp(d.getLinkBulkReactionCoeff');
disp('getLinkWallReactionCoeff:');      disp(d.getLinkWallReactionCoeff');
disp('getLinkLeakArea:');               disp(d.getLinkLeakArea');
disp('getLinkExpansionProperties:');    disp(d.getLinkExpansionProperties');
disp('getLinkLeakageRate:');            disp(d.getLinkLeakageRate');

d.unload;
