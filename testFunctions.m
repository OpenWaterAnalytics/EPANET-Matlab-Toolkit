%% EPANET-Matlab Class Test
% 
clc
clear

% Create EPANET object using the INP file
%d=epanet('Net1_Rossman2000.inp');
inpname='Net1_Rossman2000.inp';
version='epanet20012x64';
d=epanet(inpname,version);



%% Controls
Controls=d.getControls
disp('Press any key to continue...')
%pause


%% Counts
NodeCount=d.getNodeCount
NodeTankReservoirCount=d.getNodeTankReservoirCount
LinkCount=d.getLinkCount
PatternCount=d.getPatternCount
CurveCount=d.getCurveCount
ControlRulesCount=d.getControlRulesCount
NodeTankCount=d.getNodeTankCount
NodeReservoirCount=d.getNodeReservoirCount
NodeJunctionsCount=d.getNodeJunctionsCount
LinkPipeCount=d.getLinkPipeCount
LinkPumpCount=d.getLinkPumpCount
LinkValveCount=d.getLinkValveCount
disp('Press any key to continue...')
%pause

%% Errors
for e=[0:6,101:106,110,120,200,202:207,223:224, 240:241, 250:251, 301:309]
    d.getError(e)
end


%%
d.getFlowUnits
d.getLinkNameID
d.getLinkPipeNameID
d.getLinkPumpNameID
d.getLinkValveNameID
d.getLinkIndex
d.getLinkPipeIndex
d.getLinkPumpIndex
d.getLinkValveIndex
d.getLinkNodesIndex
d.getNodesConnectingLinksID
d.getLinkType
d.getLinkTypeIndex
d.getLinkDiameter
d.getLinkLength
d.getLinkRoughnessCoeff
d.getLinkMinorLossCoeff
d.getLinkInitialStatus
d.getLinkInitialSetting
d.getLinkBulkReactionCoeff
d.getLinkWallReactionCoeff
d.getLinkFlows % This is called dynamically in a loop
d.getLinkVelocity
d.getLinkHeadloss
d.getLinkStatus
d.getLinkSettings
d.getLinkPumpEnergy 

d.getNodeNameID
d.getNodeReservoirNameID
d.getNodeJunctionNameID
d.getNodeIndex
d.getNodeReservoirIndex
d.getNodeJunctionIndex
d.getNodeType
d.getNodeTypeIndex
d.getNodeElevations
d.getNodeBaseDemands
d.getNodeDemandPatternIndex
d.getNodeEmitterCoeff
d.getNodeInitialQuality
d.getNodeSourceQuality
d.getNodeSourcePatternIndex
d.getNodeSourceType
d.getNodeTankInitialLevel

d.getNodeActualDemand % This is called dynamically in a loop
d.getNodeActualDemandSensingNodes([1 2 34 25 5])  
d.getNodeHydaulicHead
d.getNodePressure
d.getNodeActualQuality
d.getNodeMassFlowRate
d.getNodeActualQualitySensingNodes([1 2 34 25 5]) 

d.getNodeTankMixiningModel
d.getNodeTankMixingModelCode
d.getNodeTankMixingModelType

d.getNodeTankMinimumWaterVolume 
d.getNodeTankVolumeCurveIndex
d.getNodeTankMinimumWaterLevel
d.getNodeTankMaximumWaterLevel
d.getNodeTankMinimumFraction
d.getNodeTankBulkReactionCoeff
d.getNodeTankIndex
d.getNodeTankNameID
d.getOptionsMaxTrials
d.getOptionsAccuracyValue
d.getOptionsQualityTolerance
d.getOptionsEmitterExponent
d.getOptionsPatternDemandMultiplier
d.getPatternID
d.getPatternIndex
d.getPatternLengths
d.getPattern
d.getPatternValue(1,12)
d.getQualityType
d.getQualityCode
d.getQualityTraceNodeIndex
d.getTimeSimulationDuration
d.getTimeHydraulicStep
d.getTimeQualityStep
d.getTimePatternStep
d.getTimePatternStart
d.getTimeReportingStep
d.getTimeReportingStart

d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.getVersion

%% To check
d.getTimeRuleControlStep % bug: It always returns Zero!
d.getTimeReportingPeriods% Check this
d.getNodeTankMixZoneVolume% bug
d.getNodeTankDiameter% bug: Produces a different diameter
d.getNodeTankInitialWaterVolume% bug: When the initial volume is zero, and then save using ENsaveinpfile, it stores a number different than zero


%% Simulate all
HTS=d.getComputedHydraulicTimeSeries % Also are included: obj.openHydraulicAnalysis;obj.initializeHydraulicAnalysis;obj.runHydraulicAnalysis;obj.nextHydraulicAnalysisStep;obj.closeHydraulicAnalysis;
QTS=d.getComputedQualityTimeSeries% Also are included: obj.openQualityAnalysis;obj.initializeQualityAnalysis;obj.runQualityAnalysis;obj.stepQualityAnalysisTimeLeft;obj.closeQualityAnalysis;
  


d.addPattern('NewPat1')
d.addPattern('NewPat2', [0.8, 1.1, 1.4, 1.1, 0.8, 0.7]); 
d.getPattern

d.getControls
d.setControl(1,1,13,1,11,150); 
d.getControls

d.getLinkDiameter
d.setLinkDiameter(2*d.getLinkDiameter);
d.getLinkDiameter

d.getLinkLength
d.setLinkLength(2*d.getLinkLength)
d.getLinkLength

d.getLinkRoughnessCoeff
d.setLinkRoughnessCoeff(2*d.getLinkRoughnessCoeff)
d.getLinkRoughnessCoeff

d.getLinkMinorLossCoeff
d.setLinkMinorLossCoeff(d.getLinkMinorLossCoeff+1.1)
d.getLinkMinorLossCoeff

d.getLinkInitialStatus
d.setLinkInitialStatus(0*d.getLinkInitialStatus)
d.getLinkInitialStatus

linkstatus=d.getLinkInitialSetting
linkstatus(end)=108;
d.setLinkInitialSetting(linkstatus)
d.getLinkInitialSetting

d.getLinkBulkReactionCoeff
d.setLinkBulkReactionCoeff(d.getLinkBulkReactionCoeff-0.055)
d.getLinkBulkReactionCoeff

d.getLinkWallReactionCoeff
d.setLinkWallReactionCoeff(-1.1*d.getLinkWallReactionCoeff)
d.getLinkWallReactionCoeff

d.getLinkStatus %dynamic
d.setLinkStatus(0*d.getLinkStatus)
d.getLinkStatus 

values = d.getLinkSettings %dynamic
values(end)=111;
d.setLinkSettings(values)
d.getLinkSettings

values = d.getNodeElevations
values(end)=720;
d.setNodeElevations(values)
d.getNodeElevations

values = d.getNodeBaseDemands
values(2)=160;
d.setNodeBaseDemands(values)
d.getNodeBaseDemands

values = d.getNodeDemandPatternIndex
values(2)=0;
d.setNodeDemandPatternIndex(values)
d.getNodeDemandPatternIndex

values = d.getNodeEmitterCoeff
values(2)=0.5;
d.setNodeEmitterCoeff(values)
d.getNodeEmitterCoeff

values = d.getNodeInitialQuality
values(2)=0.6;
d.setNodeInitialQuality(values)
d.getNodeInitialQuality

values = d.getNodeTankInitialLevel
values(end)=100; 
d.setNodeTankLevelInitial(values)  
d.getNodeTankInitialLevel

values = d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
values{end}='MIX2';
d.setNodeTankMixingModelType(values);
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
values = d.getNodeTankMixingModelType 
values{end}='FIFO';
d.setNodeTankMixingModelType(values);
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode
values = d.getNodeTankMixingModelType 
values{end}='LIFO';
d.setNodeTankMixingModelType(values);
d.getNodeTankMixingModelType 
d.getNodeTankMixingModelCode

values = d.getNodeTankDiameter 
values(end)= 60;
d.setNodeTankDiameter(values) 
d.getNodeTankDiameter 

values = d.getNodeTankMinimumWaterLevel
values(end)= 10;
d.setNodeTankMinimumWaterLevel(values)  
d.getNodeTankMinimumWaterLevel

values = d.getNodeTankMinimumWaterVolume
values(end)= 10;
d.setNodeTankMinimumWaterVolume(values) 
d.getNodeTankMinimumWaterVolume

values = d.getNodeTankMaximumWaterLevel
values(end)= 210;
d.setNodeTankMaximumWaterLevel(values) 
d.getNodeTankMaximumWaterLevel

values = d.getNodeTankMinimumFraction
values(end)= 0.5; %takes values 0-1
d.setNodeTankMinimumFraction(values) 
d.getNodeTankMinimumFraction

values = d.getNodeTankBulkReactionCoeff
values(end)= 1; 
d.setNodeTankBulkReactionCoeff(values) 
d.getNodeTankBulkReactionCoeff

d.getNodeSourceType
d.setNodeSourceType(2,'MASS')
d.getNodeSourceType
d.setNodeSourceType(2,'CONCEN')
d.getNodeSourceType
d.setNodeSourceType(2,'SETPOINT')
d.getNodeSourceType
d.setNodeSourceType(2,'FLOWPACED')
d.getNodeSourceType

values = d.getNodeSourceQuality
values(2)=0.5;
d.setNodeSourceQuality(values)
d.getNodeSourceQuality

values = d.getNodeSourcePatternIndex
values(6)=1; 
d.setNodeSourcePatternIndex(values)
d.getNodeSourcePatternIndex

d.getOptionsMaxTrials
d.setOptionsMaxTrials(45)
d.getOptionsMaxTrials

d.getOptionsAccuracyValue
d.setOptionsAccuracyValue(0.015)
d.getOptionsAccuracyValue

d.getOptionsQualityTolerance
d.setOptionsQualityTolerance(0.02)
d.getOptionsQualityTolerance

d.getOptionsEmitterExponent
d.setOptionsEmitterExponent(0.55)
d.getOptionsEmitterExponent

d.getOptionsPatternDemandMultiplier
d.setOptionsPatternDemandMultiplier(1.1)
d.getOptionsPatternDemandMultiplier

d.getTimeSimulationDuration
d.setTimeSimulationDuration(86500)
d.getTimeSimulationDuration

d.getTimeHydraulicStep
d.setTimeHydraulicStep(3500)
d.getTimeHydraulicStep

d.getTimeQualityStep
d.setTimeQualityStep(250)
d.getTimeQualityStep

d.getTimePatternStep
d.setTimePatternStep(7000)
d.getTimePatternStep

d.getTimePatternStart
d.setTimePatternStart(100)
d.getTimePatternStart

d.getTimeReportingStep
d.setTimeReportingStep(3500)
d.getTimeReportingStep

d.getTimeReportingStart
d.setTimeReportingStart(200)
d.getTimeReportingStart

d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.setTimeStatisticsType('MINIMUM')
d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.setTimeStatisticsType('MAXIMUM')
d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.setTimeStatisticsType('RANGE')
d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.setTimeStatisticsType('AVERAGE')
d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.setTimeStatisticsType('NONE')
d.getTimeStatisticsType
d.getTimeStatisticsIndex

d.getTimeRuleControlStep 
d.setTimeRuleControlStep(100) % bug: Does not change Time Rule Control Step
d.getTimeRuleControlStep

d.getPattern 
d.setPattern(1,1:0.01:2)
d.getPattern

values = d.getPattern
values(1,end)=3;
d.setPatternMatrix(values)
d.getPattern

d.getPatternValue(1,10)
d.setPatternValue(1,10,1.2)
d.getPatternValue(1,10)

d.getQualityType
d.getQualityCode
d.setQualityType('none')
d.getQualityCode
d.getQualityType
d.setQualityType('age')
d.getQualityType
d.getQualityCode
d.setQualityType('chem','mg/L')
d.getQualityType
d.getQualityCode
tankid=d.getNodeTankNameID 
d.setQualityType('trace',tankid{1})
d.getQualityType
d.getQualityCode
d.saveInputFile([pwd,'\TEST_INP_TEMP.inp']);


d.writeLineInReportFile('Line-writting testing')
%open('temp.txt'); % bug, write in status report > tmprpt.txt

d.unload
disp('Press any key to continue...')
pause


%% Report Preparation
d=epanet(inpname,version);

% Compute ranges (max - min) 
d.setTimeStatisticsType('RANGE')
d.setTimeStatisticsType('MINIMUM')
% StatisticsType('AVERAGE')
d.setTimeStatisticsType('NONE')
d.setTimeStatisticsType('MAXIMUM')


% Solve hydraulics 
d.solveCompleteHydraulics % solves internally the hydraulics (does not return something)
d.saveHydraulicsOutputReportingFile %creates a BIN file (see EPANET documentation)

% Solve quality
d.solveCompleteQuality

% Define contents of the report
d.setReportFormatReset
d.setReport('FILE TestReport.txt');
d.setReport('PAGESIZE 0')
d.setReport('NODES ALL')%/ALL/node1 node2 
d.setReport('LINKS ALL')%/ALL/link1 link2
d.setReport('PRESSURE PRECISION 1')
d.setReport('PRESSURE ABOVE 20')
d.setReport('STATUS YES')%YES/NO/FULL 
d.setReport('SUMMARY YES')%YES/NO 
d.setReport('MESSAGES YES')%YES/NO 
d.setReport('ENERGY YES')%YES/NO 
%Nodes parameters
%YES/NO/BELOW/ABOVE/PRECISION
d.setReport('ELEVATION YES')
d.setReport('DEMAND YES')
d.setReport('HEAD YES') 
d.setReport('PRESSURE YES') 
d.setReport('QUALITY YES') 
%Links parameters
%BELOW/ABOVE/PRECISION
d.setReport('LENGTH YES')
d.setReport('DIAMETER YES')
d.setReport('FLOW YES')
d.setReport('LENGTH YES')
d.setReport('VELOCITY YES')
d.setReport('HEADLOSS YES')
d.setReport('QUALITY PRECISION 1')
d.setReport('STATUS YES')
d.setReport('SETTING YES')
d.setReport('REACTION YES')
d.setReport('F-FACTOR YES')

%Write the report to file 
d.writeReport
%movefile('TestReport.txt',[pwd,'\RESULTS\','TestReport.txt']);
%open('TestReport.txt');

d.setReportFormatReset
d.setReport('FILE TestReport2.txt'); 
d.setTimeStatisticsType('AVERAGE')
d.setReport('NODES 10')
d.setReport('HEAD YES')
d.setReport('DEMAND NO')
d.setReport('PRESSURE NO')
d.setReport('QUALITY NO')
d.writeReport
%movefile('TestReport2.txt',[pwd,'\RESULTS\','TestReport2.txt']);
%open('TestReport2.txt');

d.setReportFormatReset
d.setReport('FILE TestReport3.txt'); 
d.setReport('NODES ALL')
d.setReport('LINKS ALL')
d.writeReport
%movefile('TestReport3.txt',[pwd,'\RESULTS\','TestReport3.txt']);
%open('TestReport3.txt')

d.setReportFormatReset
d.setReport('FILE TestReport4.txt'); 
d.setReport('STATUS YES') 
d.writeReport
%movefile('TestReport4.txt',[pwd,'\RESULTS\','TestReport4.txt']);
%open('TestReport4.txt')

d.setReportFormatReset
d.setReport('FILE TestReport5.txt'); 
d.setTimeStatisticsType('NONE')
d.setReport('LINKS 10')
d.setReport('LINKS 11')
d.setReport('LINKS 12')
d.setReport('FLOW YES')
d.setReport('HEADLOSS NO') %bug: It shoes Headloss even though it shouldn't
d.setReport('VELOCITY NO')
d.writeReport
%movefile('TestReport5.txt',[pwd,'\RESULTS\','TestReport5.txt']);
%open('TestReport5.txt');

d.setReportFormatReset
d.setReport('FILE TestReport6.txt'); 
d.setTimeStatisticsType('MINIMUM')
d.setReport('NODES ALL')
d.writeReport
%movefile('TestReport6.txt',[pwd,'\RESULTS\','TestReport6.txt']);
%open('TestReport6.txt'); 

d.setReportFormatReset
d.setReport('FILE TestReport7.txt'); 
d.setTimeStatisticsType('NONE')
d.setReport('LINKS ALL')
d.writeReport
%movefile('TestReport7.txt',[pwd,'\RESULTS\','TestReport7.txt']);
%open('TestReport7.txt'); 

d.unload
disp('Press any key to continue...')
pause


%% Create Hydraulics file

d=epanet(inpname,version);

d.solveCompleteHydraulics % Only call this ONLY once (see ENsolveH for more details)
d.saveHydraulicFile([pwd,'\hydraulics.hyd'])
d.useHydraulicFile([pwd,'\hydraulics.hyd'])
d.saveHydraulicsOutputReportingFile

d.unload
disp('Press any key to continue...')
pause

%% Simulation Quality
d=epanet(inpname,version);
d.setQualityType('chem','mg/L')

% Solve Hydraulics (outside the loop)
%d.solveCompleteHydraulics

% or open hydraulics files
d.useHydraulicFile([pwd,'\hydraulics.hyd'])

% Runs Quality Step-by-step
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T=[];Q=[];  
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    P=[P; d.getNodePressure];
    Q=[Q; d.getNodeActualQuality];
    T=[T; t];
    tleft = d.stepQualityAnalysisTimeLeft;
end
d.closeQualityAnalysis;
d.unload
disp('Press any key to continue...')
pause
%tstep=d.nextQualityAnalysisStep; %%% CHECK, DOES NOT SEEM TO CHANGE
%WITH SETTIMEQUALITYSTEP

%% Simulation Hydraulics
d=epanet(inpname,version);
d.setQualityType('chem','mg/L')


% Runs hydraulics Step-by-step
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1; P=[];T=[]; D=[]; H=[];F=[];
while (tstep>0)
    t=d.runHydraulicAnalysis;
    P=[P; d.getNodePressure];
    D=[D; d.getNodeActualDemand];
    H=[H; d.getNodeHydaulicHead];
    F=[F; d.getLinkFlows];
    T=[T; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis

d.unload
disp('Press any key to continue...')
pause

%% Other Functions
% OTHER FUNCTIONS
d=epanet(inpname,version);

NodeCoordinates = d.getCoordinates
d.getInputFileInfo
d.getCurveInfo
d.getLinksInfo
d.getNodesInfo
d.getControlsInfo
d.getFlowUnitsHeadlossFormula




%% Test Plot
d.plot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
disp('Press any key to continue...')
%pause


d.unload
disp('Press any key to continue...')
pause
