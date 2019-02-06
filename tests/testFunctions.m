%% EPANET-Matlab Toolkit Test Part 1
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc;
clear;
close all;

% Create EPANET object using the INP file
inpname='BWSN_Network_1.inp'; 
% Net1 Net2 Net3 BWSN_Network_1 example

% version='epanet2'; % version dev2.1
% d=epanet(inpname,version);
d=epanet(inpname);
% d=epanet(inpname, 'epanet2');

d.plot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
d.plot('line','no');
d.plot('point','no','linksindex','yes');
d.plot('linksindex','yes','fontsize',8);
d.plot('nodesindex','yes','fontsize',14);

d.getConnectivityMatrix
d.getLibFunctions

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
    d.getLinkPumpTypeCode   % returns type index of all pumps
    d.getLinkPumpType       % returns type name of all pumps: CONSTANT_HORSEPOWER, POWER_FUNCTION, CUSTOM
    d.getLinkPumpHeadCurveIndex % returns index of all pump head curve
    d.getCurveNameID        % returns all curve IDs
    d.getCurveNameID(1)     % returns specific curve ID
    d.addCurve('NewCur1')   % add new curve with ID
    indexCurve=d.addCurve('NewCur2', [1800 200; 1500 400]); % add new curve with points
    d.getCurveNameID
    d.getCurveValue(indexCurve)     % returns all points for specific curve index
    d.getCurveValue(indexCurve,2)   % returns specific point for specific curve index

    d.setCurve(3,[1900 300; 1400 200]) % Change an existing curve 
    d.getCurveValue(indexCurve) 

    d.getCurvesInfo

    len=d.getCurveLengths
    d.getCurveLengths(3)
    d.getCurveLengths('NewCur2')

    d.getCurveIndex
    d.getCurveIndex('NewCur1')

    pointindex=2
    tmppoints=d.getCurveValue(indexCurve,pointindex)
    d.setCurveValue(indexCurve,pointindex,tmppoints+100)
    d.getCurveValue(indexCurve,pointindex)


    bd1=d.getNodeBaseDemands % get an array of the base demands (some nodes may have multiple base demands for different patterns)
    bd1{1}
    node_index=5
    bd1{1}(node_index)= bd1{1}(node_index)+100
    d.setNodeBaseDemands(bd1)
    bd2=d.getNodeBaseDemands
    bd2{1}

    d.getNodeDemandCategoriesNumber
    d.getNodeDemandCategoriesNumber(end)
    d.getNodeDemandCategoriesNumber(5:end)

    % ENgetdemandpattern - Retrieves the index of a demand pattern for a specific demand category of a node
    d.getNodeDemandPatternNameID{1}
    d.getNodeDemandPatternIndex{1}

    % ENgetaveragepatternvalue - Retrieves the average value of a pattern
    d.getPatternAverageValue

    % ENgetstatistic - Retrieves hydraulic simulation statistic
    d.getStatistic

    %[int32, singlePtr, singlePtr] ENgetcoord(int32, singlePtr, singlePtr)

    d.plot;
    nodeCoords=d.getNodeCoordinates; 
    indexNode=1;
    nodeCoords{1}(indexNode)=nodeCoords{1}(indexNode)+10;%X
    nodeCoords{2}(indexNode)=nodeCoords{2}(indexNode)+20;%Y
    d.setNodeCoordinates(nodeCoords)
    d.plot;

    % Quality Info
    d.getQualityInfo

    % Others
    n=d.getComputedHydraulicTimeSeries; % EN_TANKVOLUME - ENgetnodevalue
    n.TankVolume(:,d.NodeTankIndex)

    % EN_STARTTIME  - ENgettimeparam
    d.getTimeStartTime

    % EN_HTIME - ENgettimeparam
    d.getTimeHTime 

    % EN_HALTFLAG - ENgettimeparam
    d.getTimeHaltFlag

    %find the lesser of the hydraulic time step length, or the time to next fill/empty
    d.getTimeNextEvent

    % EN_MAXVOLUME - ENgetnodevalue
    d.getNodeTankMaximumWaterVolume 

    % Curves Info
    d.getCurvesInfo

    % Pump pattern
    d.getLinkPumpPatternNameID % EN_LINKPATTERN - ENgetlinkvalue
    d.getLinkPumpPatternIndex
end
d.unload;

%% Controls
d=epanet(inpname);
Controls=d.getControls
disp('Press any key to continue...')
pause

%% Counts
NodeCount=d.getNodeCount
NodeTankReservoirCount=d.getNodeTankReservoirCount
LinkCount=d.getLinkCount
PatternCount=d.getPatternCount
CurveCount=d.getCurveCount
ControlRulesCount=d.getControlRulesCount
NodeTankCount=d.getNodeTankCount
NodeReservoirCount=d.getNodeReservoirCount
NodeJunctionsCount=d.getNodeJunctionCount
LinkPipeCount=d.getLinkPipeCount
LinkPumpCount=d.getLinkPumpCount
LinkValveCount=d.getLinkValveCount
disp('Press any key to continue...')
pause

%% Errors
for e=[0:6,101:106,110,120,200,202:207,223:224, 240:241, 250:251, 301:309]
    d.getError(e)
end


%% Get Functions
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
d.getLinkEnergy 

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
d.getNodePatternIndex
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
d.getPatternNameID
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
d.getTimeReportingPeriods% Check this
d.getNodeTankMixZoneVolume % OK
d.getNodeTankDiameter% bug: Produces a different diameter % fix in dev.2.1 OK now
d.getNodeTankInitialWaterVolume % OK


%% Simulate all
HTS=d.getComputedHydraulicTimeSeries % Also are included: obj.openHydraulicAnalysis;obj.initializeHydraulicAnalysis;obj.runHydraulicAnalysis;obj.nextHydraulicAnalysisStep;obj.closeHydraulicAnalysis;
QTS=d.getComputedQualityTimeSeries% Also are included: obj.openQualityAnalysis;obj.initializeQualityAnalysis;obj.runQualityAnalysis;obj.stepQualityAnalysisTimeLeft;obj.closeQualityAnalysis;
  
d.addPattern('NewPat1')
d.addPattern('NewPat2', [0.8, 1.1, 1.4, 1.1, 0.8, 0.7]); 
d.getPattern

try
    d.getControls(1)
    d.setControls(1, 'Link 12 OPEN AT TIME 2'); 
    d.getControls(1)
catch e
end

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

linkset=d.getLinkInitialSetting
linkset(end)=108;
if d.LinkValveCount
    linkset(d.LinkValveIndex)=0;
end
d.setLinkInitialSetting(linkset)
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

if nF==1
    values = d.getNodeBaseDemands
    values{1}(2)=160; 
    d.setNodeBaseDemands(values)
    d.getNodeBaseDemands{1}
    
    values = d.getNodeDemandPatternIndex
    values{1}(1)=2;
    d.setNodeDemandPatternIndex(values)
    d.getNodePatternIndex
else
    values = d.getNodeBaseDemands
    values(2)=160; 
    d.setNodeBaseDemands(values)
    d.getNodeBaseDemands
end

values = d.getNodePatternIndex
values(2)=0;
d.setNodeDemandPatternIndex(values)
d.getNodePatternIndex

values = d.getNodeEmitterCoeff
values(2)=0.5;
d.setNodeEmitterCoeff(values)
d.getNodeEmitterCoeff

values = d.getNodeInitialQuality
values(2)=0.6;
d.setNodeInitialQuality(values)
d.getNodeInitialQuality

if d.getNodeTankCount
    values = d.getNodeTankInitialLevel
    values(end)=100; 
    d.setNodeTankInitialLevel(values)  
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
end

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
d.setTimeRuleControlStep(100) 
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
d.setQualityType('trace',d.NodeNameID{2})
d.getQualityType
d.getQualityCode
d.saveInputFile([pwd,'/TEST_INP_TEMP.inp']);

d=epanet(inpname);
% write line in report file
% Solve hydraulics 
d.solveCompleteHydraulics % solves internally the hydraulics (does not return something)
% Solve quality
d.solveCompleteQuality
d.writeLineInReportFile('Line-writting testing!!') % Check this! at the second time is work
d.writeReport
% open([d.BinTempfile(1:end-4),'.txt'])

disp('Press any key to continue...')
pause
d.unload;

%% Report Preparation
d=epanet(inpname);
% Compute ranges (max - min) 
d.setTimeStatisticsType('RANGE')
d.getTimeStatisticsType
d.setTimeStatisticsType('MINIMUM')
d.getTimeStatisticsType
% StatisticsType('AVERAGE')
d.setTimeStatisticsType('NONE')
d.getTimeStatisticsType
d.setTimeStatisticsType('MAXIMUM')
d.getTimeStatisticsType

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
%open('TestReport2.txt');

d.setReportFormatReset
d.setReport('FILE TestReport3.txt'); 
d.setReport('NODES ALL')
d.setReport('LINKS ALL')
d.writeReport
%open('TestReport3.txt')

d.setReportFormatReset
d.setReport('FILE TestReport4.txt'); 
d.setReport('STATUS YES')   % is not appear - check
d.writeReport
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
%open('TestReport5.txt');

d.setReportFormatReset
d.setReport('FILE TestReport6.txt'); 
d.setTimeStatisticsType('MINIMUM')
d.setReport('NODES ALL')
d.writeReport
%open('TestReport6.txt'); 

d.setReportFormatReset
d.setReport('FILE TestReport7.txt'); 
d.setTimeStatisticsType('NONE')
d.setReport('LINKS ALL')
d.writeReport
%open('TestReport7.txt'); 

disp('Press any key to continue...')
pause
d.unload;


%% Create Hydraulics file
d=epanet(inpname);
d.solveCompleteHydraulics % Only call this ONLY once (see ENsolveH for more details)
d.saveHydraulicFile([pwd,'/hydraulics.hyd'])
d.useHydraulicFile([pwd,'/hydraulics.hyd'])
d.saveHydraulicsOutputReportingFile

disp('Press any key to continue...')
pause
d.unload;

%% Simulation Quality
d=epanet(inpname);
d.setQualityType('chem','mg/L')

% Solve Hydraulics (outside the loop)
%d.solveCompleteHydraulics

% or open hydraulics files
d.useHydraulicFile([pwd,'/hydraulics.hyd'])

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
disp('Press any key to continue...')
pause
%tstep=d.nextQualityAnalysisStep; %%% CHECK, DOES NOT SEEM TO CHANGE
%WITH SETTIMEQUALITYSTEP
d.unload;

%% Simulation Hydraulics
d=epanet(inpname);
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
disp('Press any key to continue...')
pause

d=epanet(inpname);
%% Test Plot
d.plot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
disp('Press any key to continue...')
pause

% Other functions
d.getNodeCoordinates
d.unload % delete txt and temp files

delete('TestR*','hydraulics*','*_INP*')
fprintf('Test finished.\n')