%% EPANET Matlab Class Test
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the right of each line number.

disp('Start environment')
fclose all;close all;
clc;
clear all;
clear class;


disp('Create EPANET Class')
% inpname='Net1_Rossman2000'; 
% inpname='ky11_Jolly2013';  
inpname='Net2_Rossman2000';  
d=epanet([inpname,'.inp']);
% d.msx([inpname,'.msx']);

d.plot;

%%EPANET
d.getControls
d.getNodeCount
d.getNodeTankReservoirCount
d.getLinkCount
d.getPatternCount
d.getCurveCount
d.getControlRulesCount
d.getNodeTankCount
d.getNodeReservoirCount
d.getNodeJunctionsCount
d.getLinkPipeCount
d.getLinkPumpCount
d.getLinkValveCount
d.getError(0) %bug at epanet lever, triggers a 251 issue
d.getError(1) 
d.getError(2) 
d.getError(3)
d.getError(4)
d.getError(5)
d.getError(6)
d.getError(101) %similar sense,Insufficient memory 
d.getError(102) %No network data to process
d.getError(103) %Hydraulics solver not initialized 
d.getError(104) %No hydraulic results available 
d.getError(105) %Water quality solver not initialized 
d.getError(106) %No results to report on 
d.getError(110) 
d.getError(120)
d.getError(200)
d.getError(202) %Illegal numeric value in function call 
d.getError(203) %Undefined node in function call 
d.getError(204) %Undefined link in function call 
d.getError(205) %Undefined time pattern in function call 
d.getError(207) %Attempt made to control a check valve 
d.getError(223)
d.getError(224)
d.getError(240) %Undefined source in function call 
d.getError(241) %Undefined control statement in function call 
d.getError(250) %Function argument has invalid format 
d.getError(251) %Illegal parameter code in function call 
d.getError(301)
d.getError(302)
d.getError(303)
d.getError(304)
d.getError(305)
d.getError(306) %Invalid hydraulics file 
d.getError(307)
d.getError(308)
d.getError(309) %Cannot write report to file 
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
d.getLinkFlows%dynamic
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

d.getNodeActualDemand%dynamic
d.getNodeActualDemandSensingNodes([1 2 34 25 5])  
d.getNodeHydaulicHead
d.getNodePressure
d.getNodeActualQuality
d.getNodeMassFlowRate
d.getNodeActualQualitySensingNodes([1 2 34 25 5]) 

d.getNodeTankInitialWaterVolume% bug  & in ENsaveinpfile
d.getNodeTankMixiningModel
d.getNodeTankMixingModelCode
d.getNodeTankMixingModelType
d.getNodeTankMixZoneVolume% bug
d.getNodeTankDiameter% bug
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
d.getTimeRuleControlStep% bug
d.getTimeStatisticsType
d.getTimeStatisticsIndex
d.getTimeReportingPeriods% :s
d.getVersion
d.getComputedHydraulicTimeSeries % obj.openHydraulicAnalysis;obj.initializeHydraulicAnalysis;obj.runHydraulicAnalysis;obj.nextHydraulicAnalysisStep;obj.closeHydraulicAnalysis;
d.getComputedQualityTimeSeries% obj.openQualityAnalysis;obj.initializeQualityAnalysis;obj.runQualityAnalysis;obj.stepQualityAnalysisTimeLeft;obj.closeQualityAnalysis;
% % d.solveCompleteHydraulics  
% % d.solveCompleteQuality      

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
d.setLinkMinorLossCoeff(d.getLinkMinorLossCoeff+1)
d.getLinkMinorLossCoeff

d.getLinkInitialStatus
d.setLinkInitialStatus(0*d.getLinkInitialStatus)
d.getLinkInitialStatus

tmp=d.getLinkInitialSetting
tmp(end)=108
d.setLinkInitialSetting(tmp)
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
d.setNodeTankDiameter(values) % bug
d.getNodeTankDiameter % bug

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
values(end)= 1; % for 2 not work..
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

d.getTimeRuleControlStep % bug
d.setTimeRuleControlStep(100) % bug
d.getTimeRuleControlStep

d.getPattern % Change this to repeat smaller length patterns
d.setPattern(1,1:0.01:2)
d.getPattern

values = d.getPattern
values(1,end)=800;
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
a=d.getNodeTankNameID 
d.setQualityType('trace',a{1})
d.getQualityType
d.getQualityCode
d.saveInputFile([pwd,'\RESULTS\','TestInpFile.inp']);
% have been applied in the end of the function
% % d.closeHydraulicAnalysis
% % d.closeQualityAnalysis
% % d.saveHydraulicFile 
% % d.useHydraulicFile 
% % d.initializeHydraulicAnalysis
% % d.initializeQualityAnalysis
% % d.tstep = nextHydraulicAnalysisStep
% % d.tstep = nextQualityAnalysisStep
% % d.openHydraulicAnalysis
% % d.openQualityAnalysis
% % d.tstep = runHydraulicAnalysis
% % d.tstep = runQualityAnalysis
% % tleft=d.stepQualityAnalysisTimeLeft

d=epanet([inpname,'.inp']);

d.writeLineInReportFile('Line-writting testing')
open('temprpt.txt'); % bug, write in status report > tmprpt.txt

% Compute ranges (max - min) 
% d.setTimeStatisticsType('RANGE')
d.setTimeStatisticsType('MAXIMUM')
% d.setTimeStatisticsType('MINIMUM')
% StatisticsType('AVERAGE')
% d.setTimeStatisticsType('NONE')

% Solve hydraulics 
d.solveCompleteHydraulics
d.saveHydraulicsOutputReportingFile

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
movefile('TestReport.txt',[pwd,'\RESULTS\','TestReport.txt']);
open('TestReport.txt');

d.setReportFormatReset
d.setReport('FILE TestReport2.txt'); 
d.setTimeStatisticsType('AVERAGE')
d.setReport('NODES 10')
d.setReport('HEAD YES')
d.setReport('DEMAND NO')
d.setReport('PRESSURE NO')
d.setReport('QUALITY NO')
d.writeReport
movefile('TestReport2.txt',[pwd,'\RESULTS\','TestReport2.txt']);
open('TestReport2.txt');

d.setReportFormatReset
d.setReport('FILE TestReport3.txt'); 
d.setReport('NODES ALL')
d.setReport('LINKS ALL')
d.writeReport
movefile('TestReport3.txt',[pwd,'\RESULTS\','TestReport3.txt']);
open('TestReport3.txt')

d.setReportFormatReset
d.setReport('FILE TestReport4.txt'); 
d.setReport('STATUS YES') 
d.writeReport
movefile('TestReport4.txt',[pwd,'\RESULTS\','TestReport4.txt']);
open('TestReport4.txt')

d.setReportFormatReset
d.setReport('FILE TestReport5.txt'); 
d.setTimeStatisticsType('NONE')
d.setReport('LINKS 10')
d.setReport('LINKS 11')
d.setReport('LINKS 12')
d.setReport('FLOW YES')
d.setReport('HEADLOSS NO') %bug
d.setReport('VELOCITY NO')
d.writeReport
movefile('TestReport5.txt',[pwd,'\RESULTS\','TestReport5.txt']);
open('TestReport5.txt');

d.setReportFormatReset
d.setReport('FILE TestReport6.txt'); 
d.setTimeStatisticsType('MINIMUM')
d.setReport('NODES ALL')
d.writeReport
movefile('TestReport6.txt',[pwd,'\RESULTS\','TestReport6.txt']);
open('TestReport6.txt'); 

d.setReportFormatReset
d.setReport('FILE TestReport7.txt'); 
d.setTimeStatisticsType('NONE')
d.setReport('LINKS ALL')
d.writeReport
movefile('TestReport7.txt',[pwd,'\RESULTS\','TestReport7.txt']);
open('TestReport7.txt'); 


d.unload 

%%%%%%%%%%%%%%%%%%%%%%Solve Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose all;
clc;
clear all;
clear class;
% Input Files
inpname='Net1_Rossman2000';
d=epanet([inpname,'.inp']);

% Simulate all times
d.solveCompleteHydraulics
d.solveCompleteQuality

d.setQualityType('chem','mg/L')
d.getQualityType

% Runs Quality Step-by-step
d.solveCompleteHydraulics
d.saveHydraulicFile([pwd,'\RESULTS\','hydraulics.hyd'])
d.useHydraulicFile([pwd,'\RESULTS\','hydraulics.hyd'])
d.saveHydraulicsOutputReportingFile
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T=[];Q=[];  
while (tleft>0)
    t=d.runQualityAnalysis;
    P=[P; d.getNodePressure];
    Q=[Q; d.getNodeActualQuality];
    T=[T; t];
    %tstep=d.nextQualityAnalysisStep; %%% CHECK, DOES NOT SEEM TO CHANGE
    %WITH SETTIMEQUALITYSTEP
    tleft = d.stepQualityAnalysisTimeLeft;
end
d.closeQualityAnalysis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose all;
clc;
clear all;
clear class;
% Input Files
inpname='Net1_Rossman2000';
d=epanet([inpname,'.inp']);

% Simulate all times
d.solveCompleteHydraulics
d.solveCompleteQuality

d.setQualityType('chem','mg/L')
d.getQualityType

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MSX
fclose all;
clc;
clear all;
clear class;
inpname='example';
% inpname='Net2_Rossman2000';

d=epanet([inpname,'.inp']);
d.msx([inpname,'.msx'])
d
d.getMsxEquationsTerms
d.getMsxEquationsPipes
d.getMsxEquationsTanks
d.getMsxSpeciesCount
d.getMsxConstantsCount
d.getMsxParametersCount
d.getMsxPatternsCount
d.getMsxSpeciesNameID
d.getMsxSpeciesType
d.getMsxSpeciesUnits
d.getMsxSpeciesATOL  
d.getMsxSpeciesRTOL
d.getMsxSpeciesIndex
d.getMsxConstantsNameID
d.getMsxConstantsValue
d.getMsxConstantsIndex
d.getMsxParametersNameID
d.getMsxParametersIndex
d.getMsxParametersTanksValue
d.getMsxParametersPipesValue
d.getMsxPatternsNameID
d.getMsxPatternsIndex
d.getMsxPatternsLengths
d.getMsxNodeInitqualValue
d.getMsxLinkInitqualValue
d.getMsxSources
d.getMsxSourceType
d.getMsxSourceLevel
d.getMsxSourcePatternIndex
d.getMsxPattern %Mass flow rate per minute of a chemical source
d.getMsxPatternValue(1,5) %Mass flow rate per minute of a chemical source
% % d.getMsxSpeciesConcentration
d.setTimeHydraulicStep(3600)
d.setTimeQualityStep(3600)
d.getMsxComputedQualityNode
d.MsxPlotConcentrationSpeciesOfNodes
d.getMsxComputedQualityLink
d.MsxPlotConcentrationSpeciesOfLinks
d.getMsxError(0)   % bug, no error
d.getMsxError(200) % bug, cannot read EPANET-MSX file
d.getMsxError(501)
d.getMsxError(502)
d.getMsxError(503)
d.getMsxError(504)
d.getMsxError(505)
d.getMsxError(506)
d.getMsxError(507)
d.getMsxError(508)
d.getMsxError(509)
d.getMsxError(510) % bug, could not open algebraic "equation" solver.
d.getMsxError(511)
d.getMsxError(512)
d.getMsxError(513)
d.getMsxError(514)
d.getMsxError(515)
d.getMsxError(516)
d.getMsxError(517)
d.getMsxError(518)
d.getMsxError(519)
d.getMsxError(520)
d.getMsxError(521)
d.getMsxError(522)
d.getMsxError(523)
d.getMsxError(524)

% Solve for hydraulics & water quality
d.MsxSolveCompleteHydraulics
d.MsxSolveCompleteQuality
% Write results to the “TestMsxReport” file
d.MsxWriteReport %a specific water quality report file is named in the [REPORT] section of the MSX input file. %BUG
copyfile([pwd,'\RESULTS\','temprpt.txt'],[pwd,'\RESULTS\','TestMsxReport.txt']);
open('TestMsxReport.txt');
%or
copyfile([pwd,'\LIBRARIES\','epanetmsx.exe'],[pwd,'\RESULTS\','epanetmsx.exe']);
copyfile([pwd,'\LIBRARIES\','epanetmsx.dll'],[pwd,'\RESULTS\','epanetmsx.dll']);
copyfile([pwd,'\LIBRARIES\','epanet2.dll'],[pwd,'\RESULTS\','epanet2.dll']);
fid = fopen('ReportMsx.bat','w');
r = sprintf('cd RESULTS \nepanetmsx %s %s %s','temp.inp','temp.msx','temp.txt'); 
fprintf(fid,'%s \n',r);fclose all;
!ReportMsx.bat
movefile('ReportMsx.bat',[pwd,'\RESULTS\','ReportMsx.bat']);
copyfile([pwd,'\RESULTS\','temp.txt'],[pwd,'\RESULTS\','TestMsxReport2.txt']);
open('TestMsxReport2.txt')
    
d.MsxAddPattern('testpat',[2 .3 .4 6 5 2 4]);
d.getMsxPatternsNameID
d.getMsxPatternsIndex
d.getMsxPatternsLengths  

v=d.getMsxSources
node = 1;
spec=1;
type = 0;
level=0.2;
pat = 1;
d.setMsxSources(node, spec, type, level, pat)
v=d.getMsxSources

d.getMsxConstantsValue     
value = [2 10 8];%index[1 2 3]
d.setMsxConstantsValue(value);
d.getMsxConstantsNameID
d.getMsxConstantsValue
d.getMsxConstantsIndex

d.getMsxParametersTanksValue
d.getMsxParametersPipesValue      

if d.getMsxParametersCount 
    % d.setMsxParametersPipesValue(pipeIndex,value) 
    d.setMsxParametersPipesValue(1,[1.5 2]) 
    d.getMsxParametersPipesValue{1}        

    a=d.getNodeTankIndex
    d.setMsxParametersTanksValue(a(1),100)  
    d.getMsxParametersTanksValue{a(1)}  
end

values = d.getMsxLinkInitqualValue
nodeIndex=1; speciesIndex=1;
values{nodeIndex}(speciesIndex)=1000;%
d.setMsxLinkInitqualValue(values)     
d.getMsxLinkInitqualValue 

linkIndex=1; speciesIndex=1;
values = d.getMsxNodeInitqualValue
values{linkIndex}(speciesIndex)=1500;%
d.setMsxNodeInitqualValue(values)
d.getMsxNodeInitqualValue   

d.setMsxPatternMatrix([.1 .2 .5 .2 1 .9]);
d.getMsxPattern

d.setMsxPatternValue(1,1,2);
d.getMsxPattern 

d.setMsxPattern(1,[1 0.5 0.8 2 1.5]);
d.getMsxPattern 
d.getMsxPatternValue(1,5) 

d.MsxSaveFile([pwd,'\RESULTS\','testMsx.msx']);                                                               
          
d.MsxSaveQualityFile([pwd,'\RESULTS\','testMsxQuality.bin'])

d.saveHydraulicsOutputReportingFile
d.saveHydraulicFile([pwd,'\RESULTS\','testMsxHydraulics.hyd'])

d.MsxUseHydraulicFile([pwd,'\RESULTS\','testMsxHydraulics.msx'])

% % MsxInitializeQualityAnalysis
% % MsxStepQualityAnalysisTimeLeft

d.MsxUnload 
d.unload

% OTHER FUNCTIONS
inpname='Net1_Rossman2000';
d=epanet([inpname,'.inp']);

NodeCoordinates = d.getCoordinates
d.getInputFileInfo
d.getCurveInfo
d.getLinksInfo
d.getNodesInfo
d.getControlsInfo
d.getFlowUnitsHeadlossFormula


%Delete s files 
a='abcdefghijklmnopqrstuvwxyz';
for i=1:length(a)
    s=sprintf('s%s*',a(i));
    delete(s)
end
for i=1:9
    s=sprintf('s%.f*',i);
    delete(s)
end
rmpath(genpath(pwd));