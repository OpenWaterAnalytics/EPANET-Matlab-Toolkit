%% Testing
fclose all;
clc;
clear all;
clear class;
% unloadlibrary('epanet2')

%  TEST - EPANET
% Input Files
% d=Epanet('Net2_Rossman2000.inp')
% d=Epanet('example.inp')
% d=Epanet('Template.inp')
d=Epanet('Net1_Rossman2000.inp')
% d=Epanet('BWSN1_Ostfeld2008.inp')
d.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Links Pipe Pump & Valves
%PIPE
d.getTimeSimulationDuration
d.setTimeSimulationDuration(100)
d.addPipe('P1','32','10')   
d.getTimeSimulationDuration
d.plot('nodes','yes','fontsize',20)
%PUMP
d.getTimeSimulationDuration
d.setTimeSimulationDuration(2102)
d.addPump('PUMP','23','32','1')   
d.getTimeSimulationDuration
d.plot('nodes','yes')
%VALVES
%PRV OR..
d.addValvePRV('P1','11','22') 
%PSV OR..
d.addValvePSV('P2','32','10') 
d.plot('nodes','yes')
%PBV 
d.addValvePBV('P3','12','23') 
d.plot('nodes','yes')
%FCV
d.addValveFCV('P4','21','9') 
d.plot('nodes','yes')
d.addValveFCV('P4','31','12') 
d.plot('nodes','yes')
%TCV
d.addValveTCV('P5','22','13') 
d.plot('nodes','yes')
%GPV
d.addValveGPV('P6','23','32') %%%%%%%%%%%%ERROR
d.plot('nodes','yes')


% SET UNITS examples
d.setFlowUnitsLPM % Net1.. GPM to LPM
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsGPM % Net1.. LPM to GPM
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsCFS % Net1.. GPM to CFS    
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsMGD % Net1.. CFS to MGD
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsIMGD % Net1..MGD to IMGD
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsAFD % Net1.. IMGD to AFD
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsLPS % Net1.. AFD to LPS
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsLPM % Net1.. LPS to LPM      
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsMLD % Net1.. LPM to MLD
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsCMD % Net1.. MLD to CMD
d.getFlowUnitsHeadlossFormula  
d.setFlowUnitsCMH % Net1.. CMD to CMH
d.getFlowUnitsHeadlossFormula  

% HeadlossFormula
d.setHeadlossDW    
d.getFlowUnitsHeadlossFormula  
d.setHeadlossCM
d.getFlowUnitsHeadlossFormula  
d.setHeadlossHW
d.getFlowUnitsHeadlossFormula  


% Remove - functions
% Nodes
d.getNodesInfo
d.removeNodeID('9') 
d.removeNodeID('2') 

% d.removeNodeID('A') %example.inp
% d.addNewPipe('PP','Source','C');
% d.plot

% d.removeNodeID('C') %example.inp
% d.plot
% d.removeNodeID('D')
% d.plot

d.getTimeSimulationDuration
d.setTimeSimulationDuration(86500)
d.removeNodeID('9') %"Net1_Rossman2000.inp"
d.plot
d.removeNodeID('22') %"Net1_Rossman2000.inp"
d.plot
d.getTimeSimulationDuration

% Links
d.setTimeSimulationDuration(10500)
d.plot('links','yes','nodes','yes')
d.removeLinkID('9') %pump of "Net1_Rossman2000.inp"
d.removeNodeID('9') %"Net1_Rossman2000.inp"
d.getTimeSimulationDuration
d.plot

% Curves
d.setTimeSimulationDuration(22500)
d.removeCurveID('1') % must be removed the pump 9, Warning: Pump 9 refers to undefined curve. Automatos svinei kai to antistoixo Pump, diagrafwntas to curve.
% Warning: Node 9 disconnected. 
d.removeNodeID('9')
d.getTimeSimulationDuration
d.getCurveInfo
d.addNewCurve('new',1995,290)
d.getCurveInfo
d.getTimeSimulationDuration


% Controls
v=d.getControlsInfo
d.setTimeSimulationDuration(500)
% LINK x status IF NODE y ABOVE/BELOW z
d.addControl('10','OPEN','10','ABOVE',100)
d.getTimeSimulationDuration

% LINK x status AT TIME t
d.addControl('10','OPEN','10.00')
% LINK x status AT CLOCKTIME c AM/PM
d.addControl('10','OPEN','12.00','AM')

v=d.getControlsInfo
d.setTimeSimulationDuration(10000)
d.removeControlLinkID(v.linksID{1});
d.getTimeSimulationDuration

v=d.getControlsInfo
d.removeControlNodeID(v.nodesID{1});
d.getTimeSimulationDuration

v=d.getControlsInfo
d.removeControlLinkID(v.linksID{1});
d.getTimeSimulationDuration


% Nodes & Link Info
d.getNodesInfo
d.getLinksInfo
d.getPumpInfo % for the specific Pump curve

% Add node & pipe
d.setTimeSimulationDuration(5000)
d.plot
[x,y]=ginput
d.addJunction('J1',x,y)
d.addPipe('P1','9','J1')
d.plot
d.getTimeSimulationDuration

d.setTimeSimulationDuration(35000)
d.plot
[x,y]=ginput
d.addReservoir('S1',x,y)
d.addPipe('P2','J1','S1')
d.plot
d.getTimeSimulationDuration

d.setTimeSimulationDuration(15000)
d.plot
[x,y]=ginput
d.addTank('T1',x,y)
d.addPipe('P3','32','T1')   
d.plot
d.getTimeSimulationDuration

d.setTimeSimulationDuration(25000)   
d.removeNodeID('T1');
d.plot
d.getTimeSimulationDuration

d.setTimeSimulationDuration(5000)
d.addPipe('MSK','23','32')
d.plot
d.getTimeSimulationDuration





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=d.getComputedHydraulicTimeSeries

d.setTimeQualityStep(3600);
Q = d.getComputedQualityTimeSeries;

d.setReportFormatReset
d.setReport('NODES ALL')
d.setReportStatus('full')
d.writeLineInReportFile('Line-writting testing')
d.getReport

%CONTROLS: EPANET cannot add new controls
d.getControls
d.setControl(1,1,13,1,11,150)
d.getControls


d.getCountNodes
d.getCountTanksReservoirs
d.getCountLinks
d.getCountPatterns
d.getCountCurves
d.getCountControls        

d.getError(101)

d.getFlowUnits 
d.getFlowUnitsHeadlossFormula %read from inp file

d.getLinkID
d.getLinkID([1 3 6])

d.getLinkIndex(d.getLinkID)
d.getLinkIndex('121')
d.getLinkIndex({'121' '31'})

d.getLinkNodes
d.getLinkType

d.getLinkDiameter
d.setLinkDiameter(2*d.getLinkDiameter);
d.getLinkDiameter

d.getLinkLength
d.setLinkLength(2*d.getLinkLength)
d.getLinkLength


d.getLinkRoughness
d.setLinkRoughness(2*d.getLinkRoughness)
d.getLinkRoughness

d.getLinkMinorLossCoeff
d.setLinkMinorLossCoeff(d.getLinkMinorLossCoeff+1)
d.getLinkMinorLossCoeff


d.getLinkInitialStatus
d.setLinkInitialStatus(0*d.getLinkInitialStatus)
d.getLinkInitialStatus

tmp=d.getLinkInitialSettings
tmp(2)=105
d.setLinkInitialSettings(tmp)
tmp(13)=0
d.setLinkInitialSettings(tmp)
d.getLinkInitialSettings

d.getLinkBulkReactionCoeff
d.setLinkBulkReactionCoeff(d.getLinkBulkReactionCoeff-0.055)
d.getLinkBulkReactionCoeff

d.getLinkWallReactionCoeff
d.setLinkWallReactionCoeff(-1.1*d.getLinkWallReactionCoeff)
d.getLinkWallReactionCoeff

d.getLinkStatus %dynamic
d.setLinkStatus(0*d.getLinkStatus)
d.getLinkStatus 


d.getLinkFlows %dynamic
d.getLinkVelocity %dynamic
d.getLinkHeadloss %dynamic
d.getLinkStatus %dynamic
d.getLinkEnergy %dynamic


d.getNodeID
d.getNodeID([1 3 5])

d.getNodeIndex(d.getNodeID)
d.getNodeIndex('12')
d.getNodeIndex('A')

d.getNodeType

d.getNodeActualDemand % dynamic
d.getNodeHydaulicHead % dynamic
d.getNodePressure % dynamic
d.getNodeActualQuality % dynamic
d.getNodeMassFlowRate % dynamic


d.getOptionTrial
d.setOptionTrial(45)
d.getOptionTrial

d.getOptionAccuracy
d.setOptionAccuracy(0.015)
d.getOptionAccuracy

d.getOptionTolerance
d.setOptionTolerance(0.02)
d.getOptionTolerance

d.getOptionEmitterExponent
d.setOptionEmitterExponent(0.55)
d.getOptionEmitterExponent

d.getOptionDemandMult
d.setOptionDemandMult(1.1)
d.getOptionDemandMult

d.addPattern('NewPat1')
d.addPattern('NewPat2', [0.8, 1.1, 1.4, 1.1, 0.8, 0.7])

d.getPatternID
d.getPatternID([2])

d.getPatternIndex
d.getPatternIndex('NewPat1')

d.getPatternLength
d.getPatternLength(1)
d.getPatternLength([1 2])
d.getPatternLength('1')
d.getPatternLength({'1' 'NewPat2'})

d.getPattern % Change this to repeat smaller length patterns
d.setPattern(1,1:0.01:2)
d.getPattern
% 
d.getPatternValue(1,10)
d.setPatternValue(1,10,1.2)
d.getPatternValue(1,10)

d.getQualityType
d.setQualityType('age')
d.getQualityType
d.setQualityType('chem','mg/L')
d.setQualityType('trace','11')
d.getQualityType

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
% 
d.getTimeReportingStep
d.setTimeReportingStep(3500)
d.getTimeReportingStep

d.getTimeReportingStart
d.setTimeReportingStart(200)
d.getTimeReportingStart



d.getTimeStatistics
d.setTimeStatistics('MINIMUM')
d.getTimeStatistics

d.getTimeReportingPeriods

d.getNodeSourceType
d.setNodeSourceType(2,'MASS')
d.getNodeSourceType

d.getTimeRuleControlStep
d.setTimeRuleControlStep(300)
d.getTimeRuleControlStep


d.saveInputFile('Netmodified.inp')
movefile('Netmodified.inp',[pwd,'\Results']);

% d.getPumpType

%%%
d.getCoordinates

d.plot
d.plot('highlightnode',{'A'})
d.plot('highlightnode',{'2'})
d.plot('highlightnode',{'A'},'Nodes','yes')
d.plot('highlightnode',{'10'},'Nodes','yes')
d.plot('highlightlink',{'2','111'})
d.plot('Nodes','yes')
d.plot('Links','yes')
d.plot('Nodes','yes','highlightnode',{'10'})
d.plot('Nodes','yes','highlightnode',{'A'})

values = d.getLinkDiameter
values(1) = 91;
d.setLinkDiameter(values)
d.getLinkDiameter

values = d.getLinkLength
values(1)=91;
d.setLinkLength(values)
d.getLinkLength


values = d.getLinkRoughness
values(1)=91;
d.setLinkRoughness(values)
d.getLinkRoughness

values = d.getLinkMinorLossCoeff
values(1:end)=0.1;
d.setLinkMinorLossCoeff(values)
d.getLinkMinorLossCoeff

values = d.getLinkInitialStatus
values(1)=0;%0 or 1
d.setLinkInitialStatus(values)
d.getLinkInitialStatus

values = d.getLinkInitialSettings
values(2)=105;
d.setLinkInitialSettings(values)
d.getLinkInitialSettings
values(13)=0;
d.setLinkInitialSettings(values)
d.getLinkInitialSettings

values = d.getLinkBulkReactionCoeff
values(2)=-0.55;
d.setLinkBulkReactionCoeff(values)
d.getLinkBulkReactionCoeff

values = d.getLinkWallReactionCoeff
values(2)=-1.1;
d.setLinkWallReactionCoeff(values)
d.getLinkWallReactionCoeff

values = d.getLinkStatus %dynamic
values(2)=1;
d.setLinkStatus(values)
d.getLinkStatus 

values = d.getLinkSettings %dynamic
values(2)=111;
d.setLinkSettings(values)
d.getLinkSettings

values = d.getNodeElevation
values(2)=720;
d.setNodeElevation(values)
d.getNodeElevation

values = d.getNodeBaseDemand
values(2)=160;
d.setNodeBaseDemand(values)
d.getNodeBaseDemand

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

values = d.getTankLevelInitial
values(11)=125;
d.setTankLevelInitial(values)
d.getTankLevelInitial


d.getNodeSourceType
d.setNodeSourceType(2,'MASS')
d.getNodeSourceType

values = d.getNodeSourceQuality
values(2)=0.5;
d.setNodeSourceQuality(values)
d.getNodeSourceQuality

values = d.getNodeSourcePatternIndex
values(6)=1; 
d.setNodeSourcePatternIndex(values)
d.getNodeSourcePatternIndex

values = d.getPattern
values(3,end)=800;
d.setPatternMatrix(values)
d.getPattern

%%%


d.getCurveInfo

d.setQualityType('chem','mg/L')
    
v=d.getComputedQualityTimeSeries('time','demand')
v=d.getComputedQualityTimeSeries('quality')
v=d.getComputedQualityTimeSeries('mass')
v=d.getComputedQualityTimeSeries('qualitySensingNodes',1:d.CountNodes)

v=d.getComputedHydraulicTimeSeries 

% d.solveCompleteHydraulics;
% d.saveHydraulicFile('hydraulics.hyd')
% 
% for i=1:d.CountNodes
%     d.setNodeSourceType(i,'SETPOINT');
% end
%    
% movefile('hydraulics.hyd',[pwd,'\Results']);
% 
% d.useHydraulicFile('hydraulics.hyd');

% TEST - MSX
% d.LoadMSX('Net2_Rossman2000.msx')
d.LoadMSX('example.msx')

% Hydraulic analysis
d.getCountSpeciesMsx     
d.getCountConstantsMsx            
d.getCountParametersMsx           
d.getCountPatternsMsx 

% Patterns
d.addPatternMsx('testpat',[2 .3 .4 6 5 2 4]);
d.getPatternLengthMsx
d.getPatternIndexMsx
d.getPatternsIDMsx                 

d.setPatternMatrixMsx([.1 .2 .5 .2 1 .9]);
d.getPatternMsx 

d.setPatternValueMsx(1,1,2);
d.getPatternMsx 

d.setPatternMsx(1,[1 0.5 0.8 2 1.5]);
d.getPatternMsx 
d.getPatternValueMsx(1,5)%1.5

% Sources
% SourceTypeMsx,SourceLevelMsx,SourcePatternIndexMsx,SourceNodeIDMsx
v=d.getSourcesMsx
node = 1;
spec=1;
type = 0;
level=0.2;
pat = 1;
d.setSourceMsx(node, spec, type, level, pat)
v=d.getSourcesMsx

% Species
d.getSpeciesIndexMsx
d.getSpeciesIDMsx
% d.getSpeciesConcentration(type, index, species)
d.getSpeciesConcentration(0,1,1)

% Constants
d.getConstantValueMsx     
value = [2 10 8];%index[1 2 3]
d.setConstantValueMsx(value);
d.getConstantValueMsx     
d.getConstantsIDMsx                                         
d.getConstantsIndexMsx 

% Parameters
d.getParametersIDMsx              
d.getParametersIndexMsx 

d.getParameterPipeValueMsx        
d.getParameterTankValueMsx        

if d.getCountParametersMsx 
    % d.setParameterPipeValueMsx(pipeIndex,value) 
    d.setParameterPipeValueMsx(1,[1.5 2]) 
    d.getParameterPipeValueMsx{1}        

    d.TankIndex
    d.setParameterTankValueMsx(d.TankIndex(1),100)  
    d.getParameterTankValueMsx{d.TankIndex(1)}  
end

% Initial Quality
values = d.getInitqualLinkValueMsx
nodeIndex=1; speciesIndex=1;
values{nodeIndex}(speciesIndex)=1000;%
d.setInitqualLinkValueMsx(values)     
d.getInitqualLinkValueMsx 

linkIndex=1; speciesIndex=1;
values = d.getInitqualNodeValueMsx
values{linkIndex}(speciesIndex)=1500;%
d.setInitqualNodeValueMsx(values)
d.getInitqualNodeValueMsx   

d.getErrorMsx(501)

d.saveMsxFile('msxsavedtest.msx');                  
% movefile('msxsavedtest.msx',[pwd,'\Results']);

d.getReportMsx   


l = d.getComputedQualityLinkMsx                                  
n = d.getComputedQualityNodeMsx                                                
          

d.unloadMsx 

d.unload


%%%% delete s-files 
    %Delete s files 
a='abcdefghijklmnoqrstuvwxyz';
for i=1:length(a)
    s=sprintf('s%s*',a(i));
    delete(s)
end
for i=1:9
    s=sprintf('s%.f*',i);
    delete(s)
end

rmpath(genpath(pwd));


% TEST - EPANET

% % Simulate all times
% d.solveCompleteHydraulics
% d.solveCompleteQuality
% 
% 
% % Runs Quality Step-by-step
% d.solveCompleteHydraulics
% d.saveHydraulicFile('hydraulics.hyd')
% d.useHydraulicFile('hydraulics.hyd')
% d.saveHydraulicsOutputReportingFile
% d.openQualityAnalysis
% d.initializeQualityAnalysis
% tleft=1; P=[];T=[]; D=[]; H=[]; Q=[]; M=[];
% while (tleft>0)
%     t=d.runQualityAnalysis;
%     P=[P; d.getNodePressure];
%     Q=[Q; d.getNodeActualQuality];
%     T=[T; t];
%     %tstep=d.nextQualityAnalysisStep; %%% CHECK, DOES NOT SEEM TO CHANGE
%     %WITH SETTIMEQUALITYSTEP
%     tleft = d.stepQualityAnalysisTimeLeft;
% end
% d.closeQualityAnalysis;
% 
% 
% % Runs hydraulics Step-by-step
% d.openHydraulicAnalysis;
% d.initializeHydraulicAnalysis;
% tstep=1; P=[];T=[]; D=[]; H=[]; Q=[]; M=[];F=[];
% while (tstep>0)
%     t=d.runHydraulicAnalysis;
%     P=[P; d.getNodePressure];
%     D=[D; d.getNodeActualDemand];
%     H=[H; d.getNodeHydaulicHead];
%     F=[F; d.getLinkFlows];
%     T=[T; t];
%     tstep=d.nextHydraulicAnalysisStep;
% end
% d.closeHydraulicAnalysis
    
