%% EPANET-Matlab Toolkit Test Part 5
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc;
clear;
close all;
clear class;

% Create EPANET object using the INP file
inpname='Net1.inp';
%Net1 Net2 Net3 net2-cl2 BWSN_Network_1
tic;d=epanet(inpname);toc
% d=epanet(inpname, 'epanet2');
d.addPattern('NewPat2', [0.8, 1.1, 1.4, 1.1, 0.8, 0.7]); 
d.BinUpdateClass; % must be run if use Bin functions

if d.Errcode
    return; 
end

%% GET, SET FLOW UNITS
% SET Flow Units 
d.getBinOptionsInfo
d.getBinNodeCoordinates

errcode=d.setBinFlowUnitsLPM % Net1.. GPM to LPM
d.getFlowUnits
tic;d.getComputedHydraulicTimeSeries
toc
tic;d.getBinComputedAllParameters
toc
% similar
errcode=d.setBinFlowUnitsGPM % Net1.. LPM to GPM
d.getFlowUnits
d.getComputedHydraulicTimeSeries  
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsCFS % Net1.. GPM to CFS 
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsMGD % Net1.. CFS to MGD
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsIMGD% Net1.MGD to IMGD
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsAFD % Net1.. IMGD to AFD
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsLPS % Net1.. AFD to LPS
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsLPM % Net1.. LPS to LPM    
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters 

errcode=d.setBinFlowUnitsMLD % Net1.. LPM to MLD
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsCMD % Net1.. MLD to CMD
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinFlowUnitsCMH % Net1.. CMD to CMH
d.getFlowUnits
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

disp('Press any key to continue...')
pause

% HeadlossFormula
errcode=d.setBinHeadlossDW
d.getBinOptionsInfo
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters 

errcode=d.setBinHeadlossCM
d.getBinOptionsInfo
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters  

errcode=d.setBinHeadlossHW
d.getBinOptionsInfo
d.getComputedHydraulicTimeSeries
d.getBinComputedAllParameters 

errcode=d.setBinFlowUnitsGPM;
d.getFlowUnits
disp('Press any key to continue...')
pause



%% GET, SET PIPES PARAMETERS
% Pipe
% case 1
tic;
Diameters=d.BinLinkPipeDiameters;
Diameters(end)=300;
errcode=d.setBinLinkPipeDiameters(Diameters);
d.getLinkDiameter

lengths=d.BinLinkPipeLengths;
lengths(end)=1111111;
errcode=d.setBinLinkPipeLengths(lengths);
d.getLinkLength

roughness=d.BinLinkPipeRoughness;
roughness(end)=222;
errcode=d.setBinLinkPipeRoughness(roughness);
d.getLinkRoughnessCoeff

minorloss=d.BinLinkPipeMinorLoss;
minorloss(end)=10;
errcode=d.setBinLinkPipeMinorLoss(minorloss);
d.getLinkMinorLossCoeff

status=d.BinLinkPipeStatus;
status{end}='closed';
errcode=d.setBinLinkPipeStatus(status);
d.getLinkInitialStatus
case1=toc
disp('Press any key to continue...')
pause


% case 2
tic
Diameters=d.BinLinkPipeDiameters;
Diameters(end)=300*10;

lengths=d.BinLinkPipeLengths;
lengths(end)=1111111*10;

roughness=d.BinLinkPipeRoughness;
roughness(end)=222*20;

errcode=d.setBinLinkPipesParameters('diameter',Diameters,'length',lengths,'roughness',roughness);
d.getLinkDiameter
d.getLinkLength
d.getLinkRoughnessCoeff

minorloss=d.BinLinkPipeMinorLoss;
minorloss(end)=0;

status=d.BinLinkPipeStatus;
status{end}='Open';

errcode=d.setBinLinkPipesParameters('minorloss',minorloss,'status',status);
d.getLinkMinorLossCoeff
d.getLinkInitialStatus
case2=toc
disp('Press any key to continue...')
pause

d.BinUpdateClass



%% GET, SET JUNCTIONS PARAMETERS
% Node parameters
d.getBinNodesInfo
% case 1
toc
elevations=d.BinNodeJunctionElevation;
elevations(end)=20;
errcode=d.setBinNodeJunctionElevation(elevations);
d.getNodeElevations

basedemands=d.BinNodeJunctionsBaseDemands;
basedemands(end)=323;
errcode=d.setBinNodeJunctionsBaseDemands(basedemands);
d.getNodeBaseDemands

patterns=d.BinNodeJunDemandPatternNameID;
% patternsid=d.BinPatternNameID;
errcode=d.addBinPattern('new',1:0.1:2);
patterns{1}='new';
errcode=d.setBinNodeJunDemandPatternNameID(patterns); %No for Reservoirs (use setBinNodeResDemandPatternNameID)
d.getPatternNameID(d.getNodePatternIndex)
case1node=toc
disp('Press any key to continue...')
pause

% case 2
tic
elevations=d.BinNodeJunctionElevation;
elevations(end)=2002;
basedemands=d.BinNodeJunctionsBaseDemands;
basedemands(end)=325;
patterns=d.getPatternNameID(d.getNodePatternIndex);
errcode=d.addBinPattern('new1',1:0.1:2);
patterns{2}='new1';
d.setBinNodeJunctionsParameters('elevation',elevations,'basedemand',basedemands,'demandpattern',patterns);
d.getNodeElevations
d.getNodeBaseDemands
d.getPatternNameID(d.getNodePatternIndex)
case2node=toc
disp('Press any key to continue...')
pause


%% GET, SET RESERVOIRS PARAMETERS
if d.NodeReservoirCount
    tic
    elevationsReservoirs=d.BinNodeReservoirElevation;
    elevationsReservoirs(end)=15;
    errcode=d.setBinNodeReservoirElevation(elevationsReservoirs)
    d.getNodeElevations

    patres=d.BinNodeResDemandPatternNameID;
    patres{end}='new';
    errcode=d.setBinNodeResDemandPatternNameID(patres);
    d.getPatternNameID(d.getNodePatternIndex)
    caseres1=toc
    disp('Press any key to continue...')
    pause

    d.BinUpdateClass
    tic
    elevationsReservoirs=d.BinNodeReservoirElevation;
    elevationsReservoirs(end)=190;
    patres=d.BinNodeResDemandPatternNameID;
    errcode=d.addBinPattern('pat2',1:0.1:2);
    patres{end}='pat2';
    errcode=d.setBinNodeReservoirParameters('elevation',elevationsReservoirs,'pattern',patres);
    d.getNodeElevations
    d.getPatternNameID(d.getNodePatternIndex)
    caseres2=toc
    disp('Press any key to continue...')
    pause
end

%% GET, SET TANKS PARAMETERS
if strcmp(inpname,'Net1.inp')
if d.BinNodeTankCount
    tic
    tankElevation=d.BinNodeTankElevation;
    tankElevation(end)= tankElevation(end)+10;
    errcode=d.setBinNodeTankElevation(tankElevation);
    d.getNodeElevations
    
    tankInitlevel=d.BinNodeTankInitialLevel;
    tankInitlevel(end)=110;
    errcode=d.setBinNodeTankInitialLevel(tankInitlevel);
    d.getNodeTankInitialLevel
    
    d=epanet(inpname);d.BinUpdateClass
    tankMinlevel=d.BinNodeTankMinimumWaterLevel;
    tankMinlevel(end)=tankMinlevel(end)+5;%5
    errcode=d.setBinNodeTankMinimumWaterLevel(tankMinlevel);
    d.getNodeTankMinimumWaterLevel	
    
    tankMaxlevel=d.BinNodeTankMaximumWaterLevel;
    tankMaxlevel(end)=300;
    errcode=d.setBinNodeTankMaximumWaterLevel(tankMaxlevel);    
    d.getNodeTankMaximumWaterLevel

    tankDiameter=d.BinNodeTankDiameter;
    tankDiameter(end)=60;
    errcode=d.setBinNodeTankDiameter(tankDiameter); %bug
    d.getNodeTankDiameter
    
    tankMinvol=d.BinNodeTankMinimumWaterVolume;
    tankMinvol(end)=10;
    errcode=d.setBinNodeTankMinimumWaterVolume(tankMinvol);
    d.getNodeTankMinimumWaterVolume
    case1tank=toc
    
    tic
    tankElevation=d.BinNodeTankElevation;
    tankElevation(end)=100;
    tankInitlevel=d.BinNodeTankInitialLevel;
    tankInitlevel(end)=550;
    tankMinlevel=d.BinNodeTankMinimumWaterLevel;
    tankMinlevel(end)=250;    
    tankMaxlevel=d.BinNodeTankMaximumWaterLevel;
    tankMaxlevel(end)=3000;
    tankDiameter=d.BinNodeTankDiameter;
    tankDiameter(end)=1000;
    tankMinvol=d.BinNodeTankMinimumWaterVolume;
    tankMinvol(end)=100;
errcode=d.setBinNodeTankParameters('elevation',tankElevation,'initlevel',tankInitlevel,'minlevel',tankMinlevel,'maxlevel',tankMaxlevel,'diameter',tankDiameter,'minvol',tankMinvol);
    d.getNodeElevations
    d.getNodeTankInitialLevel
    d.getNodeTankMinimumWaterVolume
    d.getNodeTankDiameter
    d.getNodeTankMaximumWaterLevel
    d.getNodeTankMinimumWaterLevel
    case2tank=toc     
    
    tankMixModel=d.getNodeTankMixingModelType
    tankMixModel{end}='2Comp'; % Constants for mixing models: 'MIX1','2Comp', 'FIFO','LIFO' (2Comp=MIX2)

    tankMixFraction=d.getNodeTankMinimumFraction
    tankMixFraction(end)=2;
    errcode=d.setBinNodeTankParameters('mixmodel',tankMixModel,'mixfraction',tankMixFraction); % Bug on mixfraction
    d.getNodeTankMinimumFraction
    d.getNodeTankMixingModelType
    disp('Press any key to continue...')
    pause
end



%% PLOTS
d = epanet(inpname);
sim=d.getBinComputedAllParameters

d.Binplot('nodes','yes','links','yes','highlightnode',{'10','11'},'colornode',{'r','k'},'highlightlink',{'10'},'fontsize',8);
nodeSet1={'10','11','22'};
nodeSet2={'21','23','31'};
colorNodeSet1={'r','r','r'};
colorNodeSet2={'g','g','g'};
linkSet1={'111','122','121'};
linkSet2={'110','12','113'};
colorLinkSet1={'k','k','k'};
colorLinkSet2={'y','y','y'};
d.Binplot('nodes','yes','links','yes','highlightnode',[nodeSet1 nodeSet2],'colornode',[colorNodeSet1 colorNodeSet2],...
'highlightlink',[linkSet1 linkSet2],'colorlink',[colorLinkSet1 colorLinkSet2])
d.Binplot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
d.Binplot('nodes','yes','links','yes','highlightnode',{'10','11'},'colornode',{'r','k'},'highlightlink',{'10'},'fontsize',8);
d.Binplot('links','yes','highlightnode',{'10','11'},'colornode',{'r','k'},'highlightlink',{'10'},'fontsize',8,'point','no');
disp('Press any key to continue...')
pause

%% ADD, REMOVE NODES & LINKS
d.Binplot('nodes','yes');
% addBinJunction
% arguments: newNodeID,X,Y,ToNodeID,newElevation,newBaseDemand,newDemandPattern,newPipeID,
% newLength,newDiameter,newRoughness,Code
% Code - Constants for valves: 'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'
fprintf('Add nodes and links. Select on plot.\n')
% PIPE
% Add Junction + pipe
newID='J1';
[x,y]=ginput(1);
ToNodeID='10'; 
newElevation=500; %ft
newBaseDemand=0;
newDemandPattern='1';
newPipeID='P1';
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='PIPE';
errcode=d.addBinJunction(newID,x,y,newElevation,newBaseDemand,newDemandPattern,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);

% [errcode]=addBinPipe(newLink,fromNode,toNode,newLength,newDiameter,newRoughness)
[errcode]=d.addBinPipe('P2','J1','31',1000,10,100);
d.Binplot('nodes','yes','links','yes');
d.plot;


% Add Reservoir + pipe
newID='S1';
[x,y]=ginput(1);
newElevation=500; %ft
ToNodeID='10'; 
newPipeID='P3';
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='PIPE';
[errcode]=d.addBinReservoir('S1',x,y,newElevation,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);
d.plot('nodes','yes','links','yes');

% Add Tank + pipe
[x,y]=ginput(1);
MaxLevel=20;
Diameter=50;
Initlevel=10;
newElevation=500;
initqual=0;
MinLevel=0;
MinVol=0;
newPipeID='P4';
ToNodeID='32'; 
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='PIPE';
errcode=d.addBinTank('T1',x,y,MaxLevel,Diameter,Initlevel,newElevation,initqual,MinLevel,MinVol,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);
d.plot('nodes','yes','links','yes');

% Remove Tank
errcode=d.removeBinNodeID('T1');
d.plot('nodes','yes','links','yes');

errcode=d.addBinPipe('pp1','23','32',1000,10,100);
d.plot('nodes','yes','highlightlink',{'pp1'},'fontsize',8);
disp('Press any key to continue...')
%pause

% PUMP
% Add junction + pump
d.Binplot;
newID='J2';
[x,y]=ginput(1);
ToNodeID='J1'; 
newElevation=500; %ft
newBaseDemand=0;
newDemandPattern='1';
newPumpID='PU1';
Code='PUMP';
newCurveIDofPump='C-1'; 
newCurveXvalue=1500;
newCurveYvalue=250;
newCurveType='PUMP'; % PUMP, EFFICIENCY, VOLUME, HEADLOSS            
errcode=d.addBinJunction(newID,x,y,newElevation,newBaseDemand,newDemandPattern,newPumpID,...
ToNodeID,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType,'PUMP');
d.Binplot('nodes','yes','highlightlink',{'PU1'},'fontsize',10);

% Add Reservoir + pump
newID='S2';
[x,y]=ginput(1);
newElevation=500; %ft
ToNodeID='11'; 
newPumpID='PU2';
Code='PUMP';
newCurveIDofPump='C-1n'; 
newCurveXvalue=1500;
newCurveYvalue=250;
[errcode]=d.addBinReservoir(newID,x,y,newElevation,newPumpID,...
ToNodeID,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType,Code);
d.Binplot('nodes','yes','highlightlink',{'PU2'},'fontsize',10);
% Remove reservoir S1 AND S2
[errcode]=d.removeBinNodeID('S1');
[errcode]=d.removeBinNodeID('S2');
d.plot('nodes','yes','links','yes');

% Add Tank + pump
newID='T2';
[x,y]=ginput(1);
MaxLevel=20;
Diameter=50;
Initlevel=10;
newElevation=500;
initqual=0;
MinLevel=0;
MinVol=0;
newPumpID='PU3';
ToNodeID='32'; 
Code='PUMP';
newCurveIDofPump='C-1n2'; 
newCurveXvalue=[1500 1800 2000];
newCurveYvalue=[250 200 0];
errcode=d.addBinTank(newID,x,y,MaxLevel,Diameter,Initlevel,newElevation,initqual,MinLevel,MinVol,newPumpID,...
ToNodeID,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType,Code);
d.plot('nodes','yes','links','yes');

% add pump and add curve
errcode=d.addBinPump('PUMP1','21','32','C-2',1955,250,'PUMP') % CurveType: PUMP, EFFICIENCY, VOLUME, HEADLOSS            
power=18;
errcode=d.addBinPump('PUMP2','12','23',power) % POWER PUMP        
d.plot('nodes','yes','highlightlink',{'PUMP1','PUMP2'},'fontsize',10);


% VALVE ['PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV']
% Add junction + PRV (any valve)
d.Binplot;
newID='J3';
[x,y]=ginput(1);
ToNodeID='10'; 
newElevation=500; %ft
newBaseDemand=0;
newDemandPattern='1';
newValveID='V1prv'; 
newValveDiameter=100; 
Code='PRV'; 
newValveSetting=15; 
errcode=d.addBinJunction(newID,x,y,newElevation,newBaseDemand,newDemandPattern,newValveID,...
ToNodeID,newValveDiameter,newValveSetting,Code);
d.Binplot('nodes','yes','highlightlink',{'V1prv'},'fontsize',10);

% Add Reservoir + PSV
newID='S3';
[x,y]=ginput(1);
newElevation=500; %ft
ToNodeID='11'; 
newValveID='V2psv'; 
newValveDiameter=100; 
Code='PSV'; 
newValveSetting=15; 
[errcode]=d.addBinReservoir(newID,x,y,newElevation,newValveID,...
ToNodeID,newValveDiameter,newValveSetting,Code);
d.Binplot('nodes','yes','highlightlink',{'V2psv'},'fontsize',10);
if errcode
    [errcode]=d.removeBinNodeID('S3');
end

% Add Tank + PBV
newID='T3';
[x,y]=ginput(1);
MaxLevel=20;
Diameter=50;
Initlevel=10;
newElevation=500;
initqual=0;
MinLevel=0;
MinVol=0;
newValveID='V3pbv'; 
newValveDiameter=100; 
Code='PBV'; 
newValveSetting=15; 
errcode=d.addBinTank(newID,x,y,MaxLevel,Diameter,Initlevel,newElevation,initqual,MinLevel,MinVol,newValveID,...
ToNodeID,newValveDiameter,newValveSetting,Code);
d.plot('nodes','yes','links','yes');
disp('Press any key to continue...')
pause


%% Remove NODE, LINK - functions
% Remove node
[errcode]=d.removeBinNodeID('J3');
[errcode]=d.removeBinNodeID('T3');
d.plot('nodes','yes','links','yes');

% Remove Link
d.plot('highlightlink',{'121'})
errcode=d.removeBinLinkID('121')
d.plot('nodes','yes','highlightlink',{'PUMP1','PU1'},'fontsize',10);
errcode=d.removeBinLinkID('PUMP1')
errcode=d.removeBinLinkID('PU1');d.Binplot('nodes','yes','links','yes');
%Warning: Node J2 disconnected. 

%VALVES
%PRV 
newValveID='V1'; 
toNode='J2'; 
fromNode='22';
newValveDiameter=110;
newValveSetting=10; 
errcode=d.addBinValvePRV(newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
d.plot('nodes','yes','links','yes','highlightlink',{'V1'},'highlightnode',{'J1'})
errcode=d.removeBinNodeID('J1');
d.plot('nodes','yes','links','yes');
%PSV 
errcode=d.addBinValvePSV('V2','32','10',newValveDiameter,newValveSetting);
d.plot('nodes','yes','highlightlink',{'V2','V1'},'fontsize',8);
%PBV
errcode=d.addBinValvePBV('V3','12','23',newValveDiameter,newValveSetting);
d.plot('nodes','yes','highlightlink',{'V3'},'fontsize',8);
% FCV
errcode=d.addBinValveFCV('V4','21','10',newValveDiameter,newValveSetting);
d.plot('nodes','yes','highlightlink',{'V4'},'fontsize',8);
%TCV
errcode=d.addBinValveTCV('V5','32','13',newValveDiameter,newValveSetting);
d.plot('nodes','yes','highlightlink',{'V5'},'fontsize',8);
d.getBinComputedAllParameters  %"Run was unsuccessful."
d.getComputedHydraulicTimeSeries
% %GPV
errcode=d.addBinValveGPV('V6','11','22',newValveDiameter,newValveSetting); %%%%%%%%%%%%ERROR
if errcode
    d.Binplot('nodes','yes','highlightlink',{'V6'},'fontsize',8);
    errcode=d.removeBinLinkID('V6');
end
d.plot;
disp('Press any key to continue...')
pause
close all;
end



%% GET, REMOVE RULES CONTROLS
% section [RULES]
if strcmp(inpname,'BWSN_Network_1.inp')
    p=d.getBinRulesControlsInfo
    errcode=d.removeBinRulesControlLinkID(p.BinRulesControlLinksID{4}{3}); % PUMP-170
    errcode=d.removeBinRulesControlNodeID(p.BinRulesControlNodesID{1}{2}); % TANK-130
    v=d.getBinRulesControlsInfo
    disp('Press any key to continue...')
    pause
end



%% GET, SET LINK STATUS 
% SECTION [STATUS]
LinksInfo=d.getBinLinksInfo;
if d.BinLinkPumpCount
    ps=LinksInfo.BinLinkPumpStatus;
    idpumps=LinksInfo.BinLinkPumpStatusNameID;
    ps{1}='Closed';
    errcode=d.setBinLinkPumpStatus(ps) 
    pppnew=d.getBinLinksInfo
    pppnew.BinLinkInitialStatus
    d.getLinkInitialStatus
    disp('Press any key to continue...')
    pause
end




%% GET, SET REACTIONS
% SECTION [REACTIONS]
pp=LinksInfo.BinLinkWallReactionCoeff;
pp2=LinksInfo.BinLinkBulkReactionCoeff;
errcode=d.setBinLinkReactionCoeff('wall',pp,'bulk',pp2);
d.getLinkBulkReactionCoeff
d.getLinkWallReactionCoeff

errcode=d.setBinLinkGlobalWallReactionCoeff(1.99)
errcode=d.setBinLinkGlobalBulkReactionCoeff(2)
disp('Press any key to continue...')
pause



%% GET, SET Initial Quality 
% SECTION [QUALITY]
d.BinUpdateClass
p=d.BinNodeInitialQuality;
p(end)=1;
errcode=d.setBinNodeInitialQuality(p)
isequal(d.getBinNodesInfo.BinNodeInitialQuality,d.getNodeInitialQuality)
disp('Press any key to continue...')
pause



%% GET, SET OPTIONS
% SECTION [OPTIONS]
d.getBinOptionsInfo
d.BinQualityType
errcode=d.setBinQualityChem
d.getQualityType
errcode=d.setBinQualityAge
d.getQualityType
errcode=d.setBinQualityNone
d.getQualityType
traceID=d.BinQualityTraceNodeID;
n=d.BinNodeNameID;
errcode=d.setBinQualityTrace(n{2})
d.getQualityType
d.getQualityTraceNodeIndex

species = 'Chlorine'; 
units = 'mg/L';
errcode = d.setBinQualType(species,units)
d.getQualityInfo.QualityChemName
d.getQualityInfo.QualityChemUnits

disp('Press any key to continue...')
pause



%% GET, SET TIMES
% SECTION [TIMES]
d.getBinTimesInfo
d.BinTimeSimulationDuration
d.BinTimeHydraulicStep
d.BinTimeQualityStep
d.BinTimePatternStep
d.BinTimePatternStart
d.BinTimeReportingStep
d.BinTimeReportingStart
d.BinTimeStatisticsIndex
d.BinTimeStatistics

errcode=d.setBinTimeSimulationDuration(48*3600);
errcode=d.setBinTimeHydraulicStep(1000);
errcode=d.setBinTimeQualityStep(1000);
errcode=d.setBinTimePatternStep(1000);
errcode=d.setBinTimePatternStart(1);
errcode=d.setBinTimeReportingStep(1000);
errcode=d.setBinTimeReportingStart(1);
errcode=d.setBinTimeStatisticsAverage;
errcode=d.setBinTimeStatisticsMaximum;
errcode=d.setBinTimeStatisticsMinimum;
errcode=d.setBinTimeStatisticsNone;
errcode=d.setBinTimeStatisticsRange;

P=d.getBinTimesInfo
isequal(d.getTimeSimulationDuration,P.BinTimeSimulationDuration)
isequal(d.getTimeHydraulicStep,P.BinTimeHydraulicStep)
isequal(d.getTimeQualityStep,P.BinTimeQualityStep)
isequal(d.getTimePatternStep,P.BinTimePatternStep)
isequal(d.getTimePatternStart,P.BinTimePatternStart)
isequal(d.getTimeReportingStep,P.BinTimeReportingStep)
isequal(d.getTimeReportingStart,P.BinTimeReportingStart)
isequal(d.getTimeStatisticsIndex,P.BinTimeStatisticsIndex)

disp('Press any key to continue...')
pause



%% ADD, GET, REMOVE CURVES
d.getBinCurvesInfo.BinCurveCount
d.getCurveCount
errcode=d.addBinCurvePump('C-1n3',[1500 1800 2000],[250 200 0])%Flow-Head
errcode=d.addBinCurveEfficiency('C-2n',[1500 1800 2000],[250 200 0])%Flow-Efficiency
errcode=d.addBinCurveVolume('C33',[1500 1800 2000],[250 200 0])%Heigh-Volume
errcode=d.addBinCurveHeadloss('C44',[1500 1800 2000],[250 200 0])%Flow-Headloss
d.getBinCurvesInfo
errcode=d.removeBinCurveID('C-1');
errcode=d.removeBinLinkID('PU3');
d.getBinCurvesInfo
% warning Flow & Heigh
errcode=d.addBinCurvePump('C-11',[2000 1500 1800],[250 200 0])
errcode=d.addBinCurveEfficiency('C-22',[1500 2000 1800],[250 200 0])%Flow-Efficiency
errcode=d.addBinCurveVolume('C333',[1500 2000 1800],[250 200 0])%Heigh-Volume
errcode=d.addBinCurveHeadloss('C244',[1500 2000 1800],[250 200 0])%Flow-Headloss
d.getBinCurvesInfo.BinCurveCount
d.getCurveCount
disp('Press any key to continue...')
pause



%% GET, SET Source Info
d.BinNodeSourcePatternIndex
d.BinNodeSourceQuality
d.BinNodeSourceTypeIndex
d.BinNodeSourceType
d.BinNodeSourcePatternNameID

v=d.getBinNodeSourceInfo
v.BinNodeSourceQuality
v.BinNodeSourceType
v.BinNodeSourcePatternNameID
v.BinNodeSourceQuality(2)=20;
b=d.getBinPatternsInfo
v.BinNodeSourcePatternNameID(2)=b.BinPatternNameID(1);
v.BinNodeSourceType{2}='SETPOINT'; % {'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}
errcode=d.setBinNodeSourceQuality(v)

v=d.getBinNodeSourceInfo
v.BinNodeSourceQuality
v.BinNodeSourceType
v.BinNodeSourcePatternNameID
v.BinNodeSourceQuality(3)=20;
b=d.getBinPatternsInfo
v.BinNodeSourcePatternNameID(3)=b.BinPatternNameID(1);
v.BinNodeSourceType{3}='CONCEN'; % {'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}
errcode=d.setBinNodeSourceQuality(v)
v=d.getBinNodeSourceInfo

d.getNodeSourceType
d.getNodeSourceQuality
d.getNodeSourcePatternIndex
disp('Press any key to continue...')
pause



%% GET, SET VALVES
if strcmp(inpname,'BWSN_Network_1.inp')
    if d.BinLinkValveCount
        index=1;
        m=d.BinLinkValveMinorLoss; 
        m(index)=0.2;
        dv=d.BinLinkValveDiameters;
        dv(index)=180;
        t=d.BinLinkValveType;
        t{index}='TCV';
        s=d.BinLinkValveSetting;
        s(index)=10;
        statuss=d.BinLinkValveStatus;
        statuss{index}='closed';
        errcode=d.setBinLinkValvesParameters('minorloss',m,'diameter',dv,'type',t,'setting',s,'status',statuss);
        LinkValveIndex=d.getLinkValveIndex;
        LinkType=d.getLinkType;
        LinkType(LinkValveIndex(index))
        d.getLinkMinorLossCoeff(LinkValveIndex)
        d.getLinkDiameter(LinkValveIndex)
        d.getLinkInitialSetting(LinkValveIndex)  %if write in [STATUS] section the values of setting in valves is diff! Check!
        d.getLinkInitialStatus(LinkValveIndex)
        disp('Press any key to continue...')
        pause
    end
end




%% ADD, SET, GET PATTERNS 
% add pattern
errcode=d.addBinPattern('new',1:0.1:2);
errcode=d.addBinPattern('new2',1:0.1:2);
d.getPatternNameID

% set pattern
bb=d.getBinPatternsInfo
idp=bb.BinPatternNameID;
values=bb.BinPatternValue{2};
values(1:end)=1;
errcode=d.setBinPattern(idp{1},values);
d.getBinPatternsInfo
errcode=d.setBinPattern(idp{2},values);
d.getPattern
d.getPatternNameID

disp('Press any key to continue...')
pause



%% GET Functions Info
d.getBinControlsInfo
d.getBinNodesInfo
d.getBinNodeSourceInfo
d.getBinCurvesInfo
d.getBinOptionsInfo
d.getBinTimesInfo
d.getBinPatternsInfo
d.getBinLinksInfo
d.BinUpdateClass

disp('Press any key to continue...')
pause

if strcmp(inpname,'Net1.inp')
%% GET, ADD, REMOVE CONTROLS
    v=d.getBinControlsInfo;
    errcode=d.removeBinControlLinkID(v.BinControlLinksID{1});
    % LINK x status AT TIME t
    v=d.getBinControlsInfo
    errcode=d.addBinControl('10','OPEN','10.00');
    % LINK x status AT CLOCKTIME c AM/PM
    errcode=d.addBinControl('10','OPEN','12.00','AM');
    v=d.getBinControlsInfo
    if ~errcode
        errcode=d.removeBinControlLinkID(v.BinControlLinksID{1});
    end
    v=d.getBinControlsInfo
    % LINK x status IF NODE y ABOVE/BELOW z
    errcode=d.addBinControl('12','OPEN','12','ABOVE',100)
    v=d.getBinControlsInfo
    if ~errcode
        errcode=d.removeBinControlNodeID(v.BinControlNodesID{1});
    end
    v=d.getBinControlsInfo
    disp('Press any key to continue...')
    pause
end


%% Other Functions
d = epanet(inpname);

d.getBinNodeNameID
d.getNodeNameID
d.getBinLinkNameID
d.getLinkNameID
d.getBinNodeCoordinates

a=d.getBinComputedAllParameters;
d.getBinElevationEachNode
d.getBinLengthEachLink
d.getBinDiameterEachLink
d.getBinComputedPumpIndexListLinks
d.getBinComputedPumpUtilization
d.getBinComputedAverageEfficiency
d.getBinComputedAverageKwattsOrMillionGallons
d.getBinComputedAverageKwatts
d.getBinComputedPeakKwatts
d.getBinComputedAverageCostPerDay
d.getBinComputedNodeDemand
d.getBinComputedNodeHead
d.getBinComputedNodePressure
d.getBinComputedNodeQuality
d.getBinComputedLinkFlow
d.getBinComputedLinkVelocity
d.getBinComputedLinkHeadloss
d.getBinComputedLinkQuality
d.getBinComputedLinkStatus
d.getBinComputedLinkSetting
d.getBinComputedLinkReactionRate
d.getBinComputedLinkFrictionFactor
d.getBinComputedAverageBulkReactionRate
d.getBinComputedAverageWallReactionRate
d.getBinComputedAverageTankReactionRate
d.getBinComputedAverageSourceInflow
          
disp('Press any key to continue...')
pause


%% Run epanet2d.exe 
% Computes quality and hydraulics
d.getBinComputedAllParameters

disp('Press any key to continue...')
pause


%% OTHER PROPERTIES
d.BinUpdateClass

d.InputFile
d.BinTempfile                

d.BinControlLinksID
d.BinControlNodesID
d.BinControlRulesCount
d.BinControlsInfo

d.BinCountInitialQualitylines
d.BinCountPatternlines
d.BinCountReactionlines
d.BinCountStatuslines

d.BinCurveAllLines
d.BinCurveCount
d.BinCurveNameID
d.BinCurveTypes
d.BinCurveXvalue
d.BinCurveYvalue


d.BinLinkCount
d.BinLinkFlowUnits
d.BinLinkNameID

d.BinLinkDiameters
d.BinLinkInitialStatus
d.BinLinkInitialStatusNameID
d.BinLinkLengths
d.BinLinkRoughnessCoeff       
d.BinLinkSettings             
d.BinLinkType  

d.BinLinkToNode  
d.BinLinkFromNode

d.BinLinkBulkReactionCoeff
d.BinLinkWallReactionCoeff   
d.BinLinkGlobalBulkReactionCoeff
d.BinLinkGlobalWallReactionCoeff

d.BinLinkPipeCount   
d.BinLinkPipeIndex            
d.BinLinkPipeNameID           

d.BinLinkPipeDiameters        
d.BinLinkPipeLengths          
d.BinLinkPipeMinorLoss        
d.BinLinkPipeRoughness        
d.BinLinkPipeStatus   

d.BinLinkPumpCount           
d.BinLinkPumpCurveNameID      
d.BinLinkPumpIndex          
d.BinLinkPumpNameID          
d.BinLinkPumpNameIDPower   

d.BinLinkPumpPatterns         
d.BinLinkPumpPower            
d.BinLinkPumpStatus           
d.BinLinkPumpStatusNameID  

d.BinLinkValveCount  
d.BinLinkValveIndex          
d.BinLinkValveNameID         

d.BinLinkValveDiameters     
d.BinLinkValveMinorLoss     
d.BinLinkValveSetting      
d.BinLinkValveStatus          
d.BinLinkValveStatusNameID    
d.BinLinkValveType  

d.BinNodeCount               
d.BinNodeNameID           
d.BinNodeBaseDemands         
d.BinNodeCoordinates          
d.BinNodeJunDemandPatternNameID  
d.BinNodeElevations           
d.BinNodeInitialQuality   
d.BinNodeType               

d.BinNodeJunctionCount 
d.BinNodeJunctionIndex        
d.BinNodeJunctionNameID 
d.BinNodeJunctionElevation    
d.BinNodeJunctionsBaseDemands
d.BinNodeJunctionsBaseDemandsID
d.BinNodePressureUnits    

d.BinNodeResDemandPatternNameID
d.BinNodeReservoirCount       
d.BinNodeReservoirElevation   
d.BinNodeReservoirIndex       
d.BinNodeReservoirNameID 

d.BinNodeSourcePatternIndex   
d.BinNodeSourcePatternNameID
d.BinNodeSourceQuality        
d.BinNodeSourceType        
d.BinNodeSourceTypeIndex   

d.BinNodeTankCount         
d.BinNodeTankDiameter         
d.BinNodeTankElevation        
d.BinNodeTankIndex           
d.BinNodeTankInitialLevel     
d.BinNodeTankMaximumWaterLevel
d.BinNodeTankMinimumFraction  
d.BinNodeTankMinimumWaterLevel
d.BinNodeTankMinimumWaterVolume
d.BinNodeTankMixID           
d.BinNodeTankMixModel        
d.BinNodeTankNameID           
d.BinNodeTankReservoirCount   

d.BinOptionsAccuracyValue     
d.BinOptionsDiffusivity       
d.BinOptionsEmitterExponent   
d.BinOptionsHeadloss         
d.BinOptionsMaxTrials         
d.BinOptionsPattern           
d.BinOptionsPatternDemandMultiplier
d.BinOptionsQualityTolerance  
d.BinOptionsSpecificGravity   
d.BinOptionsUnbalanced        
d.BinOptionsViscosity    

d.BinPatternCount            
d.BinPatternLengths           
d.BinPatternMatrix            
d.BinPatternNameID            
d.BinPatternValue     

d.BinQualityCode              
d.BinQualityTraceNodeID    
d.BinQualityTraceNodeIndex   
d.BinQualityType            
d.BinQualityUnits      

d.BinRulesControlLinksID      
d.BinRulesControlNodesID    
d.BinRulesControlsInfo      
d.BinRulesCount  

d.BinTimeHydraulicStep       
d.BinTimePatternStart         
d.BinTimePatternStep          
d.BinTimeQualityStep         
d.BinTimeReportingStart       
d.BinTimeReportingStep       
d.BinTimeSimulationDuration   
d.BinTimeStatistics          
d.BinTimeStatisticsIndex    

d.BinUnits_US_Customary     
d.BinUnits_SI_Metric   

d.BinUnits                    
d.BinUnits.BinLinkFlowUnits
d.BinUnits.BinQualityUnits
d.BinUnits.BinNodePressureUnits
d.BinUnits.BinPatternDemandsUnits
d.BinUnits.BinLinkPipeDiameterUnits
d.BinUnits.BinNodeTankDiameterUnits
d.BinUnits.BinEnergyEfficiencyUnits
d.BinUnits.BinNodeElevationUnits
d.BinUnits.BinNodeDemandUnits
d.BinUnits.BinNodeEmitterCoefficientUnits
d.BinUnits.BinEnergyUnits
d.BinUnits.BinLinkFrictionFactorUnits
d.BinUnits.BinNodeHeadUnits
d.BinUnits.BinLinkLengthsUnits
d.BinUnits.BinLinkMinorLossCoeffUnits
d.BinUnits.BinLinkPumpPowerUnits
d.BinUnits.BinQualityReactionCoeffBulkUnits
d.BinUnits.BinQualityReactionCoeffWallUnits
d.BinUnits.BinLinkPipeRoughnessCoeffUnits
d.BinUnits.BinQualitySourceMassInjectionUnits
d.BinUnits.BinLinkVelocityUnits
d.BinUnits.BinNodeTankVolumeUnits
d.BinUnits.BinQualityWaterAgeUnits
        
d.BinClose
d.unload
fprintf('Test finished.\n')