%% EPANET-Matlab Toolkit Test Part 7
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc; clear; close all;
addpath(genpath(pwd));
  
% Create EPANET object using the INP file
inpname='Net1.inp';  
% Net1 Net2 Net3 BWSN_Network_1

d=epanet(inpname,'loadfile'); %Using loadfile: properties not calculated
% d=epanet(inpname);

if d.LinkPumpCount
    [HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;
    indexCurve=d.addCurve('NewCur2', [1800 200; 1500 400]);
    fromNode = d.getNodeNameID{2};
    toNode = d.getNodeNameID{6};
    index=d.addLinkPump('Pump',fromNode,toNode);
    d.setLinkPumpHeadCurveIndex(index,indexCurve);
    [HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;
    d.deleteLink(index);
end

%% Add/delete Node
d.plot('nodes','yes');
disp('Add a junction')
d.getNodeNameID
d.getNodeJunctionCount
[x,y]=ginput(1);
index = d.addNodeJunction('Junction');     
d.setNodeCoordinates(index,[x,y]);
d.getNodeNameID
d.getNodeJunctionCount
d.plot('nodes','yes');

disp('Add a reservoir')
d.getNodeReservoirCount
[x,y]=ginput(1);
index = d.addNodeReservoir('Reservoir');
d.setNodeCoordinates(index,[x,y]);
d.getNodeNameID
d.getNodeReservoirCount
d.plot('nodes','yes');

disp('Add a tank')
d.getNodeTankCount
[x,y]=ginput(1);
index = d.addNodeTank('Tank');
d.setNodeCoordinates(index,[x,y]);
d.getNodeNameID
d.getNodeTankCount
d.plot('nodes','yes');

disp('Delete a node')
d.getNodeCount
nodeIndex = d.getNodeIndex('Junction');
d.deleteNode(nodeIndex);
d.getNodeNameID
d.getNodeCount
d.plot('nodes','yes');


          
%% Add/delete Link
disp('Add a cv pipe')
d.getLinkNameID
d.getLinkPipeCount
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
index=d.addLinkPipeCV('CVPipe',fromNode,toNode);
d.getLinkNameID
d.getLinkPipeCount
d.getLinkType
d.plot('links','yes');
d.deleteLink(index);
d.plot('links','yes');

disp('Add a pipe')
d.getLinkNameID
d.getLinkPipeCount
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
index=d.addLinkPipe('Pipe',fromNode,toNode);
d.getLinkNameID
d.getLinkPipeCount
d.deleteLink(index);

disp('Add a pump')
d.getLinkPumpCount
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
index=d.addLinkPump('Pump',fromNode,toNode);
d.getLinkNameID
d.getLinkPumpCount
d.getLinkType
d.deleteLink(index);

%similar
disp('Add a valve')
index=d.addLinkValvePRV('PRV-V',fromNode,toNode); 
index1=d.addLinkValveFCV('FCV',fromNode,toNode); 
index2=d.addLinkValveGPV('GPV',fromNode,toNode); 
index3=d.addLinkValvePBV('PBV',fromNode,toNode); 
index4=d.addLinkValvePSV('PSV',fromNode,toNode); 
index5=d.addLinkValveTCV('TCV',fromNode,toNode); 
d.getLinkType

d.unload
fprintf('Test finished.\n')