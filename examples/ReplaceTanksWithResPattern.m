function ReplaceTanksWithResPattern()
%% Choose Network
fclose all;clear class;
close all;clear all;clc;
addpath(genpath(pwd));
[inpname] = enterNetwork();

%% Load network and simulate
d=epanet(inpname);
allParameters=d.getBinComputedAllParameters;
H = allParameters.BinnodeHead;

%% Find tank patterns:
tankInd = d.NodeTankIndex;
if isempty(tankInd)
    disp('No tanks in network. Please Choose other network.');
    d.unload; return;
end
tankPat = H(:,tankInd);

%% Create patterns for new reservoirs from tanks
tankPatInd=0;
for i=tankInd
    tankPatInd = tankPatInd +1;
    patInd(tankPatInd)=d.addPattern(['P-',d.NodeNameID{i}],tankPat(tankPatInd));
end
d.BinUpdateClass;

%% Replace tanks with reservoirs
tankInd = d.NodeTankIndex;
for i = 1:d.NodeTankCount
    newID{i}=['Res-',d.getNodeNameID{tankInd(1)}];
    d = replaceTankWithReservoirBin(d,tankInd(1),newID{i});
    tankInd = d.getNodeTankIndex;
end
j=1;
indicesRes =d.getNodeIndex(newID); 
for i=indicesRes
    d.setNodeDemandPatternIndex(i,patInd(j));j=j+1;
end

%% Save new input file
filename=['New_',inpname];
d.saveInputFile(filename);
d.unload
end


%%%%%%%%%%%%%%%%% Other functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inpname] = enterNetwork()
clc
dirName = [pwd,'\networks\*.inp'];
Allinpnames = dir(dirName);
disp('Choose Benchmark Water Network:')
for i=1:length(Allinpnames)
    disp([num2str(i),'. ', Allinpnames(i).name])
end
x = input(sprintf('\nEnter Network Number: '));
inpname=Allinpnames(x).name;    
end

function d = replaceTankWithReservoirBin(d,indexTank,newID)

    tankID = d.getNodeNameID(indexTank);
    Coordinates = d.getNodeCoordinates;
    X = Coordinates{1}(indexTank);
    Y = Coordinates{2}(indexTank);

    % Vriskoume ta pipes p en enwmena panw tou gia na ta 3ana valoume
    LinkNodesIndex=d.getLinkNodesIndex;
    pipesIndexEnwmenaStoTank=find(sum(d.getLinkNodesIndex==indexTank,2))';

    head = d.getNodeElevations(indexTank);
    status=d.getBinLinksInfo.BinLinkPipeStatus; %connected pipes status
    for i=1:length(pipesIndexEnwmenaStoTank)
        ind=find(LinkNodesIndex(pipesIndexEnwmenaStoTank(i),:)~=indexTank,1);
        nodeIndex(i)=LinkNodesIndex(pipesIndexEnwmenaStoTank(i),ind);
        toNodeID(i) = d.getNodeNameID(nodeIndex(i));
        linklengths(i)=d.getLinkLength(pipesIndexEnwmenaStoTank(i));
        linkDiameters(i)=d.getLinkDiameter(pipesIndexEnwmenaStoTank(i));
        linkRoughness(i)=d.getLinkRoughnessCoeff(pipesIndexEnwmenaStoTank(i));
        linkType(i)=d.getLinkType(pipesIndexEnwmenaStoTank(i));
        linkID(i)=d.getLinkNameID(pipesIndexEnwmenaStoTank(i));
    end
    
    errcode=d.removeBinNodeID(tankID);

    newElevation=1;
    ToNodeID=toNodeID{1};
    newPipeID=linkID{1};
    newLength=linklengths(1);
    newDiameter=linkDiameters(1); 
    newRoughness=linkRoughness(1);
    Code=linkType{1};
    [errcode]=d.addBinReservoir(newID,X,Y,newElevation,newPipeID,...
    ToNodeID,newLength,newDiameter,newRoughness,Code);

    for i=2:length(pipesIndexEnwmenaStoTank)
        [errcode]=d.addBinPipe(linkID{i},newID,toNodeID{i},linklengths(i),linkDiameters(i),linkRoughness(i));
    end
    errcode=d.setBinLinkPipeStatus(status);

end
