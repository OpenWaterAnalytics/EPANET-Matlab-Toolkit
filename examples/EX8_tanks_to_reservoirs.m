%% Replace tanks with reservoirs using bin functions
% This example contains:
%   Choose Network.
%   Load network and simulate.
%   Find tank patterns.
%   Create patterns for new reservoirs from tanks.
%   Replace tanks with reservoirs.
%   Unload library.

%%
function EX8_tanks_to_reservoirs()
%Clear 
clear; close('all'); clc;
start_toolkit;

% Choose Network.
inpname = 'Net1.inp';

% Load network and simulate.
d = epanet('Net1.inp');
res = d.getComputedTimeSeries;
H = res.Head;

% Find tank patterns.
tank_index = d.getNodeTankIndex;
if isempty(tank_index)
    disp('No tanks in network. Please Choose other network.');
    d.unload; return;
end
tankPat = H(:, tank_index);

% Create patterns for new reservoirs from tanks.
tankPatInd=0;
for i = tank_index
    tankPatInd = tankPatInd +1;
    patInd(tankPatInd) = d.addPattern(['P-',d.getNodeNameID{i}],tankPat(tankPatInd));
end
d.BinUpdateClass;

% Replace tanks with reservoirs.
i = 1;
new_id = {};
for tank_i = tank_index
    new_id{i} = ['Res-', d.getNodeNameID{tank_i}];
    replace_tanks_with_reservoirs(d, tank_i, new_id{i});
    tank_i = d.getNodeTankIndex;
    i = i+1;
end
j=1;
indices_res = d.getNodeIndex(new_id); 
for i=indices_res
    d.setNodeDemandPatternIndex(i, patInd(j));j=j+1;
end

d.plot;
% Save new input file
%filename=['new_',inpname];
%d.saveInputFile(filename);

% Unload library.
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

function d = replace_tanks_with_reservoirs(d, tank_index, new_id)

    tank_id = d.getNodeNameID(tank_index);
    Coordinates = d.getNodeCoordinates;
    X = Coordinates{1}(tank_index);
    Y = Coordinates{2}(tank_index);

    % Find pipes who connected with specific tank
    LinkNodesIndex = d.getLinkNodesIndex;
    pipes_conn_on_tanks = find(sum(d.getLinkNodesIndex == tank_index, 2))';

    head = d.getNodeElevations(tank_index);
    status = d.getLinkStatus;
    for i=1:length(pipes_conn_on_tanks)
        ind = find(LinkNodesIndex(pipes_conn_on_tanks(i),:)~=tank_index,1);
        nodeIndex(i) = LinkNodesIndex(pipes_conn_on_tanks(i),ind);
        toNodeID(i) = d.getNodeNameID(nodeIndex(i));
        linklengths(i) = d.getLinkLength(pipes_conn_on_tanks(i));
        linkDiameters(i) = d.getLinkDiameter(pipes_conn_on_tanks(i));
        linkRoughness(i) = d.getLinkRoughnessCoeff(pipes_conn_on_tanks(i));
        linkType(i) = d.getLinkType(pipes_conn_on_tanks(i));
        linkID(i) = d.getLinkNameID(pipes_conn_on_tanks(i));
    end
    
    d.removeBinNodeID(tank_id);

    newElevation = head;
    ToNodeID=toNodeID{1};
    newPipeID=linkID{1};
    newLength=linklengths(1);
    newDiameter=linkDiameters(1); 
    newRoughness=linkRoughness(1);
    Code=linkType{1};
    
    node_index = d.addBinNodeReservoir(new_id, [X,Y], newElevation, '',0, ... 
    {Code, {newPipeID}, new_id, ToNodeID, newLength, newDiameter, newRoughness});

    for i=2:length(pipes_conn_on_tanks)
        d.addBinLinkPipe(linkID{i},new_id,toNodeID{i},linklengths(i),linkDiameters(i),linkRoughness(i));
    end
    
    d.setLinkStatus(status);
    
end
