%% Add CV pipe with bin Functions addBinJunction
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

disp('Add CV PIPE with Bin Functions addBinJunction')

% Plot newtork with bin functions
d.Binplot('nodes','yes');
% addBinJunction
% arguments: newNodeID,X,Y,ToNodeID,newElevation,newBaseDemand,newDemandPattern,newPipeID,
% newLength,newDiameter,newRoughness,Code
% Add Junction + CV pipe
newID='J1';
[x,y]=ginput(1);
ToNodeID='10'; 
newElevation=500; %ft
newBaseDemand=0;
newDemandPattern='1';
newPipeID='CV-P2';
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='CVPIPE';%'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'
errcode=d.addBinJunction(newID,x,y,newElevation,newBaseDemand,newDemandPattern,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);

% Plot normal with new pipe
d.plot('links','yes');

% Unload library
d.unload

