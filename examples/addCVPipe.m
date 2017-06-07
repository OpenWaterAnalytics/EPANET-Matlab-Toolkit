%% Example 1: Add CV PIPE 
d=epanet('Net1.inp');
disp('Add a CV pipe')
d.plot;
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
index=d.addLinkPipeCV('CVPipe',fromNode,toNode);
d.plot('links','yes');
d.unload;

%% Example 2: Add CV PIPE with Bin Functions
d=epanet('Net1.inp');
d.Binplot();
[errcode]=d.addBinCVPipe('CV-P1',fromNode,toNode,1000,10,100);
d.plot('links','yes');
d.unload;

%% Example 3: Add CV PIPE with Bin Functions addBinJunction
d=epanet('Net1.inp');
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
d.plot('links','yes');
d.unload;