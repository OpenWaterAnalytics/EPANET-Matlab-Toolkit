%% Example: Replace pump with pipe
clc;clear all;
inpname='Net1_Rossman2000.inp';
tic;d=epanet(inpname);toc
d.BinUpdateClass; 

d.plot('nodes','yes','links','yes');


errcode=d.removeBinNodeID(d.LinkPumpNameID{1});

newID='9';  % add reservoir because is removed automatically.
% [x,y]=ginput(1);
x=9.9786;
y=70.2286;
newElevation=800; %ft
ToNodeID='10'; 
newPipeID='newPipe';
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='PIPE';
[errcode]=d.addBinReservoir(newID,x,y,newElevation,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);
if ~errcode
    warning('Replaced pump with pipe is successful.');
end
d.plot('nodes','yes','links','yes');

d.unload
delete('*_temp*',[d.inputfile(1:end-4),'.txt'])