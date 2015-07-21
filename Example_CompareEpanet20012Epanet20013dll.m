%% Example: Compare epanet20012 - epanet20013 versions
clc; clear all; close all;
%tankdiameter

inpname='Net2_Rossman2000.inp';
if strcmp(computer('arch'),'win32')           
    epanetversion1='epanet20012x86'; %epanet20012x86 epanet20012x64
    epanetversion2='epanet20013patchx86'; %epanet20013patchx86 epanet20013patchx64
elseif strcmp(computer('arch'),'win64')
    epanetversion1='epanet20012x64'; %epanet20012x86 epanet20012x64
    epanetversion2='epanet20013patchx64'; %epanet20013patchx86 epanet20013patchx64
end

d=epanet(inpname,epanetversion1); 
a=d.getNodeTankDiameter;
d.setNodeTankDiameter(a);
res12=d.getComputedHydraulicTimeSeries;

d=epanet(inpname,epanetversion2); 
a=d.getNodeTankDiameter;
d.setNodeTankDiameter(a);
res13=d.getComputedHydraulicTimeSeries;

figure;
plot(res12.Time/3600,res12.TankVolume(:,d.NodeTankIndex));
hold on;
plot(res13.Time/3600,res13.TankVolume(:,d.NodeTankIndex),'r');
legend({'epanet20012','epanet20013'});
s=sprintf('Tank Volume (%s)',d.NodeTankVolumeUnits);
ylabel(s);
xlabel('Time (hrs)');
s=sprintf('Tank Volume at tank %s',char(d.NodeNameID(d.NodeTankIndex(1))));
title(s);

figure;
index=10;
plot(res12.Time/3600,res12.Pressure(:,index));
hold on;
plot(res13.Time/3600,res13.Pressure(:,index),'r');
legend({'epanet20012','epanet20013'});
s=sprintf('Pressure (%s)',d.NodePressureUnits);
ylabel(s);
xlabel('Time (hrs)');
s=sprintf('Pressure at node %s',char(d.NodeNameID(index)));
title(s);


figure;
index=10;
plot(res12.Time/3600,res12.Flow(:,index));
hold on;
plot(res13.Time/3600,res13.Flow(:,index),'r');
legend({'epanet20012','epanet20013'});
s=sprintf('Flow (%s)',char(d.LinkFlowUnits));
ylabel(s);  
xlabel('Time (hrs)');
s=sprintf('Flow at link %s',char(d.LinkNameID(index)));
title(s);


figure;
index=10;
plot(res12.Time/3600,res12.Head(:,index));
hold on;
plot(res13.Time/3600,res13.Head(:,index),'r');
legend({'epanet20012','epanet20013'});
s=sprintf('Head (%s)',d.NodeHeadUnits);
ylabel(s);
xlabel('Time (hrs)');
s=sprintf('Head at node %s',char(d.NodeNameID(index)));
title(s);


% d=epanet(inpname,epanetversion1); 
% a=d.getNodeTankDiameter;
% d.setNodeTankDiameter(a);
% d.getComputedHydraulicTimeSeries;
% res12=d.getComputedQualityTimeSeries
% 
% d=epanet(inpname,epanetversion2); 
% a=d.getNodeTankDiameter;
% d.setNodeTankDiameter(a);
% d.getComputedHydraulicTimeSeries;
% res13=d.getComputedQualityTimeSeries
% 
% figure;
% index=10;
% plot(res12.Time/3600,res12.Quality(:,index));
% hold on;
% plot(res13.Time/3600,res13.Quality(:,index),'r');
% legend({'epanet20012','epanet20013'});
% s=sprintf('Head (%s)',d.NodeHeadUnits);
% ylabel(s);
% xlabel('Time (hrs)');
% s=sprintf('Head at node %s',char(d.NodeNameID(index)));
% title(s);
