%% Example: Calculate the quality at node 2

%% with EPANET libraries epanet20012
clc; clear all; close all;
network='Net2_Rossman2000.inp';
if strcmp(computer('arch'),'win32')           
    d=epanet(network,'epanet20012x86');
elseif strcmp(computer('arch'),'win64')
    d=epanet(network,'epanet20012x64');
end
nodeindex=2;step=300;
h=ones(1,d.NodeCount)*0.5;
h(d.getNodeIndex('10'))=0.1;
d.setNodeInitialQuality(h);
m=zeros(1,d.NodeCount);
d.setNodeSourcePatternIndex(m);
m(1)=0.8;
d.setNodeSourceQuality(m);
l=d.getLinkBulkReactionCoeff;
l=-ones(1,length(l))*0.3;
d.setLinkBulkReactionCoeff(l)
d.setLinkWallReactionCoeff(-ones(1,length(l)));
d.setTimeQualityStep(step);
d.saveInputFile(d.Bintempfile);
d.solveCompleteHydraulics
res=d.getComputedQualityTimeSeries;
%d.plot('nodes','yes','highlightnode',{'2'})
figure;
plot(res.Time,res.Quality(:,nodeindex),'b','LineWidth',2);    hold on; 
xlabel('Time(hours)');
ylabel('Fluride(mg/L');
str=sprintf('Fluride for node ID "%s"',char(d.NodeNameID(nodeindex)));
title(str);grid on;

% epanet20013
if strcmp(computer('arch'),'win32')           
    d=epanet(network,'epanet20013patchx86');
elseif strcmp(computer('arch'),'win64')
    d=epanet(network,'epanet20013patchx64');
end
nodeindex=2;step=300;
h=ones(1,d.NodeCount)*0.5;
h(d.getNodeIndex('10'))=0.1;
d.setNodeInitialQuality(h);
m=zeros(1,d.NodeCount);
d.setNodeSourcePatternIndex(m);
m(1)=0.8;
d.setNodeSourceQuality(m);
l=d.getLinkBulkReactionCoeff;
l=-ones(1,length(l))*0.3;
d.setLinkBulkReactionCoeff(l)
d.setLinkWallReactionCoeff(-ones(1,length(l)));
d.setTimeQualityStep(step);
d.saveInputFile(d.Bintempfile);
d.solveCompleteHydraulics
res=d.getComputedQualityTimeSeries;
%d.plot('nodes','yes','highlightnode',{'2'})
plot(res.Time,res.Quality(:,nodeindex),'g','LineWidth',2);    hold on; 
xlabel('Time(hours)');
ylabel('Fluride(mg/L');
str=sprintf('Fluride for node ID "%s"',char(d.NodeNameID(nodeindex)));
title(str);grid on;


%% without EPANET libraries
d=epanet(network,'bin');
d.setBinNodeInitialQuality(h);
d.setBinTimeQualityStep(step);
d.setBinLinkGlobalBulkReactionCoeff(-0.3);
d.setBinLinkGlobalWallReactionCoeff(-1);
v=d.getBinNodeSourceInfo;
v.BinNodeSourceQuality(1)=0.8;
v.BinNodeSourcePatternNameID(1)='';
v.BinNodeSourceType{1}='CONCEN'; % {'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}
errcode=d.setBinNodeSourceQuality(v);
mn=d.getBinTimesInfo;
times=0:mn.BinTimeQualityStep:mn.BinTimeSimulationDuration;
resbin=d.getBinComputedNodeQuality';
plot(times,resbin(:,nodeindex),'r');

%% with MSX libraries
if strcmp(computer('arch'),'win32')           
    d=epanet(network,'epanet20012x86');
elseif strcmp(computer('arch'),'win64')
    d=epanet(network,'epanet20012x64');
end
d.msx([network(1:end-4),'.msx']);
nnn=d.getMsxComputedQualityNode;
% cmap=hsv(5);
% for i=1:d.getMsxSpeciesCount
plot(nnn.Time,nnn.Quality{nodeindex}{1},'k','LineWidth',2);
%     hold on; 
% end
legend('epanet20012','epanet20013','WithOutLib','WithMSX')
delete('*_temp*',[d.BinInputFile(1:end-9),'.txt'])
