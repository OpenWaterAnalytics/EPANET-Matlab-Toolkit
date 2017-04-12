clear all; close all; clc;
addpath(genpath(pwd));
d=epanet('Net1.inp');tic;
Results = d.getComputedTimeSeries %Using ENepanet, create and read binary file
tocResults = toc;tic;
Hydraulics = d.getComputedHydraulicTimeSeries % Using the functions(ENopenH, ENinit, ENrunH, ENgetnodevalue/&ENgetlinkvalue, ENnextH, ENcloseH)
tocHydraulics = toc;tic;
Quality = d.getComputedQualityTimeSeries % ENopenQ, ENinitQ, ENrunQ, ENgetnodevalue/&ENgetlinkvalue, ENstepQ, ENcloseQ
tocQuality = toc;tic;

pipeindex=4;
nodeindex=6;

%For step by step Hydraulic Analysis
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1; P=[];T=[]; D=[]; H=[];F=[];
while (tstep>0)
    t=d.runHydraulicAnalysis;
    P=[P; d.getNodePressure];
    D=[D; d.getNodeActualDemand];
    H=[H; d.getNodeHydaulicHead];
    F=[F; d.getLinkFlows];
    T=[T; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis

%For step by step Quality Analysis
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T=[];Q=[];  
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    P=[P; d.getNodePressure];
    Q=[Q; d.getNodeActualQuality];
    T=[T; t];
    tleft = d.stepQualityAnalysisTimeLeft;
end
d.closeQualityAnalysis;
d.unload

%Run time d.getComputedTimeSeries
disp(['Run Time of function d.getComputedTimeSeries: ', num2str(tocResults), '(sec)'])
disp(['Run Time of function d.getComputedHydraulicTimeSeries: ', num2str(tocHydraulics), '(sec)'])
disp(['Run Time of function d.getComputedQualityTimeSeries: ', num2str(tocQuality), '(sec)'])

figure;
plot(Results.Flow(:,pipeindex));
title('d.getComputedTimeSeries (Ignore events)');
xlabel('Time(hours)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])

figure;
plot(Hydraulics.Flow(:,pipeindex));
title('d.getComputedHydraulicTimeSeries');
xlabel('Time(hours)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])

h=figure;
plot(0);axis off
whitebg('w');

figure;
plot(F(:,pipeindex));
title('step by step Hydraulic Analysis');
xlabel('Time(hours)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])

figure;
plot(Results.NodeQuality(:,nodeindex));
title('d.getComputedTimeSeries (Ignore events)');
xlabel('Time(hours)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Node ID "',d.NodeNameID{nodeindex},'"'])

figure;
plot(Quality.NodeQuality(:,nodeindex));
title('d.getComputedQualityTimeSeries');
xlabel('Time(hours)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Link ID "',d.NodeNameID{nodeindex},'"'])

figure;
plot(Q(:,nodeindex));
title('step by step Quality Analysis');
xlabel('Time(hours)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Link ID "',d.NodeNameID{nodeindex},'"'])

