d=epanet('Net1.inp');tic;

hours = 100;

d.setTimeSimulationDuration(hours*3600);
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
tstep=1;P=[];T_H=[];D=[];H=[];F=[];
while (tstep>0)
    t=d.runHydraulicAnalysis;
    P=[P; d.getNodePressure];
    D=[D; d.getNodeActualDemand];
    H=[H; d.getNodeHydaulicHead];
    F=[F; d.getLinkFlows];
    T_H=[T_H; t];
    tstep=d.nextHydraulicAnalysisStep;
end
d.closeHydraulicAnalysis

%For step by step Quality Analysis
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; P=[];T_Q=[];Q=[];  
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    Q=[Q; d.getNodeActualQuality];
    T_Q=[T_Q; t];
    tleft = d.stepQualityAnalysisTimeLeft;
end
d.closeQualityAnalysis;
d.unload

%Run time d.getComputedTimeSeries
fprintf('\nSimulation duration: %d hours\n', hours);
disp(['Run Time of function d.getComputedTimeSeries: ', num2str(tocResults), '(sec)'])
disp(['Run Time of function d.getComputedHydraulicTimeSeries: ', num2str(tocHydraulics), '(sec)'])
disp(['Run Time of function d.getComputedQualityTimeSeries: ', num2str(tocQuality), '(sec)'])

figure;
plot(Results.Time, Results.Flow(:,pipeindex));
title('d.getComputedTimeSeries (Ignore events)');
xlabel('Time (sec)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

figure;
plot(Hydraulics.Time, Hydraulics.Flow(:,pipeindex));
title('d.getComputedHydraulicTimeSeries');
xlabel('Time (sec)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

h=figure;
plot(0);axis off
whitebg('w');

figure;
plot(T_H, F(:,pipeindex));
title('step by step Hydraulic Analysis');
xlabel('Time (sec)'); ylabel(['Flow (',d.LinkFlowUnits{1},') - Link ID "',d.LinkNameID{pipeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

figure;
plot(Results.Time, Results.NodeQuality(:,nodeindex));
title('d.getComputedTimeSeries (Ignore events)');
xlabel('Time (sec)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Node ID "',d.NodeNameID{nodeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

figure;
plot(Quality.Time, Quality.NodeQuality(:,nodeindex));
title('d.getComputedQualityTimeSeries');
xlabel('Time (sec)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Link ID "',d.NodeNameID{nodeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

figure;
plot(T_Q, Q(:,nodeindex));
title('step by step Quality Analysis');
xlabel('Time (sec)'); ylabel(['Node Quality (',d.QualityChemUnits,') - Link ID "',d.NodeNameID{nodeindex},'"'])
set(gca,'XTickLabel',num2str(get(gca,'XTick').'));

