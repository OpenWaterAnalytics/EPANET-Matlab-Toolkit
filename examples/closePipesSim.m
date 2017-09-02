%% Closing pipes during simulation

d=epanet('Net1.inp');

d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
tstep=1; F=[];
index=2;

Status = [0 0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 1 1]';

i=1;
while (tstep>0)
    t=d.runHydraulicAnalysis;
    d.setLinkStatus(index ,Status(i)); i=i+1;
    F=[F; d.getLinkFlows];
    tstep=d.nextHydraulicAnalysisStep;
end

d.closeHydraulicAnalysis;
Flows = F(:,index);
T = table(Flows,Status);

fprintf('\nFlows and status for node index 2:\n')
T
d.unload;