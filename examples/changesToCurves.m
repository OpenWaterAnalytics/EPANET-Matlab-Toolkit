% Changes to curves..
% |getComputedTimeSeries|Computed Hydraulic & Quality Time Series via ENepanet binary file|

% example 1
d = epanet('Net1.inp');
Results = d.getComputedTimeSeries;
figure;subplot(2,1,1);
indexNode = d.getNodeIndex('10');
plot(Results.Time, Results.Pressure(:,indexNode));
title('Before set curve (without events)');

headCurve = d.getLinkPumpHeadCurveIndex;
% d.getCurvesInfo
d.setCurve(headCurve,[2000 250]); 

Results = d.getComputedTimeSeries;  
subplot(2,1,2);
plot(Results.Time, Results.Pressure(:,indexNode));
title('After set curve (without events)');
suptitle('Example 1');
d.unload; 

% example 2
d = epanet('Net1.inp');

Hyd = d.getComputedHydraulicTimeSeries;  
figure;subplot(2,1,1);
indexNode = d.getNodeIndex('10');
plot(Hyd.Time, Hyd.Pressure(:,indexNode));
title('Before set curve        ');

headCurve = d.getLinkPumpHeadCurveIndex;
% d.getCurvesInfo
d.setCurve(headCurve,[2000 250]); 

% save and open "new" file
d.saveInputFile(d.BinTempfile);
d.loadEPANETFile(d.BinTempfile);

Hyd = d.getComputedHydraulicTimeSeries;  
subplot(2,1,2);
plot(Hyd.Time, Hyd.Pressure(:,indexNode));
title('After set curve        ');
suptitle('Example 2');
d.unload;
