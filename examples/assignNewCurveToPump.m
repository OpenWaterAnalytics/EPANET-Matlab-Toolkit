%% Assing a new curves to a specific pump

d=epanet('Net1.inp');
fprintf('\n');

indexCurve=d.addCurve('NewCurve', [1800 200; 1500 400]);

pumpIndex = d.getLinkPumpIndex(1);

[HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;

disp(['Head Curve Index: ' num2str(HeadCurveIndex)] );
disp(['On pump index: ' num2str(PumpIndex)] );

d.setLinkPumpHeadCurveIndex(pumpIndex,indexCurve);
fprintf(['\nAssign new curve to pump: ' num2str(PumpIndex),'\n\n'] );

[HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;
disp(['New Head Curve Index: ' num2str(HeadCurveIndex)] );
disp(['On pump index: ' num2str(PumpIndex)] );

fprintf('\n');
d.unload
