%% Assing a new curves to a specific pump
%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

fprintf('\n');

% Add new curve in the network
indexCurve=d.addCurve('NewCurve', [1800 200; 1500 400]);

% Get pump index
pumpIndex = d.getLinkPumpIndex(1);

% Get head curve index
[HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;

disp(['Head Curve Index: ' num2str(HeadCurveIndex)] );
disp(['On pump index: ' num2str(PumpIndex)] );

% Assing new curve index on the specific pump 
d.setLinkPumpHeadCurveIndex(pumpIndex,indexCurve);
fprintf(['\nAssign new curve to pump: ' num2str(PumpIndex),'\n\n'] );
[HeadCurveIndex,PumpIndex] = d.getLinkPumpHeadCurveIndex;
disp(['New Head Curve Index: ' num2str(HeadCurveIndex)] );
disp(['On pump index: ' num2str(PumpIndex)] );

fprintf('\n');

% Unload library
d.unload
