start_toolkit
d = epanet('ky10.inp');

d.runsCompleteSimulation;

%% all
disp('getNodeTankInitialLevel');            d.getNodeTankInitialLevel
disp('getNodeTankInitialWaterVolume');      d.getNodeTankInitialWaterVolume
disp('getNodeTankMinimumWaterVolume');      d.getNodeTankMinimumWaterVolume
disp('getNodeTankDiameter');                d.getNodeTankDiameter
disp('getNodeTankMixZoneVolume');           d.getNodeTankMixZoneVolume
disp('getNodeTankMaximumWaterVolume');      d.getNodeTankMaximumWaterVolume
disp('getNodeTankVolumeCurveIndex');        d.getNodeTankVolumeCurveIndex
disp('getNodeTankMinimumWaterLevel');       d.getNodeTankMinimumWaterLevel
disp('getNodeTankMaximumWaterLevel');       d.getNodeTankMaximumWaterLevel
disp('getNodeTankMixingFraction');          d.getNodeTankMixingFraction
disp('getNodeTankBulkReactionCoeff');       d.getNodeTankBulkReactionCoeff
disp('getNodeTankVolume');                  d.getNodeTankVolume
disp('getNodeTankCanOverFlow');             d.getNodeTankCanOverFlow
disp('getNodeTankMixingModelCode');         d.getNodeTankMixingModelCode
disp('getNodeTankMixingModelType');         d.getNodeTankMixingModelType
disp('getNodeTankMixingModelCode');         d.getNodeTankMixingModelCode

%% single
disp('getNodeTankInitialLevel(1)');            d.getNodeTankInitialLevel(1)
disp('getNodeTankInitialWaterVolume(1)');      d.getNodeTankInitialWaterVolume(1)
disp('getNodeTankMinimumWaterVolume(1)');      d.getNodeTankMinimumWaterVolume(1)
disp('getNodeTankDiameter(1)');                d.getNodeTankDiameter(1)
disp('getNodeTankMixZoneVolume(1)');           d.getNodeTankMixZoneVolume(1)
disp('getNodeTankMaximumWaterVolume(1)');      d.getNodeTankMaximumWaterVolume(1)
disp('getNodeTankVolumeCurveIndex(1)');        d.getNodeTankVolumeCurveIndex(1)
disp('getNodeTankMinimumWaterLevel(1)');       d.getNodeTankMinimumWaterLevel(1)
disp('getNodeTankMaximumWaterLevel(1)');       d.getNodeTankMaximumWaterLevel(1)
disp('getNodeTankMixingFraction(1)');          d.getNodeTankMixingFraction(1)
disp('getNodeTankBulkReactionCoeff(1)');       d.getNodeTankBulkReactionCoeff(1)
disp('getNodeTankVolume(1)');                  d.getNodeTankVolume(1)
disp('getNodeTankCanOverFlow(1)');             d.getNodeTankCanOverFlow(1)
disp('getNodeTankMixingModelCode(1)');         d.getNodeTankMixingModelCode(1)
disp('getNodeTankMixingModelType(1)');         d.getNodeTankMixingModelType(1)
disp('getNodeTankMixingModelCode(1)');         d.getNodeTankMixingModelCode(1)
%% range
disp('getNodeTankInitialLevel(:)');            d.getNodeTankInitialLevel(1:4)
disp('getNodeTankInitialWaterVolume(:)');      d.getNodeTankInitialWaterVolume(1:4)
disp('getNodeTankMinimumWaterVolume(:)');      d.getNodeTankMinimumWaterVolume(1:4)
disp('getNodeTankDiameter(:)');                d.getNodeTankDiameter(1:4)
disp('getNodeTankMixZoneVolume(:)');           d.getNodeTankMixZoneVolume(1:4)
disp('getNodeTankMaximumWaterVolume(:)');      d.getNodeTankMaximumWaterVolume(1:4)
disp('getNodeTankVolumeCurveIndex(:)');        d.getNodeTankVolumeCurveIndex(1:4)
disp('getNodeTankMinimumWaterLevel(:)');       d.getNodeTankMinimumWaterLevel(1:4)
disp('getNodeTankMaximumWaterLevel(:)');       d.getNodeTankMaximumWaterLevel(1:4)
disp('getNodeTankMixingFraction(:)');          d.getNodeTankMixingFraction(1:4)
disp('getNodeTankBulkReactionCoeff(:)');       d.getNodeTankBulkReactionCoeff(1:4)
disp('getNodeTankVolume(:)');                  d.getNodeTankVolume(1:4)
disp('getNodeTankCanOverFlow(:)');             d.getNodeTankCanOverFlow(1:4)
disp('getNodeTankMixingModelCode(:)');         d.getNodeTankMixingModelCode(1:4)
disp('getNodeTankMixingModelType(:)');         d.getNodeTankMixingModelType(1:4)
disp('getNodeTankMixingModelCode(:)');          d.getNodeTankMixingModelCode(1:4)
d.unload;

d = epanet('Net3.inp');

%% All pumps
disp('Running getLinkPumpPower');         d.getLinkPumpPower'
disp('Running getLinkPumpHCurve');        d.getLinkPumpHCurve'
disp('Running getLinkPumpECurve');        d.getLinkPumpECurve'
disp('Running getLinkPumpECost');         d.getLinkPumpECost'
disp('Running getLinkPumpEPat');          d.getLinkPumpEPat'
disp('Running getLinkPumpPatternIndex');  d.getLinkPumpPatternIndex'

%% Single pump (1)
disp('Running getLinkPumpPower(1)');         d.getLinkPumpPower(1)
disp('Running getLinkPumpHCurve(1)');        d.getLinkPumpHCurve(1)
disp('Running getLinkPumpECurve(1)');        d.getLinkPumpECurve(1)
disp('Running getLinkPumpECost(1)');         d.getLinkPumpECost(1)
disp('Running getLinkPumpEPat(1)');          d.getLinkPumpEPat(1)
disp('Running getLinkPumpPatternIndex(1)');  d.getLinkPumpPatternIndex(1)

%% Range of pumps
disp('Running getLinkPumpPower(:)');         d.getLinkPumpPower(1:2)
disp('Running getLinkPumpHCurve(:)');        d.getLinkPumpHCurve(1:2)
disp('Running getLinkPumpECurve(:)');        d.getLinkPumpECurve(1:2)
disp('Running getLinkPumpECost(:)');         d.getLinkPumpECost(1:2)
disp('Running getLinkPumpEPat(:)');          d.getLinkPumpEPat(1:2)
disp('Running getLinkPumpPatternIndex(:)');  d.getLinkPumpPatternIndex(1:2)


d.unload;