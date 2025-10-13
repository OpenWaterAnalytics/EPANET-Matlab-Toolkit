start_toolkit
d = epanet('ky10.inp');
d.getNodeTankCount
firstTank = d.getNodeTankIndex(1)

%% all
disp('getNodeTankInitialLevel');             d.getNodeTankInitialLevel
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
disp('getNodeTankInitialLevel(923)');            d.getNodeTankInitialLevel(firstTank)
disp('getNodeTankInitialWaterVolume(923)');      d.getNodeTankInitialWaterVolume(firstTank)
disp('getNodeTankMinimumWaterVolume(923)');      d.getNodeTankMinimumWaterVolume(firstTank)
disp('getNodeTankDiameter(923)');                d.getNodeTankDiameter(firstTank)
disp('getNodeTankMixZoneVolume(923)');           d.getNodeTankMixZoneVolume(firstTank)
disp('getNodeTankMaximumWaterVolume(923)');      d.getNodeTankMaximumWaterVolume(firstTank)
disp('getNodeTankVolumeCurveIndex(923)');        d.getNodeTankVolumeCurveIndex(firstTank)
disp('getNodeTankMinimumWaterLevel(923)');       d.getNodeTankMinimumWaterLevel(firstTank)
disp('getNodeTankMaximumWaterLevel(923)');       d.getNodeTankMaximumWaterLevel(firstTank)
disp('getNodeTankMixingFraction(923)');          d.getNodeTankMixingFraction(firstTank)
disp('getNodeTankBulkReactionCoeff(923)');       d.getNodeTankBulkReactionCoeff(firstTank)
disp('getNodeTankVolume(923)');                  d.getNodeTankVolume(firstTank)
disp('getNodeTankCanOverFlow(923)');             d.getNodeTankCanOverFlow(firstTank)
disp('getNodeTankMixingModelCode(923)');         d.getNodeTankMixingModelCode(firstTank)
disp('getNodeTankMixingModelType(923)');         d.getNodeTankMixingModelType(firstTank)
disp('getNodeTankMixingModelCode(923)');          d.getNodeTankMixingModelCode(firstTank)
%% range
disp('getNodeTankInitialLevel(:)');            d.getNodeTankInitialLevel(firstTank:firstTank+3)
disp('getNodeTankInitialWaterVolume(:)');      d.getNodeTankInitialWaterVolume(firstTank:firstTank+3)
disp('getNodeTankMinimumWaterVolume(:)');      d.getNodeTankMinimumWaterVolume(firstTank:firstTank+3)
disp('getNodeTankDiameter(:)');                d.getNodeTankDiameter(firstTank:firstTank+3)
disp('getNodeTankMixZoneVolume(:)');           d.getNodeTankMixZoneVolume(firstTank:firstTank+3)
disp('getNodeTankMaximumWaterVolume(:)');      d.getNodeTankMaximumWaterVolume(firstTank:firstTank+3)
disp('getNodeTankVolumeCurveIndex(:)');        d.getNodeTankVolumeCurveIndex(firstTank:firstTank+3)
disp('getNodeTankMinimumWaterLevel(:)');       d.getNodeTankMinimumWaterLevel(firstTank:firstTank+3)
disp('getNodeTankMaximumWaterLevel(:)');       d.getNodeTankMaximumWaterLevel(firstTank:firstTank+3)
disp('getNodeTankMixingFraction(:)');          d.getNodeTankMixingFraction(firstTank:firstTank+3)
disp('getNodeTankBulkReactionCoeff(:)');       d.getNodeTankBulkReactionCoeff(firstTank:firstTank+3)
disp('getNodeTankVolume(:)');                  d.getNodeTankVolume(firstTank:firstTank+3)
disp('getNodeTankCanOverFlow(:)');             d.getNodeTankCanOverFlow(firstTank:firstTank+3)
disp('getNodeTankMixingModelCode(:)');         d.getNodeTankMixingModelCode(firstTank:firstTank+3)
disp('getNodeTankMixingModelType(:)');         d.getNodeTankMixingModelType(firstTank:firstTank+3)
disp('getNodeTankMixingModelCode(:)');          d.getNodeTankMixingModelCode(firstTank:firstTank+3)
d.unload;

d = epanet('Net3.inp');
firstPump = d.getLinkPumpIndex(1)

%% All pumps
disp('Running getLinkPumpPower');         d.getLinkPumpPower'
disp('Running getLinkPumpHCurve');        d.getLinkPumpHCurve'
disp('Running getLinkPumpECurve');        d.getLinkPumpECurve'
disp('Running getLinkPumpECost');         d.getLinkPumpECost'
disp('Running getLinkPumpEPat');          d.getLinkPumpEPat'
disp('Running getLinkPumpPatternIndex');  d.getLinkPumpPatternIndex'

%% Single pump (118)
disp('Running getLinkPumpPower(118)');         d.getLinkPumpPower(firstPump)
disp('Running getLinkPumpHCurve(118)');        d.getLinkPumpHCurve(firstPump)
disp('Running getLinkPumpECurve(118)');        d.getLinkPumpECurve(firstPump)
disp('Running getLinkPumpECost(118)');         d.getLinkPumpECost(firstPump)
disp('Running getLinkPumpEPat(118)');          d.getLinkPumpEPat(firstPump)
disp('Running getLinkPumpPatternIndex(118)');  d.getLinkPumpPatternIndex(firstPump)

%% Range of pumps
disp('Running getLinkPumpPower(:)');         d.getLinkPumpPower(firstPump:firstPump+1)
disp('Running getLinkPumpHCurve(:)');        d.getLinkPumpHCurve(firstPump:firstPump+1)
disp('Running getLinkPumpECurve(:)');        d.getLinkPumpECurve(firstPump:firstPump+1)
disp('Running getLinkPumpECost(:)');         d.getLinkPumpECost(firstPump:firstPump+1)
disp('Running getLinkPumpEPat(:)');          d.getLinkPumpEPat(firstPump:firstPump+1)
disp('Running getLinkPumpPatternIndex(:)');  d.getLinkPumpPatternIndex(firstPump:firstPump+1)


d.unload;