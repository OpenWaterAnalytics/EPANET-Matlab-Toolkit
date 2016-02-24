EPANET-Matlab-Class [EpaClass 2.1]
==================================

The `EPANET-Matlab Class` is an open-source software which operates within the Matlab environment, for providing a programming interface for the latest version of [EPANET](https://github.com/OpenWaterAnalytics/epanet), a hydraulic and quality modeling software created by the US EPA, with Matlab, a  high-level technical computing software. The goal of the Matlab Class is to serve as a common programming framework for research and development in the growing field of smart water networks.

The `EPANET-Matlab Class` features easy to use commands/wrappers for viewing, modifying, simulating and plotting results produced by the EPANET libraries.  

# Requirements #
* [Matlab](http://www.mathworks.com/)
* [EPANET](https://github.com/OpenWaterAnalytics/epanet) The OpenWaterAnalytics (OWA) EPANET version is recommended instead of the one provided by the EPA.gov website, as a number of bugs have been fixed in the OWA version. 

# How to use the class #
Download the folder `Download ZIP`, set the run path in Matlab within the saved folder, and run `RunTests.m`. This will execute all the commands which have been implemented in the Class.

Example:

d=epanet(networks/Net1_Rossman2000.inp)

d.getNodeCount

d.getNodeElevations


# How to fix/report bugs #
To fix a bug `Fork` the `EPANET-Matlab Class`, `Edit` the code and make the appropriate change, and then `Pull` it so that we evaluate it. 

Keep in mind that some bugs may exist in the `EPANET` libraries, in case you are not receiving the expected results.

# Licenses #
* `EPANET`: Public Domain
* `EPANET-MSX`: GNU Lesser General Public License
* `EPANET-Matlab Class`: EUPL 

# Acknowledgements #
* Marios Kyriakou
* Demetrios Eliades

The `EPANET-Matlab Class` is based/inspired on the [EPANET-Matlab Toolkit](http://www.mathworks.com/matlabcentral/fileexchange/25100-epanet-matlab-toolkit) as well as the OpenWaterAnalytics [EPANET-Matlab Wrappers](https://github.com/OpenWaterAnalytics/epanet-matlab)

# List of Matlab Class Functions #
|Function|Description|
|--------|-----------|
|String getError(Integer)| Returns the description of an error code| 
|LoadFile
|addCurve
|addPattern
|addlistener
|closeHydraulicAnalysis
|closeNetwork
|closeQualityAnalysis
|getComputedHydraulicTimeSeries
|getComputedQualityTimeSeries
|getConnectivityMatrix
|getControlRulesCount
|getControls
|getCurveCount
|getCurveIndex
|getCurveLengths
|getCurveNameID
|getCurveValue
|getCurveXY
|getENfunctionsImpemented
|getError
|getFlowUnits
|getHeadCurveIndex
|getLibFunctions
|getLinkBulkReactionCoeff
|getLinkCount
|getLinkDiameter
|getLinkFlows
|getLinkHeadloss
|getLinkIndex
|getLinkInitialSetting
|getLinkInitialStatus
|getLinkLength
|getLinkMinorLossCoeff
|getLinkNameID
|getLinkNodesIndex
|getLinkPipeCount
|getLinkPipeIndex
|getLinkPipeNameID
|getLinkPumpCount
|getLinkPumpEnergy
|getLinkPumpIndex
|getLinkPumpNameID
|getLinkPumpPatternIndex
|getLinkPumpPatternNameID
|getLinkPumpType
|getLinkPumpTypeCode
|getLinkQuality
|getLinkRoughnessCoeff
|getLinkSettings
|getLinkStatus
|getLinkType
|getLinkTypeIndex
|getLinkValveCount
|getLinkValveIndex
|getLinkValveNameID
|getLinkVelocity
|getLinkWallReactionCoeff
|getNodeActualDemand
|getNodeActualDemandSensingNodes
|getNodeActualQuality
|getNodeActualQualitySensingNodes
|getNodeBaseDemands
|getNodeCoordinates
|getNodeCount
|getNodeDemandPatternIndex
|getNodeDemandPatternNameID
|getNodeDemandPatternsIndex
|getNodeElevations
|getNodeEmitterCoeff
|getNodeHydaulicHead
|getNodeIndex
|getNodeInitialQuality
|getNodeJunctionCount
|getNodeJunctionIndex
|getNodeJunctionNameID
|getNodeMassFlowRate
|getNodeNameID
|getNodeNumDemandCategories
|getNodePressure
|getNodeReservoirCount
|getNodeReservoirIndex
|getNodeReservoirNameID
|getNodeSourcePatternIndex
|getNodeSourceQuality
|getNodeSourceType
|getNodeTankBulkReactionCoeff
|getNodeTankCount
|getNodeTankDiameter
|getNodeTankIndex
|getNodeTankInitialLevel
|getNodeTankInitialWaterVolume
|getNodeTankMaxVolume
|getNodeTankMaximumWaterLevel
|getNodeTankMinimumFraction
|getNodeTankMinimumWaterLevel
|getNodeTankMinimumWaterVolume
|getNodeTankMixZoneVolume
|getNodeTankMixingModelCode
|getNodeTankMixingModelType
|getNodeTankMixiningModel
|getNodeTankNameID
|getNodeTankReservoirCount
|getNodeTankVolume
|getNodeTankVolumeCurveIndex
|getNodeType
|getNodeTypeIndex
|getNodesConnectingLinksID
|getOptionsAccuracyValue
|getOptionsEmitterExponent
|getOptionsMaxTrials
|getOptionsPatternDemandMultiplier
|getOptionsQualityTolerance
|getPattern
|getPatternAveragePatternValue
|getPatternCount
|getPatternIndex
|getPatternLengths
|getPatternNameID
|getPatternValue
|getQualityCode
|getQualityInfo
|getQualityTraceNodeIndex
|getQualityType
|getStatistic
|getTimeHTime
|getTimeHaltFlag
|getTimeHydraulicStep
|getTimeNextEvent
|getTimePatternStart
|getTimePatternStep
|getTimeQualityStep
|getTimeReportingPeriods
|getTimeReportingStart
|getTimeReportingStep
|getTimeRuleControlStep
|getTimeSimulationDuration
|getTimeStartTime
|getTimeStatisticsIndex
|getTimeStatisticsType
|getVersion
initializeHydraulicAnalysis
initializeQualityAnalysis
nextHydraulicAnalysisStep
nextQualityAnalysisStep
openHydraulicAnalysis
openQualityAnalysis
runHydraulicAnalysis
runQualityAnalysis
saveHydraulicFile
saveHydraulicsOutputReportingFile
saveInputFile
plot
readInpFile
|setControl
|setCurve
|setCurveValue
|setLinkBulkReactionCoeff
|setLinkDiameter
|setLinkInitialSetting
|setLinkInitialStatus
|setLinkLength
|setLinkMinorLossCoeff
|setLinkRoughnessCoeff
|setLinkSettings
|setLinkStatus
|setLinkWallReactionCoeff
|setNodeBaseDemands
|setNodeCoordinates
|setNodeDemandPatternIndex
|setNodeElevations
|setNodeEmitterCoeff
|setNodeInitialQuality
|setNodeSourcePatternIndex
|setNodeSourceQuality
|setNodeSourceType
|setNodeTankBulkReactionCoeff
|setNodeTankDiameter
|setNodeTankInitialLevel
|setNodeTankMaximumWaterLevel
|setNodeTankMinimumFraction
|setNodeTankMinimumWaterLevel
|setNodeTankMinimumWaterVolume
|setNodeTankMixingModelType
|setOptionsAccuracyValue
|setOptionsEmitterExponent
|setOptionsMaxTrials
|setOptionsPatternDemandMultiplier
|setOptionsQualityTolerance
|setPattern
|setPatternMatrix
|setPatternValue
|setQualityType
|setReport
|setReportFormatReset
|setReportStatus
|setTimeHTime
|setTimeHaltFlag
|setTimeHydraulicStep
|setTimePatternStart
|setTimePatternStep
|setTimeQualityStep
|setTimeReportingStart
|setTimeReportingStep
|setTimeRuleControlStep
|setTimeSimulationDuration
|setTimeStatisticsType
solveCompleteHydraulics
solveCompleteQuality
stepQualityAnalysisTimeLeft
unload
useHydraulicFile
writeLineInReportFile
writeReport

msx
MsxAddPattern
MsxInitializeQualityAnalysis
MsxPlotConcentrationSpeciesOfLinks
MsxPlotConcentrationSpeciesOfNodes
MsxSaveFile
MsxSaveQualityFile
MsxSolveCompleteHydraulics
MsxSolveCompleteQuality
MsxStepQualityAnalysisTimeLeft
MsxUnload
MsxUseHydraulicFile
MsxWriteReport
|getMsxAreaUnits
|getMsxAtol
|getMsxCompiler
|getMsxComputedQualityLink
|getMsxComputedQualityNode
|getMsxConstantsCount
|getMsxConstantsIndex
|getMsxConstantsNameID
|getMsxConstantsValue
|getMsxCoupling
|getMsxEquationsPipes
|getMsxEquationsTanks
|getMsxEquationsTerms
|getMsxError
|getMsxLinkInitqualValue
|getMsxNodeInitqualValue
|getMsxParametersCount
|getMsxParametersIndex
|getMsxParametersNameID
|getMsxParametersPipesValue
|getMsxParametersTanksValue
|getMsxPattern
|getMsxPatternValue
|getMsxPatternsCount
|getMsxPatternsIndex
|getMsxPatternsLengths
|getMsxPatternsNameID
|getMsxRateUnits
|getMsxRtol
|getMsxSolver
|getMsxSourceLevel
|getMsxSourceNodeNameID
|getMsxSourcePatternIndex
|getMsxSourceType
|getMsxSources
|getMsxSpeciesATOL
|getMsxSpeciesConcentration
|getMsxSpeciesCount
|getMsxSpeciesIndex
|getMsxSpeciesNameID
|getMsxSpeciesRTOL
|getMsxSpeciesType
|getMsxSpeciesUnits
|getMsxTimeStep
|setMsxAreaUnitsCM2
|setMsxAreaUnitsFT2
|setMsxAreaUnitsM2
|setMsxAtol
|setMsxCompilerGC
|setMsxCompilerNONE
|setMsxCompilerVC
|setMsxConstantsValue
|setMsxCouplingFULL
|setMsxCouplingNONE
|setMsxLinkInitqualValue
|setMsxNodeInitqualValue
|setMsxParametersPipesValue
|setMsxParametersTanksValue
|setMsxPattern
|setMsxPatternMatrix
|setMsxPatternValue
|setMsxRateUnitsDAY
|setMsxRateUnitsHR
|setMsxRateUnitsMIN
|setMsxRateUnitsSEC
|setMsxRtol
|setMsxSolverEUL
|setMsxSolverRK5
|setMsxSolverROS2
|setMsxSources
|setMsxTimeStep

BinClose
BinUpdateClass
Binplot
addBinControl
addBinCurveEfficiency
addBinCurveHeadloss
addBinCurvePump
addBinCurveVolume
addBinJunction
addBinPattern
addBinPipe
addBinPump
addBinReservoir
addBinTank
addBinValveFCV
addBinValveGPV
addBinValvePBV
addBinValvePRV
addBinValvePSV
addBinValveTCV
remAddBinCurvesID
removeBinControlLinkID
removeBinControlNodeID
removeBinCurveID
removeBinLinkID
removeBinNodeID
removeBinRulesControlLinkID
removeBinRulesControlNodeID
saveBinInpFile
|getBinComputedAllParameters
|getBinComputedAverageBulkReactionRate
|getBinComputedAverageCostPerDay
|getBinComputedAverageEfficiency
|getBinComputedAverageKwatts
|getBinComputedAverageKwattsOrMillionGallons
|getBinComputedAverageSourceInflow
|getBinComputedAvera|getankReactionRate
|getBinComputedAverageWallReactionRate
|getBinComputedLinkFlow
|getBinComputedLinkFrictionFactor
|getBinComputedLinkHeadloss
|getBinComputedLinkQuality
|getBinComputedLinkReactionRate
|getBinComputedLinkSetting
|getBinComputedLinkStatus
|getBinComputedLinkVelocity
|getBinComputedNodeDemand
|getBinComputedNodeHead
|getBinComputedNodePressure
|getBinComputedNodeQuality
|getBinComputedPeakKwatts
|getBinComputedPumpIndexListLinks
|getBinComputedPumpUtilization
|getBinControlsInfo
|getBinCoordinatesSection
|getBinCurvesInfo
|getBinDiameterEachLink
|getBinElevationEachNode
|getBinLengthEachLink
|getBinLinkIndex
|getBinLinkNameID
|getBinLinksInfo
|getBinNodeCoordinates
|getBinNodeIndex
|getBinNodeNameID
|getBinNodeSourceInfo
|getBinNodesInfo
|getBinNumberReportingPeriods
|getBinOptionsInfo
|getBinPatternIndex
|getBinPatternsInfo
|getBinRulesControlsInfo
|getBinRulesSection
|getBinSimulationDuration
|getBinTimesInfo
|getBinUnits
|setBinFlowUnitsAFD
|setBinFlowUnitsCFS
|setBinFlowUnitsCMD
|setBinFlowUnitsCMH
|setBinFlowUnitsGPM
|setBinFlowUnitsIMGD
|setBinFlowUnitsLPM
|setBinFlowUnitsLPS
|setBinFlowUnitsMGD
|setBinFlowUnitsMLD
|setBinHeadlossCM
|setBinHeadlossDW
|setBinHeadlossHW
|setBinLinkGlobalBulkReactionCoeff
|setBinLinkGlobalWallReactionCoeff
|setBinLinkPipeDiameters
|setBinLinkPipeLengths
|setBinLinkPipeMinorLoss
|setBinLinkPipeRoughness
|setBinLinkPipeStatus
|setBinLinkPipesParameters
|setBinLinkPumpStatus
|setBinLinkReactionCoeff
|setBinLinkValvesParameters
|setBinNodeDemandPatternNameID
|setBinNodeInitialQuality
|setBinNodeJunctionElevation
|setBinNodeJunctionsBaseDemands
|setBinNodeJunctionsParameters
|setBinNodeResDemandPatternNameID
|setBinNodeReservoirElevation
|setBinNodeReservoirParameters
|setBinNodeSourceQuality
|setBinNodeTankDiameter
|setBinNodeTankElevation
|setBinNodeTankInitLevel
|setBinNodeTankMaxLevel
|setBinNodeTankMinLevel
|setBinNodeTankMinVol
|setBinNodeTankParameters
|setBinPattern
|setBinQualityAge
|setBinQualityChem
|setBinQualityNone
|setBinQualityTrace
|setBinTimeHydraulicStep
|setBinTimePatternStart
|setBinTimePatternStep
|setBinTimeQualityStep
|setBinTimeReportingStart
|setBinTimeReportingStep
|setBinTimeSimulationDuration
|setBinTimeStatisticsAverage
|setBinTimeStatisticsMaximum
|setBinTimeStatisticsMinimum
|setBinTimeStatisticsNone
|setBinTimeStatisticsRange


# List of EPANET 2.1 Functions Supported #

|Function|Description|
|--------|-----------|
|ENgetpumptype|Retrieves the type of a pump for specific link index|
|ENgetheadcurveindex|Retrieves index of a head curve for specific link index|
|ENsetcurvevalue|Sets x,y point for a specific point and curve|
|ENsetcurve|Sets x,y values for a specific curve|
|ENaddcurve|Adds a new curve appended to the end of the existing curves|
|ENgetcurvevalue|Retrieves x,y point for a specific point number and curve|
|ENgetcurvelen|Retrieves number of points in a curve|
|ENgetcurveid|Retrieves ID of a curve with specific index|
|ENgetcurveindex|Retrieves index of curve with specific ID|
|ENsetcoord|Sets coordinate x, y for a node|
|ENgetcoord|Retrieves coordinate x, y for a node|
|ENgetstatistic|Retrieves hydraulic simulation statistic|
|ENgetnumdemands|Retrieves the number of demand categories for a node|
|ENgetbasedemand|Retrieves the nodes base demand for a category|
|ENgetdemandpattern|Retrieves the index of a demand pattern for a specific demand category of a node|
|ENsetbasedemand|Sets the nodes base demand for a category|
|ENgetaveragepatternvalue|Retrieves the average value of a pattern|
