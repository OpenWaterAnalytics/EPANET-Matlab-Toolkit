EPANET-Matlab-Toolkit
==================================

The `EPANET-Matlab Toolkit` is an open-source software which operates within the Matlab environment, for providing a programming interface for the latest version of [EPANET](https://github.com/OpenWaterAnalytics/epanet), a hydraulic and quality modeling software created by the US EPA, with Matlab, a  high-level technical computing software. The goal of the Matlab Class is to serve as a common programming framework for research and development in the growing field of smart water networks.

The `EPANET-Matlab Toolkit` features easy to use commands/wrappers for viewing, modifying, simulating and plotting results produced by the EPANET libraries.  

# Requirements #
* [Matlab](http://www.mathworks.com/)
* [EPANET](https://github.com/OpenWaterAnalytics/epanet) The OpenWaterAnalytics (OWA) EPANET version is recommended instead of the one provided by the EPA.gov website, as a number of bugs have been fixed in the OWA version. 

# How to use the class #
Download the folder `Download ZIP`, set the run path in Matlab within the saved folder, and run `RunTests.m`. This will execute all the commands which have been implemented in the Class.

Example:

d=epanet('networks/Net1_Rossman2000.inp')

d.getNodeCount

d.getNodeElevations


# How to fix/report bugs #
To fix a bug `Fork` the `EPANET-Matlab Toolkit`, `Edit` the code and make the appropriate change, and then `Pull` it so that we evaluate it. 

Keep in mind that some bugs may exist in the `EPANET` libraries, in case you are not receiving the expected results.

# Licenses #
* `EPANET`: Public Domain
* `EPANET-MSX`: GNU Lesser General Public License
* `EPANET-Matlab Toolkit`: EUPL 

# Acknowledgements #
* Marios Kyriakou
* Demetrios Eliades

The `EPANET-Matlab Toolkit` is based/inspired on the [EPANET-Matlab Toolkit](http://www.mathworks.com/matlabcentral/fileexchange/25100-epanet-matlab-toolkit) as well as the OpenWaterAnalytics [EPANET-Matlab Wrappers](https://github.com/OpenWaterAnalytics/epanet-matlab)

# List of Matlab Class Functions #
|Function|Description|
|---------|---------|
|Struct epanet(String)| Load Input file and open the EPANET Toolkit system|
|unload()|Unload library and close the EPANET Toolkit system|
|Integer LoadFile(String)| Open the EPANET Toolkit system|
|String getError(Integer)| Returns the description of an error code| 
|Struct getComputedHydraulicTimeSeries()|Computed Hydraulic Time Series|
|Struct getComputedQualityTimeSeries()|Computed Quality Time Series|
|Array getConnectivityMatrix|Return connectivity matrix of the network|
|getControlRulesCount|Retrieves the number of control rules|
|getControls|Retrieves the controls|
|getCurveCount|Retrieves the number of curves|
|getCurveIndex|Retrieves index of curve with specific ID|
|getCurveLengths|Retrieves number of points in a curve|
|getCurveNameID|Retrieves curve id|
|getCurveValue|Retrieves (x,y) values of specific curve index|
|getCurveXY|Retrieves (x,y) values of all curves|
|getENfunctionsImpemented|Retrieves the epanet functions that have been developed|
|getFlowUnits|Retrieves the units used to express all flow rates|
|getHeadCurveIndex|Retrieves index of a head curve for specific link index|
|getLibFunctions|Retrieves the functions of DLL|
|getLinkBulkReactionCoeff|Retrieves the value of all link bulk reaction coefficients|
|getLinkCount|Retrieves the number of links|
|getLinkDiameter|Retrieves the value of all link diameters|
|getLinkFlows|Retrieves the value of all computed link flow rates|
|getLinkHeadloss|Retrieves the value of all computed link headloss|
|getLinkIndex|Retrieves the indices of all links, or the indices of an ID set of links|
|getLinkInitialSetting|Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves|
|getLinkInitialStatus|Retrieves the value of all link initial status|
|getLinkLength|Retrieves the value of all link lengths|
|getLinkMinorLossCoeff|Retrieves the value of all link minor loss coefficients|
|getLinkNameID|Retrieves the ID label(s) of all links, or the IDs of an index set of links|
|getLinkNodesIndex|Retrieves the indexes of the from/to nodes of all links|
|getLinkPipeCount|Retrieves the number of pipes|
|getLinkPipeIndex|Retrieves the indices of pipes|
|getLinkPipeNameID|Retrieves the pipe IDs|
|getLinkPumpCount|Retrieves the number of pumps|
|getLinkPumpEnergy|Retrieves the value of all computed energy in kwatts|
|getLinkPumpIndex|Retrieves the indices of pumps|
|getLinkPumpNameID|Retrieves the pump IDs|
|getLinkPumpPatternIndex|Retrieves the pump pattern indices|
|getLinkPumpPatternNameID|Retrieves the pump pattern IDs|
|getLinkPumpType|Retrieves the type of a pump for specific link index|
|getLinkPumpTypeCode|Retrieves the type code of a pump for specific link index|
|getLinkQuality|Retrieves the quality of links|
|getLinkRoughnessCoeff|Retrieves the value of all link roughness|
|getLinkSettings|Retrieves the value of all computed link roughness for pipes or actual speed for pumps or actual setting for valves|
|getLinkStatus|Retrieves the value of all computed link status (0 = closed, 1 = open)|
|getLinkType|Retrieves the link-type for all links|
|getLinkTypeIndex|Retrieves the link-type code for all links.|
|getLinkValveCount|Retrieves the number of valves|
|getLinkValveIndex|Retrieves the indices of valves|
|getLinkValveNameID|Retrieves the valve IDs|
|getLinkVelocity|Retrieves the value of all computed link velocities|
|getLinkWallReactionCoeff|Retrieves the value of all link wall reaction coefficients|
|getNodeActualDemand|Retrieves the computed value of all actual demands|
|getNodeActualDemandSensingNodes|Retrieves the computed demand values at some sensing nodes|
|getNodeActualQuality|Retrieves the computed values of the actual quality for all nodes|
|getNodeActualQualitySensingNodes|Retrieves the computed quality values at some sensing nodes|
|getNodeBaseDemands|Retrieves the value of all node base demands|
|getNodeCoordinates|Retrieves coordinate x, y, and x, y vertices for a node|
|getNodeCount|Retrieves the number of nodes|
|getNodeDemandPatternIndex|Retrieves the value of all node demand pattern indices|
|getNodeDemandPatternNameID|Retrieves the value of all node demand pattern IDs|
|getNodeDemandPatternsIndex|Retrieves the value of all node demand pattern indices dev-2.1|
|getNodeElevations|Retrieves the value of all node elevations|
|getNodeEmitterCoeff|Retrieves the value of all node emmitter coefficients|
|getNodeHydaulicHead|Retrieves the computed values of all hydraulic heads|
|getNodeIndex|Retrieves the indices of all nodes or some nodes with a specified ID|
|getNodeInitialQuality|Retrieves the value of all node initial quality|
|getNodeJunctionCount|Retrieves the number of junctions|
|getNodeJunctionIndex|Retrieves the junctions indices|
|getNodeJunctionNameID|Retrieves the junctions IDs|
|getNodeMassFlowRate|Retrieves the computed mass flow rates per minute of chemical sources|
|getNodeNameID|Retrieves the ID label of all nodes or some nodes with a specified index|
|getNodeNumDemandCategories|Retrieves the number of demand categories for a node|
|getNodePressure|Retrieves the computed values of all node pressures|
|getNodeReservoirCount|Retrieves the number of reservoirs|
|getNodeReservoirIndex|Retrieves the indices of reservoirs|
|getNodeReservoirNameID|Retrieves the reservoirs IDs|
|getNodeSourcePatternIndex|Retrieves the value of all node source pattern index|
|getNodeSourceQuality|Retrieves the value of all nodes source quality|
|getNodeSourceType|Retrieves the value of all node source type|
|getNodeTankBulkReactionCoeff|Retrieves the tank bulk rate coefficient|
|getNodeTankCount|Retrieves the number of tanks|
|getNodeTankDiameter|Retrieves the tank diameters|
|getNodeTankIndex|Retrieves the indices of tanks|
|getNodeTankInitialLevel|Retrieves the value of all tank initial water levels|
|getNodeTankInitialWaterVolume|Retrieves the tank initial volume|
|getNodeTankMaxVolume|Retrieves maximum water volume|
|getNodeTankMaximumWaterLevel|Retrieves the tank maximum water level|
|getNodeTankMinimumFraction|Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank|
|getNodeTankMinimumWaterLevel|Retrieves the tank minimum water level|
|getNodeTankMinimumWaterVolume|Retrieves the tank minimum volume|
|getNodeTankMixZoneVolume
|getNodeTankMixingModelCode
|getNodeTankMixingModelType
|getNodeTankMixiningModel
|getNodeTankNameID
|getNodeTankReservoirCount
|getNodeTankVolumer
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
|addCurve
|addPattern
|closeHydraulicAnalysis
|closeNetwork
|closeQualityAnalysis
|initializeHydraulicAnalysis
|initializeQualityAnalysis
|nextHydraulicAnalysisStep
|nextQualityAnalysisStep
|openHydraulicAnalysis
|openQualityAnalysis
|runHydraulicAnalysis
|runQualityAnalysis
|saveHydraulicFile
|saveHydraulicsOutputReportingFile
|saveInputFile
|plot
|readInpFile
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
|solveCompleteHydraulics
|solveCompleteQuality
|stepQualityAnalysisTimeLeft
|useHydraulicFile
|writeLineInReportFile
|writeReport
|<b> MSX Functions </b>
|msx
|MsxAddPattern
|MsxInitializeQualityAnalysis
|MsxPlotConcentrationSpeciesOfLinks
|MsxPlotConcentrationSpeciesOfNodes
|MsxSaveFile
|MsxSaveQualityFile
|MsxSolveCompleteHydraulics
|MsxSolveCompleteQuality
|MsxStepQualityAnalysisTimeLeft
|MsxUnload
|MsxUseHydraulicFile
|MsxWriteReport
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
|<b>Bin Functions</b>
|BinClose
|BinUpdateClass
|Binplot
|addBinControl
|addBinCurveEfficiency
|addBinCurveHeadloss
|addBinCurvePump
|addBinCurveVolume
|addBinJunction
|addBinPattern
|addBinPipe
|addBinPump
|addBinReservoir
|addBinTank
|addBinValveFCV
|addBinValveGPV
|addBinValvePBV
|addBinValvePRV
|addBinValvePSV
|addBinValveTCV
|remAddBinCurvesID
|removeBinControlLinkID
|removeBinControlNodeID
|removeBinCurveID
|removeBinLinkID
|removeBinNodeID
|removeBinRulesControlLinkID
|removeBinRulesControlNodeID
|saveBinInpFile
|getBinComputedAllParameters
|getBinComputedAverageBulkReactionRate
|getBinComputedAverageCostPerDay
|getBinComputedAverageEfficiency
|getBinComputedAverageKwatts
|getBinComputedAverageKwattsOrMillionGallons
|getBinComputedAverageSourceInflow                                                                                             
|getBinComputedAverageTankReactionRate
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
