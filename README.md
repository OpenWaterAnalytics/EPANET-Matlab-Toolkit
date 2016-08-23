EPANET-MATLAB-Toolkit
==================================

The `EPANET-Matlab Toolkit` is an open-source software, originally developed by the [KIOS Research Center for Intelligent Systems and Networks of the University of Cyprus](http://www.kios.ucy.ac.cy/) which operates within the Matlab environment, for providing a programming interface for the latest version of [EPANET](https://github.com/OpenWaterAnalytics/epanet), a hydraulic and quality modeling software created by the US EPA, with Matlab, a  high-level technical computing software. The goal of the EPANET Matlab Toolkit is to serve as a common programming framework for research and development in the growing field of smart water networks.

The `EPANET-Matlab Toolkit` features easy to use commands/wrappers for viewing, modifying, simulating and plotting results produced by the EPANET libraries.  

For support, please use the OWA community forum : http://community.wateranalytics.org/

# Requirements #
* [Matlab](http://www.mathworks.com/)
* [EPANET 2.1](https://github.com/OpenWaterAnalytics/epanet) 

# How to use the Toolkit #
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

# Contributors #
* Marios Kyriakou, [KIOS Research Center for Intelligent Systems and Networks, University of Cyprus](http://www.kios.ucy.ac.cy/)
* Demetrios Eliades, [KIOS Research Center for Intelligent Systems and Networks, University of Cyprus](http://www.kios.ucy.ac.cy/c)

The `EPANET-Matlab Toolkit` is based/inspired on the [EPANET-Matlab Toolkit](http://www.mathworks.com/matlabcentral/fileexchange/25100-epanet-matlab-toolkit) as well as the OpenWaterAnalytics [EPANET-Matlab Wrappers](https://github.com/OpenWaterAnalytics/epanet-matlab)
# List of Matlab Class Functions #
|Function|Description|
|---------|---------|
|epanet| Load Input file and open the EPANET Toolkit system|
|unload|Unload library and close the EPANET Toolkit system|
|loadEPANETFile| Open the EPANET Toolkit system|
|getError| Returns the description of an error code| 
|getComputedHydraulicTimeSeries|Computed Hydraulic Time Series|
|getComputedQualityTimeSeries|Computed Quality Time Series|
|getConnectivityMatrix|Return connectivity matrix of the network|
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
|getNodePatternIndex|Retrieves the value of all node pattern indices|
|getNodeDemandPatternIndex|Retrieves the value of all node demand pattern indices|
|getNodeDemandPatternNameID|Retrieves the value of all node demand pattern IDs|
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
|getNodeTankMixZoneVolume|Retrieves the tank mixing zone volume|
|getNodeTankMixingModelCode|Retrieves the tank mixing model code|
|getNodeTankMixingModelType|Retrieves the tank mixing model type (mix1, mix2, fifo, lifo)|
|getNodeTankMixiningModel|Retrieves the tank mixing model|
|getNodeTankNameID|Retrieves the tanks IDs|
|getNodeTankReservoirCount|Retrieves the number of tanks|
|getNodeTankVolume|Retrieves the tank volume|
|getNodeTankVolumeCurveIndex|Retrieves the tank volume curve index|
|getNodeType|Retrieves the node-type for all nodes|
|getNodeTypeIndex|Retrieves the node code-index for all nodes|
|getNodesConnectingLinksID|Retrieves the id of the from/to nodes of all links|
|getOptionsAccuracyValue|Retrieve the analysis convergence criterion (0.001)|
|getOptionsEmitterExponent|Retrieve power exponent for the emmitters (0.5)|
|getOptionsMaxTrials|Retrieve maximum number of analysis trials|
|getOptionsPatternDemandMultiplier|Retrieve the demand multiplier (x1)|
|getOptionsQualityTolerance|Retrieve the water quality analysis tolerance|
|getPattern|Retrieves the multiplier factor for all patterns and all times|
|getPatternAveragePatternValue|Retrieves the average value of a pattern|
|getPatternCount|Retrieves the number of patterns|
|getPatternIndex|Retrieves the index of all or some time patterns IDs|
|getPatternLengths|Retrieves the number of time periods in all or some patterns|
|getPatternNameID|Retrieves the patterns IDs|
|getPatternValue|Retrieves the multiplier factor for a certain pattern and time|
|getQualityCode|Retrieves the code of water quality analysis type|
|getQualityInfo|Retrieves the quality info - bug in ENgetqualinfo|
|getQualityTraceNodeIndex|Retrieves the trace node index of water quality analysis type|
|getQualityType|Retrieves the type of water quality analysis type|
|getStatistic|Retrieves hydraulic simulation statistic|
|getTimeHTime|Retrieves the number of htime|
|getTimeHaltFlag|Retrieves the number of  halt flag|
|getTimeHydraulicStep|Retrieves the value of the hydraulic time step|
|getTimeNextEvent|Retrieves the number of next event|
|getTimePatternStart|Retrieves the value of pattern start time|
|getTimePatternStep|Retrieves the value of the pattern time step|
|getTimeQualityStep|Retrieves the value of the water quality time step|
|getTimeReportingPeriods|Retrieves the number of reporting periods saved to the binary|
|getTimeReportingStart|Retrie ves the value of the reporting start time|
|getTimeReportingStep|Retrieves the value of the reporting time step|
|getTimeRuleControlStep|Retrieves the time step for evaluating rule-based controls|
|getTimeSimulationDuration|Retrieves the value of simulation duration|
|getTimeStartTime|Retrieves the number of start time|
|getTimeStatisticsType|Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')|
|getTimeStatisticsIndex|Retrieves the type of time series post-processing|
|getVersion|Retrieve the current EPANET version of DLL|
|getNodesInfo|Retrieves nodes info e.g. elevations, demand pattern indices, emitter coeff. , initial quality, source quality, source pattern indices, source type code, type indices|
|getLinksInfo|Retrieves links info e.g. diameters, lengths, roughness coeff. , minor loss coeff. , initial status, initial settings, bulk reaction coeff. , wall reaction coeff. , nodes connecting link indices, type indices|
|addCurve|Adds a new curve appended to the end of the existing curves|
|addPattern|Adds a new time pattern to the network|
|closeHydraulicAnalysis|Closes the hydraulic analysis system, freeing all allocated memory|
|closeNetwork|Closes down the Toolkit system|
|closeQualityAnalysis|Closes the water quality analysis system, freeing all allocated memory|
|initializeHydraulicAnalysis|Initializes storage tank levels, link status and settings, and the simulation clock time prior to running a hydraulic analysis|
|initializeQualityAnalysis|Initializes water quality and the simulation clock time prior to running a water quality analysis|
|nextHydraulicAnalysisStep|Determines the length of time until the next hydraulic event occurs in an extended period simulation|
|nextQualityAnalysisStep|Advances the water quality simulation to the start of the next hydraulic time period|
|openHydraulicAnalysis|Opens the hydraulics analysis system|
|openQualityAnalysis|Opens the water quality analysis system|
|runHydraulicAnalysis|Runs a single period hydraulic analysis, retrieving the current simulation clock time t|
|runQualityAnalysis|Makes available the hydraulic and water quality results that occur at the start of the next time period of a water quality analysis, where the start of the period is returned in t|
|saveHydraulicFile|Saves the current contents of the binary hydraulics file to a file|
|saveHydraulicsOutputReportingFile|Transfers results of a hydraulic simulation from the binary Hydraulics file to the binary Output file, where results are only reported at uniform reporting intervals|
|saveInputFile|Writes all current network input data to a file using the format of an EPANET input file|
|plot|Plot the network input file|
|setControl|Sets the parameters of a simple control statement|
|setCurve|Sets x,y values for a specific curve|
|setCurveValue|Retrieves x,y point for a specific point number and curve|
|setLinkBulkReactionCoeff|Sets the values of bulk reactions|
|setLinkDiameter|Sets the values of diameters|
|setLinkInitialSetting|Sets the values of initial settings|
|setLinkInitialStatus|Sets the values of initial status|
|setLinkLength|Sets the values of lengths|
|setLinkMinorLossCoeff|Sets the values of minor loss coeff.|
|setLinkRoughnessCoeff|Sets the values of roughness coeff.|
|setLinkSettings|Sets the values of settings|
|setLinkStatus|Sets the values of status|
|setLinkWallReactionCoeff|Sets the values of wall reactions|
|setNodeBaseDemands|Sets the values of demands|
|setNodeCoordinates|Sets node coordinates|
|setNodeDemandPatternIndex|Sets the values of demand pattern indices|
|setNodeElevations|Sets the values of elevations|
|setNodeEmitterCoeff|Sets the values of emitter coeff.|
|setNodeInitialQuality|Sets the values of initial qualities|
|setNodeSourcePatternIndex|Sets the values of source pattern indices|
|setNodeSourceQuality|Sets the values of source qualities|
|setNodeSourceType|Sets the values of source types|
|setNodeTankBulkReactionCoeff|Sets the values of tank bulk reaction coeff.|
|setNodeTankDiameter|Sets the values of tanks diameter|
|setNodeTankInitialLevel|Sets the values of tanks initial level|
|setNodeTankMaximumWaterLevel|Sets the values of tanks maximum water level|
|setNodeTankMinimumWaterLevel|Sets the values of tanks minimum water level|
|setNodeTankMinimumFraction|Sets the values of tanks mix fraction|
|setNodeTankMinimumWaterVolume|Sets the values of tanks minimum water volume|
|setNodeTankMixingModelType|Sets the values of tanks model|
|setOptionsAccuracyValue|Sets the value of accurancy|
|setOptionsEmitterExponent|Sets the value of emitter exponent|
|setOptionsMaxTrials|Sets the value of max trials|
|setOptionsPatternDemandMultiplier|Sets the value of pattern demand multiplier|
|setOptionsQualityTolerance|Sets the value of tolerance|
|setPattern|Sets all of the multiplier factors for a specific time pattern|
|setPatternMatrix|Sets all of the multiplier factors for all patterns|
|setPatternValue|Sets the multiplier factor for a specific period within a time pattern|
|setQualityType|Sets the type of water quality analysis called for|
|setReport|Issues a report formatting command. Formatting commands are the same as used in the [REPORT] section of the EPANET Input file|
|setReportFormatReset|Clears any report formatting commands that either appeared in the [REPORT] section of the EPANET Input file or were issued with the ENsetreport function|
|setReportStatus|Sets the level of hydraulic status reporting|
|setTimeHTime|Sets the htime|
|setTimeHaltFlag|Sets the halt flag|
|setTimeHydraulicStep|Sets the hydraulic step|
|setTimePatternStart|Sets the pattern start|
|setTimePatternStep|Sets the pattern step|
|setTimeQualityStep|Sets the quality step|
|setTimeReportingStart|Sets the reporting start|
|setTimeReportingStep|Sets the reporting step|
|setTimeRuleControlStep|Sets the rule control step|
|setTimeSimulationDuration|Sets the simulation duration|
|setTimeStatisticsType|Sets the statistic type|
|solveCompleteHydraulics|Runs a complete hydraulic simulation with results for all time periods written to the binary Hydraulics file|
|solveCompleteQuality|Runs a complete water quality simulation with results at uniform reporting intervals written to EPANET's binary Output file|
|stepQualityAnalysisTimeLeft|Advances the water quality simulation one water quality time step. The time remaining in the overall simulation is returned in tleft|
|useHydraulicFile|Uses the contents of the specified file as the current binary hydraulics file|
|writeLineInReportFile|Writes a line of text to the EPANET report file|
|writeReport|Writes a formatted text report on simulation results to the Report file|
|<b> MSX Functions </b>
|loadMSXFile|Opens the EPANET-MSX toolkit system|
|addMSXPattern|Adds a new, empty MSX source time pattern to the project|
|initializeMSXQualityAnalysis|Initializes the MSX system before solving for water quality results in step-wise fashion|
|saveMSXFile|Saves the data associated with the current MSX project into a new MSX input file|
|saveMSXQualityFile|Saves water quality results computed for each node, link and reporting time period to a named binary file|
|solveMSXCompleteHydraulics|Solves for system hydraulics over the entire simulation period saving results to an internal scratch file|
|solveMSXCompleteQuality|Solves for water quality over the entire simulation period and saves the results to an internal scratch file|
|stepMSXQualityAnalysisTimeLeft|Advances the water quality solution through a single water quality time step when performing a step-wise simulation|
|writeMSXFile|Write a new MSX file|
|writeMSXReport|Writes water quality simulations results as instructed by the MSX input file to a text file|
|useMSXHydraulicFile|Uses a previously saved EPANET hydraulics file as the source of hydraulic information|
|plotMSXConcentrationSpeciesOfLinks|Plots the concentration species of links|
|plotMSXConcentrationSpeciesOfNodes|Plots the concentration species of nodes|
|runMSXexe|Writes water quality simulations results as instructed by the MSX input file to a text file using the epanetmsx.exe|
|unloadMSX|Closes the EPANET-MSX toolkit system|
|getMSXAtol|Retrieves the absolute concentration tolerance|
|getMSXRtol|Retrieves the relative concentration tolerance|
|getMSXComputedQualityLink|Retrieves the concentration of a chemical species at a specific link of the network at the current simulation time step|
|getMSXComputedQualityNode|Retrieves the concentration of a chemical species at a specific node of the network at the current simulation time step.|
|getMSXConstantsCount|Retrieves the number of constants|
|getMSXConstantsIndex|Retrieves the internal index number of constants (given its ID name)|
|getMSXConstantsNameID|Retrieves the ID name of constants (given its internal index number)|
|getMSXConstantsValue|Retrieves the value of a particular reaction constant|
|getMSXError|Returns the text for an error message given its error code|
|getMSXLinkInitqualValue|Retrieves the initial concentration of chemical species assigned to links of the pipe network|
|getMSXNodeInitqualValue|Retrieves the initial concentration of chemical species assigned to nodes|
|getMSXParametersCount|Retrieves the number of parameters|
|getMSXParametersIndex|Retrieves the indices of parameters|
|getMSXParametersNameID|Retrieves the ID name of parameters|
|getMSXParametersPipesValue|Retrieves the value of reaction parameters for pipes|
|getMSXParametersTanksValue|Retrieves the value of reaction parameters for tanks|
|getMSXPattern|Retrieves the multiplier factor for all patterns and all times|
|getMSXPatternValue|Retrieves the multiplier at a specific time period for a given source time pattern|
|getMSXPatternsCount|Retrieves the number of patterns|
|getMSXPatternsIndex|Retrieves the indices of patterns|
|getMSXPatternsLengths|Retrieves the number of time periods in all or some patterns|
|getMSXPatternsNameID|Retrieves the patterns IDs|
|getMSXSourceLevel|Retrieves the value of all nodes source level|
|getMSXSourceNodeNameID|Retrieves the indices of parameters|
|getMSXSourcePatternIndex|Retrieves the value of all node source pattern index|
|getMSXSourceType|Retrieves the value of all node source type|
|getMSXSources|Retrieves the source info|
|getMSXSpeciesATOL|Retrieves the atol|
|getMSXSpeciesRTOL|Retrieves the rtol|
|getMSXSpeciesConcentration|Retrieves the concentration of chemical species for nodes and links|
|getMSXSpeciesCount|Retrieves the number of species|
|getMSXSpeciesIndex|Retrieves the indices of species|
|getMSXSpeciesNameID|Retrieves the species IDs|
|getMSXSpeciesType|Retrieves the value of all species|
|getMSXSpeciesUnits|Retrieves the species units|
|getMSXTimeStep|Retrieves the time step|
|getMSXRateUnits|Retrieves the rate/time units (SEC/MIN/HR/DAY)|
|getMSXAreaUnits|Retrieves the area units (FT2/M2/CM2)|
|getMSXCompiler|Retrieves the compiler (NONE/VC/GC)|
|getMSXCoupling|Retrieves the coupling (FULL/NONE)|
|getMSXEquationsPipes|Retrieves the species dynamics in pipes|
|getMSXEquationsTanks|Retrieves the species dynamics in tanks|
|getMSXEquationsTerms|Retrieves the species dynamics in terms|
|getMSXSolver|Retrieves the solver (EUL/RK5/ROS2)|
|setMSXAreaUnitsCM2|Sets area units to CM2|
|setMSXAreaUnitsFT2|Sets area units to FT2|
|setMSXAreaUnitsM2|Sets area units to M2|
|setMSXAtol|Sets the value of Atol|
|setMSXRtol|Sets the value of Rtol|
|setMSXCompilerGC|Sets compilet to GC|
|setMSXCompilerNONE|Sets compiler to None|
|setMSXCompilerVC|Sets compiler to VC|
|setMSXConstantsValue|Assigns a new value to a specific reaction constant|
|setMSXCouplingFULL|Sets coupling option to FULL|
|setMSXCouplingNONE|Sets coupling option to NONE|
|setMSXLinkInitqualValue|Assigns an initial concentration of chemical species to links|
|setMSXNodeInitqualValue|Assigns an initial concentration of chemical species to nodes|
|setMSXParametersPipesValue|Assigns a value to a particular reaction parameter for given pipes|
|setMSXParametersTanksValue|Assigns a value to a particular reaction parameter for given tanks|
|setMSXPattern|Sets all of the multiplier factors for a specific time pattern|
|setMSXPatternMatrix|Sets all of the multiplier factors for all patterns|
|setMSXPatternValue|Assigns a new value to the multiplier for a specific time period in a given MSX source time pattern|
|setMSXRateUnitsDAY|Sets rate units to DAY|
|setMSXRateUnitsHR|Sets rate units to HR|
|setMSXRateUnitsMIN|Sets rate units to MIN|
|setMSXRateUnitsSEC|Sets rate units to SEC|
|setMSXSolverEUL|Sets solver to EUL (standard Euler integrator)|
|setMSXSolverRK5|Sets solver to RK5 (Runge-Kutta 5th order integrator)|
|setMSXSolverROS2|Sets solver to ROS2 (2nd order Rosenbrock integrator)|
|setMSXSources|Sets the attributes of an external source of a particular chemical species to a specific node of the pipe network|
|setMSXTimeStep|Sets time step|
|<b>Bin Functions</b>
|BinClose|Close binary files and delete|
|BinUpdateClass|Run all bin functions and update the results|
|Binplot|Plot the network input file|
|addBinControl|Adds a new control to the network|
|addBinCurveEfficiency|Adds a new curve efficiency to the network|
|addBinCurveHeadloss|Adds a new curve headloss to the network|
|addBinCurvePump|Adds a new curve pump to the network|
|addBinCurveVolume|Adds a new curve volume to the network|
|addBinJunction|Adds a new junction to the network|
|addBinPattern|Adds a new time pattern to the network|
|addBinPipe|Adds a new pipe to the network|
|addBinPump|Adds a new pump to the network|
|addBinReservoir|Adds a new reservoir to the network|
|addBinTank|Adds a new tank to the network|
|addBinValveFCV|Adds a new valve FCV to the network|
|addBinValveGPV|Adds a new valve GPV to the network|
|addBinValvePBV|Adds a new valve PBV to the network|
|addBinValvePRV|Adds a new valve PRV to the network|
|addBinValvePSV|Adds a new valve PSV to the network|
|addBinValveTCV|Adds a new valve TCV to the network|
|removeBinControlLinkID|Removes a specific control based on link ID|
|removeBinControlNodeID|Removes a specific control based on node ID|
|removeBinCurveID|Removes a specific curve based on ID|
|removeBinLinkID|Removes a specific link based on ID|
|removeBinNodeID|Removes a specific node based on ID|
|removeBinRulesControlLinkID|Removes a specific rule based on link ID|
|removeBinRulesControlNodeID|Removes a specific rule based on node ID|
|saveBinInpFile|Writes all current network input data to a file using the format of an EPANET input file|
|getBinComputedAllParameters|Computes hydraulic and quality time series|
|getBinComputedAverageBulkReactionRate|Computes the average bulk reaction rate|
|getBinComputedAverageCostPerDay|Computes the average cost per day|
|getBinComputedAverageEfficiency|Computes the average efficiency|
|getBinComputedAverageKwatts|Computes the average Kwatts|
|getBinComputedAverageKwattsOrMillionGallons|Computes the average Kwatts or million gallons|
|getBinComputedAverageSourceInflow|Computes the average source inflow|                                                                                             
|getBinComputedAverageTankReactionRate|Computes the average tank reaction rate|
|getBinComputedAverageWallReactionRate|Computes average wall reaction rate|
|getBinComputedLinkFlow|Computes the flow of links|
|getBinComputedLinkFrictionFactor|Computes the link friction factor|
|getBinComputedLinkHeadloss|Computes the headloss of links|
|getBinComputedLinkQuality|Computes the quality of links|
|getBinComputedLinkReactionRate|Computes the reaction rate of links|
|getBinComputedLinkSetting|Computes the setting of links|
|getBinComputedLinkStatus|Computes the status of links|
|getBinComputedLinkVelocity|Computes the velocity of links|
|getBinComputedNodeDemand|Computes the demand of nodes|
|getBinComputedNodeHead|Computes the head of nodes|
|getBinComputedNodePressure|Computes the pressure of nodes|
|getBinComputedNodeQuality|Computes the quality of nodes|
|getBinComputedPeakKwatts|Computes the peak Kwatts|
|getBinComputedPumpIndexListLinks|Retrieves the pump indices|
|getBinComputedPumpUtilization|Computes the pump utilization|
|getBinDiameterEachLink|Retrieves the diameter of each link|
|getBinLengthEachLink|Retrieves the length of each link|
|getBinLinkIndex|Retrieves the indices of all links|
|getBinLinkNameID|Retrieves the ID label(s) of all links|
|getBinElevationEachNode|Retrieves the elevation of each node|
|getBinNodeCoordinates|Retrieves coordinate x, y, and x, y vertices for a node|
|getBinNodeIndex|Retrieves the indices of all nodes|
|getBinNodeNameID|Retrieves the ID label(s) of all nodes|
|getBinNumberReportingPeriods|Retrieves the number of reporting periods|
|getBinControlsInfo|Retrieves the controls info|
|getBinCurvesInfo|Retrieves the curves info|
|getBinLinksInfo|Retrieves the links info|
|getBinNodesInfo|Retrieves the nodes info|
|getBinNodeSourceInfo|Retrieves the sources info|
|getBinOptionsInfo|Retrieves the options info|
|getBinPatternsInfo|Retrieves the patterns info|
|getBinRulesControlsInfo|Retrieves the controls info|
|getBinTimesInfo|Retrieves the times info|
|getBinPatternIndex|Retrieves the indices of all patterns|
|getBinSimulationDuration|Retrieves the value of simulation duration|
|getBinUnits|Retrieves the units used to express all flow rates|
|setBinFlowUnitsAFD|Sets flow units to AFD|
|setBinFlowUnitsCFS|Sets flow units to CFS|
|setBinFlowUnitsCMD|Sets flow units to CMD|
|setBinFlowUnitsCMH|Sets flow units to CMH|
|setBinFlowUnitsGPM|Sets flow units to GPM|
|setBinFlowUnitsIMGD|Sets flow units to IMGD|
|setBinFlowUnitsLPM|Sets flow units to LPM|
|setBinFlowUnitsLPS|Sets flow units to LPS|
|setBinFlowUnitsMGD|Sets flow units to MGD|
|setBinFlowUnitsMLD|Sets flow units to MLD|
|setBinHeadlossCM|Sets headloss to C-M|
|setBinHeadlossDW|Sets headloss to D-W|
|setBinHeadlossHW|Sets headloss to H-W|
|setBinLinkGlobalBulkReactionCoeff|Sets the global bulk reaction rate coeff.|
|setBinLinkGlobalWallReactionCoeff|Sets the global wall reaction rate coeff.|
|setBinLinkPipeDiameters|Sets the values of pipe diameters|
|setBinLinkPipeLengths|Sets the values of pipe lengths|
|setBinLinkPipeMinorLoss|Sets the values of pipe minor losses|
|setBinLinkPipeRoughness|Sets the values of pipe roughness|
|setBinLinkPipeStatus|Sets the values of pipe status|
|setBinLinkPipesParameters|Sets the values of pipe parameters (diameters, lengths, minor losses, roughness, status)|
|setBinLinkPumpStatus|Sets the values of pump status|
|setBinLinkReactionCoeff|Sets the values of bulk and wall reaction coeff.|
|setBinLinkValvesParameters|Sets the values of valve parameters (diameters, types, settings, minor losses)|
|setBinNodeDemandPatternNameID|Sets the names of demand pattern IDs|
|setBinNodeInitialQuality|Sets the values of initial qualities|
|setBinNodeJunctionElevation|Sets the values of elevations for junctions|
|setBinNodeJunctionsBaseDemands|Sets the values of base demands|
|setBinNodeJunctionsParameters|Sets the values of junction parameters (elevations, base demands, demand patterns)|
|setBinNodeResDemandPatternNameID|Sets the names of demand pattern IDs for reservoirs|
|setBinNodeReservoirElevation|Sets the values of elevations for reservoirs|
|setBinNodeReservoirParameters|Sets the values of reservoir parameters (elevations, patterns)|
|setBinNodeSourceQuality|Sets the values of source qualities|
|setBinNodeTankDiameter|Sets the values of tanks diameter|
|setBinNodeTankElevation|Sets the values of tanks elevation|
|setBinNodeTankInitLevel|Sets the values of tanks initial level|
|setBinNodeTankMaxLevel|Sets the values of tanks maximum water level|
|setBinNodeTankMinLevel|Sets the values of tanks minimum water level|
|setBinNodeTankMinVol|Sets the values of tanks minimum water volume|
|setBinNodeTankParameters|Sets the values of reservoir parameters (elevations, initialLevels, minLevels, maxLevels, diameters, minVolume, mixfraction)|
|setBinPattern|Sets all of the multiplier factors for a specific time pattern|
|setBinQualityAge|Sets the type of water quality analysis to Age|
|setBinQualityChem|Sets the type of water quality analysis to Chem|
|setBinQualityNone|Sets the type of water quality analysis to None|
|setBinQualityTrace|Sets the type of water quality analysis to Trace|
|setBinTimeHydraulicStep|Sets the hydraulic step|
|setBinTimePatternStart|Sets the pattern start|
|setBinTimePatternStep|Sets the pattern step|
|setBinTimeQualityStep|Sets the quality step|
|setBinTimeReportingStart|Sets the reporting start|
|setBinTimeReportingStep|Sets the reporting step|
|setBinTimeSimulationDuration|Sets the simulation duration|
|setBinTimeStatisticsAverage|Sets the statistic type to Average|
|setBinTimeStatisticsMaximum|Sets the statistic type to Maximum|
|setBinTimeStatisticsMinimum|Sets the statistic type to Minimum|
|setBinTimeStatisticsNone|Sets the statistic type to None|
|setBinTimeStatisticsRange|Sets the statistic type to Range|


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
