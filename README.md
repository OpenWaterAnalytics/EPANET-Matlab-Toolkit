<a href="http://www.kios.ucy.ac.cy"><img src="http://www.kios.ucy.ac.cy/templates/favourite/images/kios_logo_hover.png"/><a>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.831493.svg)](https://doi.org/10.5281/zenodo.831493)


EPANET-MATLAB-Toolkit
==================================

The `EPANET-Matlab Toolkit` is an open-source software, originally developed by the [KIOS Research Center for Intelligent Systems and Networks of the University of Cyprus](http://www.kios.ucy.ac.cy/) which operates within the Matlab environment, for providing a programming interface for the latest version of [EPANET](https://github.com/OpenWaterAnalytics/epanet), a hydraulic and quality modeling software created by the US EPA, with Matlab, a  high-level technical computing software. The goal of the EPANET Matlab Toolkit is to serve as a common programming framework for research and development in the growing field of smart water networks.

The `EPANET-Matlab Toolkit` features easy to use commands/wrappers for viewing, modifying, simulating and plotting results produced by the EPANET libraries.  

For support, please use the OWA community forum : http://community.wateranalytics.org/

## Table of Contents

- [How to cite](#how-to-cite)
- [Requirements](#requirements)
- [How to install necessary compilers](#How-to-install-necessary-compilers)
- [How to use the Toolkit](#How-to-use-the-Toolkit)
- [How to fix/report bugs](#How-to-fixreport-bugs)
- [Licenses](#Licenses)
- [Contributors](#Contributors)
- [List of Matlab Class Functions](#List-of-Matlab-Class-Functions)
- [List of EPANET 2.1 Functions Supported](#List-of-EPANET-21-Functions-Supported)
- [List of EPANET 2.012 Functions Supported](#List-of-EPANET-2012-Functions-Supported)
- [List of EPANET 2.2 Functions Supported](#List-of-EPANET-22-Functions-Supported)

## How to cite 

D.G. Eliades, M. Kyriakou, S. Vrachimis and M.M. Polycarpou, "EPANET-MATLAB Toolkit: An Open-Source Software for Interfacing EPANET with MATLAB", in *Proc. 14th International Conference on Computing and Control for the Water Industry (CCWI)*, The Netherlands, Nov 2016, p.8. (doi:10.5281/zenodo.831493)

```
@INPROCEEDINGS{Eliades2016, 
author={Eliades, Demetrios G. and Kyriakou, Marios and Vrachimis, Stelios and Polycarpou, Marios M.}, 
title={EPANET-MATLAB Toolkit: An Open-Source Software for Interfacing EPANET with MATLAB}, 
booktitle={Proc. 14th International Conference on Computing and Control for the Water Industry (CCWI)}, 
year={2016},
pages={8},
address = {The Netherlands},
month={Nov},
DOI={10.5281/zenodo.831493}}
```

&uparrow; [Back to top](#table-of-contents)

## Requirements

* [Matlab](http://www.mathworks.com/)
* [EPANET 2.1](https://github.com/OpenWaterAnalytics/epanet) 

&uparrow; [Back to top](#table-of-contents)

## How to install necessary compilers

In order to use the EPANET-MATLAB-Toolkit the <b> MinGW-w64 compiler </b> must be install: <p>
<a href="http://www.youtube.com/watch?feature=player_embedded&v=R_RABL3_6EY
" target="_blank"><img src="http://img.youtube.com/vi/R_RABL3_6EY/0.jpg" 
alt="How to install MinGW-w64 compiler #Matlab" width="240" height="180" border="5" /></a>

In case you have a version of matlab older than 2015b install the <b> Windows SDK compiler </b>: <p>
<a href="http://www.youtube.com/watch?feature=player_embedded&v=hc3OkDypd24
" target="_blank"><img src="http://img.youtube.com/vi/hc3OkDypd24/0.jpg" 
alt="How to install MinGW-w64 compiler #Matlab" width="240" height="180" border="5" /></a>

&uparrow; [Back to top](#table-of-contents)

## How to use the Toolkit

Detailed examples on how to use the toolkit can be found in the [publication](https://zenodo.org/record/831493#.W9B69PZRXIV) of the Toolkit , together with the [code](https://zenodo.org/record/437778). A presentation of its use is also provided [here](https://github.com/KIOS-Research/CCWI2016/blob/master/CCWI2016/Presentation/Eliades_CCWI2016.ppt).

To start, you need to download the folder from GitHub (e.g., `Download ZIP`), set the run path in Matlab within the saved folder, and run `RunTests.m`. This will execute all the commands which have been implemented in the Class.

Minimum Example:

d=epanet('Net1.inp')

d.getNodeCount

d.getNodeElevations

Help Functions:

help d.plot

&uparrow; [Back to top](#table-of-contents)

## How to fix/report bugs

To fix a bug `Fork` the `EPANET-Matlab Toolkit`, `Edit` the code and make the appropriate change, and then `Pull` it so that we evaluate it. 

Keep in mind that some bugs may exist in the `EPANET` libraries, in case you are not receiving the expected results.

&uparrow; [Back to top](#table-of-contents)

## Licenses

* `EPANET`: Public Domain
* `EPANET-MSX`: GNU Lesser General Public License
* `EPANET-Matlab Toolkit`: EUPL 

&uparrow; [Back to top](#table-of-contents)

## Contributors

* Marios Kyriakou, [KIOS Research Center for Intelligent Systems and Networks, University of Cyprus](http://www.kios.ucy.ac.cy/)
* Demetrios Eliades, [KIOS Research Center for Intelligent Systems and Networks, University of Cyprus](http://www.kios.ucy.ac.cy/)
* Stelios Vrachimis, [KIOS Research Center for Intelligent Systems and Networks, University of Cyprus](http://www.kios.ucy.ac.cy/)

The `EPANET-Matlab Toolkit` is based/inspired on the [EPANET-Matlab Toolkit](http://www.mathworks.com/matlabcentral/fileexchange/25100-epanet-matlab-toolkit) as well as the OpenWaterAnalytics [EPANET-Matlab Wrappers](https://github.com/OpenWaterAnalytics/epanet-matlab)

&uparrow; [Back to top](#table-of-contents)

## List of Matlab Class Functions

|Function|Description|
|---------|---------|
|epanet| Load Input file and open the EPANET Toolkit system|
|unload|Unload library and close the EPANET Toolkit system|
|loadEPANETFile| Open the EPANET Toolkit system|
|getError| Returns the description of an error code| 
|getComputedHydraulicTimeSeries|Computed Hydraulic Time Series|
|getComputedQualityTimeSeries|Computed Quality Time Series|
|getComputedTimeSeries|Computed Hydraulic & Quality Time Series using the bimary file who created from executable|
|getComputedTimeSeries_ENepanet|Computed Hydraulic & Quality Time Series via ENepanet binary file|
|getConnectivityMatrix|Return connectivity matrix of the network|
|getCounts|Retrieves the number of network components|
|getControlRulesCount|Retrieves the number of control rules|
|getControls|Retrieves the controls|
|getCurveCount|Retrieves the number of curves|
|getCurveIndex|Retrieves index of curve with specific ID|
|getCurveLengths|Retrieves number of points in a curve|
|getCurveNameID|Retrieves curve id|
|getCurveType|Retrieves the curve-type (VOLUME, PUMP, EFFICIENCY, HEADLOSS, GENERAL)|
|getCurveTypeIndex|Retrieves the curve-type index for all curves|
|getCurveValue|Retrieves (x,y) values of specific curve index|
|getCurveXY|Retrieves (x,y) values of all curves|
|getDemandModel|Retrieves the type of demand model in use and its parameters|
|getENfunctionsImpemented|Retrieves the epanet functions that have been developed|
|getFlowUnits|Retrieves the units used to express all flow rates|
|getLibFunctions|Retrieves the functions of DLL|
|getLinkActualQuality|Current computed link quality (read only)|
|getLinkBulkReactionCoeff|Bulk chemical reaction coefficient|
|getLinkComment|Retrieves the comment string assigned to the link object|
|getLinkCount|Retrieves the number of links|
|getLinkDiameter|Retrieves the value of all link diameters|
|getLinkEnergy|Current computed pump energy usage (read only)|
|getLinkFlows|Current computed flow rate (read only)|
|getLinkHeadloss|Current computed head loss (read only)|
|getLinkIndex|Retrieves the indices of all links, or the indices of an ID set of links|
|getLinkInitialSetting|Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves|
|getLinkInitialStatus|Retrieves the value of all link initial status|
|getLinkLength|Retrieves the value of all link lengths|
|getLinkMinorLossCoeff|Retrieves the value of all link minor loss coefficients|
|getLinkNameID|Retrieves the ID label(s) of all links, or the IDs of an index set of links|
|getLinkNodesIndex-getNodesConnectingLinksIndex|Retrieves the indexes of the from/to nodes of all links|
|getLinkPipeCount|Retrieves the number of pipes|
|getLinkPipeIndex|Retrieves the indices of pipes|
|getLinkPipeNameID|Retrieves the pipe IDs|
|getLinkPumpCount|Retrieves the number of pumps|
|getLinkPumpEfficiency|Retrieves the value of all computed efficiency|
|getLinkPumpHeadCurveIndex|Retrieves index of a head curve for specific link index|
|getLinkPumpIndex|Retrieves the indices of pumps|
|getLinkPumpNameID|Retrieves the pump IDs|
|getLinkPumpPatternIndex|Pump speed time pattern index|
|getLinkPumpPatternNameID|Retrieves the pump pattern IDs|
|getLinkPumpPower|Pump constant power rating|
|getLinkPumpHCurve|Pump head v. flow curve index|
|getLinkPumpECurve|Pump efficiency v. flow curve index|
|getLinkPumpECost|Pump average energy price|
|getLinkPumpEPat|Pump energy price time pattern index|
|getLinkPumpType|Retrieves the type of a pump for specific link index|
|getLinkPumpTypeCode|Retrieves the type code of a pump for specific link index|
|getLinkPumpState|Current computed pump state (read only) (see @ref EN_PumpStateType)|
|getLinkPumpSwitches|Calculates the number of pump switches|
|getLinkRoughnessCoeff|Retrieves the value of all link roughness|
|getLinkSettings|Retrieves the value of all computed link roughness for pipes or actual speed for pumps or actual setting for valves|
|getLinkStatus|Current link status (see @ref EN_LinkStatusType)|
|getLinkType|Retrieves the link-type for all links|
|getLinkTypeIndex|Retrieves the link-type code for all links|
|getLinkValveCount|Retrieves the number of valves|
|getLinkValveIndex|Retrieves the indices of valves|
|getLinkValveNameID|Retrieves the valve IDs|
|getLinkVelocity|Current computed flow velocity (read only)|
|getLinkVertices|Retrieves the coordinate's of a vertex point assigned to a link.|
|getLinkVerticesCount|Retrieves the number of internal vertex points assigned to a link|
|getLinkWallReactionCoeff|Pipe wall chemical reaction coefficient|
|getNodeActualDemand|Retrieves the computed value of all actual demands|
|getNodeActualDemandSensingNodes|Retrieves the computed demand values at some sensing nodes|
|getNodeActualQuality|Retrieves the computed values of the actual quality for all nodes|
|getNodeActualQualitySensingNodes|Retrieves the computed quality values at some sensing nodes|
|getNodeBaseDemands|Retrieves the value of all node base demands|
|getNodeComment|Retrieves the comment string assigned to the node object|
|getNodeCoordinates|Retrieves coordinate x, y, and x, y vertices for a node|
|getNodeCount|Retrieves the number of nodes|
|getNodePatternIndex|Retrieves the value of all node pattern indices|
|getNodeDemandDeficit|Retrieves the amount that full demand is reduced under PDA. (EPANET Version 2.2)|
|getNodeDemandPatternIndex|Retrieves the value of all node demand pattern indices|
|getNodeDemandPatternNameID|Retrieves the value of all node demand pattern IDs|
|getNodeElevations|Retrieves the value of all node elevations|
|getNodeEmitterCoeff|Retrieves the value of all node emmitter coefficients|
|getNodeHydaulicHead|Retrieves the computed values of all hydraulic heads|
|getNodeIndex|Retrieves the indices of all nodes or some nodes with a specified ID|
|getNodeInitialQuality|Retrieves the value of all node initial quality|
|getNodeJunctionCount|Retrieves the number of junctions|
|getNodeJunctionDemandName|Gets the name of a node's demand category|
|getNodeJunctionIndex|Retrieves the junctions indices|
|getNodeJunctionNameID|Retrieves the junctions IDs|
|getNodeJunctionDemandIndex|Retrieves the demand index of the junctions. (EPANET Version 2.2)|
|getNodeMassFlowRate|Retrieves the computed mass flow rates per minute of chemical sources|
|getNodeNameID|Retrieves the ID label of all nodes or some nodes with a specified index|
|getNodeDemandCategoriesNumber|Retrieves the number of demand categories for a node|
|getNodePressure|Retrieves the computed values of all node pressures|
|getNodeReservoirCount|Retrieves the number of reservoirs|
|getNodeReservoirIndex|Retrieves the indices of reservoirs|
|getNodeReservoirNameID|Retrieves the reservoirs IDs|
|getNodeSourcePatternIndex|Retrieves the value of all node source pattern index|
|getNodeSourceQuality|Retrieves the value of all nodes source quality|
|getNodeSourceType|Retrieves the value of all node source type|
|getNodeTankData|Retrieves a group of properties for a tank. (EPANET Version 2.2)|
|getNodeTankBulkReactionCoeff|Retrieves the tank bulk rate coefficient|
|getNodeTankCanOverFlow|Retrieves the tank can overflow (= 1) or not (= 0)|
|getNodeTankCount|Retrieves the number of tanks|
|getNodeTankDiameter|Retrieves the tank diameters|
|getNodeTankIndex|Retrieves the indices of tanks|
|getNodeTankInitialLevel|Retrieves the value of all tank initial water levels|
|getNodeTankInitialWaterVolume|Retrieves the tank initial volume|
|getNodeTankMaximumWaterVolume|Retrieves maximum water volume|
|getNodeTankMaximumWaterLevel|Retrieves the tank maximum water level|
|getNodeTankMixingFraction|Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank|
|getNodeTankMinimumWaterLevel|Retrieves the tank minimum water level|
|getNodeTankMinimumWaterVolume|Retrieves the tank minimum volume|
|getNodeTankMixZoneVolume|Retrieves the tank mixing zone volume|
|getNodeTankMixingModelCode|Retrieves the tank mixing model code|
|getNodeTankMixingModelType|Retrieves the tank mixing model type (mix1, mix2, fifo, lifo)|
|getNodeTankNameID|Retrieves the tanks IDs|
|getNodeTankReservoirCount|Retrieves the number of tanks|
|getNodeTankVolume|Retrieves the tank volume|
|getNodeTankVolumeCurveIndex|Retrieves the tank volume curve index|
|getNodeType|Retrieves the node-type for all nodes|
|getNodeTypeIndex|Retrieves the node code-index for all nodes|
|getNodesConnectingLinksID|Retrieves the id of the from/to nodes of all links|
|getOptionsAccuracyValue|Retrieve the analysis convergence criterion (0.001)|
|getOptionsDemandCharge|Retrieve energy price pattern|
|getOptionsSpecificGravity|Retrieves the specific gravity (EPANET Version 2.2)|
|getOptionsSpecificViscosity|Retrieves the specific viscosity (EPANET Version 2.2)|
|getOptionsExtraTrials|Retrieves the extra trials allowed if hydraulics don't converge (EPANET Version 2.2)|
|getOptionsCheckFrequency|Retrieves the frequency of hydraulic status checks (EPANET Version 2.2)|
|getOptionsMaximumCheck|Retrieves the maximum trials for status checking. (EPANET Version 2.2)|
|getOptionsEmitterExponent|Retrieve power exponent for the emmitters (0.5)|
|getOptionsFlowChange|Retrieve flow change|
|getOptionsGlobalEffic|Retrieve global efficiency pumps|
|getOptionsGlobalPrice|Retrieve global average energy price per kW-Hour|
|getOptionsGlobalPattern|Retrieve global pattern|
|getOptionsHeadError|Retrieve the head error|
|getOptionsHeadLossFormula|Retrieve headloss formula code (Hazen-Williams, Darcy-Weisbach or Chezy-Manning)|
|getOptionsMaxTrials|Retrieve maximum number of analysis trials|
|getOptionsPatternDemandMultiplier|Retrieve the demand multiplier (x1)|
|getOptionsQualityTolerance|Retrieve the water quality analysis tolerance|
|getOptionsDampLimit|Retrieves the accuracy level where solution damping begins. (EPANET Version 2.2)|
|getOptionsSpecificDiffusivity|Retrieves the specific diffusivity (relative to chlorine at 20 deg C). (EPANET Version 2.2)|
|getOptionsPipeBulkReactionOrder|Retrieves the bulk water reaction order for pipes. (EPANET Version 2.2)|
|getOptionsPipeWallReactionOrder|Retrieves the wall reaction order for pipes (either 0 or 1). (EPANET Version 2.2)|
|getOptionsTankBulkReactionOrder|Retrieves the bulk water reaction order for tanks. (EPANET Version 2.2)|
|getOptionsLimitingConcentration|Retrieves the limiting concentration for growth reactions. (EPANET Version 2.2)|
|getPattern|Retrieves the multiplier factor for all patterns and all times|
|getPatternAveragePatternValue|Retrieves the average value of a pattern|
|getPatternComment|Retrieves the comment string assigned to the pattern object|
|getPatternCount|Retrieves the number of patterns|
|getPatternIndex|Retrieves the index of all or some time patterns IDs|
|getPatternLengths|Retrieves the number of time periods in all or some patterns|
|getPatternNameID|Retrieves the patterns IDs|
|getPatternValue|Retrieves the multiplier factor for a certain pattern and time|
|getQualityCode|Retrieves the code of water quality analysis type|
|getQualityInfo|Retrieves the quality info - bug in ENgetqualinfo|
|getQualityTraceNodeIndex|Retrieves the trace node index of water quality analysis type|
|getQualityType|Retrieves the type of water quality analysis type|
|getRules|Retrieves the rule - based control statements. (EPANET Version 2.2)|
|getRuleCount|Retrieves the number of rules. (EPANET Version 2.2)|
|getRuleID|Retrieves the ID name of a rule-based control given its index. (EPANET Version 2.2)|
|getRuleInfo|Retrieves summary information about a rule-based control given it's index. (EPANET Version 2.2)|
|getStatistic|Retrieves hydraulic simulation statistic|
|getTimeHTime|Retrieves the number of htime|
|getTimeQTime|Retrieves the number of qtime|
|getTimeHaltFlag|Retrieves the number of  halt flag|
|getTimeHydraulicStep|Retrieves the value of the hydraulic time step|
|getTimeNextEvent|Retrieves the number of next event|
|getTimeNextEventTank|Retrieves the index of tank with shortest time to become empty or full|
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
|getTitle|Retrieves the title lines of the project|
|getUnits|Retrieves the Units of Measurement|
|getVersion|Retrieve the current EPANET version of DLL|
|getNodesInfo|Retrieves nodes info e.g. elevations, demand pattern indices, emitter coeff. , initial quality, source quality, source pattern indices, source type code, type indices|
|getLinksInfo|Retrieves links info e.g. diameters, lengths, roughness coeff. , minor loss coeff. , initial status, initial settings, bulk reaction coeff. , wall reaction coeff. , nodes connecting link indices, type indices|
|addControls|Adds a new simple control. (EPANET Version 2.2)|
|addCurve|Adds a new curve appended to the end of the existing curves|
|addPattern|Adds a new time pattern to the network|
|addNodeJunction|Adds a new junction|
|addNodeJunctionDemand|Adds a new demand to a junction given the junction index, base demand, demand time pattern and demand name category. (EPANET Version 2.2)|
|addNodeReservoir|Adds a new reservoir|
|addNodeTank|Adds a new tank|
|addLinkPipeCV|Adds a new CV pipe|
|addLinkPipe|Adds a new pipe|
|addLinkPump|Adds a new pump|
|addLinkValvePRV|Adds a new PRV valve|
|addLinkValvePSV|Adds a new PSV valve|
|addLinkValvePBV|Adds a new PBV valve|
|addLinkValveFCV|Adds a new FCV valve|
|addLinkValveTCV|Adds a new TCV valve|
|addLinkValveGPV|Adds a new GPV valve|
|addRules|Adds a new rule-based control to a project. (EPANET Version 2.2)|
|deleteCurve|Deletes a data curve from the project|
|deleteLink|Deletes a link|
|deleteNode|Deletes a node|
|deletePattern|Deletes a time pattern from a project|
|deleteControls|Deletes an existing simple control. (EPANET Version 2.2)|
|deleteRules|Deletes an existing rule-based control given it's index. (EPANET Version 2.2)|
|clearReport|Clears the contents of a project's report file. (EPANET Version 2.2)|
|copyReport|Copies the current contents of a project's report file to another file. (EPANET Version 2.2)|
|closeHydraulicAnalysis|Closes the hydraulic analysis system, freeing all allocated memory|
|closeNetwork|Closes down the Toolkit system|
|closeQualityAnalysis|Closes the water quality analysis system, freeing all allocated memory|
|runsCompleteSimulation|Runs a complete hydraulic and water simulation to create binary & report files with default name net_temp.bin or you can use argument to run via ENepanet|
|initializeHydraulicAnalysis|Initializes storage tank levels, link status and settings, and the simulation clock time prior to running a hydraulic analysis|
|initializeQualityAnalysis|Initializes water quality and the simulation clock time prior to running a water quality analysis|
|nextHydraulicAnalysisStep|Determines the length of time until the next hydraulic event occurs in an extended period simulation|
|nextQualityAnalysisStep|Advances the water quality simulation to the start of the next hydraulic time period|
|openAnyInp|Open as on matlab editor any EPANET input file|
|openCurrentInp|Open EPANET input file who is loaded|
|openHydraulicAnalysis|Opens the hydraulics analysis system|
|openQualityAnalysis|Opens the water quality analysis system|
|runHydraulicAnalysis|Runs a single period hydraulic analysis, retrieving the current simulation clock time t|
|runQualityAnalysis|Makes available the hydraulic and water quality results that occur at the start of the next time period of a water quality analysis, where the start of the period is returned in t|
|saveHydraulicFile|Saves the current contents of the binary hydraulics file to a file|
|saveHydraulicsOutputReportingFile|Transfers results of a hydraulic simulation from the binary Hydraulics file to the binary Output file, where results are only reported at uniform reporting intervals|
|saveInputFile|Writes all current network input data to a file using the format of an EPANET input file|
|plot|Plot the network input file|
|setControls|Sets the parameters of a simple control statement|
|setCurve|Sets x,y values for a specific curve|
|setCurveNameID|Sets the name ID of a curve given it's index and the new ID. (EPANET Version 2.2)|
|setCurveValue|Retrieves x,y point for a specific point number and curve|
|setDemandModel|Sets the type of demand model to use and its parameters|
|setFlowUnitsAFD|Sets flow units to AFD|
|setFlowUnitsCFS|Sets flow units to CFS|
|setFlowUnitsCMD|Sets flow units to CMD|
|setFlowUnitsCMH|Sets flow units to CMH|
|setFlowUnitsGPM|Sets flow units to GPM|
|setFlowUnitsIMGD|Sets flow units to IMGD|
|setFlowUnitsLPM|Sets flow units to LPM|
|setFlowUnitsLPS|Sets flow units to LPS|
|setFlowUnitsMGD|Sets flow units to MGD|
|setFlowUnitsMLD|Sets flow units to MLD|
|setLinkBulkReactionCoeff|Sets the values of bulk reactions|
|setLinkComment|Sets the comment string assigned to the link object|
|setLinkDiameter|Sets the values of diameters|
|setLinkPipeData|Sets a group of properties for a pipe. (EPANET Version 2.2)|
|setLinkPumpHeadCurveIndex|Sets the curves index for pumps index|
|setLinkPumpPatternIndex|Sets the pump speed time pattern index. (EPANET Version 2.2)|
|setLinkPumpPower|Sets the power for pumps. (EPANET Version 2.2)|
|setLinkPumpHCurve|Sets the pump head v. flow curve index. (EPANET Version 2.2)|
|setLinkPumpECurve|Sets the pump efficiency v. flow curve index. (EPANET Version 2.2)|
|setLinkPumpECost|Sets the pump average energy price. (EPANET Version 2.2)|
|setLinkPumpEPat|Sets the pump energy price time pattern index. (EPANET Version 2.2)|
|setLinkInitialSetting|Sets the values of initial settings|
|setLinkInitialStatus|Sets the values of initial status|
|setLinkLength|Sets the values of lengths|
|setLinkMinorLossCoeff|Sets the values of minor loss coeff.|
|setLinkNameID|Sets the ID name for links|
|setLinkNodesIndex|Sets the indexes of a link's start- and end-nodes. (EPANET Version 2.2)|
|setLinkRoughnessCoeff|Sets the values of roughness coeff.|
|setLinkSettings|Sets the values of settings|
|setLinkStatus|Sets the values of status|
|setLinkTypePipe|Set the link type pipe for a specified link|
|setLinkTypePipeCV|Set the link type cvpipe for a specified link|
|setLinkTypePump|Set the link type pump for a specified link|
|setLinkTypeValveFCV|Set the link type valve FCV for a specified link|
|setLinkTypeValveGPV|Set the link type valve PCV for a specified link|
|setLinkTypeValvePBV|Set the link type valve PBV for a specified link|
|setLinkTypeValvePRV|Set the link type valve PRV for a specified link|
|setLinkTypeValvePSV|Set the link type valve PSV for a specified link|
|setLinkTypeValveTCV|Set the link type valve TCV for a specified link|
|setLinkVertices|Assigns a set of internal vertex points to a link|
|setLinkWallReactionCoeff|Sets the values of wall reactions|
|setNodeBaseDemands|Sets the values of demands|
|setNodeComment|Sets the comment string assigned to the node object|
|setNodeCoordinates|Sets node coordinates|
|setNodeDemandPatternIndex|Sets the values of demand pattern indices|
|setNodeElevations|Sets the values of elevations|
|setNodeEmitterCoeff|Sets the values of emitter coeff.|
|setNodeInitialQuality|Sets the values of initial qualities|
|setNodeJunctionDemandName|Assigns a name to a node's demand category|
|setNodeNameID|Sets the ID name for nodes|
|setNodeSourcePatternIndex|Sets the values of source pattern indices|
|setNodeSourceQuality|Sets the values of source qualities|
|setNodeSourceType|Sets the values of source types|
|setNodeTankData| Sets a group of properties for a tank. (EPANET Version 2.2)|
|setNodeTankBulkReactionCoeff|Sets the values of tank bulk reaction coeff.|
|setNodeTankDiameter|Sets the values of tanks diameter|
|setNodeTankCanOverFlow|Sets the value of tank can overflow (= 1) or not (= 0)|
|setNodeTankInitialLevel|Sets the values of tanks initial level|
|setNodeTankMaximumWaterLevel|Sets the values of tanks maximum water level|
|setNodeTankMinimumWaterLevel|Sets the values of tanks minimum water level|
|setNodeTankMixingFraction|Sets the values of tanks mix fraction|
|setNodeTankMinimumWaterVolume|Sets the values of tanks minimum water volume|
|setNodeTankMixingModelType|Sets the values of tanks model|
|setOptionsAccuracyValue|Sets the value of accurancy|
|setOptionsGlobalEffic|Sets the value of global pump efficiency(percent) (EPANET Version 2.2)|
|setOptionsGlobalPrice|Sets the value of global energy price per KWH (EPANET Version 2.2)|
|setOptionsGlobalPattern|Sets the index of a global energy price pattern (EPANET Version 2.2)|
|setOptionsDemandCharge|Sets the energy demand charge per max. KW usage (EPANET Version 2.2)|
|setOptionsSpecificGravity|Sets the specific gravity (EPANET Version 2.2)|
|setOptionsSpecificViscosity|Sets the specific viscosity (EPANET Version 2.2)|
|setOptionsExtraTrials|Sets the extra trials allowed if hydraulics don't converge (EPANET Version 2.2)|
|setOptionsMaximumCheck|Sets the maximum trials for status checking. (EPANET Version 2.2)|
|setOptionsCheckFrequency|Sets the frequency of hydraulic status checks (EPANET Version 2.2)|
|setOptionsDampLimit|Sets the accuracy level where solution damping begins. (EPANET Version 2.2)|
|setOptionsSpecificDiffusivity|Sets the specific diffusivity (relative to chlorine at 20 deg C). (EPANET Version 2.2)|
|setOptionsPipeBulkReactionOrder|Sets the bulk water reaction order for pipes. (EPANET Version 2.2)|
|setOptionsPipeWallReactionOrder|Sets the wall reaction order for pipes (either 0 or 1). (EPANET Version 2.2)|
|setOptionsTankBulkReactionOrder|Sets the bulk water reaction order for tanks. (EPANET Version 2.2)|
|setOptionsLimitingConcentration|Sets the limiting concentration for growth reactions. (EPANET Version 2.2)|
|setOptionsEmitterExponent|Sets the value of emitter exponent|
|setOptionsMaxTrials|Sets the value of max trials|
|setOptionsPatternDemandMultiplier|Sets the value of pattern demand multiplier|
|setOptionsQualityTolerance|Sets the value of tolerance|
|setPattern|Sets all of the multiplier factors for a specific time pattern|
|setPatternComment|Sets the comment string assigned to the pattern object|
|setPatternNameID|Sets the name ID of a time pattern given it's index and the new ID. (EPANET Version 2.2)|
|setPatternMatrix|Sets all of the multiplier factors for all patterns|
|setPatternValue|Sets the multiplier factor for a specific period within a time pattern|
|setQualityType|Sets the type of water quality analysis called for|
|setReport|Issues a report formatting command. Formatting commands are the same as used in the [REPORT] section of the EPANET Input file|
|setReportFormatReset|Clears any report formatting commands that either appeared in the [REPORT] section of the EPANET Input file or were issued with the ENsetreport function|
|setReportStatus|Sets the level of hydraulic status reporting|
|setRules|Sets a rule - based control. (EPANET Version 2.2)|
|setRuleElseAction|Sets rule - based control else actions. (EPANET Version 2.2)|
|setRulePremise|Sets the premise of a rule - based control. (EPANET Version 2.2)|
|setRulePremiseObejctNameID|Sets the ID of an object in a premise of a rule-based control. (EPANET Version 2.2)|
|setRulePremiseStatus|Sets the status being compared to in a premise of a rule-based control. (EPANET Version 2.2)|
|setRulePremiseValue|Sets the value being compared to in a premise of a rule-based control. (EPANET Version 2.2)|
|setRulePriority|Sets rule - based control priority. (EPANET Version 2.2)|
|setRuleThenAction|Sets rule - based control then actions. (EPANET Version 2.2)|
|setTimeHydraulicStep|Sets the hydraulic step|
|setTimePatternStart|Sets the pattern start|
|setTimePatternStep|Sets the pattern step|
|setTimeQualityStep|Sets the quality step|
|setTimeReportingStart|Sets the reporting start|
|setTimeReportingStep|Sets the reporting step|
|setTimeRuleControlStep|Sets the rule control step|
|setTimeSimulationDuration|Sets the simulation duration|
|setTimeStatisticsType|Sets the statistic type|
|setTitle|Sets the title lines of the project|
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
|writeMSXReportExe|Writes water quality simulations results as instructed by the MSX input file to a specific name text file|
|useMSXHydraulicFile|Uses a previously saved EPANET hydraulics file as the source of hydraulic information|
|plotMSXConcentrationSpeciesOfLinks|Plots the concentration species of links|
|plotMSXConcentrationSpeciesOfNodes|Plots the concentration species of nodes|
|runMSXexe|Writes water quality simulations results as instructed by the MSX input file to a text file using the epanetmsx.exe|
|unloadMSX|Closes the EPANET-MSX toolkit system|
|getMSXAtol|Retrieves the absolute concentration tolerance|
|getMSXRtol|Retrieves the relative concentration tolerance|
|getMSXComputedQualitySpecie|Retrieves the quality values for specific specie (e.g getMSXComputedQualitySpecie('CL2'))|
|getMSXComputedQualityLink|Retrieves the concentration of a chemical species at a specific link of the network at the current simulation time step|
|getMSXComputedQualityNode|Retrieves the concentration of a chemical species at a specific node of the network at the current simulation time step.|
|getMSXConstantsCount|Retrieves the number of constants|
|getMSXConstantsIndex|Retrieves the internal index number of constants (given its ID name)|
|getMSXConstantsNameID|Retrieves the ID name of constants (given its internal index number)|
|getMSXConstantsValue|Retrieves the value of a particular reaction constant|
|getMSXError|Returns the text for an error message given its error code|
|getMSXLinkInitqualValue|Retrieves the initial concentration of chemical species assigned to links of the pipe network|
|getMSXNodeInitqualValue|Retrieves the initial concentration of chemical species assigned to nodes|
|getMSXOptions|Retrieves all the msx option parameters|
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
|getMSXSourceNodeNameID|Retrieves the ID label of all nodes|
|getMSXSourcePatternIndex|Retrieves the value of all node source pattern index|
|getMSXSourceType|Retrieves the value of all node source type|
|getMSXSources|Retrieves the source info|
|getMSXSpeciesATOL|Retrieves the atol|
|getMSXSpeciesRTOL|Retrieves the rtol|
|getMSXSpeciesConcentration|Retrieves the concentration of chemical species for nodes and links|
|getMSXSpeciesCount|Retrieves the number of species|
|getMSXSpeciesIndex|Retrieves the indices of species|
|getMSXSpeciesNameID|Retrieves the species IDs|
|getMSXSpeciesType|Retrieves the type of all species (BULK/WALL)|
|getMSXSpeciesUnits|Retrieves the species mass units|
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
|addBinNodeJunction|Adds a new junction to the network|
|addBinPattern|Adds a new time pattern to the network|
|addBinPipe|Adds a new pipe to the network|
|addBinLinkPipe|Adds a new pipe to the network|
|addBinPump|Adds a new pump to the network|
|addBinLinkPump|Adds a new pump to the network|
|addBinReservoir|Adds a new reservoir to the network|
|addBinNodeReservoir|Adds a new reservoir to the network|
|addBinTank|Adds a new tank to the network|
|addBinNodeTank|Adds a new tank to the network|
|addBinValveFCV|Adds a new valve FCV to the network|
|addBinValveGPV|Adds a new valve GPV to the network|
|addBinValvePBV|Adds a new valve PBV to the network|
|addBinValvePRV|Adds a new valve PRV to the network|
|addBinValvePSV|Adds a new valve PSV to the network|
|addBinValveTCV|Adds a new valve TCV to the network|
|addBinLinkValve|Adds a new valve to the network|
|addBinLinkVertices|Adds interior vertex points to network links|
|removeBinControlLinkID|Removes a specific control based on link ID|
|removeBinControlNodeID|Removes a specific control based on node ID|
|removeBinCurveID|Removes a specific curve based on ID|
|removeBinLinkID|Removes a specific link based on ID|
|removeBinNodeID|Removes a specific node based on ID|
|removeBinRulesControlLinkID|Removes a specific rule based on link ID|
|removeBinRulesControlNodeID|Removes a specific rule based on node ID|
|deleteBinLinkVertices|Deletes interior vertex points of network links|
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
|getBinLimitingPotential|Retrieves limiting potential|
|getBinNodesInfo|Retrieves the nodes info|
|getBinNodeSourceInfo|Retrieves the sources info|
|getBinOptionsInfo|Retrieves the options info|
|getBinPatternsInfo|Retrieves the patterns info|
|getBinRulesControlsInfo|Retrieves the controls info|
|getBinTimesInfo|Retrieves the times info|
|getBinPatternIndex|Retrieves the indices of all patterns|
|getBinSimulationDuration|Retrieves the value of simulation duration|
|getBinSections|Retrieves some basic sections from inp file|
|getBinUnits|Retrieves the units used to express all flow rates|
|getBinLinkVertices|Retrieves the link vertices|
|getBinLinkVerticesCount|Retrieves the number of vertices|
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
|setBinLimitingPotential|Sets limiting potential|
|setBinLinkPipeDiameters|Sets the values of pipe diameters|
|setBinLinkPipeLengths|Sets the values of pipe lengths|
|setBinLinkPipeMinorLoss|Sets the values of pipe minor losses|
|setBinLinkPipeRoughness|Sets the values of pipe roughness|
|setBinLinkPipeStatus|Sets the values of pipe status|
|setBinLinkPipesParameters|Sets the values of pipe parameters (diameters, lengths, minor losses, roughness, status)|
|setBinLinkPumpStatus|Sets the values of pump status|
|setBinLinkReactionCoeff|Sets the values of bulk and wall reaction coeff.|
|setBinLinkValvesParameters|Sets the values of valve parameters (diameters, types, settings, minor losses)|
|setBinNodeJunDemandPatternNameID|Sets the names of demand pattern IDs for junctions|
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
|setBinNodeTankInitialLevel|Sets the values of tanks initial level|
|setBinNodeTankMaximumWaterLevel|Sets the values of tanks maximum water level|
|setBinNodeTankMinimumWaterLevel|Sets the values of tanks minimum water level|
|setBinNodeTankMinimumWaterVolume|Sets the values of tanks minimum water volume|
|setBinNodeTankParameters|Sets the values of reservoir parameters (elevations, initialLevels, minLevels, maxLevels, diameters, minVolume, mixfraction)|
|setBinPattern|Sets all of the multiplier factors for a specific time pattern|
|setBinQualityAge|Sets the type of water quality analysis to Age|
|setBinQualityChem|Sets the type of water quality analysis to Chem|
|setBinQualityNone|Sets the type of water quality analysis to None|
|setBinQualityTrace|Sets the type of water quality analysis to Trace|
|setBinQualType|Sets the type of water quality analysis to any chem e.g. chlorine|
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
|setBinLinkVertices|Sets interior vertex points of network links|

&uparrow; [Back to top](#table-of-contents)

## List of EPANET 2.012 Functions Supported

|Function|Description|
|--------|-----------|
|ENaddpattern|Adds a new time pattern to the network|
|ENclose|Closes down the Toolkit system (including all files being processed)|
|ENcloseH|Closes the hydraulic analysis system, freeing all allocated memory|
|ENcloseQ|Closes the water quality analysis system, freeing all allocated memory|
|ENepanet|Runs a complete EPANET simulation|
|ENgetcount|Retrieves the number of network components of a specified type|
|ENgetcontrol|Retrieves the parameters of a simple control statement|
|ENgeterror|Retrieves the text of the message associated with a particular error or warning code|
|ENgetflowunits|Retrieves a code number indicating the units used to express all flow rates|
|ENgetlinkid|Retrieves the ID label of a link with a specified index|
|ENgetlinkindex|Retrieves the index of a link with a specified ID|
|ENgetlinknodes|Retrieves the indexes of the end nodes of a specified link|
|ENgetlinktype|Retrieves the link-type code for a specific link|
|ENgetlinkvalue|Retrieves the value of a specific link parameter|
|ENgetnodeid|Retrieves the ID label of a node with a specified index|
|ENgetnodeindex|Retrieves the index of a node with a specified ID|
|ENgetnodetype|Retrieves the node-type code for a specific node|
|ENgetnodevalue|Retrieves the value of a specific link parameter|
|ENgetoption|Retrieves the value of a particular analysis option|
|ENgetpatternid|Retrieves the ID label of a particular time pattern|
|ENgetpatternindex|Retrieves the index of a particular time pattern|
|ENgetpatternlen|Retrieves the number of time periods in a specific time pattern|
|ENgetpatternvalue|Retrieves the multiplier factor for a specific time period in a time pattern|
|ENgetqualtype|Retrieves the type of water quality analysis called for|
|ENgettimeparam|Retrieves the value of a specific analysis time parameter|
|ENgetversion|Retrieves the version number|
|ENinitH|Initializes hydraulic analysis|
|ENinitQ|Initializes water quality analysis|
|ENnextH|Determine time (in seconds) until next hydraulic event|
|ENnextQ|Advances WQ simulation to next hydraulic event|
|ENopen|Opens EPANET input file & reads in network data|
|ENopenH|Sets up data structures for hydraulic analysis|
|ENopenQ|Sets up data structures for WQ analysis|
|ENreport|Writes simulation report to the report file|
|ENresetreport|Resets report options to default values|
|ENrunH|Run a hydraulic solution period|
|ENrunQ|Retrieves hydraulic & WQ results at time t|
|ENsaveH|Saves hydraulic results to binary file|
|ENsavehydfile|Copies binary hydraulics file to disk|
|ENsaveinpfile|Saves current data to "INP" formatted text file|
|ENsetcontrol|Specify parameters to define a simple control|
|ENsetlinkvalue|Set a proprty value for a link|
|ENsetnodevalue|Set a property value for a node|
|ENsetoption|Set a value for an anlysis option|
|ENsetpattern|Set multipliers for a specific pattern|
|ENsetpatternvalue|Set the multiplier for a specific pattern at a specific period|
|ENsetqualtype|Sets the type of water quality analysis called|
|ENsetreport|Processes a reporting format command|
|ENsetstatusreport|Sets the level of hydraulic status reporting|
|ENsettimeparam|Set the value for a time parameter|
|ENsolveH|Solves the network hydraulics for all time periods|
|ENsolveQ|Solves for network water quality in all time periods|
|ENstepQ|Advances WQ simulation by a single WQ time step|
|ENusehydfile|Opens previously saved binary hydraulics file|
|ENwriteline|Writes line of text to the report file|

&uparrow; [Back to top](#table-of-contents)

## List of EPANET 2.1 Functions Supported

|Function|Description|
|--------|-----------|
|ENaddcurve|Adds a new curve appended to the end of the existing curves|
|ENgetaveragepatternvalue|Retrieves the average value of a pattern|
|ENgetbasedemand|Retrieves the nodes base demand for a category|
|ENgetcoord|Retrieves coordinate x, y for a node|
|ENgetcurve|Retrieves a curve's properties|
|ENgetcurveid|Retrieves ID of a curve with specific index|
|ENgetcurveindex|Retrieves index of curve with specific ID|
|ENgetcurvelen|Retrieves number of points in a curve|
|ENgetcurvevalue|Retrieves x,y point for a specific point number and curve|
|ENgetdemandpattern|Retrieves the index of a demand pattern for a specific demand category of a node|
|ENgetheadcurveindex|Retrieves index of a head curve for specific link index|
|ENgetnumdemands|Retrieves the number of demand categories for a node|
|ENgetpumptype|Retrieves the type of a pump for specific link index|
|ENgetqualinfo|Retrieves quality analysis information (type, chemical name, units, trace node ID)|
|ENgetstatistic|Retrieves hydraulic simulation statistic|
|ENsetbasedemand|Sets the nodes base demand for a category|
|ENsetcoord|Sets coordinate x, y for a node|
|ENsetcurve|Sets x,y values for a specific curve|
|ENsetcurvevalue|Sets x,y point for a specific point and curve|

&uparrow; [Back to top](#table-of-contents)

## List of EPANET 2.2 Functions Supported

|Function|Description|
|--------|-----------|
|ENaddcontrol|Specify parameters to add a new simple control|
|ENaddlink|Adds a new link|
|ENaddnode|Adds a new node|
|ENaddrule|Adds a new rule-based control to a project|
|ENadddemand|Appends a new demand to a junction node demands list|
|ENclearreport|Clears the contents of a project's report file|
|ENcopyreport|Copies the current contents of a project's report file to another file|
|ENdeletelink|Deletes a link|
|ENdeletenode|Deletes a node|
|ENsetcurveid|Changes the ID name of a data curve given its index|
|ENsetpatternid|Changes the ID name of a time pattern given its index|
|ENsetdemandpattern|Sets the index of the demand pattern assigned to a node for a category index|
|ENsetheadcurveindex|Sets the curve index for a specified pump index|
|ENgetcurvetype|Retrieves the type of a curve|
|ENgetdemandindex|Retrieves the index of a node's named demand category|
|ENgetpremise|Gets the properties of a premise in a rule-based control|
|ENgetelseaction|Gets the properties of an ELSE action in a rule-based control|
|ENgetruleid|Gets the ID name of a rule-based control given its index|
|ENgetrule|Retrieves summary information about a rule-based control|
|ENgetthenaction|Gets the properties of a THEN action in a rule-based control|
|ENsetflowunits|Sets the flow units|
|ENgetdemandmodel|Retrieves the type of demand model in use and its parameters|
|ENsetdemandmodel|Sets the type of demand model to use and its parameters|
|ENsetelseaction|Sets the properties of an ELSE action in a rule-based control|
|ENsetnodeid|Change the ID name for a node|
|ENsetlinkid|Change the ID name for a link|
|ENsetpipedata|Sets a group of properties for a pipe link|
|ENsetpremise|Sets the properties of a premise in a rule-based control|
|ENsetpremiseindex|Sets the index of an object in a premise of a rule-based control|
|ENsetpremisestatus|Sets the status being compared to in a premise of a rule-based control|
|ENsetpremisevalue|Sets the value in a premise of a rule-based control|
|ENsetrulepriority|Sets the priority of a rule-based control|
|ENsettankdata|Sets a group of properties for a tank node|
|ENsetthenaction|Sets the properties of a THEN action in a rule-based control|
|ENgettitle|Retrieves the title lines of the project|
|ENsettitle|Sets the title lines of the project|
|ENsetlinknodes|Sets the indexes of a link's start- and end-nodes|
|ENsetlinktype|Changes the type of a particular link (e.g. pipe to pump)|
|ENgetdemandname|Gets the name of a node's demand category|
|ENsetdemandname|Assigns a name to a node's demand category|
|ENgetcomment|Retrieves the comment string assigned to the object (NODE, LINK, TIMEPAT or CURVE)|
|ENsetcomment|Sets the comment string assigned to the object (NODE, LINK, TIMEPAT or CURVE)|
|ENdeletepattern|Deletes a time pattern from a project|
|ENdeletecurve|Deletes a data curve from the project|
|ENdeletecontrol|Deletes an existing simple control|
|ENdeleterule|Deletes an existing rule-based control|
|ENsetjuncdata|Sets a group of properties for a junction node|
|ENgetvertex|Retrieves the coordinate's of a vertex point assigned to a link|
|ENgetvertexcount|Retrieves the number of internal vertex points assigned to a link|
|ENsetvertices|Assigns a set of internal vertex points to a link|
|ENgetresultindex|Retrieves the order in which a node's or link's results were saved to an output file|

&uparrow; [Back to top](#table-of-contents)
