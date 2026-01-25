### EPANET Matlab Toolkit (EMT) v2.3.3.0

It introduces full support for new EPANET 2.3 API calls, pressure-driven analysis (PDA), leakage modeling, extended curve and valve types, new units, enhanced control/rule handling, and expanded project I/O capabilities.  

- Restores initial fully-open valve status handling
- Fixes PCV update behavior via API and controls
- Fixes `EN_setpipedata` minor loss coefficient assignment
- Fixes `EN_setlinkvalue` Open/Closed status assignment

New test file: [testFunctions2_3](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/dev/tests/testFunctions2_3.m)

Leakage example: [EX28_Leakage](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/dev/examples/m-files/EX28_Leakage.m)

#### List of New EPANET 2.3 Functions Supported
| Function              | Description
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `apiENgetcontrolenabled` | Get enabled/disabled flag for a simple control by index. Returns boolean/int (1 = enabled, 0 = disabled).                                                  |
| `apiENsetcontrolenabled` | Enable or disable a specific simple control by index.                                                                                                      |
| `apiENgetruleenabled`    | Get enabled/disabled flag for a rule-based control by rule index.                                                                                          |
| `apiENsetruleenabled`    | Enable or disable a rule by rule index.                                                                                                                    |
| `apiENgetlinkvalues`     | Bulk-retrieve property values for all links for a given link property code. Returns an array/vector.                                                       |
| `apiENgetnodevalues`     | Bulk-retrieve property values for all nodes for a given node property code. Returns an array/vector.                                                       |
| `apiENsetcurvetype`      | Set the type of a curve object (volume, pump, efficiency, headloss, generic, valve).                                                                       |
| `apiENsetvertex`         | Set the coordinates of a link’s intermediate vertex (polyline point) by vertex index.                                                                      |
| `apiENtimetonextevent`   | Return the type of event that terminates the current time step (e.g., hydraulic step, WQ step, tank level event, control event) and the time to the event. |
| `apiENloadpatternfile`   | Load time patterns from an external file into the current project under a specific pattern ID.                                                             |
| `apiENopenX`             | Open an input file even if it has certain formatting errors (lenient parsing mode).                                                                        |


#### EMT Functions added/updated

##### Valve Links
| Function               | Description                                                                   |
| ---------------------- | ----------------------------------------------------------------------------- |
| `addLinkValvePCV`      | Create a new **Pressure Control Valve (PCV)** link between two nodes.         |
| `setLinkTypeValvePCV`  | Convert an existing link to a **PCV (Pressure Control Valve)** type.          |
| `getLinkValveCurveGPV` | Retrieve the control curve associated with a **GPV (General Purpose Valve)**. |
| `setLinkValveCurveGPV` | Assign or update the control curve for a **GPV** valve.                       |
| `getLinkValveCurvePCV` | Retrieve the control curve associated with a **PCV** valve.                   |
| `setLinkValveCurvePCV` | Assign or update the control curve for a **PCV** valve.                       |

##### Curves
| Function                 | Description                                    |
| ------------------------ | ---------------------------------------------- |
| `setCurveType`           | Set the type of a curve object.                |
| `setCurveTypeVolume`     | Set the type of a curve object to Volume.          |
| `setCurveTypePump`       | Set the type of a curve object to Pump.            |
| `setCurveTypeEfficiency` | Set the type of a curve object to Efficiency.     |
| `setCurveTypeHeadloss`   | Set the type of a curve object to Headloss.        |
| `setCurveTypeGeneral`    | Set the type of a curve object to General. |
| `setCurveTypeValveCurve` | Set the type of a curve object to Curve.   |


##### Controls & Rules
| Function           | Description                                                             |
| ------------------ | ---------------------------------------------------------------------- |
| `getControlState`  | Retrieve simple control enabled flag. |
| `getRuleEnabled`   | Retrieve rule enabled flag.             |
| `setRuleEnabled`   | Enable or disable a rule-based control.                                 |
| `getLinkInControl` | Identify links referenced by control rules.                             |
| `getNodeInControl` | Identify nodes referenced by control rules.                             |


##### Events
| Function             | Description                                                                                                                   |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `getTimetoNextEvent` | Return the type of event that terminates the current time step. |

##### Patterns
| Function                        | Description                                |
| ------------------------------- |--------------------------------------------|
| `loadPatternFile`               | Load a time pattern file into the project. |
| `getPatternAverageDefaultValue` | Report average value used when a pattern is missing.                                          |


##### Options & Units
| Function                  | Description                                |
| ------------------------- |--------------------------------------------|
| `getOptionsPressureUnits`   | Retrieve current pressure units. |
| `setOptionsPressureUnits` | Set pressure units (EN_PSI / EN_KPA / EN_METERS).                                          |
| `setOptionsPressureUnitsMeters` | Set pressure units to meters.                                          |
| `setOptionsPressureUnitsPSI` | Set pressure units to PSI.                                          |
| `setOptionsPressureUnitsKPA` | Set pressure units to kPa.                                          |
| `getOptionsStatusReport` | Get current status report level.                                          |
| `setOptionsStatusReport` | Set status report level (EN_NO_REPORT / EN_NORMAL_REPORT / EN_FULL_REPORT).                                          |
| `setOptionsStatusReportNo` | Set report level to no report.                                          |
| `setOptionsStatusReportNormal` | Set report level to normal.                                          |
| `setOptionsStatusReportFull` | Set report level to full.                                          |
| `getOptionsDemandPattern` | Get default demand pattern behavior.                                          |
| `setOptionsDemandPattern` | Set default demand pattern behavior.                                          |
| `getOptionsEmitterBackFlow` | Get emitter backflow setting.                                          |
| `setOptionsEmitterBackFlowAllowed` | Allow emitter backflow.                                          |
| `setOptionsEmitterBackFlowDisallowed` | Disallow emitter backflow.                                          |
| `setFlowUnitsCMS` | Set flow units to CMS.                                          |

##### Leakage & Demand
| Function                        | Description                                     |
| ------------------------------- |-------------------------------------------------|
| `getLinkLeakArea`               | Get leakage area for a link.                    |
| `setLinkLeakArea` | Set leakage area for a link.                    |
| `getLinkExpansionProperties` | Get expansion properties for a link.            |
| `setLinkExpansionProperties` | Set expansion properties for a link.            |
| `getLinkLeakageRate` | Get leakage rate for a link (read only).        |
| `getNodeLeakageFlow` | Get leakage flow at a node (read only).         |
| `getNodeEmitterFlow` | Get emitter flow at a node (read only).         |
| `getConsumerDemandRequested` | Retrieve consumer demand requested (read only). |
| `getConsumerDemandDelivered` | Retrieve consumer demand delivered (read only). |

##### Project I/O
| Function                        | Description                                |
| ------------------------------- |--------------------------------------------|
| `openX`               | Open an input file even if it has formatting errors (lenient parsing mode). |


##### Statistics
| Function                       | Description                                                |
| ------------------------------ |------------------------------------------------------------|
| `getStatisticIterations`       | Retrieves the number of iterations taken in the simulation . |
| `getStatisticRelativeError`    | Retrieves the relative error statistic from the simulation. |
| `getStatisticDeficientNodes`   | Retrieve number of deficient nodes.                        |
| `getStatisticDemandReduction`  | Retrieve demand reduction statistics.                      |
| `getStatisticTotalLeakageLoss` | Retrieve total leakage loss value.                   |


##### Constants / Enums

`d.ToolkitConstants.EN_`

```
EN_NODE_INCONTROL = 28 
EN_EMITTERFLOW = 29 
EN_LEAKAGEFLOW = 30 
EN_DEMANDFLOW = 31 
EN_FULLDEMAND = 32 
EN_LINK_INCONTROL = 23 
EN_GPV_CURVE = 24 
EN_PCV_CURVE = 25 
EN_LEAK_AREA = 26 
EN_LEAK_EXPAN = 27 
EN_LINK_LEAKAGE = 28
```
###### Link & unit types
```
EN_PCV = 9 
EN_CMS = 10
```
###### Pressure unit type
```
EN_PSI = 0 
EN_KPA = 1 
EN_METERS = 2
```
###### Demand model
```
EN_DDA = 0 # demand-driven 
EN_PDA = 1 # pressure-driven
```
###### Option codes
```
EN_DEMANDPATTERN = 23 
EN_EMITBACKFLOW = 24 
EN_PRESS_UNITS = 25 
EN_STATUS_REPORT = 26
```
###### Pump curve type
```
EN_CONST_HP = 0 # Constant horsepower 
EN_POWER_FUNC= 1 # Power function 
EN_CUSTOM = 2 # User-defined 
EN_NOCURVE = 3 # No pump curve
```
###### Curve classes
```
EN_VOLUME_CURVE = 0 
EN_PUMP_CURVE = 1 
EN_EFFIC_CURVE = 2 
EN_HLOSS_CURVE = 3 
EN_GENERIC_CURVE= 4 
EN_VALVE_CURVE = 5
EN_UNCONDITIONAL = 0 
EN_CONDITIONAL = 1
```

###### Status report level
```
EN_NO_REPORT = 0 
EN_NORMAL_REPORT= 1 
EN_FULL_REPORT = 2
```

###### Statistics (added leakage loss)
``` EN_LEAKAGELOSS = 7 ```

###### Timestep end events
``` 
EN_STEP_REPORT = 0 
EN_STEP_HYD = 1 
EN_STEP_WQ = 2 
EN_STEP_TANKEVENT = 3 
EN_STEP_CONTROLEVENT= 4
EN_MISSING = -1.0E10 
EN_SET_CLOSED= -1.0E10 
EN_SET_OPEN = 1.0E10
```

### EPANET Matlab Toolkit (EMT) v2.2.6.1

- Add the example [Fasted Parallel computation](./EX27_Fasted_Parallel_computation.m) 
- Update the function `getComputedTimeSeries_ENepanet` to work in parallel
- Add the function `readEpanetBinaryFile`
- Add the function `getNodeJunctionBaseDemands` #230 
- Add the function `getNodeJunctionActualDemand` #230 

### EPANET Matlab Toolkit (EMT) v2.2.5

- Add the function `getComputedAnalysisTimeSeries` (Computed Hydraulic and Quality analysis)
- Update Movie example (use epanet class function instead of properties)
- Load epanet file with specific DLL without fill properties: 
`d = epanet(inpname, 'epanet2', 'loadfile')`; 
- Update the functions `getAdjacencyMatrix`, `getFlowDirections`
  help `d.getAdjacencyMatrix` 
  help `d.getFlowDirections`
- Add the function `getMSXComputedLinkQualitySpecie` (Returns the link quality for specific specie)
- Add the function `getMSXComputedNodeQualitySpecie` (Returns the node quality for specific specie)
- Update in the function `setMSXTimeStep`
- Update the function `loadMSXFile`, e.g.  
    ```
  d.loadMSXFile('net2-cl2.msx', 'epanetmsx');
    d.loadMSXFile('net2-cl2.msx', 'epanetmsx', 'loadfile');
    d.loadMSXFile('net2-cl2.msx', 'loadfile');
  ```

### EPANET Matlab Toolkit (EMT) v2.2.4

- Update the LICENSE file to the last version of EUPL v. 1.2
- Minor fix in the function `setNodeBaseDemands` for the demand category (Thanks Roya @RPM-2022).
- Add the function `getMSXComputedTimeSeries`.
- Minor fix in the functions `plotMSXSpeciesLinkConcentration`, `plotMSXSpeciesNodeConcentration`
- Show a message when you set the wrong species index in the functions (getMSXComputed..)
- Add the library epanetmsx_thunk64 and msxepanet.m to help deploy an app with epanet msx.
- Add `loadMSXEPANETFile` and `loadMSXlibrary` functions. Used for parallel simulations.
- Minor fix in the function getLinkVolumes

### EPANET Matlab Toolkit (EMT) v2.2.3

- Create [CITATION.cff](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/CITATION.md)
- Create [networkviewer.mlapp](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/gui/networkviewer.mlapp)
- Minor fix in the function [getAdjacencyMatrix](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/8c1af5736dc445b754b6e7d9ac6547000b824e46)
- Add the function [getLinkVolumes](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/f7bc19a7c9fc9988510caae1ee7506b6385266ac)
- Update the function [getNodeLinks](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/0d4adad2d7b4f3555d32ad01a60612eb1337ae86) to return all the links to which all nodes are connected to
- Add the function [plotDiGraph](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/9c35084091202747d92f18ae5b050bc3061171ab)
- Add the function [getFlowDirections](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/e9a44771cf39bb4e09d616a82e8264028732dc80)
- Fix issue #222, a bug in highlightnode and legendposition
- Some [updates ](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/commit/d8bd5b4cec1b4a1846b86ce5f0121bfefde3db95)in the movie example 
- Fix issue #223, Error Recognizing setControl Function
- Merge release notes in one file

### EPANET Matlab Toolkit (EMT) v2.2.2

- Performance updates (getComputed..)
- Update [EX24_Parallel_computations.m](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/EX24_Parallel_computations.m)
- Cleanup and some fixes
- When running `d = epanet('Net1.inp')` creates a temporary input file using the Matlab function `copyfile` instead of `ENsaveinpfile` the first time.
- Update [Movie Examples](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/tree/master/examples/movie-example)
- We have prepared a tutorial as part of the CCWI-WDSA 2022 conference. You can find the files in the `tutorial` folder at the following [GitHub repository](https://github.com/KIOS-Research/CCWI2022-EMT-Tutorial).
  - Basic functionality
  - Run first EPANET analysis
  - Run first MSX analysis
  - How to include uncertainties in simulations
  - How to analyze network graph
  - How to place pressure sensors
  - How to create leakage events
  - How to detect leakage events
  - How to place quality sensors
  - How to create contamination events (using MSX)
  - How to detect contamination events

### EPANET Matlab Toolkit (EMT) v2.2.1

* The epanet class initialization function can now run without an input:
  i.e., d = epanet, so that the user can create a network from scratch. Check the following example: [Toolkit_EX4_Network_Building](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/Toolkit_EX4_Network_Building.mlx)


* Use of EN_ functions, using a project handle (ph). This creates the possibility 
  of working with more than one networks (EPANET Input Files), i.e., many EPANET classes simultaneously!
  Check the following example: [EX6_load_two_inp_files](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/EX6_load_two_inp_files.mlx)


* The user is now able to use EN and MSX functions explicitily. Api has been added in front of EN and MSX functions (e.g., apiENopen and apiMSΧopen). Check the new examples featuring only apiEN and apiMSX functions:
  * [Toolkit_api_EX1_using_EN_functions](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/Toolkit_api_EX1_using_EN_functions.mlx)
  * [Toolkit_api_EX2_using_MSX_functions](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/Toolkit_api_EX2_using_MSX_functions.mlx)


* Help text for all apiEN and apiMSX functions.
```
>> help d.apiENgetcount
--- help for epanet.apiENgetcount ---

  Retrieves the number of objects of a given type in a project.
 
  apiENgetcount(countcode, LibEPANET, ph)
 
  Parameters:
  countcode number of objects of the specified type.
  LibEPANET epanet library DLL name
  ph        epanet project handle.
 
  Returns:
  an error code.
```


* Help text for all MSX functions. 

* Fix bugs, cleanup, and updates.

* Transitioned Examples to Live Code File Format (.mlx). 

[List of all the examples](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/tree/master/examples#readme)


Updated functions to accept more input data, e.g.,               getMSXComputedQualitySpecie:

```
>> help d.getMSXComputedQualitySpecie
--- help for epanet/getMSXComputedQualitySpecie ---

  Returns the node/link quality for specific specie.
 
  Example 1:    
    d = epanet('net2-cl2.inp');
    d.loadMSXFile('net2-cl2.msx');
    MSX_comp = d.getMSXComputedQualitySpecie('CL2')
    MSX_comp.NodeQuality % row: time, col: node index
    MSX_comp.LinkQuality % row: time, col: link index
    MSX_comp.Time
 
  Example 2:
    d = epanet('example.inp');            
    d.loadMSXFile('example.msx');
    MSX_comp = d.getMSXComputedQualitySpecie  % Computes quality for all
                                               the species in 3D arrays.
    MSX_comp.NodeQuality(:,:,1) % Gets node quality for the first specie.
    MSX_comp.LinkQuality(:,:,5) % Gets link quality for the fith specie.
  
  Example 3:
    d = epanet('example.inp');            
    d.loadMSXFile('example.msx');
    MSX_comp = d.getMSXComputedQualitySpecie({'AStot', 'AS5s'}) % Computes quality for 'AStot', 'AS5s'
                                                                species in 3D arrays.
    MSX_comp.NodeQuality(:,:,1) % Gets node quality for the first specie.
    MSX_comp.NodeQuality(:,:,2) % Gets node quality for the second specie.
 
  See also getMSXComputedQualityNode, getMSXComputedQualityLink.
```
More updated functions: getNodeTankNameID, getLinkValveNameID, getLinkPumpNameID, getNodeTankIndex, getLinkValveIndex, getLinkPumpIndex.

## List of New EPANET 2.2.1 Matlab Class Functions
|Function|Description|
|---------|---------|
|appRotateNetwork|Rotates the network by theta degrees counter-clockwise|
|appShiftNetwork|Shifts the network in x and y directions|
|createProject|Creates an EPANET project|
|deletePatternAll|Deletes all time patterns from a project|
|deleteProject|Deletes an EPANET project|
|getEN_functionsImpemented|Retrieves the EPANET EN_ functions that have been developed|
|getGraph|Retrieves the graph of the current EPANET network|
|getNodeLinks|Retrieves the links which a specific node is connected to|
|getLinkQuality|Retrieves the value of link quality|
|plotGraph|Plots the graph of the current EPANET network|
|reverseLinkNodeIndices|Reverses the node indices that connect a link|
|runProject|Runs a complete EPANET simulation|
|setNodeTypeJunction|Transforms a node to junction|
|setNodeTypeReservoir|Transforms a node to reservoir|
|setNodeTypeTank|Transforms a node to tank|
|splitPipe|Splits a pipe, creating two new pipes and adds a junction in between them|
|toJson|Creates json text variable|
|toJsonFile|Creates a .json file and adds the input values in json format|

## List of New EPANET 2.2 Functions Supported
|Function|Description|
|---------|---------|
|apiEN_createproject|Creates a new EPANET project|
|apiEN_deleteproject|Deletes the EPANET project|
|apiEN_runproject|Runs a complete EPANET simulation|
