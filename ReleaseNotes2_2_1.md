## New functionalities
* The epanet class initialization function can now run without an input:
  i.e., d = epanet, so that the user can create a network from scratch. Check the following example: [Toolkit_EX4_Network_Building](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/Toolkit_EX4_Network_Building.m)


* Use of EN_ functions, using a project handle (ph). This creates the possibility 
  of working with more than one networks (EPANET Input Files), i.e., many EPANET classes simultaneously!
  Check the following example: [EX6_load_two_inp_files](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/master/examples/EX6_load_two_inp_files.m )


* The user is now able to use EN and MSX functions explicitily. Api has been added in front of EN and MSX functions (e.g., apiENopen and apiMSΧopen). Check the new examples featuring only apiEN and apiMSX functions:
  * [Toolkit_api_EX1_using_EN_functions](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/dev/examples/Toolkit_api_EX1_using_EN_functions.m)
  * [Toolkit_api_EX2_using_MSX_functions](https://github.com/OpenWaterAnalytics/EPANET-Matlab-Toolkit/blob/dev/examples/Toolkit_api_EX2_using_MSX_functions.m)


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
