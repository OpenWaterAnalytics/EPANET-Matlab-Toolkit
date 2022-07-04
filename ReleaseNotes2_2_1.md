## New functionalities
* The epanet function can now run without an input:
  e.g., d = epanet
* Use of EN_ functions, using a project handle (ph). This creates the possibility of 
  of working with more than one nets, i.e., many epanet classes simultaneously: e.g., d = epanet(netname, 'ph');
  Check examples/EX6_load_two_inp_files.m.   
* Added api in front of EN and MSX functions (e.g., apiENopen and apiMSΧopen)
* New examples using only apiEN and apiMSX functions:
  * Toolkit_api_EX1_using_EN_functions
  * Toolkit_api_EX2_using_MSX_functions
* Help text for all apiEN and apiMSX functions
* Help text for all MSX functions 


### Updated functions (accept more input data)
  * getNodeTankNameID
  * getLinkValveNameID
  * getLinkPumpNameID  
  * getNodeTankIndex 
  * getLinkValveIndex
  * getLinkPumpIndex
  * getMSXComputedQualitySpecie 



## List of New EPANET 2.2.1 Matlab Class Functions
|Function|Description|
|---------|---------|
|appRotateNetwork|Rotates the network by theta degrees counter-clockwise|
|appShiftNetwork|Shifts the network in x and y directions|
|createProject|Creates an epanet project|
|deleteProject|Deletes an epanet project|
|getEN_functionsImpemented|Retrieves the epanet EN_ functions that have been developed|
|getGraph|Retrieves the graph of the current epanet network|
|getLinkQuality|Retrieves the value of link quality|
|plotGraph|Plots the graph of the current epanet network|
|runProject|Runs a complete EPANET simulation|
|setNodeTypeJunction|Transforms a node to junction|
|setNodeTypeReservoir|Transforms a node to reservoir|
|setNodeTypeTank|Transforms a node to tank|
|splitPipe|Splits a pipe, creating two new pipes and adds a junction in between them|
|toJson|Creates json text variable|
|toJsonFile|Creates a .json file and adds the input values in json format|
