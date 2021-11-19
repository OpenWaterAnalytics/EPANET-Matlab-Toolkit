## New functionalities
* The epanet function can now run without an input:
  e.g., d = epanet
* Added api in front of EN and MSX functions (e.g., apiENopen and apiMSΧopen)
* New examples using only apiEN and apiMSX functions:
  * Toolkit_api_EX1_using_EN_functions
  * Toolkit_api_EX2_using_MSX_functions
* Help text for all apiEN and apiMSX functions 
## List of EPANET 2.2.1 Matlab Class Functions
|Function|Description|
|---------|---------|
|appShiftNetwork|Shifts the network in x and y directions|
|appRotateNetwork|Rotates the network by theta degrees counter-clockwise|
|getLinkQuality|Retrieves the value of link quality|
|setNodeTypeJunction|Transforms a node to junction|
|setNodeTypeReservoir|Transforms a node to reservoir|
|setNodeTypeTank|Transforms a node to tank|
|splitPipe|Splits a pipe, creating two new pipes and adds a junction in between them|

