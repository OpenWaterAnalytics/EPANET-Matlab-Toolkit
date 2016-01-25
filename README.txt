EPANET-Matlab-Class dev-2.1

32 bit 
....... 
testFunctions.m [OK]
testMsxFunctions.m [OK]
testBinFunctions.m [OK]
testBinWithoutLib.m [OK]

64bit
.......
testFunctions.m [OK]
testMsxFunctions.m [OK]
testBinFunctions.m [OK]
testBinWithoutLib.m [OK]


----------------------------------------
New functions in 2.1:
Function                    Description
-----------------------------------------

ENgetpumptype               Retrieves the type of a pump for specific link index
ENgetheadcurveindex         Retrieves index of a head curve for specific link index
ENgetcoord                  Retrieves coordinate x, y for a node
ENsetcoord                  Sets coordinate x, y for a node
ENgetstatistic              Retrieves hydraulic simulation statistic
ENgetnumdemands             Retrieves the number of demand categories for a node
ENgetbasedemand             Retrieves the node's base demand for a category
ENgetdemandpattern          Retrieves the index of a demand pattern for a specific demand category of a node
ENgetaveragepatternvalue 	Retrieves the average value of a pattern
ENsetbasedemand             Sets the node's base demand for a category
ENgetcurveid                Retrieves ID of a curve with specific index
ENgetcurvelen               Retrieves number of points in a curve
ENgetcurveindex             Retrieves index of curve with specific ID
ENgetcurvevalue             Retrieves x,y point for a specific point number and curve
ENsetcurve                  Sets x,y values for a specific curve
ENaddcurve                  Adds a new curve appended to the end of the existing curves
ENsetcurvevalue             Sets x,y point for a specific point and curve

