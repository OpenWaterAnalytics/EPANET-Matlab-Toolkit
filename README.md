EPANET-Matlab-Class
===================
The `EPANET-Matlab Class` is an open-source software which operates within the Matlab environment, for providing a programming interface for the latest version of [EPANET](https://github.com/OpenWaterAnalytics/epanet), a hydraulic and quality modeling software created by the US EPA, with Matlab, a  high-level technical computing software. The goal of the Matlab Class is to serve as a common programming framework for research and development in the growing field of smart water networks.

The `EPANET-Matlab Class` features easy to use commands/wrappers for viewing, modifying, simulating and plotting results produced by the EPANET libraries.  

# Requirements #
* [Matlab](http://www.mathworks.com/)
* [EPANET](https://github.com/OpenWaterAnalytics/epanet) The OpenWaterAnalytics (OWA) EPANET version is recommended instead of the one provided by the EPA.gov website, as a number of bugs have been fixed in the OWA version. 

# How to use the class #
Download the folder `Download ZIP`, set the run path in Matlab within the saved folder, and run `TEST.m`. This will execute all the commands which have been implemented in the Class.

Example:

d=epanet('Net1.inp')

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


