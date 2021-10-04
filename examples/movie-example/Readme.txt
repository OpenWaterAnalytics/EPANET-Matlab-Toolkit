This directory contains various Matlab routines that add
additional capabilities to the standard Epanet and Epanet-MSX
toolkits. Documentation on usage of these routines are included
in the code preambles.

The utilities available include:

NetworkFrame

Plots a static frame of a network graph described in Epanet
input format, with nodes and/or links colored using specified
values. On first call for a given network, the line and node
objects are rendered, and the object handles stored so that
subsequent calls are more efficient (e.g., for generating movie
frames recording simulation results).

NetworkMovie

Uses NetworkFrame() to plot frames of the network graph with
nodes and links colored by the columns of matrices V and L. The
resulting animation is stored as an *.avi movie.

getQualityData
Epanet or Epanet-MSX Matlab application to run a simulation and return 
quality results

main_create_quality_movie

Example to illustrate the use of NetworkMovie and NetworkFrame
to animate either a specified bulk or wall specie, or both, for
an Epanet/MSX simulation, and store the result as an *.avi
movie

main_create_hydraulic_movie

Example to illustrate the use of NetworkMovie and NetworkFrame
to animate pressure and flows and store the result as an *.avi
movie


Original code: 
--------------
https://github.com/OpenWaterAnalytics/epanet-matlab/tree/master/Utils

---

The original code modified to work with EPANET-Matlab-toolkit