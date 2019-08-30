<a href="http://www.kios.ucy.ac.cy"><img src="http://www.kios.ucy.ac.cy/templates/favourite/images/kios_logo_hover.png"/><a>


EPANET-MATLAB-Toolkit
==================================

## Table of Contents

- [How to cite](#how-to-cite)
- [Example titles](#example-titles)
- [Example descriptions](#example-descriptions)

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

## Example titles

- Example 1: <b> [EX1_Plot_network_topology] </b>, Plots the network topology.
- Example 2: <b> [EX2_Hydraulic_analysis] </b>, Runs the hydraulic analysis of a network.
- Example 3: <b> [EX3_Quality_analysis] </b>, Runs the quality analysis of a network.
- Example 4: <b> [EX4_Plot_time_series] </b>, Visualises/plots time series for node pressures, water velocity and water flow.
- Example 5: <b> [EX5_Plot_values_parameters] </b>, Plots all flow, pressure, etc values on network map for a specific time step.
- Example 6: <b> [EX6_load_two_inp_files] </b>, Loads 2 different networks.
- Example 7: <b> [EX7_set_pump_curves] </b>, Sets pump curves.
- Example 8: <b> [EX8_tanks_to_reservoirs] </b>, Replaces tanks with reservoirs.
- Example 9: <b> [EX9_compare_simulations] </b>, Compares hydraulics and quality analysis functions.
- Example 10: <b> [EX10_close_pipes_during_sim] </b>, Closes pipes during simulation.
- Example 11: <b> [EX11_assing_new_curve_pump] </b>, Assings a new curve to a specific pump.
- Example 12: <b> [EX12_add_multiple_controlpatterns_bin] </b>, Defines pattern time for valves setting via bin functions.
- Example 13a: <b> [EX13a_add_cvpipe_junction] </b>, Adds CV pipe and junction in a network with the bin Function addBinJunction.
- Example 13b: <b> [EX13b_add_cvpipe_bin] </b>, Adds CV pipe in a network with Bin Function addBinCVPipe.
- Example 13c: <b> [EX13c_add_cvpipe] </b>, Adds a CV Pipe in a network.
- Example 14: <b> [EX14_hydraulic_and_quality_analysis] </b>, Runs the Hydraulic and Quality analysis of a network.
- Example 15: <b> [EX15_write_msx_file] </b>, Writes MSX File e.g.Net2.msx.
- Example 16: <b> [EX16_create_multiple_scenarios] </b>, Creates scenarios based on node count of inp file EPANET.
- Example 17a: <b> [EX17a_add_multiple_controls_pipestatus_bin] </b>, Defines pattern time for pipe status via bin & normal functions.
- Example 17b: <b> [EX17b_add_multiple_controls_pipestatus] </b>, Defines pattern time for pipe status.
- Example 18: <b> [EX18_change_status_pipes] </b>, Changes status randomly at pipes during simulation.
- Example 19: <b> [EX19_rotate_network] </b>, Rotates EPANET Inp File.
- Example 20a: <b> [EX20a_external_controls] </b>, a) External Controls, Changes control status of pump STEP BY STEP.
- Example 20b: <b> [EX20b_external_controls] </b>, b) External Controls, Adds controls in hydraulic analysis STEP-BY-STEP.
- Example 21: <b> [EX21_Pressure_driven_analysis_option] </b>, Tests the Pressure Driven Analysis option.
- Example 22: <b> [EX22_Overflow_option_for_tanks] </b>, Tests the overflow option for tanks.

&uparrow; [Back to top](#table-of-contents)

## Example descriptions

- Example 1: <b> EX1_Plot_network_topology </b>

Plots the network topology.

This example contains:
 Load a network.
 Plot network topology.
 Plot node IDs.
 Plot links IDs.
 Plot link indices.
 Plot node indices.
 Highlight specific links.
 Highlight specific nodes.
 Plot only links.
 Plot only nodes.
 Highlight multiple links with different colors.
 Highlight multiple nodes with different colors.
 Hide legend.
 Extent network on figure.
 Legend Position.
 Plot at specific axes.
 Plot for matlab apps (GUIs).
 Plot colors.
 Unload library.
 
 
- Example 2: <b> EX2_Hydraulic_analysis </b>

Runs the hydraulic analysis of a network.

This example contains:
 Load a network.
 Set simulation time duration.
 Hydraulic analysis using ENepanet binary file.
 Hydraulic analysis using epanet2d.exe binary file.
 Hydraulic analysis.
 Hydraulic analysis step-by-step.
 Unload library.
 
 
 - Example 3: <b> EX3_Quality_analysis </b>
 
Runs the quality analysis of a network.

This example contains: 
 Load a network.
 Compute Quality without MSX.
 Compute Quality step by step.
 Load EPANET-MSX files.
 Compute Quality with MSX (specify type).
 Get quality of specific nodes.
 Get quality of specific links.
 Get species names.
 Get quality for specific species type (nodes and links).
 Plot concentration for specific node and all species.
 Plot concentration for specific link and all species.
 Unload libraries.


- Example 4: <b> EX4_Plot_time_series </b>

Visualises/plots time series for node pressures, water velocity and water flow.

This example contains:
 Load a network.
 Hydraulic analysis using ENepanet binary file.
 Change time-stamps from seconds to hours.
 Plot node pressures for specific nodes.
 Plot water velocity for specific links.
 Plot water flow for specific links.
 Unload library.
 
 
- Example 5: <b> EX5_Plot_values_parameters </b>

Plots all flow, pressure, etc values on network map for a specific time step.

This function contains:
 Load a network.
 Hydraulic analysis using ENepanet binary file.
 Get node and link info.
 Set custom offset for text.
 Get plot axes.
 Plot flow values on map.
 Plot pressure values on map.
 Unload library.
 
 
- Example 6: <b> EX6_load_two_inp_files </b>

Loads 2 different networks.

This example contains:
 Using 2 dll files. Load 2 Input files.
 Load networks.
 Disp elevations for the two networks.
 Plot networks.
 Unload libraries.
 Delete files who created.
 
 
- Example 7: <b> EX7_set_pump_curves </b>

Sets pump curves.

This example contains:
 Load a network.
 Computed Hydraulic & Quality Time Series via ENepanet binary file.
 Plot pressures.
 Get head curve.
 Set new head curve values.
 Computed hydraulics.
 Unload library.
 
 
- Example 8: <b> EX8_tanks_to_reservoirs </b>

Replaces tanks with reservoirs.

This example contains:
 Choose Network.
 Load network and simulate.
 Find tank patterns.
 Create patterns for new reservoirs from tanks.
 Replace tanks with reservoirs.
 Unload library.
 
 
- Example 9: <b> EX9_compare_simulations </b>

Compares hydraulics and quality analysis functions.
 
This example contains:
 Load a network.
 Set simulation duration.
 Test hydraulics and quality analysis functions.
 Step by step hydraulic analysis.
 Step by step quality analysis.
 Unload library.
 Run time d.getComputedTimeSeries.
 
 
- Example 10: <b> EX10_close_pipes_during_sim </b>

Closes pipes during simulation.

This example contains:
 Load a network.
 Link index for change the status.
 Run step by step hydraulic analysis.
 Get flows for the specific link index.
 Unload library.
 
 
- Example 11: <b> EX11_assing_new_curve_pump </b>

Assings a new curve to a specific pump.

This example contains:
 Load a network.
 Add new curve in the network.
 Get pump index.
 Get head curve index.
 Assing new curve index on the specific pump.
 Unload library.
 
 
- Example 12: <b> EX12_add_multiple_controlpatterns_bin </b>

Defines pattern time for valves setting via bin functions.

This example contains:
 Add paths and load inp file.
 Get times.
 Create random settings.
 Settings for valve 'VALVE-173'.
 Check pressures at node after valve.
 Run simulation with executable.
 
 
- Example 13a: <b> EX13a_add_cvpipe_junction </b>

Adds CV pipe with the bin Function addBinJunction.

This example contains:
 Load a network.
 Plot newtork with bin functions.
 Add Junction + CV pipe.
 Plot normal with new pipe.
 Unload library.
 
 
- Example 13b: <b> EX13b_add_cvpipe_bin </b>

Adds CV pipe in a network with Bin Function addBinCVPipe.

This example contains:
 Load a network.
 Get from and to nodes for add cv pipe.
 Plot network.
 Plot network with changes.
 Unload library.
 
 
- Example 13c: <b> EX13c_add_cvpipe </b>

Adds a CV Pipe in a network.

This example contains:
 Load a network.
 Plot network components.
 Get from and to nodes for add cv pipe.
 Plot new network with changes.
 Unload library.
 
 
- Example 14: <b> EX14_hydraulic_and_quality_analysis </b>

Runs the Hydraulic and Quality analysis of a network.
 
This example contains:
 Load a network.
 Hydraulic and Quality analysis STEP-BY-STEP.
 Display nodes pressures, links flows, nodes actual qualities, links actual qualities.
 Unload library. 
 
 
- Example 15: <b> EX14_hydraulic_and_quality_analysis </b>

Writes MSX File e.g.Net2.msx.

This example contains:
 Load a network.
 Set Input Arguments:
    Filename
    Section Title
    Section Options
    Section Species
    Section Coefficients
    Section Terms
    Section Pipes
    Section Tanks
    Section Sources
    Section Quality Global
    Section Quality
    Section Parameters
    Section Patterns
 Write MSX File.
 Load MSX File.
 Compute.
 Unload libraries.
 
 
- Example 16: <b> EX16_create_multiple_scenarios </b>

Creates scenarios based on node count of inp file EPANET.

This example contains:
 Load a network.
 Create scenario parameters.
 Create Multiple Scenarios.
 Set quality type.
 Run Scenarios.
 Get computed quality time series.
 Plot Quality VS Time.
 Unload library.
 
 
- Example 17a: <b> EX17a_add_multiple_controls_pipestatus_bin </b>
 
Defines pattern time for pipe status via bin & normal functions.
 
This example contains:
 Add paths and load inp file.
 Get times.
 Get Pipe infos.
 Find the max length of link pipe names.
 Times.
 Create random status for pipes every hour.
 Add all controls.
 Unload libraries.
 
 
- Example 17b: <b> EX17b_add_multiple_controls_pipestatus </b>

Defines pattern time for pipe status.

This example contains:
 Add paths and load inp file.
 Get times.
 Get Pipe infos.
 Find the max length of link pipe names.
 Times.
 Create random status for pipes every hour.
 Unload library.
 
 
- Example 18: <b> EX18_change_status_pipes </b>

Changes status randomly at pipes during simulation.

This example contains:
 Load a network.
 Get pipe count.
 Get pipe indices.
 Run step by step hydraulic analysis.
 Set status random 0/1 for pipes.
 Plot flows.
 Unload library.
 
 
- Example 19: <b> EX19_rotate_network </b>
 
Rotates EPANET Inp File.

This example contains:
 Load network and paths.
 Plot nework initial.
 Rotate degrees theta.
 Define the x- and y-data for the original line we would like to rotate.
 Create a matrix of these points, which will be useful in future calculations.
 Choose a point which will be the center of rotation.
 Create a matrix which will be used later in calculations.
 Define a 60 degree counter-clockwise rotation matrix.
 Do the rotation:
   Shift points in the plane so that the center of rotation is at the origin.
   Apply the rotation about the origin.
   Shift again so the origin goes back to the desired center of rotation.
 Pick out the vectors of rotated x- and y-data.
 Plot rotated network.
 Save rotated Inp file.
 
 
- Example 20a: <b> EX20a_external_controls </b>
 
a) External Controls, Changes control status of pump STEP BY STEP.
 
This example contains:
 Load network.
 Delete Controls.
 Hydraulic analysis STEP-BY-STEP.
 CONTROLS.
 Add new controls in live.
 Unload library.
 
 
- Example 20b: <b> EX20b_external_controls </b>

b) External Controls, Adds controls in hydraulic analysis STEP-BY-STEP.

This example contains:
 Load network.
 Delete Controls.
 Hydraulic analysis STEP-BY-STEP.
 CONTROLS.
 Add new controls in live.
 Delete controls.
 Unload library.
 
 
- Example 21: <b> EX21_Pressure_driven_analysis_option </b>

Tests the Pressure Driven Analysis option.

This example contains:
 Set Demand Multiplier to 10 to cause negative pressures.
 Run single period analysis.
 Solve hydraulics with default DDA option which will return with neg. pressure warning code.
 Check that 4 demand nodes have negative pressures.
 Switch to PDA with pressure limits of 20 - 100 psi.
 Solve hydraulics again.
 Check that 6 nodes had demand reductions totaling 32.66%.
 Check that Junction 12 had full demand.
 Check that Junction 21 had a demand deficit of 413.67.
 Unload library.
 
 
- Example 22: <b> EX22_Overflow_option_for_tanks </b>

Tests the overflow option for tanks.

This example contains:
 Get index of the tank and its inlet/outlet pipe.
 Set initial & maximum level to 130.
 Set duration to 1 hr.
 Solve hydraulics with default of no tank spillage allowed.
 Check that tank remains full.
 Check that there is no spillage.
 Check that inflow link is closed.
 Turn tank overflow option on.
 Solve hydraulics again.
 Check that tank remains full.
 Check that there is spillage equal to tank inflow
 (inflow has neg. sign since tank is start node of inflow pipe).
 Save project to file and then close it.
 Re-open saved file & run it.
 Check that tank spillage has same value as before.
 Unload library.
 
 &uparrow; [Back to top](#table-of-contents)