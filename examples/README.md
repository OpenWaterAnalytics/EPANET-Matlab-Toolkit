<a href="http://www.kios.ucy.ac.cy"><img src="http://www.kios.ucy.ac.cy/templates/favourite/images/kios_logo_hover.png"/><a>


EPANET-MATLAB-Toolkit - Examples
==================================

## Table of Contents

- [Examples](#examples)
- [Descriptions](#descriptions)

## Examples
- [Example 1:](#example-1---plot-network-topology) Plot network topology. ``EX1_Plot_network_topology``
- [Example 2:](#example-2---hydraulic-analysis) Hydraulic analysis. `` EX2_Hydraulic_analysis``
- [Example 3:](#example-3---quality-analysis) Quality analysis. ``EX3_Quality_analysis``
- [Example 4:](#example-4---plot-time-series) Plot time series. ``EX4_Plot_time_series``
- [Example 5:](#example-5---plot-values-parameters) Plot values parameters. ``EX5_Plot_values_parameters``
- [Example 6:](#example-6---load-two-inp-files) Load two inp files. ``EX6_load_two_inp_files``
- [Example 7:](#example-7---set-pump-curves) Set pump curves. ``EX7_set_pump_curves``
- [Example 8:](#example-8---tanks-to-reservoirs) Tanks to reservoirs. ``EX8_tanks_to_reservoirs``
- [Example 9:](#example-9---compare-simulations) Compare Simulations. ``EX9_compare_simulations``
- [Example 10:](#example-10---close-pipes-during-simulation) Close pipes during simulation. ``EX10_close_pipes_during_sim``
- [Example 11:](#example-11---assign-a-new-curve-to-a-pump) Assign a new curve to a pump. ``EX11_assing_new_curve_pump``
- [Example 12:](#example-12---add-multiple-control-patterns) Add multiple control patterns. ``EX12_add_multiple_controlpatterns_bin``
- [Example 13a:](#example-13a---add-pipe-and-junction-via-bin-function) Add pipe and junction via bin function. ``EX13a_add_cvpipe_junction``
- [Example 13b:](#example-13b---add-pipe-via-bin-function) Add pipe via bin function. ``EX13b_add_cvpipe_bin`` 
- [Example 13c:](#example-13c---add-pipe-via-normal-function) Add pipe via normal function. ``EX13c_add_cvpipe`` 
- [Example 14:](#example-14---hydraulic-and-quality-analysis) Hydraulic and Quality analysis. ``EX14_hydraulic_and_quality_analysis``
- [Example 15:](#example-15---write-msx-file) Write MSX file. ``EX15_write_msx_file``
- [Example 16:](#example-16---creating-scenarios) Creating scenarios. ``EX16_create_multiple_scenarios``
- [Example 17a:](#example-17a---add-controls-via-bin-functions) Add controls via bin functions. ``EX17a_add_multiple_controls_pipestatus_bin``
- [Example 17b:](#example-17b---add-controls-via-normal-functions) Add controls via normal functions. ``EX17b_add_multiple_controls_pipestatus``
- [Example 18:](#example-18---change-status-of-pipe-during-simulation) Change status of pipe during simulation. ``EX18_change_status_pipes``
- [Example 19:](#example-19---rotate-network) Rotate network. ``EX19_rotate_network``
- [Example 20a:](#example-20a---change-of-external-controls) Change of external Controls. ``EX20a_external_controls``
- [Example 20b:](#example-20b---addition-of-external-controls) Addition of external Controls. ``EX20b_external_controls``
- [Example 21:](#example-21---pressure-driven-analysis) Pressure driven analysis. ``EX21_Pressure_driven_analysis_option``
- [Example 22:](#example-22---overflow-of-tanks) Overflow of tanks. ``EX22_Overflow_option_for_tanks``
- [Example 23:](#example-23---change-connection-of-links) Change connection of the links. ``EX23_Change_connection_links``
- [Example 24:](#example-24---delete-all-patterns) Delete all patterns. ``EX24_delete_all_patterns``

&uparrow; [Back to top](#table-of-contents)

## Descriptions

#### Example 1 - Plot network topology

Plots the network topology.

This example contains:
- Load a network.
- Plot network topology.
- Plot node IDs.
- Plot links IDs.
- Plot link indices.
- Plot node indices.
- Highlight specific links.
- Highlight specific nodes.
- Plot only links.
- Plot only nodes.
- Highlight multiple links with different colors.
- Highlight multiple nodes with different colors.
- Hide legend.
- Extent network on figure.
- Legend Position.
- Plot at specific axes.
- Plot for matlab apps (GUIs).
- Plot colors.
- Unload library.

&uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 2 - Hydraulic analysis

Runs the hydraulic analysis of a network.

This example contains:
- Load a network.
- Set simulation time duration.
- Hydraulic analysis using ENepanet binary file.
- Hydraulic analysis using epanet2d.exe binary file.
- Hydraulic analysis.
- Hydraulic analysis step-by-step.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 3 - Quality analysis
 
Runs the quality analysis of a network.

This example contains: 
- Load a network.
- Compute Quality without MSX.
- Compute Quality step by step.
- Load EPANET-MSX files.
- Compute Quality with MSX (specify type).
- Get quality of specific nodes.
- Get quality of specific links.
- Get species names.
- Get quality for specific species type (nodes and links).
- Plot concentration for specific node and all species.
- Plot concentration for specific link and all species.
- Unload libraries.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 4 - Plot time series

Visualises/plots time series for node pressures, water velocity and water flow.

This example contains:
- Load a network.
- Hydraulic analysis using ENepanet binary file.
- Change time-stamps from seconds to hours.
- Plot node pressures for specific nodes.
- Plot water velocity for specific links.
- Plot water flow for specific links.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 5 - Plot values parameters

Plots all flow, pressure, etc values on network map for a specific time step.

This function contains:
- Load a network.
- Hydraulic analysis using ENepanet binary file.
- Get node and link info.
- Set custom offset for text.
- Get plot axes.
- Plot flow values on map.
- Plot pressure values on map.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 6 - Load two inp files

Loads 2 different networks.

This example contains:
- Using 2 dll files. Load 2 Input files.
- Load networks.
- Disp elevations for the two networks.
- Plot networks.
- Unload libraries.
- Delete files who created.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 7 - Set pump curves

Sets pump curves.

This example contains:
- Load a network.
- Computed Hydraulic & Quality Time Series via ENepanet binary file.
- Plot pressures.
- Get head curve.
- Set new head curve values.
- Computed hydraulics.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 8 - Tanks to reservoirs

Replaces tanks with reservoirs.

This example contains:
- Choose Network.
- Load network and simulate.
- Find tank patterns.
- Create patterns for new reservoirs from tanks.
- Replace tanks with reservoirs.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 9 - Compare Simulations

Compares hydraulics and quality analysis functions.
 
This example contains:
- Load a network.
- Set simulation duration.
- Test hydraulics and quality analysis functions.
- Step by step hydraulic analysis.
- Step by step quality analysis.
- Unload library.
- Run time d.getComputedTimeSeries.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 10 - Close pipes during simulation

Closes pipes during simulation.

This example contains:
- Load a network.
- Link index for change the status.
- Run step by step hydraulic analysis.
- Get flows for the specific link index.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 11 - Assign a new curve to a pump

Assings a new curve to a specific pump.

This example contains:
- Load a network.
- Add new curve in the network.
- Get pump index.
- Get head curve index.
- Assing new curve index on the specific pump.
- Unload library.
  
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 12 - Add multiple control patterns

Defines pattern time for valves setting via bin functions.

This example contains:
- Add paths and load inp file.
- Get times.
- Create random settings.
- Settings for valve 'VALVE-173'.
- Check pressures at node after valve.
- Run simulation with executable.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 13a - Add pipe and junction via bin function

Adds CV pipe with the bin Function addBinJunction.

This example contains:
- Load a network.
- Plot newtork with bin functions.
- Add Junction + CV pipe.
- Plot normal with new pipe.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 13b - Add pipe via bin function

Adds CV pipe in a network with Bin Function addBinCVPipe.

This example contains:
- Load a network.
- Get from and to nodes for add cv pipe.
- Plot network.
- Plot network with changes.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 13c - Add pipe via normal function

Adds a CV Pipe in a network.

This example contains:
- Load a network.
- Plot network components.
- Get from and to nodes for add cv pipe.
- Plot new network with changes.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 14 - Hydraulic and Quality analysis

Runs the Hydraulic and Quality analysis of a network.
 
This example contains:
- Load a network.
- Hydraulic and Quality analysis STEP-BY-STEP.
- Display nodes pressures, links flows, nodes actual qualities, links actual qualities.
- Unload library. 

 &uparrow; [Back to top](#table-of-contents)
 <br />
  
#### Example 15 - Write MSX file
Writes MSX File e.g.Net2.msx.

This example contains:
- Load a network.
- Set Input Arguments:
  - Filename
  - Section Title
  - Section Options
  - Section Species
  - Section Coefficients
  - Section Terms
  - Section Pipes
  - Section Tanks
  - Section Sources
  - Section Quality Global
  - Section Quality
  - Section Parameters
  - Section Patterns
- Write MSX File.
- Load MSX File.
- Compute.
- Unload libraries.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 16 - Creating scenarios

Creates scenarios based on node count of inp file EPANET.

This example contains:
- Load a network.
- Create scenario parameters.
- Create Multiple Scenarios.
- Set quality type.
- Run Scenarios.
- Get computed quality time series.
- Plot Quality VS Time.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 17a - Add controls via bin functions
 
Defines pattern time for pipe status via bin & normal functions.
 
This example contains:
- Add paths and load inp file.
- Get times.
- Get Pipe infos.
- Find the max length of link pipe names.
- Times.
- Create random status for pipes every hour.
- Add all controls.
- Unload libraries.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 17b - Add controls via normal functions

Defines pattern time for pipe status.

This example contains:
- Add paths and load inp file.
- Get times.
- Get Pipe infos.
- Find the max length of link pipe names.
- Times.
- Create random status for pipes every hour.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 18 - Change status of pipe during simulation

Changes status randomly at pipes during simulation.

This example contains:
 Load a network.
 Get pipe count.
 Get pipe indices.
 Run step by step hydraulic analysis.
 Set status random 0/1 for pipes.
 Plot flows.
 Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 19 - Rotate network
 
Rotates EPANET Inp File.

This example contains:
- Load network and paths.
- Plot nework initial.
- Rotate degrees theta.
- Define the x- and y-data for the original line we would like to rotate.
- Create a matrix of these points, which will be useful in future calculations.
- Choose a point which will be the center of rotation.
- Create a matrix which will be used later in calculations.
- Define a 60 degree counter-clockwise rotation matrix.
- Do the rotation:
   - Shift points in the plane so that the center of rotation is at the origin.
   - Apply the rotation about the origin.
   - Shift again so the origin goes back to the desired center of rotation.
- Pick out the vectors of rotated x- and y-data.
- Plot rotated network.
- Save rotated Inp file.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 20a - Change of external controls
 
a) External Controls, Changes control status of pump STEP BY STEP.
 
This example contains:
- Load network.
- Delete Controls.
- Hydraulic analysis STEP-BY-STEP.
- CONTROLS.
- Add new controls in live.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 20b - Addition of external controls

b) External Controls, Adds controls in hydraulic analysis STEP-BY-STEP.

This example contains:
- Load network.
- Delete Controls.
- Hydraulic analysis STEP-BY-STEP.
- CONTROLS.
- Add new controls in live.
- Delete controls.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 21 - Pressure driven analysis

Tests the Pressure Driven Analysis option.

This example contains:
- Set Demand Multiplier to 10 to cause negative pressures.
- Run single period analysis.
- Solve hydraulics with default DDA option which will return with neg. pressure warning code.
- Check that 4 demand nodes have negative pressures.
- Switch to PDA with pressure limits of 20 - 100 psi.
- Solve hydraulics again.
- Check that 6 nodes had demand reductions totaling 32.66%.
- Check that Junction 12 had full demand.
- Check that Junction 21 had a demand deficit of 413.67.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)
 <br />
 
#### Example 22 - Overflow of tanks

Tests the overflow option for tanks.

This example contains:
- Get index of the tank and its inlet/outlet pipe.
- Set initial & maximum level to 130.
- Set duration to 1 hr.
- Solve hydraulics with default of no tank spillage allowed.
- Check that tank remains full.
- Check that there is no spillage.
- Check that inflow link is closed.
- Turn tank overflow option on.
- Solve hydraulics again.
- Check that tank remains full.
- Check that there is spillage equal to tank inflow
 (inflow has neg. sign since tank is start node of inflow pipe).
- Save project to file and then close it.
- Re-open saved file & run it.
- Check that tank spillage has same value as before.
- Unload library.
 
 &uparrow; [Back to top](#table-of-contents)

#### Example 23 - Change connection of links

- This example contains:
- Load network.
- Get all link node indices.
- Set link node indices for specific link index.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)

#### Example 24 - Delete all patterns

- This example contains:
- Load network.
- Delete all patterns.
- Save new file without patterns.
- Unload library.

 &uparrow; [Back to top](#table-of-contents)

