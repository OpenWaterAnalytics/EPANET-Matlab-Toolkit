%% Write MSX File e.g.Net2.msx

%% 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net2.inp');

% Input Arguments
msx={};
msx.FILENAME = 'Net2tmp.msx';

% section Title
msx.TITLE = 'Example: Net2tmp MSX TEST';

% section Options
msx.AREA_UNITS = 'FT2'; %AREA_UNITS FT2/M2/CM2
msx.RATE_UNITS = 'DAY'; %TIME_UNITS SEC/MIN/HR/DAY
msx.SOLVER = 'EUL'; %SOLVER EUL/RK5/ROS2
msx.COUPLING = 'NONE'; %COUPLING FULL/NONE
msx.COMPILER = 'NONE'; %COMPILER NONE/VC/GC
msx.TIMESTEP = 300; %TIMESTEP in seconds
msx.ATOL = 0.001;  %ATOL value
msx.RTOL = 0.001;  %RTOL value

% section Species
% <type> <specieID> <units> (<atol> <rtol>)
msx.SPECIES={'BULK CL2 MG 0.01 0.001'}; %type [BULK/WALL] [specieID] [units UG/MG] [atol] [rtol]

% section Coefficients 
% CONSTANT name value % PARAMETER name value
msx.COEFFICIENTS = {'PARAMETER Kb 0.3', 'PARAMETER Kw 1'}; %[name] [value]

% section Terms
% <termID> <expression>
msx.TERMS = {'Kf 1.5826e-4 * RE^0.88 / D'}; % [termID] [expression]

% section Pipes
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.PIPES = {'RATE CL2 -Kb*CL2-(4/D)*Kw*Kf/(Kw+Kf)*CL2'}; % [type] [specieID] [expression]

% section Tanks
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.TANKS = {'RATE CL2 -Kb*CL2'}; % [type] [specieID] [expression]

% section Sources
% <type> <nodeID> <specieID> <strength> (<patternID>)
msx.SOURCES = {'CONC 1 CL2 0.8 '}; %[CONC/MASS/FLOW/SETPOINT] [nodeID] [specieID] [strength] [patternID]

% section Quality Global
% GLOBAL <specieID> <value>
msx.GLOBAL = {'Global CL2 0.5'}; % [specieID] [value]

% others
% NODE <nodeID> <bulkSpecieID> <value>
% LINK <linkID> <wallSpecieID> <value>
msx.QUALITY = {'NODE 26 CL2 0.1'}; %[NODE/LINK] [ID] [bulkSpecieID/wallSpecieID] [value]

% section Parameters
% PIPE <pipeID> <paramID> <value>
% TANK <tankID> <paramID> <value>
msx.PARAMETERS = {''};

% section Patterns
% <patternID> <multiplier> <multiplier> 
msx.PATERNS = {''}; % [patternID] [multiplier]     

% Write MSX File
d.writeMSXFile(msx);

% Load MSX File
d.loadMSXFile(msx.FILENAME)

% Compute
d.getMSXComputedQualityNode

% Unload libraries
d.unloadMSX;
d.unload;
