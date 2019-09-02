%% Writes MSX File e.g.Net2.msx.
% This example contains:
%   Load a network.
%   Set Input Arguments:
%       Filename
%       Section Title
%       Section Options
%       Section Species
%       Section Coefficients
%       Section Terms
%       Section Pipes
%       Section Tanks
%       Section Sources
%       Section Quality Global
%       Section Quality
%       Section Parameters
%       Section Patterns
%   Write MSX File.
%   Load MSX File.
%   Compute.
%   Unload libraries.

%% 
% Clear
clear; close('all'); clc;
start_toolkit;

% Load a network.
d = epanet('Net2.inp');

% Set Input Arguments:
% Filename
msx={};
msx.FILENAME = 'Net2tmp.msx';

% Section Title
msx.TITLE = 'Example: Net2tmp MSX TEST';

% Section Options
msx.AREA_UNITS = 'FT2'; %AREA_UNITS FT2/M2/CM2
msx.RATE_UNITS = 'DAY'; %TIME_UNITS SEC/MIN/HR/DAY
msx.SOLVER = 'EUL'; %SOLVER EUL/RK5/ROS2
msx.COUPLING = 'NONE'; %COUPLING FULL/NONE
msx.COMPILER = 'NONE'; %COMPILER NONE/VC/GC
msx.TIMESTEP = 300; %TIMESTEP in seconds
msx.ATOL = 0.001;  %ATOL value
msx.RTOL = 0.001;  %RTOL value

% Section Species
% <type> <specieID> <units> (<atol> <rtol>)
msx.SPECIES={'BULK CL2 MG 0.01 0.001'}; %type [BULK/WALL] [specieID] [units UG/MG] [atol] [rtol]

% Section Coefficients
% CONSTANT name value % PARAMETER name value
msx.COEFFICIENTS = {'PARAMETER Kb 0.3', 'PARAMETER Kw 1'}; %[name] [value]

% Section Terms
% <termID> <expression>
msx.TERMS = {'Kf 1.5826e-4 * RE^0.88 / D'}; % [termID] [expression]

% Section Pipes
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.PIPES = {'RATE CL2 -Kb*CL2-(4/D)*Kw*Kf/(Kw+Kf)*CL2'}; % [type] [specieID] [expression]

% Section Tanks
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.TANKS = {'RATE CL2 -Kb*CL2'}; % [type] [specieID] [expression]

% Section Sources
% <type> <nodeID> <specieID> <strength> (<patternID>)
msx.SOURCES = {'CONC 1 CL2 0.8 '}; %[CONC/MASS/FLOW/SETPOINT] [nodeID] [specieID] [strength] [patternID]

% Section Quality Global
% GLOBAL <specieID> <value>
msx.GLOBAL = {'Global CL2 0.5'}; % [specieID] [value]

% Section Quality
% NODE <nodeID> <bulkSpecieID> <value>
% LINK <linkID> <wallSpecieID> <value>
msx.QUALITY = {'NODE 26 CL2 0.1'}; %[NODE/LINK] [ID] [bulkSpecieID/wallSpecieID] [value]

% Section Parameters
% PIPE <pipeID> <paramID> <value>
% TANK <tankID> <paramID> <value>
msx.PARAMETERS = {''};

% Section Patterns
% <patternID> <multiplier> <multiplier> 
msx.PATERNS = {''}; % [patternID] [multiplier]     

% Write MSX File.
d.writeMSXFile(msx);

% Load MSX File.
d.loadMSXFile(msx.FILENAME)

% Compute.
d.getMSXComputedQualityNode

% Unload libraries.
d.unloadMSX;
d.unload;
