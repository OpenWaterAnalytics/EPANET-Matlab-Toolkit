%% Write MSX File

%% 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Input Arguments
msx={};
msx.msxFile = 'Net1.msx';
% section Title
msx.titleDescription{1} = 'Example: Net1 MSX TEST';
% section Options
msx.options{1}='FT2'; %AREA_UNITS FT2/M2/CM2
msx.options{2}='DAY'; %TIME_UNITS SEC/MIN/HR/DAY
msx.options{3}='EUL'; %SOLVER EUL/RK5/ROS2
msx.options{4}='NONE'; %COUPLING FULL/NONE
msx.options{5}='NONE'; %COMPILER NONE/VC/GC
msx.options{6}=3600; %TIMESTEP in seconds
msx.options{7}=0.01;  %ATOL value
msx.options{8}=0.001;  %RTOL value
% section Species
% <type> <specieID> <units> (<atol> <rtol>)
msx.species{1}={'BULK'}; %type BULK/WALL
msx.species{2}={'CL2'}; %specieID
msx.species{3}={'MG'}; %units UG/MG
msx.species{4}={0.01}; %atol
msx.species{5}={0.001}; %rtol

% section Coefficients 
% CONSTANT name value % PARAMETER name value
msx.coefficients{1}={'PARAMETER', 'PARAMETER'}; 
msx.coefficients{2}={'Kb', 'Kw'}; 
msx.coefficients{3}={0.3, 1}; 

% section Terms
% <termID> <expression>
msx.terms{1}={'Kf'}; % termID
msx.terms{2}={'1.5826e-4 * RE^0.88 / D'}; % expression

% section Pipes
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.pipes{1} ={'RATE'}; %type
msx.pipes{2} ={'CL2'}; %specieID
msx.pipes{3} ={'-Kb*CL2 - (4/D)*Kw*Kf/(Kw+Kf)*CL2'}; %expression

% section Tanks
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.tanks{1} ={'RATE'}; %type
msx.tanks{2} ={'CL2'}; %specieID
msx.tanks{3} ={'-Kb*CL2'}; %expression

% section Sources
% <type> <nodeID> <specieID> <strength> (<patternID>)
msx.sources{1}={''}; %CONC/MASS/FLOW/SETPOINT
msx.sources{2}={''}; %nodeID
msx.sources{3}={''}; %specieID
msx.sources{4}={''}; %strength
msx.sources{5}={''}; %patternID

% section Quality Global
% GLOBAL <specieID> <value>
msx.global{1} = {''};
msx.global{2} = {''};%specieID
msx.global{3} = {''};%value
% others
% NODE <nodeID> <bulkSpecieID> <value>
% LINK <linkID> <wallSpecieID> <value>
msx.quality{1} = {''}; %NODE/LINK
msx.quality{2} = {''}; %ID
msx.quality{3} = {''}; %bulkSpecieID/wallSpecieID
msx.quality{4} = {''}; %value

% section Parameters
% PIPE <pipeID> <paramID> <value>
% TANK <tankID> <paramID> <value>
msx.parameters{1} = {''};
msx.parameters{2} = {''};
msx.parameters{3} = {''};
msx.parameters{4} = {''};

% section Patterns
% <patternID> <multiplier> <multiplier> 
msx.patterns{1} = {''}; %patternID
msx.patterns{2} = {''}; %multiplier     

% Unload libraries
d.writeMSXFile(msx);
d.unload;
