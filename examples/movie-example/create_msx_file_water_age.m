function msx = create_msx_file_water_age(input)
%CREATE_MSX_FILE_WATER_AGE - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output] = create_msx_file_water_age(input)
%
% Inputs:
%    input1 - Description
%
% Outputs:
%    output1 - Description
%
% Example: 
%    Line 1 of example
%
% Other m-files required: simulate_chlorine_residuals
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author        : Demetrios G. Eliades, Marios Kyriakou
% Work address  : KIOS Research Center, University of Cyprus
% email         : eldemet@ucy.ac.cy
% Website       : http://www.kios.ucy.ac.cy
% Last revision : September 2016

%------------- BEGIN CODE --------------

% Input Arguments
msx={};
msx.FILENAME = 'water_age.msx';
% section Title
msx.TITLE{1} = 'Example: Water Age.';
% section Options
msx.AREA_UNITS='FT2'; %AREA_UNITS FT2/M2/CM2
msx.TIME_UNITS='DAY'; %TIME_UNITS SEC/MIN/HR/DAY
msx.SOLVER='EUL'; %SOLVER EUL/RK5/ROS2
msx.COUPLING='NONE'; %COUPLING FULL/NONE
msx.COMPILER='NONE'; %COMPILER NONE/VC/GC
msx.TIMESTEP=3600; %TIMESTEP in seconds
msx.ATOL=0.01;  %ATOL value
msx.RTOL=0.001;  %RTOL value
% section Species
% <type> <specieID> <units> (<atol> <rtol>)
msx.SPECIES{1}={'BULK AGE HR 0.01 0.001'}; %rtol

% section Pipes
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.PIPES{1} ={'RATE AGE 24'}; %expression

% section Tanks
% EQUIL <specieID> <expression>
% RATE <specieID> <expression>
% FORMULA <specieID> <expression>
msx.TANKS{1} ={'RATE AGE 24'}; %expression

%------------- END OF CODE --------------

  

