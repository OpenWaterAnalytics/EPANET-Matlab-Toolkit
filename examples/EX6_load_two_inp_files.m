%% Load 2 different networks
%Clear 
clear; close('all'); clc;
start_toolkit;

% Using 2 dll files. Load 2 Input files.
arch = computer('arch');
pwdepanet = fileparts(which('epanet.m'));
if strcmpi(arch,'win64')
    LibEPANETpath = [pwdepanet,'/64bit/'];
elseif strcmpi(arch,'win32')
    LibEPANETpath = [pwdepanet,'/32bit/'];
end        
temp_lib_folder = [LibEPANETpath,'temp_lib/'];
mkdir(temp_lib_folder)
epanet2_tmp = [temp_lib_folder,'epanet2tmp.dll'];
dll_pathLib = [LibEPANETpath, 'epanet2.dll'];
header_pathLib = [LibEPANETpath, 'epanet2.h'];
try
    copyfile(dll_pathLib, epanet2_tmp);
    copyfile(header_pathLib, [temp_lib_folder,'epanet2tmp.h']);
catch
end

% Load networks
d1 = epanet('Net1.inp', epanet2_tmp);
warning('off'); % ignore warnings
d2 = epanet('Net2.inp');
warning('on');

% Disp elevations for the two networks
disp('Net1 - Elevations:')
disp('------------------')
disp(d1.getNodeElevations)
disp('Net2 - Elevations:')
disp('------------------')
disp(d2.getNodeElevations)

% Plot networks
d1.plot;
d2.plot;

% Unload libraries
d1.unload;
d2.unload;

% Delete files who created
delete(epanet2_tmp);
delete([temp_lib_folder,'epanet2tmp.h']);
rmdir(temp_lib_folder);