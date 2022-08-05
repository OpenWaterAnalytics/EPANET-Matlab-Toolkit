%% MSX no inp, msx files test 
clear;clc;
try unloadlibrary epanetmsx; catch, end
try unloadlibrary legacymsx; catch, end

%% Load MSXCore
lib_core = [fileparts(which('64bit/epanetmsx.dll')), '\epanetmsx'];
h_file_core = which('64bit/epanetmsx.h');
loadlibrary(lib_core, h_file_core);

%% Load Legacy
lib_legacy = [fileparts(which('64bit/legacymsx.dll')), '\legacymsx'];
h_file_legacy = which('64bit/legacyToolkit.h');
loadlibrary(lib_legacy, h_file_legacy);

%% Open MSX project
MSX = libpointer('voidPtr');
Errcode = calllib('epanetmsx', 'MSX_open', MSX);
setdatatype(MSX, 'ProjectPtr')

%% Builing the network from example.inp
Errcode = calllib('epanetmsx', 'MSX_setFlowFlag', MSX, 8);
Errcode = calllib('epanetmsx', 'MSX_setTimeParameter', MSX, 0, 80*3600);
Errcode = calllib('epanetmsx', 'MSX_setTimeParameter', MSX, 1, 1*3600);
Errcode = calllib('epanetmsx', 'MSX_setTimeParameter', MSX, 2, 8*3600);
Errcode = calllib('epanetmsx', 'MSX_setTimeParameter', MSX, 5, 8*3600);
Errcode = calllib('epanetmsx', 'MSX_setTimeParameter', MSX, 6, 0);

%% Add nodes
Errcode = calllib('epanetmsx', 'MSX_addNode', MSX, 'a');
Errcode = calllib('epanetmsx', 'MSX_addNode', MSX, 'b');
Errcode = calllib('epanetmsx', 'MSX_addNode', MSX, 'c');
Errcode = calllib('epanetmsx', 'MSX_addNode', MSX, 'e');
Errcode = calllib('epanetmsx', 'MSX_addReservoir', MSX, 'source', 0, 0, 0);

%% Add links
Errcode = calllib('epanetmsx', 'MSX_addLink', MSX, '1', 'source', 'a', 1000, 200, 100);
Errcode = calllib('epanetmsx', 'MSX_addLink', MSX, '2', 'a', 'b', 800, 150, 100);
Errcode = calllib('epanetmsx', 'MSX_addLink', MSX, '3', 'a', 'c', 1200, 200, 100);
Errcode = calllib('epanetmsx', 'MSX_addLink', MSX, '4', 'b', 'c', 1000, 150, 100);
Errcode = calllib('epanetmsx', 'MSX_addLink', MSX, '5', 'c', 'e', 2000, 150, 100);

%% Add Options
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 0, 'M2');
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 1, 'HR');
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 2, 'RK5');
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 4, '28800');
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 5, '0.001');
Errcode = calllib('epanetmsx', 'MSX_addOption', MSX, 6, '0.0001');

%% Add Species
Errcode = calllib('epanetmsx', 'MSX_addSpecies', MSX, 'AS3', 0, 1, 0.0, 0.0);
Errcode = calllib('epanetmsx', 'MSX_addSpecies', MSX, 'AS5', 0, 1, 0.0, 0.0);
Errcode = calllib('epanetmsx', 'MSX_addSpecies', MSX, 'AStot', 0, 1, 0.0, 0.0);
Errcode = calllib('epanetmsx', 'MSX_addSpecies', MSX, 'AS5s', 1, 1, 0.0, 0.0);
Errcode = calllib('epanetmsx', 'MSX_addSpecies', MSX, 'NH2CL', 0, 0, 0.0, 0.0);

%% Add Coefficents
Errcode = calllib('epanetmsx', 'MSX_addCoefficeint', MSX, 6, 'Ka', 10.0);
Errcode = calllib('epanetmsx', 'MSX_addCoefficeint', MSX, 6, 'Kb', 0.1);
Errcode = calllib('epanetmsx', 'MSX_addCoefficeint', MSX, 6, 'K1', 5.0);
Errcode = calllib('epanetmsx', 'MSX_addCoefficeint', MSX, 6, 'K2', 1.0);
Errcode = calllib('epanetmsx', 'MSX_addCoefficeint', MSX, 6, 'Smax', 50);

%% Add terms
Errcode = calllib('epanetmsx', 'MSX_addTerm', MSX, 'Ks', 'K1/K2');

%% Add Expressions
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 1, 1, 'AS3', '-Ka*AS3*NH2CL');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 1, 1, 'AS5', 'Ka*AS3*NH2CL-Av*(K1*(Smax-AS5s)*AS5-K2*AS5s)');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 1, 1, 'NH2CL', '-Kb*NH2CL');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 1, 3, 'AS5s', 'Ks*Smax*AS5/(1+Ks*AS5)-AS5s');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 1, 2, 'AStot', 'AS3 + AS5');

Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 2, 1, 'AS3', '-Ka*AS3*NH2CL');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 2, 1, 'AS5', 'Ka*AS3*NH2CL');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 2, 1, 'NH2CL', '-Kb*NH2CL');
Errcode = calllib('epanetmsx', 'MSX_addExpression', MSX, 2, 2, 'AStot', 'AS3+AS5');

%% Add Quality
Errcode = calllib('epanetmsx', 'MSX_addQuality', MSX, 'NODE', 'AS3', 10.0, 'source');
Errcode = calllib('epanetmsx', 'MSX_addQuality', MSX, 'NODE', 'NH2CL', 2.5, 'source');

%% Setup Report
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'NODE', 'c', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'NODE', 'e', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'LINK', '4', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'LINK', '5', 2);   
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'SPECIE', 'AStot', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'SPECIE', 'AS3', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'SPECIE', 'AS5', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'SPECIE', 'AS5s', 2);
Errcode = calllib('epanetmsx', 'MSX_setReport', MSX, 'SPECIE', 'NH2CL', 2);

%% Finish Setup
Errcode = calllib('epanetmsx', 'MSX_init', MSX);

%% Run
demands = [0.040220, 0.033353, 0.053953, 0.022562, -0.150088];
heads = [327.371979, 327.172974, 327.164185, 326.991211, 328.083984];
flows =  [0.150088, 0.039916, 0.069952, 0.006563, 0.022562];
Errcode = calllib('epanetmsx', 'MSX_setHydraulics', MSX, demands, heads, flows);

t = 0;
tleft = 1;
oldHour = -1;
newHour = 0;
count = 1;
[Errcode, ~, node_count] = calllib('epanetmsx', 'MSX_getcount', MSX, 0, 0);
while (tleft>0 && Errcode == 0)
    
    [Errcode] = calllib('legacymsx', 'MSXsaveResults', MSX);
    [Errcode, ~, t, tleft] = calllib('epanetmsx', 'MSX_step', MSX, t, tleft);
    
    % Get Quality
    TYPE = 0; SPECIES = 2; % 'AS5'
    for i=1:node_count
        [Errcode, ~, qual(count,i)] = calllib('epanetmsx', 'MSX_getQualityByIndex', MSX, TYPE, i, SPECIES, 0);
    end
    time(count) = t;
    count = count+1;
end

Errcode = calllib('legacymsx', 'MSXsaveFinalResults', MSX);
Errcode = calllib('legacymsx', 'MSXreport', MSX, 'example_report');

%% Close - Crashes when running MSX_close. 
% Errcode = calllib('epanetmsx', 'MSX_close', MSX);
