% Parallel computations
% You need to star the "Parallel Pool" (MATLAB down-left corner for drop-down)
clc; clear; close all; clear class;
start_toolkit;
%% First example
% Run the simulation and store the pressure results using parallel processing
rng(2022)
tic;
% Number of simulations
Nsim = 100;
% 5% max uncertainty in base demands
eta_bar = 0.05;
% initialize matrix to save MCS pressures
Pmcs = cell(Nsim,1);
inpname = 'Net2';
clear destinp;
dirscenarios = 'scenarios';
try
    mkdir(dirscenarios);
    addpath(dirscenarios);
catch
end
parfor i = 1:Nsim
    % Create temporary file for each scenario
    destinp{i} = [dirscenarios, '/', inpname, '_', num2str(i), '.inp'];
    copyfile(which([inpname, '.inp']), destinp{i});  
    % Load EMT, data and functions
    G(i) = epanet(destinp{i}, 'loadfile-ph');
    % Get nominal base demands
    base_demands = G(i).getNodeBaseDemands{1};
    % Compute new base demands
    delta_bd = (2*rand(1,length(base_demands))-1).*eta_bar.*base_demands;
    new_base_demands = base_demands + delta_bd;
    % Set base demands
    G(i).setNodeBaseDemands(new_base_demands);
    % Compute pressure and node quality
    Pmcs{i} = G(i).getComputedHydraulicTimeSeries('Pressure').Pressure;
    disp(['Epoch ',int2str(i)])
end
toc;
tic;
% Stop parallel
delete(gcp('nocreate'));
% Unload networks
for i = 1:Nsim
    G(i).unload();
end
try 
    % Remove all scenarios
    rmpath(dirscenarios)
    rmdir(dirscenarios, 's')
catch
end
toc;

%% Second Example
tic
d = epanet('Net1.inp'); 
iterations = 100;

parfor i = 1:iterations
    d.loadlibrary;
    d.loadEPANETFile(d.TempInpFile); 

    % set parameters
    elevations = d.getNodeElevations-d.getNodeElevations*rand(1)*.5;
    d.setNodeElevations(elevations*10);
    
    % Computed Hydraulics
    H{i} = d.getComputedHydraulicTimeSeries;
end
d.unload;
toc



