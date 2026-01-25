% Parallel computations
% You need to start the "Parallel Pool" (MATLAB down-left corner drop-down)
clc; clear; close all; clear classes;
start_toolkit;

%% First example
% Run the simulation and store the pressure results using parallel processing
rng(2022)
tic

% Number of simulations
Nsim = 100;

% 5% max uncertainty in base demands
eta_bar = 0.05;

% Initialize cell array to save MCS pressures
Pmcs = cell(Nsim, 1);

inpname = 'Net2';
dirscenarios = 'scenarios';

% Preallocate scenario filenames for parfor
destinp = cell(Nsim, 1);

try
    mkdir(dirscenarios);
    addpath(dirscenarios);
catch
end

parfor i = 1:Nsim
    % Create temporary file for each scenario
    destinp{i} = [dirscenarios, '/', inpname, '_', num2str(i), '.inp'];
    copyfile(which([inpname, '.inp']), destinp{i});

    % Load EMT, data, and functions
    Gi = epanet(destinp{i}, 'loadfile-ph');

    % Get nominal base demands
    base_demands = Gi.getNodeBaseDemands{1};

    % Compute new base demands
    delta_bd = (2 * rand(1, length(base_demands)) - 1) .* eta_bar .* base_demands;
    new_base_demands = base_demands + delta_bd;

    % Set base demands
    Gi.setNodeBaseDemands(new_base_demands);

    % Compute pressure and node quality
    Pmcs{i} = Gi.getComputedHydraulicTimeSeries('Pressure').Pressure;

    % Store object for unloading later
    G(i) = Gi;

    disp(['Epoch ', int2str(i)])

end

toc

tic

% Stop parallel pool
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

toc

%% Second example
tic

d = epanet('Net1.inp');
iterations = 100;

% Preallocate output cell array
H = cell(iterations, 1);

parfor i = 1:iterations
    di = epanet('Net1.inp');
    di.loadlibrary;
    di.loadEPANETFile(di.TempInpFile);

    % Set parameters
    elevations = di.getNodeElevations - di.getNodeElevations * rand(1) * 0.5;
    di.setNodeElevations(elevations * 10);

    % Compute hydraulics
    H{i} = di.getComputedHydraulicTimeSeries;

    di.unload;

end

d.unload;
toc
