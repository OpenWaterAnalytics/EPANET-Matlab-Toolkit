% Fasted Parallel computations using the function getComputedTimeSeries_ENepanet
% You need to star the "Parallel Pool" (MATLAB down-left corner for drop-down)
clc; clear; close all; clear class; fclose all;

% Load epanet matlab paths
start_toolkit;

% Shut down and delete the current parallel pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

% Restart the parallel pool
parpool;

% Example 1
% Run the simulation and store the pressure results using parallel processing
rng(2022)
tic;
% Number of simulations
Nsim = 3000;
% 5% max uncertainty in base demands
eta_bar = 0.05;
% initialize matrix to save results
Pmcs = cell(Nsim,1);
inpname = 'Net2.inp';

Ginit = epanet('Net2.inp', 'loadfile-ph');

% Get nominal base demands
base_demands = Ginit.getNodeBaseDemands{1};

tmpinpfile = {};
tmprptfile = {};
tmpbinfile = {};
% Create INP scenarios
parfor i = 1:Nsim
    % Load EMT, data and functions
    G = Ginit;
    G.loadlibrary;
    G.loadEPANETFile(G.TempInpFile);
    % Compute new base demands
    delta_bd = (2*rand(1,length(base_demands))-1).*eta_bar.*base_demands;
    new_base_demands = base_demands + delta_bd;
    % Set base demands
    G.setNodeBaseDemands(new_base_demands);
    % Define names of temporary INP and BIN files
    tmpinpfile{i} = ['@#inp', num2str(i),'.inp'];
    tmpbinfile{i} = ['@#bin', num2str(i),'.bin'];   
    tmprptfile{i} = ['@#rpt', num2str(i),'.rpt'];   
    % Computed time series
    Pmcs{i} =  G.getComputedTimeSeries_ENepanet(tmpinpfile{i}, tmpbinfile{i}, tmprptfile{i});
    disp(['Epoch ',int2str(i)]);
    % Delete temporary files
    delete(tmpinpfile{i}, tmpbinfile{i}, tmprptfile{i});
end
toc;
Ginit.unload();

