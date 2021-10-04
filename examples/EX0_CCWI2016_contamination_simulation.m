%%  Case Study for CCWI2016
% Below, we provide a template solution on how to design 
% an Arsenite contamination simulator based on EPANET andEPANET-MSX using the EPANET-MATLAB Toolkit. 
% More about example: https://zenodo.org/record/831493#.XvnadWgzaUl
% https://github.com/KIOS-Research/CCWI2016

close all; clear; close('all'); clc; rng(1)

%% Load EPANET Network and MSX
G = epanet('BWSN_Network_1.inp'); % Load EPANET Input file
G.loadMSXFile('Arsenite.msx'); % Load MSX file

% Sensor locations
sensor_id = {'JUNCTION-17', 'JUNCTION-83', 'JUNCTION-122', 'JUNCTION-31', 'JUNCTION-45'};
sensor_index = G.getNodeIndex(sensor_id);

%% Simulation Setup
t_d = 5; % days
G.setTimeSimulationDuration(t_d*24*60*60); % Set simulation duration of 5 days

%% Get Network data
demand_pattern = G.getPattern;
roughness_coeff = G.getLinkRoughnessCoeff;
node_id = G.getNodeNameID;

G.setMSXTimeStep(3600)

%% Scenarios
Ns = 100; % Number of scenarios to simulate
u_p = 0.20; % pattern uncertainty
u_r = 0.20; % roughness coefficient uncertainty
max_inj_conc = 2.0;
inj_start_time = 2*48; % after day 2 (Dt = 30min)
inj_duration = 24; % maximum duration of 12 hours
inj_sc=[randi(G.NodeCount,Ns,1), max_inj_conc*rand(Ns,1), randi(48,Ns,1)+inj_start_time, randi(inj_duration,Ns,1)]; % Injection location, magnitude, start time, duration 

%% Run epochs
tic
for i = 1:Ns
    disp(['Iteration ', int2str(i)])
    % Randomize hydraulics
    r_p = -u_p + 2*u_p.*rand(size(demand_pattern,1),size(demand_pattern,2));
    new_demand_pattern = demand_pattern + demand_pattern.*r_p;
    G.setPatternMatrix(new_demand_pattern); % Set new patterns
    r_r = -u_r + 2*u_r.*rand(size(roughness_coeff,1),size(roughness_coeff,2));
    new_roughness_coeff = roughness_coeff + roughness_coeff.*r_r;
    G.setLinkRoughnessCoeff(new_roughness_coeff); % Set new roughness coefficients
    % Randomize injection
    G.setMSXSources(node_id(inj_sc(i,1)), 'AsIII', 'Setpoint', inj_sc(i,2), 'AS3PAT'); % Specify Arsenite injection source
    as3_pat = zeros(1, t_d*48);
    as3_pat(inj_sc(i,3):(inj_sc(i,3)+inj_sc(i,4))) = 1; % 
    G.setMSXPattern('AS3PAT',as3_pat); % Set pattern for injection
    Q{i} = G.getMSXComputedQualityNode(sensor_index); % Solve hydraulics and MSX quality dynamics
    G.setMSXSources(node_id(inj_sc(i,1)), 'AsIII', 'Setpoint', 0, 'AS3PAT'); % Reset injection source
end
toc
G.unloadMSX

%% Plot results
for i = 1:Ns
    for j = 1:length(sensor_index)
       subplot(5,1,j)
       plot(Q{i}.Time/24/60/60, Q{i}.Quality{j}(:,1),'-','Color',[0,0.7,0.9]); hold on; grid on
    end
end
for i = 1:length(sensor_index)
   subplot(5,1,i)
   title(sensor_id{i})
   ylabel('Cl_2 (mg/L)')
   xlabel('Time (days)')
end
