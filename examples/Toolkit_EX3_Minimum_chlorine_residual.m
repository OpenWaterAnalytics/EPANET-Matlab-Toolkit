%% This example illustrates how the Toolkit could be used to determine the lowest dose of 
% chlorine applied at the entrance to a distribution system needed to ensure that a minimum 
% residual is met throughout the system. We assume that the EPANET input file contains the proper 
% set of kinetic coefficients that describe the rate at which chlorine will decay in the system being 
% studied. In the example code, the ID label of the source node is contained in SourceID, the minimum 
% residual target is given by Ctarget, and the target is only checked after a start-up duration of 5 
% days (432,000 seconds). To keep the code more readable, no error checking is made on the results 
% returned from the Toolkit function calls.

%% https://github.com/OpenWaterAnalytics/EPANET/wiki/Example-3

%Clear 
clear; close('all'); clc;
start_toolkit;

% Load a network
d = epanet('Net1.inp');

% Set Ctarget
Ctarget = 0.5;

% Source node id
SourceID = '2';

% Set simulation duration 6 days
d.setTimeSimulationDuration(6*24*3600);

% Obtain a hydraulic solution 
d.solveCompleteHydraulics;

% Get the number of nodes
nnodes = d.getNodeCount;

% Get source node's index
sourceindex = d.getNodeIndex(SourceID);

% Setup system to analyze for chlorine (in case it was not done in the input file.)
d.setQualityType('Chlorine', 'mg/L', '');

% Open the water quality solver
d.openQualityAnalysis;

% Begin the search for the source concentration
csource = 0.0;
violation = 0;
while (~violation && (csource <= 4.0))

    % Update source concentration to next level
    csource = csource + 0.1; 
    d.setNodeSourceQuality(sourceindex, csource);
    
    % Run WQ simulation checking for target violations
    d.initializeQualityAnalysis;
    tstep = 1;
    while (~violation && (tstep > 0))
      
        t = d.runQualityAnalysis;
        if (t > 432000)
            for i=1:nnodes
                c = d.getNodeActualQuality(i);
                if (c < Ctarget) 
                    violation = 1;
                    break;
                end
                
            end
        end
        tstep = d.nextQualityAnalysisStep;
    end
end

csource

% Unload library
d.unload

