function [V,L,T] = getQualityData(bulkSpecieID,wallSpecieID,msxFname,d)
%% Synopsis
%   [V,L,T,errcode] =
%   getQualityData(bulkSpecieID,wallSpecieID,msxFname,inpFname) Returns
%   vertex (V) and link (L) data from running the Epanet or Epanet-MSX
%   simulation defined by the input files msxFname and inpFname. For MSX
%   simulations the data in V is data defined on all network vertices
%   (nodes) for the bulk specie ID bulkSpecieID, and L is data defined on
%   all network links (pipes) for the wall specie ID wallSpecieID. For an
%   ordinary Epanet simulation, msxFname would be a blank character string
%   (msxFname='') or an empty matrix (msxFname=[]).  In this case both
%   bulkSpecieID and wallSpecieID are ignored, since Epanet does not know
%   about bulk or wall species; accordingly, for an ordinary Epanet
%   simulation, the return value of V contains the single-specie simulation
%   water quality data as defined in inpFname, while L=[].  For both MSX
%   and ordinary Epanet simulations, the ordering in V and/or L is node
%   and/or link index order. The output vector T contains the simulation
%   times (in seconds) associated with each frame.
%
%   If defined, the output matrices V, L, and T are (nnodex x nframes),
%   (nlinks x nframes), and (nframes x 1) respectively.  The value of
%   nframes is determined differently for Epanet and Epanet-MSX
%   simulations:
%
%   Epanet simulations: nframes is defined by the values of REPORT START
%   and REPORT STEP in the Epanet input file.  Specifically, frames are
%   stored only when Epanet simulation time t satisfies: (t >= REPORT_START
%   && mod(t,REPORT_STEP) == 0).
%
%   Epanet-MSX simulations: nframes is defined by the values of REPORT
%   START in the Epanet input file the value of TIMESTEP in the Epanet-MSX
%   input file.  Specifically, frames are stored only when Epanet-MSX
%   simulation time t satisfies: (t >= REPORT_START && mod(t,TIMESTEP) ==
%   0).
% Jim Uber
%
% modified by Marios Kyriakou 25/09/2016

%% Hydraulic and Water Quality Analyses
T = [];         % The output times
V = [];         % The vertex data
L = [];         % The link data
if ~isempty(msxFname)
    % MSX simulation
    % Get the simulation time parameters
    reportStart = d.getTimeReportingStart;

    % Hydraulic analysis
    d.solveMSXCompleteHydraulics;

    % Specie information
    bulkSpecieIndex=0;
    wallSpecieIndex=0;
    species_types = d.getMSXSpeciesType;
    if ~isempty(bulkSpecieID)
        bulkSpecieIndex = d.getMSXSpeciesIndex(bulkSpecieID);
        bulkSpecieType = species_types(bulkSpecieIndex);
        if ~sum(strcmp(bulkSpecieType,species_types))
            disp 'Specified bulk specie ID is incorrect type'
            return
        end
    end
    if ~isempty(wallSpecieID)
        wallSpecieIndex = d.getMSXSpeciesIndex(wallSpecieID);
        wallSpecieType = species_types(wallSpecieIndex);
        if ~sum(strcmp(wallSpecieType,species_types))
            disp 'Specified wall specie ID is incorrect type'
            return
        end
    end

    tleft = 1;      % Time left in simulation
    iframe = 1;     % frame index
    d.initializeMSXQualityAnalysis(1);
    while (tleft > 0)
        [t, tleft]=d.stepMSXQualityAnalysisTimeLeft;
        if t >= reportStart
            % Retrieve results at time t
            if bulkSpecieIndex
                V = [V zeros(d.NodeCount,1)];
                for in=1:d.NodeCount
                    V(in,iframe) = d.getMSXSpeciesConcentration(d.MSXTYPENODE,in,bulkSpecieIndex);
                end
            end
            if wallSpecieIndex
                L = [L zeros(d.LinkCount,1)];
                for in=1:d.LinkCount
                    L(in,iframe) = d.getMSXSpeciesConcentration(d.MSXTYPELINK,in,wallSpecieIndex);
                end
            end
            T = [T; t];
            iframe = iframe + 1;
        end
    end
else
    % Ordinary Epanet simulation
    % Get the simulation time parameters
    reportStart = d.getTimeReportingStart;
    reportStep = d.getTimeReportingStep;

    % Hydraulic analysis
    d.solveCompleteHydraulics;
    
    tleft = 1;      % Time left in simulation
    iframe = 1;     % frame index
    d.openQualityAnalysis;
    d.initializeQualityAnalysis;
    while (tleft > 0)
        t = d.runQualityAnalysis();
        if t >= reportStart && ~mod(t,reportStep)
            % Retrieve results at time t
            V = [V zeros(d.NodeCount,1)];
            for in=1:d.NodeCount
                V(in,iframe) = d.getNodeActualQuality(in);
            end
            T = [T; t];
            iframe = iframe + 1;
        end
        tleft = d.nextQualityAnalysisStep;
    end
    d.closeQualityAnalysis;
end


