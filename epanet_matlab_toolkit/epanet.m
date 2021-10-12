classdef epanet <handle
    %EPANET-Matlab Toolkit version dev2.2: A Matlab Class for EPANET and EPANET-MSX
    %libraries
    %
    %
    %   How to run:
    %   d=epanet('Net1.inp');
    %
    %   To select a different DLL version:
    %   d=epanet('Net1.inp', 'epanet2');
    %
    %   EPANET is software that models water distribution piping systems
    %   developed by the US EPA and provided under a public domain licence.
    %   This Matlab Class serves as an interface between Matlab and
    %   EPANET/EPANET-MSX, to assist researchers and the industry when
    %   solving problems related with water distribution systems.
    %
    %   EPANET and EPANET-MSX were developed by the Water Supply and Water
    %   Resources Division of the U.S. Environmental Protection Agency's
    %   National Risk Management Research Laboratory. EPANET is under the
    %   Public Domain and EPANET-MSX under GNU LGPL.
    %
    %   The latest EPANET files can downloaded at:
    %   https://github.com/OpenWaterAnalytics/epanet
    %
    %   Inspired by:
    %   EPANET-Matlab Wrappers (Jim Uber)
    %   EPANET-Matlab Toolkit (Demetrios Eliades)
    %   getwdsdata (Philip Jonkergouw)
    %
    %   EPANET-Matlab Class Licence:
    %
    %   Copyright 2016 KIOS Research Center for Intelligent Systems and
    %   Networks, University of Cyprus (www.kios.org.cy)
    %
    %   Licensed under the EUPL, Version 1.1 or - as soon they will be
    %   approved by the European Commission - subsequent libepanets of the
    %   EUPL (the "Licence"); You may not use this work except in
    %   compliance with the Licence. You may obtain a copy of the Licence
    %   at:
    %
    %   http://ec.europa.eu/idabc/eupl
    %
    %   Unless required by applicable law or agreed to in writing, software
    %   distributed under the Licence is distributed on an "AS IS" basis,
    %   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
    %   implied. See the Licence for the specific language governing
    %   permissions and limitations under the Licence.
    properties
        ControlLevelValues;          % The control level values
        ControlLinkIndex;            % Set of control links index
        ControlNodeIndex;            % Set of control nodes index
        ControlRules;                % Retrieves the parameters of all control statements
        ControlRulesCount;           % Number of controls
        Controls;                    % Controls info
        ControlSettings;             % Settings for the controls
        ControlTypes;                % Set of control types
        ControlTypesIndex;           % Index of the control types
        CurveCount;                  % Number of curves
        CurveIndex;                  % Index of curves
        CurvesInfo;                  % Curves info
        DemandModelCode;             % Demand model code DDA - 0, PDA - 1
        DemandModelPmin;             % Demand model Pmin - Pressure below which there is no demand
        DemandModelPreq;             % Demand model Preq - Pressure required to deliver full demand
        DemandModelPexp;             % Demand model Pexp - Pressure exponent in demand function
        DemandModelType;             % Demand model type DDA, PDA
        EnergyEfficiencyUnits;       % Units for efficiency
        EnergyUnits;                 % Units for energy
        Errcode;                     % Code for the EPANET error message
        InputFile;                   % Name of the input file
        Iterations;                  % Iterations to reach solution
        LibEPANET;                   % EPANET library dll
        LibEPANETpath;               % EPANET library dll path
        libFunctions;                % EPANET functions in dll
        LinkBulkReactionCoeff;       % Bulk reaction coefficient of each link
        LinkCount;                   % Number of links
        LinkDiameter;                % Diameter of each link
        LinkFlowUnits;               % Units of flow
        LinkFrictionFactorUnits;     % Units for friction factor
        LinkIndex;                   % Index of links
        LinkInitialSetting;          % Initial settings of links
        LinkInitialStatus;           % Initial status of links
        LinkLength;                  % Length of links
        LinkLengthsUnits;            % Units of length
        LinkMinorLossCoeff;          % Minor loss coefficient of links
        LinkMinorLossCoeffUnits;     % Minor loss coefficient units
        LinkNameID;                  % Name ID of links
        LinkPipeCount;               % Number of pipes
        LinkPipeDiameterUnits;       % Units for pipe diameters
        LinkPipeIndex;               % Index of pipe links
        LinkPipeNameID;              % Name ID of pipe links
        LinkPipeRoughnessCoeffUnits; % Pipe roughness coefficient units
        LinkPumpCount;               % Number of pumps
        LinkPumpHeadCurveIndex;      % Head curve indices
        LinkPumpIndex;               % Index of pumps
        LinkPumpNameID;              % Name ID of pumps
        LinkPumpPatternIndex;        % Index of pump pattern
        LinkPumpPatternNameID;       % ID of pump pattern
        LinkPumpPower;               % Power value
        LinkPumpPowerUnits;          % Units of power
        LinkPumpType;                % Pump type e.g constant horsepower, power function, user-defined custom curv
        LinkPumpTypeCode;            % Pump index/code
        LinkRoughnessCoeff;          % Roughness coefficient of links
        LinkType;                    % ID of link type
        LinkTypeIndex;               % Index of link type
        LinkValveCount;              % Number of valves
        LinkValveIndex;              % Index of valves
        LinkValveNameID;             % ID name of valves
        LinkVelocityUnits;           % Units for velocity
        LinkWallReactionCoeff;       % Wall reaction coefficient of links
        NodeBaseDemands;             % Base demands of nodes
        NodeCoordinates;             % Coordinates for each node (long/lat & intermediate pipe coordinates)
        NodeCount;                   % Number of nodes
        NodeDemandPatternIndex;      % Index of demand patterns
        NodeDemandPatternNameID;     % ID of demand patterns
        NodeDemandUnits;             % Units for demand
        NodeElevations;              % Elevation of nodes
        NodeElevationUnits;          % Units for elevation
        NodeEmitterCoeff;            % Emmitter Coefficient of nodes
        NodeEmitterCoefficientUnits; % Units for emitter coefficient
        NodeHeadUnits;               % Nodal head units
        NodeIndex;                   % Index of nodes
        NodeInitialQuality;          % Initial quality of nodes
        NodeJunctionCount;           % Number of junctions
        NodeJunctionIndex;           % Index of node junctions
        NodeJunctionNameID;          % Name ID of node junctions
        NodeNameID;                  % Name ID of all nodes
        NodeDemandCategoriesNumber;  % Number of demand categories for nodes
        NodePatternIndex;            % Node demand pattern indices
        NodePressureUnits;           % Units for Pressure
        NodeReservoirCount;          % Number of reservoirs
        NodeReservoirIndex;          % Index of reservoirs
        NodeReservoirNameID;         % Name ID of reservoirs
        NodesConnectingLinksID;      % Name IDs of nodes which connect links
        NodesConnectingLinksIndex;   % Indices of nodes which connect links
        NodeSourcePatternIndex;      % Index of pattern for node sources
        NodeSourceQuality;           % Quality of node sources
        NodeSourceTypeIndex;         % Index of source type
        NodeTankBulkReactionCoeff;   % Bulk reaction coefficients in tanks
        NodeTankCount;               % Number of tanks
        NodeTankDiameter;            % Diameters of tanks
        NodeTankDiameterUnits;       % Units for tank diameters
        NodeTankIndex;               % Indices of Tanks
        NodeTankInitialLevel;        % Initial water level in tanks
        NodeTankInitialWaterVolume;  % Initial water volume in tanks
        NodeTankMaximumWaterLevel;   % Maximum water level in tanks
        NodeTankMaximumWaterVolume;  % Maximum water volume
        NodeTankMinimumFraction;     % Fraction of the total tank volume devoted to the inlet/outlet compartment
        NodeTankMinimumWaterLevel;   % Minimum water level
        NodeTankMinimumWaterVolume;  % Minimum water volume
        NodeTankMixingModelCode;     % Code of mixing model (MIXED:0, 2COMP:1, FIFO:2, LIFO:3)
        NodeTankMixingModelType;     % Type of mixing model (MIXED, 2COMP, FIFO, or LIFO)
        NodeTankMixZoneVolume;       % Mixing zone volume
        NodeTankNameID;              % Name ID of Tanks
        NodeTankReservoirCount;      % Number of tanks and reservoirs
        NodeTankVolumeCurveIndex;    % Index of curve for tank volumes
        NodeTankVolumeUnits;         % Units for volume
        NodeType;                    % ID of node type
        NodeTypeIndex;               % Index of nodetype
        OptionsAccuracyValue;        % Convergence value (0.001 is default)
        OptionsEmitterExponent;      % Exponent of pressure at an emmiter node (0.5 is default)
        OptionsHeadLossFormula;      % Headloss formula (Hazen-Williams, Darcy-Weisbach or Chezy-Manning)
        OptionsHydraulics;           % Save or Use hydraulic soltion. *** Not implemented ***
        OptionsMaxTrials;            % Maximum number of trials (40 is default)
        OptionsPattern;              % *** Not implemented *** % but get with BinOptionsPattern
        OptionsPatternDemandMultiplier; % Multiply demand values (1 is default)
        OptionsQualityTolerance;     % Tolerance for water quality (0.01 is default)
        OptionsSpecificGravity;      % *** Not implemented *** % but get with BinOptionsSpecificGravity
        OptionsUnbalanced;           % *** Not implemented *** % but get with BinOptionsUnbalanced
        OptionsViscosity;            % *** Not implemented *** % but get with BinOptionsViscosity
        OptionsHeadError;
        OptionsFlowChange;
        Pattern;                     % Get all patterns - matrix
        PatternAverageValue;         % Average value of patterns
        PatternCount;                % Number of patterns
        PatternDemandsUnits;         % Units for demands
        PatternIndex;                % Indices of the patterns
        PatternLengths;              % Length of the patterns
        PatternNameID;               % ID of the patterns
        QualityChemName;             % Quality Chem Name
        QualityChemUnits;            % Quality Chem Units
        QualityCode;                 % Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        QualityReactionCoeffBulkUnits;   % Bulk reaction coefficient units
        QualityReactionCoeffWallUnits;   % Wall reaction coefficient units
        QualitySourceMassInjectionUnits; % Units for source mass injection
        QualityTraceNodeIndex;       % Index of trace node (0 if QualityCode<3)
        QualityType;                 % Water quality analysis type (None/Chemical/Age/Trace)
        QualityUnits;                % Units for quality concentration.
        QualityWaterAgeUnits;        % Units for water age
        RelativeError;               % Relative error - hydraulic simulation statistic
        %         RulePremises;
        %         RuleTrueActions;
        %         RuleFalseActions;
        %         RulePriority;
        TempInpFile;                 % Name of the temporary input file
        TimeHaltFlag;                % Number of halt flag
        TimeHTime;                   % Number of htime
        TimeHydraulicStep;           % Hydraulic time step
        TimeNextEvent;               % Find the lesser of the hydraulic time step length, or the time to next fill/empty
        TimePatternStart;            % Pattern start time
        TimePatternStep;             % Pattern Step
        TimeQualityStep;             % Quality Step
        TimeReportingPeriods;        % Reporting periods
        TimeReportingStart;          % Start time for reporting
        TimeReportingStep;           % Reporting time step
        TimeRuleControlStep;         % Time step for evaluating rule-based controls
        TimeSimulationDuration;      % Simulation duration
        TimeStartTime;               % Number of start time
        TimeStatisticsIndex;         % Index of time series post-processing type ('NONE':0, 'AVERAGE':1, 'MINIMUM':2, 'MAXIMUM':3, 'RANGE':4)
        TimeStatisticsType;          % Type of time series post-processing ('NONE', 'AVERAGE', 'MINIMUM', 'MAXIMUM', 'RANGE')
        ToolkitConstants;            % Contains all parameters from epanet2.h
        Units_SI_Metric;             % Equal with 1 if is SI-Metric
        Units_US_Customary;          % Equal with 1 if is US-Customary
        Version;                     % EPANET version

        % Parameters used with EPANET MSX
        MSXLibEPANET;                % MSX EPANET library dll
        MSXLibEPANETPath;            % MSX EPANET library path
        MSXConstantsNameID;          % ID name of constants
        MSXConstantsValue;           % Value of constants
        MSXConstantsCount;           % Number of constants
        MSXConstantsIndex;           % Index of constants
        MSXParametersCount;          % Number of parameters
        MSXPatternsCount;            % Number of msx patterns
        MSXSpeciesCount;             % Number of species
        MSXLinkInitqualValue;        % Initial concentration of chemical species assigned to links of the pipe network
        MSXNodeInitqualValue;        % Initial concentration of chemical species assigned to nodes
        MSXFile;                     % Name of the msx file
        MSXTempFile;                 % Name of the temp msx file
        MSXParametersNameID;         % ID name of parameters
        MSXParametersIndex;          % Index name of parameters
        MSXParametersPipesValue;     % Value of reaction parameters for pipes
        MSXParametersTanksValue;     % Value of reaction parameters for tanks
        MSXPatternsNameID;           % ID name of msx patterns
        MSXPatternsIndex;            % Index of msx patterns
        MSXPatternsLengths;          % Number of time periods in all patterns
        MSXPattern;                  % Get all msx patterns
        MSXEquationsPipes;           % Species dynamics in pipes
        MSXSources;                  % Sources info
        MSXSourceLevel;              % Value of all nodes source level
        MSXSourceNodeNameID;         % ID label of all nodes
        MSXSourcePatternIndex;       % Value of all node source pattern index
        MSXSourceType;               % Value of all node source type 'CONCEN', 'MASS', 'SETPOINT', 'FLOWPACED'
        MSXSourceTypeCode;           % Code of source type
        MSXSpeciesATOL;              % Absolute tolerance used to determine when two concentration levels of a species are the same
        MSXSpeciesIndex;             % Index name of species
        MSXSpeciesNameID;            % ID name of species
        MSXSpeciesRTOL;              % Relative accuracy level on a species concentration used to adjust time steps in the RK5 and ROS2 integration methods
        MSXSpeciesType;              % Type of all species, bulk or wall
        MSXSpeciesUnits;             % Species mass units
        MSXEquationsTanks;           % Species dynamics in tanks
        MSXEquationsTerms;           % Species dynamics in terms

        % Parameters used when the Binary mode is used
        Bin;                         % Check if use Bin functions (use saveInputFile if is 1)
        BinControlLinksID;           % Set of control links ID
        BinControlNodesID;           % Set of control nodes ID
        BinControlRulesCount;        % Number of controls
        BinControlsInfo;             % Controls info
        BinCountInitialQualitylines; % Count lines used by quality section
        BinCountPatternlines;        % Count lines used by pattern section
        BinCountReactionlines;       % Count lines used by reaction section
        BinCountStatuslines;         % Count lines used by status section
        BinCurveAllLines;            % Curves info from section
        BinCurveCount;               % Number of curves
        BinCurveNameID;              % ID name of curves
        BinCurveTypes;               % Type of curves
        BinCurveXvalue;              % X-value of curves
        BinCurveYvalue;              % Y-value of curves
        BinLinkBulkReactionCoeff;    % Bulk reaction coefficient of each link
        BinLinkCount;                % Number of links
        BinLinkDiameters;            % Diameter of each link
        BinLinkFlowUnits;            % Units of flow
        BinLinkFromNode;             % IDs of node at start of link
        BinLinkGlobalBulkReactionCoeff;% Global Bulk Reaction Coeff
        BinLinkGlobalWallReactionCoeff;% Global Wall Reaction Coeff
        BinLinkInitialStatus;        % Initial status of links
        BinLinkInitialStatusNameID;  % Name ID of links where is status refers in BinLinkInitialStatus
        BinLinkLengths;              % Length of links
        BinLinkNameID;               % Name ID of links
        BinLinkPipeCount;            % Number of pipes
        BinLinkPipeDiameters;        % Diameter of each pipe
        BinLinkPipeIndex;            % Index of pipe links
        BinLinkPipeLengths;          % Length of pipes
        BinLinkPipeMinorLoss;        % Minor loss coefficient of pipes
        BinLinkPipeNameID;           % Name ID of pipes
        BinLinkPipeRoughness;        % Roughness coefficient of pipes
        BinLinkPipeStatus;           % Initial status of pipes
        BinLinkPumpCount;            % Number of pumps
        BinLinkPumpCurveNameID;      % Curve Name ID used from pumps
        BinLinkPumpIndex;            % Index of pumps
        BinLinkPumpNameID;           % Name ID of pumps
        BinLinkPumpNameIDPower;      % Name ID of pumps with Power
        BinLinkPumpPatterns;         % Patterns used from pumps
        BinLinkPumpPower;            % Power of pumps
        BinLinkPumpStatus;           % Initial status of pumps
        BinLinkPumpStatusNameID;     % Name ID of pumps where is status refers in BinLinkPumpStatus
        BinLinkRoughnessCoeff;       % Roughness coefficient of links
        BinLinkSettings;             % Initial settings of links
        BinLinkToNode;               % IDs of node at start of link
        BinLinkType;                 % ID of link type
        BinLinkValveCount;           % Number of valves
        BinLinkValveDiameters;       % Diameter of each valve
        BinLinkValveIndex;           % Index of valve links
        BinLinkValveMinorLoss;       % Minor loss coefficient of valves
        BinLinkValveNameID;          % Name ID of valves
        BinLinkValveSetting;         % Initial settings of valves
        BinLinkValveStatus;          % Initial status of valves
        BinLinkValveStatusNameID;    % Name ID of valves where is status refers in BinLinkValveStatus
        BinLinkValveType;            % Valve type, 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'
        BinLinkWallReactionCoeff;    % Wall reaction coefficient of links
        BinNodeBaseDemands;          % Base demands of nodes
        BinNodeCoordinates;          % Coordinates for each node (long/lat & intermediate pipe coordinates) - vertices
        BinNodeCount;                % Number of nodes
        BinNodeJunDemandPatternNameID;% ID of demand patterns
        BinNodeElevations;           % Elevation of nodes
        BinNodeInitialQuality;       % Initial quality of nodes
        BinNodeJunctionCount;        % Number of junctions
        BinNodeJunctionElevation;    % Elevation of junctions
        BinNodeJunctionIndex;        % Index of junctions
        BinNodeJunctionNameID;       % Name ID of junctions
        BinNodeJunctionsBaseDemands; % Base demands of junctions
        BinNodeJunctionsBaseDemandsID;% Name ID of junctions where is basedemand refers in BinNodeJunctionsBaseDemands
        BinNodeNameID;               % Name ID of nodes
        BinNodePressureUnits;        % Units for Pressure
        BinNodeResDemandPatternNameID;% ID of demand patterns for reservoirs
        BinNodeReservoirCount;       % Number of reservoirs
        BinNodeReservoirElevation;   % Elevation of reservoirs
        BinNodeReservoirIndex;       % Index of reservoirs
        BinNodeReservoirNameID;      % Name ID of reservoirs
        BinNodeSourcePatternIndex;   % Index of pattern for node sources
        BinNodeSourcePatternNameID;  % ID of pattern for node sources
        BinNodeSourceQuality;        % Quality of node sources
        BinNodeSourceType;           % Source Types
        BinNodeSourceTypeIndex;      % Index of source type
        BinNodeTankCount;            % Number of tanks
        BinNodeTankDiameter;         % Diameters of tanks
        BinNodeTankElevation;        % Elevations of tanks
        BinNodeTankIndex;            % Index of tanks
        BinNodeTankInitialLevel;     % Initial water level in tanks
        BinNodeTankMaximumWaterLevel;% Maximum water level in tanks
        BinNodeTankMinimumFraction;  % Fraction of the total tank volume devoted to the inlet/outlet compartment
        BinNodeTankMinimumWaterLevel;% Minimum water level
        BinNodeTankMinimumWaterVolume;% Minimum water volume
        BinNodeTankMixID;            % Name ID of tanks mix
        BinNodeTankMixModel;         % Mix Model Type
        BinNodeTankNameID;           % Name ID of tanks
        BinNodeTankReservoirCount;   % Number of reservoirs
        BinNodeType;                 % ID of node type
        BinOptionsAccuracyValue;     % Convergence value (0.001 is default)
        BinOptionsDiffusivity;       % Diffusivity value (1 is default)
        BinOptionsEmitterExponent;   % Exponent of pressure at an emmiter node (0.5 is default)
        BinOptionsHeadloss;          % Headloss formula (Hazen-Williams, Darcy-Weisbach or Chezy-Manning)
        BinOptionsMaxTrials;         % Maximum number of trials (40 is default)
        BinOptionsPattern;           % Pattern to be applied to all junctions where no demand pattern was specified
        BinOptionsPatternDemandMultiplier;% Multiply demand values (1 is default)
        BinOptionsQualityTolerance;  % Tolerance for water quality (0.01 is default)
        BinOptionsSpecificGravity;   % Ratio of the density of the fluid being modeled to that of water at 4 deg. C (unitless)
        BinOptionsUnbalanced;        % Determines what happens if a hydraulic solution cannot be reached (STOP is default) STOP/CONTINUE
        BinOptionsViscosity;         % Kinematic viscosity of the fluid being modeled relative to that of water at 20 deg. C (1.0 centistoke)(1 is default)
        BinPatternCount;             % Number of patterns
        BinPatternLengths;           % Length of the patterns
        BinPatternMatrix;            % Get all patterns - matrix
        BinPatternNameID;            % ID of the patterns
        BinPatternValue;             % Patterns values for all
        BinQualityCode;              % Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        BinQualityTraceNodeID;       % ID of trace node (0 if QualityCode<3)
        BinQualityTraceNodeIndex;    % Index of trace node (0 if QualityCode<3)
        BinQualityType;              % Water quality analysis type (None:0/Chemical:1/Age:2/Trace:3)
        BinQualityUnits;             % Units for quality concentration.
        BinRulesControlLinksID;      % Set of rule links ID
        BinRulesControlNodesID;      % Set of rule nodes ID
        BinRulesControlsInfo;        % Rules info from section
        BinRulesCount;               % Number of rules
        BinUnits_SI_Metric;          % Equal with 1 if is SI-Metric
        BinTempfile;                 % Name of the temp input file
        BinTimeHydraulicStep;        % Hydraulic time step
        BinTimePatternStart;         % Pattern start time
        BinTimePatternStep;          % Pattern Step
        BinTimeQualityStep;          % Quality Step
        BinTimeReportingStart;       % Start time for reporting
        BinTimeReportingStep;        % Reporting time step
        BinTimeSimulationDuration;   % Simulation duration
        BinTimeStatistics;           % Type of time series post-processing ('NONE', 'AVERAGE', 'MINIMUM', 'MAXIMUM', 'RANGE')
        BinTimeStatisticsIndex;      % Index of time series post-processing type ('NONE':0, 'AVERAGE':1, 'MINIMUM':2, 'MAXIMUM':3, 'RANGE':4)
        BinUnits;                    % Units of all parameters
        BinUnits_US_Customary;       % Equal with 1 if is US-Customary
        CMDCODE;                     % Code=1 Hide, Code=0 Show (messages at command window)
    end
    properties (Constant = true)
        classversion='v2.2.0'; % 05/10/2021

        LOGOP={'IF', 'AND', 'OR'} % Constants for rule-based controls: 'IF', 'AND', 'OR' % EPANET Version 2.2
        RULEOBJECT={'NODE', 'LINK', 'SYSTEM'}; % Constants for rule-based controls: 'NODE','LINK','SYSTEM' % EPANET Version 2.2
        RULEVARIABLE={'DEMAND', 'HEAD', 'GRADE', 'LEVEL', 'PRESSURE', 'FLOW', 'STATUS', ... % Constants for rule-based controls: 'DEMAND', 'HEAD', 'GRADE' etc. % EPANET Version 2.2
        'SETTING', 'POWER', 'TIME', 'CLOCKTIME', 'FILLTIME', 'DRAINTIME'};
        RULEOPERATOR={'=', '~=', '<=', '>=', '<', '>', 'IS', 'NOT', 'BELOW', 'ABOVE'}; % Constants for rule-based controls: '=', '~=', '<=' etc. % EPANET Version 2.2
        RULESTATUS={'OPEN', 'CLOSED', 'ACTIVE'}; % Constants for rule-based controls: 'OPEN', 'CLOSED', 'ACTIVE' % EPANET Version 2.2
        RULEPREMISECHECK={'NODE', 'JUNCTION', 'RESERVOIR', 'TANK', 'LINK', 'PIPE', 'PUMP', 'VALVE', 'SYSTEM'}; % Constants for rule-based controls: 'NODE', 'JUNCTION' etc. % EPANET Version 2.2
        TYPECONTROL={'LOWLEVEL', 'HIGHLEVEL', 'TIMER', 'TIMEOFDAY'}; % Constants for control: 'LOWLEVEL', 'HILEVEL', 'TIMER', 'TIMEOFDAY'
        TYPECURVE={'VOLUME', 'PUMP', 'EFFICIENCY', 'HEADLOSS', 'GENERAL'}; % Constants for pump curves: 'PUMP', 'EFFICIENCY', 'VOLUME', 'HEADLOSS' % EPANET Version 2.2
        TYPELINK={'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'}; % Constants for links: 'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV', 'VALVE'
        TYPEMIXMODEL={'MIX1', 'MIX2', 'FIFO', 'LIFO'}; % Constants for mixing models: 'MIX1', 'MIX2', 'FIFO', 'LIFO'
        TYPENODE={'JUNCTION', 'RESERVOIR', 'TANK'}; % Contants for nodes: 'JUNCTION', 'RESERVOIR', 'TANK'
        TYPEPUMP={'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'}; % Constants for pumps: 'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'
        TYPEQUALITY={'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'}; % Constants for quality: 'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'
        TYPEREPORT={'YES', 'NO', 'FULL'}; % Constants for report: 'YES', 'NO', 'FULL'
        TYPESOURCE={'CONCEN', 'MASS', 'SETPOINT', 'FLOWPACED'}; % Constants for sources: 'CONCEN', 'MASS', 'SETPOINT', 'FLOWPACED'
        TYPESTATS={'NONE', 'AVERAGE', 'MINIMUM', 'MAXIMUM', 'RANGE'}; % Constants for statistics: 'NONE', 'AVERAGE', 'MINIMUM', 'MAXIMUM', 'RANGE'
        TYPEUNITS={'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'}; % Constants for units: 'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'
        TYPEHEADLOSS={'HW', 'DW', 'CM'}; % Constants of headloss types: HW: Hazen-Williams, DW: Darcy-Weisbach, CM: Chezy-Manning
        TYPESTATUS = {'CLOSED', 'OPEN'}; % Link status
        TYPEPUMPSTATE = {'XHEAD', '', 'CLOSED', 'OPEN', '', 'XFLOW'}; % Link PUMP status
        %   d.TYPEPUMPSTATE(res.State + 1)
    %   EN_PUMP_XHEAD   = 0,  //!< Pump closed - cannot supply head
    %   EN_PUMP_CLOSED  = 2,  //!< Pump closed
    %   EN_PUMP_OPEN    = 3,  //!< Pump open
    %   EN_PUMP_XFLOW   = 5   //!< Pump open - cannot supply flow
        TYPEBINSTATUS = {'CLOSED (MAX. HEAD EXCEEDED)', 'TEMPORARILY CLOSED', 'CLOSED',...
            'OPEN', 'ACTIVE(PARTIALY OPEN)', 'OPEN (MAX. FLOW EXCEEDED',...
            'OPEN (PRESSURE SETTING NOT MET)'};
        %   0 = closed (max. head exceeded)
        %   1 = temporarily closed
        %   2 = closed
        %   3 = open
        %   4 = active (partially open)
        %   5 = open (max. flow exceeded)
        %   6 = open (flow setting not met)
        %   7 = open (pressure setting not met)
        DEMANDMODEL = {'DDA', 'PDA'}; % Demand model types. DDA #0 Demand driven analysis, PDA #1 Pressure driven analysis.
        MSXTYPEAREAUNITS={'FT2', 'M2', 'CM2'}; % sets the units used to express pipe wall surface area
        MSXTYPERATEUNITS={'SEC', 'MIN', 'HR', 'DAY'}; % is the units in which all reaction rate terms are expressed
        MSXTYPESOLVER={'EUL', 'RK5', 'ROS2'}; % is the choice of numerical integration method used to solve the reaction system
        MSXTYPECOUPLING={'FULL', 'NONE'}; % is the choice of numerical integration method used to solve the reaction system
        MSXTYPECOMPILER={'NONE', 'VC', 'GC'}; % is the choice of numerical integration method used to solve the reaction system
        MSXTYPESOURCE={'NOSOURCE', 'CONCEN', 'MASS', 'SETPOINT', 'FLOWPACED'}; % Constants for sources: 'NO SOURCE''CONCEN', 'MASS', 'SETPOINT', 'FLOWPACED'
        MSXTYPENODE=0; % for a node
        MSXTYPELINK=1; % for a link
    end
    properties (Access = private)
        solve = 1;
    end
    methods (Access = private)
        function Errcode = setFlowUnits(obj, unitcode, typecode, varargin)
            if ~typecode
                [Errcode]=Options(obj, unitcode);
                if isempty(varargin{1})
                    return
                else
                    varargin = char(varargin{1});
                end
                obj.saveBinInpFile(varargin);

            elseif typecode
                [Errcode] = ENsetflowunits(obj.LibEPANET, unitcode);
                if nargin==4
                    if isempty(varargin{1})
                        return
                    else
                        varargin = char(varargin{1});
                    end
                    obj.saveInputFile(varargin);
                end
            end
        end
        function value = get_link_info(obj, constant, varargin)
            [indices, value] = getLinkIndices(obj, varargin);
            j=1;
            for i=indices
                [obj.Errcode, value(j)] = obj.ENgetlinkvalue(i, constant, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function set_Node_Link(obj, param, fun, propertie, value, varargin)
            if strcmp(param, 'tank')
                count = obj.getNodeTankCount;
                indices = obj.getNodeTankIndex;
            elseif strcmp(param, 'pump')
                count = obj.getLinkPumpCount;
                indices = obj.getLinkPumpIndex;
            end
            if isempty(varargin{1})
                j=1;
                for i=1:count
                    [obj.Errcode] = eval([fun, '(indices(i), propertie, value(j), obj.LibEPANET)']);
                    error(obj.getError(obj.Errcode));
                    if ~isscalar(value)
                        j=j+1;
                    end
                end
            else
                varargin{1} = varargin{1}{1};
                if ~ismember(value, indices)
                    value = indices(value);
                end
                j=1;
                for i=1:length(value)
                    [obj.Errcode] = eval([fun, '(value(i), propertie, varargin{1}(j), obj.LibEPANET)']);
                    error(obj.getError(obj.Errcode));
                    if ~isscalar(varargin{1})
                        j=j+1;
                    end
                end
            end
        end
        function value = get_node_link(obj, param, fun, propertie, varargin)
            if strcmp(param, 'tank')
                indices = obj.getNodeTankIndex; % obj.getNodeIndex;
            elseif strcmp(param, 'pump')
                indices = obj.getLinkPumpIndex; % obj.getLinkIndex;
            end
            if isempty(varargin{1})
                value = zeros(1, length(indices));
                j = 1;
                for i=indices
                    [obj.Errcode, value(j)] = eval([fun, '(i, propertie, obj.LibEPANET)']);
                    j=j+1;
                end
            else
                varargin{1} = varargin{1}{1};
                if ~ismember(varargin{1}, indices)
                    varargin{1} = indices(varargin{1});
                end
                value = zeros(1, length(varargin{1}));
                j = 1;
                for i=1:length(varargin{1})
                    [obj.Errcode, value(j)] = eval([fun, '(varargin{1}(i), propertie, obj.LibEPANET)']);
                    if ~isscalar(varargin{1})
                        j=j+1;
                    end
                end
            end
        end
        function set_node_demand_pattern(obj, fun, propertie, value, extra)

            categ = 1;
            if length(extra) == 2
                indices = value;
                categ = extra{1};
                param = extra{2};
            elseif length(extra) == 1
                indices = value;
                param = extra{1};
            elseif isempty(extra)
                indices = getNodeJunctionIndices(obj,[]);
                param = value;
                if iscell(param)
                    categ = length(param);
                end
            end
            for c=1:categ
                if isempty(extra) && iscell(value)
                    param = value{c};
                end
                j = 1;
                for i = indices
                    if ~ismember(i, obj.getNodeReservoirIndex) && sum(strcmp(obj.libFunctions, fun))
                        if c > obj.getNodeDemandCategoriesNumber(i)
                            obj.addNodeJunctionDemand(i, param(j));
                        else
                            [obj.Errcode] = eval([fun, '(i, c, param(j), obj.LibEPANET)']);
                        end
                    else
                        [obj.Errcode] = ENsetnodevalue(i, propertie, param(j), obj.LibEPANET);
                    end
                    j=j+1;
                end
            end
        end
        function value = get_node_tank_mixining_model(obj, varargin)
            obj.NodeTankMixingModelCode = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MIXMODEL, varargin);
            obj.NodeTankMixingModelType = obj.TYPEMIXMODEL(obj.NodeTankMixingModelCode + 1);
            value={obj.NodeTankMixingModelCode obj.NodeTankMixingModelType};
        end
        function index = check_if_numeric(obj, id)
            if isnumeric(id)
                index = id;
            else
                index = obj.getLinkIndex(id);
            end
        end

    end
    methods(Static)
        function [Errcode] = ENwriteline(line, LibEPANET)
            [Errcode]=calllib(LibEPANET, 'ENwriteline', line);
        end
        function [Errcode] = ENaddpattern(patid, LibEPANET)
            Errcode=calllib(LibEPANET, 'ENaddpattern', patid);
        end
        function [Errcode] = ENclose(LibEPANET)
            [Errcode]=calllib(LibEPANET, 'ENclose');
        end
        function [Errcode] = ENcloseH(LibEPANET)
            [Errcode]=calllib(LibEPANET, 'ENcloseH');
        end
        function [Errcode] = ENsetlinkid(index, newid, LibEPANET)
            % EPANET Version 2.2
            [Errcode]=calllib(LibEPANET, 'ENsetlinkid', index, newid);
        end
        function [Errcode, value] = ENgetnumdemands(index, LibEPANET)
            %epanet20100
            [Errcode, value]=calllib(LibEPANET, 'ENgetnumdemands', index, 0);
        end
        function [Errcode, value] = ENgetdemandpattern(index, numdemands, LibEPANET)
            %epanet20100
            [Errcode, value]=calllib(LibEPANET, 'ENgetdemandpattern', index, numdemands, 0);
        end
        function [Errcode] = ENcloseQ(LibEPANET)
            [Errcode]=calllib(LibEPANET, 'ENcloseQ');
        end
        function [Errcode, demandIndex] = ENgetdemandindex(nodeindex, demandName, LibEPANET)
            % EPANET Version 2.2
            [Errcode, ~, demandIndex]=calllib(LibEPANET, 'ENgetdemandindex', nodeindex, demandName, 0);
        end
        function [Errcode, value] = ENgetstatistic(code, LibEPANET)
            %epanet20100
            [Errcode, value]=calllib(LibEPANET, ' ENgetstatistic', code, 0);
        end
        function [Errcode, ctype, lindex, setting, nindex, level] = ENgetcontrol(cindex, LibEPANET)
            [Errcode, ctype, lindex, setting, nindex, level]=calllib(LibEPANET, 'ENgetcontrol', cindex, 0, 0, 0, 0, 0);
        end
        function [Errcode, count] = ENgetcount(countcode, LibEPANET)
            [Errcode, count]=calllib(LibEPANET, 'ENgetcount', countcode, 0);
        end
        function [errmsg, e] = ENgeterror(Errcode, LibEPANET)
            e=0; errmsg='';
            if Errcode[e, errmsg] = calllib(LibEPANET, 'ENgeterror', Errcode, char(32*ones(1, 79)), 79);
            if e, [e, errmsg] = calllib(LibEPANET, 'ENgeterror', e, char(32*ones(1, 79)), 79); end
        end end
        function [Errcode, flowunitsindex] = ENgetflowunits(LibEPANET)
            [Errcode, flowunitsindex]=calllib(LibEPANET, 'ENgetflowunits', 0);
        end
        function [Errcode, id] = ENgetlinkid(index, LibEPANET)
            id=char(32*ones(1, 31));
            [Errcode, id]=calllib(LibEPANET, 'ENgetlinkid', index, id);
        end
        function [Errcode, index] = ENgetlinkindex(id, LibEPANET)
            [Errcode, ~, index]=calllib(LibEPANET, 'ENgetlinkindex', id, 0);
        end
        function [Errcode, from, to] = ENgetlinknodes(index, LibEPANET)
            [Errcode, from, to]=calllib(LibEPANET, 'ENgetlinknodes', index, 0, 0);
            from = double(from);
            to = double(to);
        end
        function [Errcode, typecode] = ENgetlinktype(index, LibEPANET)
            [Errcode, typecode]=calllib(LibEPANET, 'ENgetlinktype', index, 0);
            typecode = double(typecode);
            if ~isnumeric(typecode), typecode = getTypeLink(typecode); end
        end
        function [Errcode, value] = ENgetlinkvalue(index, paramcode, LibEPANET)
            [Errcode, value]=calllib(LibEPANET, 'ENgetlinkvalue', index, paramcode, 0);
            value = double(value);
        end
        function [Errcode, id] = ENgetnodeid(index, LibEPANET)
            id=char(32*ones(1, 31));
            [Errcode, id]=calllib(LibEPANET, 'ENgetnodeid', index, id);
        end
        function [Errcode] = ENsetnodeid(index, newid, LibEPANET)
            % EPANET Version 2.2
            [Errcode]=calllib(LibEPANET, 'ENsetnodeid', index, newid);
        end
        function [Errcode, index] = ENgetnodeindex(id, LibEPANET)
            [Errcode, ~, index]=calllib(LibEPANET, 'ENgetnodeindex', id, 0);
        end
        function [Errcode, type] =ENgetnodetype(index, LibEPANET)
            [Errcode, type]=calllib(LibEPANET, 'ENgetnodetype', index, 0);
        end
        function [Errcode, value] = ENgetoption(optioncode, LibEPANET)
            [Errcode, value]=calllib(LibEPANET, 'ENgetoption', optioncode, 0);
        end
        function [Errcode, id] = ENgetpatternid(index, LibEPANET)
            id=char(32*ones(1, 31));
            [Errcode, id]=calllib(LibEPANET, 'ENgetpatternid', index, id);
        end
        function [Errcode, id] = ENgetcurveid(index, LibEPANET)
            % EPANET Version dev2.1
            id=char(32*ones(1, 31));
            [Errcode, id]=calllib(LibEPANET, 'ENgetcurveid', index, id);
        end
        function [Errcode] = ENsetcurveid(index, id, LibEPANET)
            % EPANET Version 2.2
            [Errcode,~]=calllib(LibEPANET, 'ENsetcurveid', index,id);
        end
        function [Errcode] = ENsetpatternid(index, id, LibEPANET)
            % EPANET Version 2.2
            [Errcode,~]=calllib(LibEPANET, 'ENsetpatternid', index,id);
        end
        function [Errcode, type] =ENgetcurvetype(index, LibEPANET)
            % EPANET Version 2.2
            [Errcode, type]=calllib(LibEPANET,ENgetcurvetype', index, 0);
        end
        function [Errcode, index] = ENgetpatternindex(id, LibEPANET)
            [Errcode, ~, index]=calllib(LibEPANET, 'ENgetpatternindex', id, 0);
        end
    end
    methods
        function obj = epanet(varargin)
            %Constructor of the EPANET Class
            try unloadlibrary('epanet2');catch; end
            try unloadlibrary('epanetmsx');catch; end
            % DLLs
            arch = computer('arch');
            pwdepanet = fileparts(which(mfilename));
            if strcmpi(arch, 'win64')% if no DLL is given, select one automatically
                obj.LibEPANETpath = [pwdepanet, '/64bit/'];
            elseif strcmpi(arch, 'win32')
                obj.LibEPANETpath = [pwdepanet, '/32bit/'];
            end
            if isunix
                obj.LibEPANETpath = [pwdepanet, '/glnx/'];
                obj.LibEPANET = 'libepanet';
            end
            if ismac
                obj.LibEPANETpath = [pwdepanet, '/mac/'];
                obj.LibEPANET = 'libepanet';
            end
            if ~isdeployed
                obj.InputFile=which(varargin{1}); % Get name of INP file
            else
                obj.InputFile=varargin{1};
            end
            % Bin functions
            if nargin>1
                if strcmpi(varargin{2}, 'BIN')
                    obj.LibEPANET = '';
                    obj.BinTempfile=[obj.InputFile(1:end-4), '_temp.inp'];
                    copyfile(obj.InputFile, obj.BinTempfile);
                    obj.InputFile=obj.BinTempfile;
                    obj.Bin=0;
                    if nargin==3, if strcmpi(varargin{3}, 'LOADFILE'); return; end; end
                    obj = BinUpdateClass(obj);
                    obj.saveBinInpFile;
                    return;
                end
            end
            obj.Bin=1;
            [~, inp]=fileparts(obj.InputFile);
            if isempty(inp)
                if nargin==2 && strcmpi(varargin{2}, 'CREATE')
                    % skip
                else
                    error(['File "', varargin{1}, '" is not a valid.']);
                end
            end
            if nargin==2 && ~strcmpi(varargin{2}, 'loadfile') && ~strcmpi(varargin{2}, 'CREATE')% e.g. d = epanet('Net1.inp', 'epanet2');
                [pwdDLL, obj.LibEPANET] = fileparts(varargin{2}); % Get DLL LibEPANET (e.g. epanet20012x86 for 32-bit)
                if isempty(pwdDLL)
                    pwdDLL = pwd;
                end
                obj.LibEPANETpath = [pwdDLL, '\'];
                try ENLoadLibrary(obj.LibEPANETpath, obj.LibEPANET, 0);
                catch
                   obj.Errcode=-1;
                   error(['File "', obj.LibEPANET, '" is not a valid win application.']);
                end
            elseif ~isunix
                obj.LibEPANET = 'epanet2';
            end
            if ~isdeployed
                if ~exist(obj.InputFile, 'file')
                    if nargin==2 && strcmpi(varargin{2}, 'CREATE')
                        % skip
                    else
                        error(['File "', varargin{1}, '" does not exist in folder. (e.g. addpath(genpath(pwd));)']);
                    end
                end
            end
            %Load EPANET Library
            ENLoadLibrary(obj.LibEPANETpath, obj.LibEPANET);
            disp([' (EMT version {', obj.classversion,'}).'])
            %Load parameters
            obj.ToolkitConstants = obj.getToolkitConstants;

            %For the getComputedQualityTimeSeries
            obj.solve = 0;
            %Open the file
            if ~isempty(obj.InputFile)
                if nargin==2 && strcmpi(varargin{2}, 'CREATE')
                    warning(['Network name "', inp ,'.inp" already exists.'])
                end
                obj.Errcode=ENopen(obj.InputFile, [obj.InputFile(1:end-4), '.txt'], '', obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                obj.InputFile = varargin{1};
                % initializes an EPANET project that isn't opened with an input file
                obj.initializeEPANET(obj.ToolkitConstants.EN_GPM, obj.ToolkitConstants.EN_HW);
                warning('Initializes the EPANET project!');
            end

            %Save the temporary input file
            obj.BinTempfile=[obj.InputFile(1:end-4), '_temp.inp'];
            obj.saveInputFile(obj.BinTempfile); %create a new INP file (Working Copy) using the SAVE command of EPANET
            obj.closeNetwork;  %ENclose; %Close input file
            %Load temporary file
            rptfile = [obj.InputFile(1:end-4), '_temp.txt'];
            binfile = [obj.InputFile(1:end-4), '_temp.bin'];
            obj.Errcode=ENopen(obj.BinTempfile, rptfile, binfile, obj.LibEPANET);
            if obj.Errcode
                error(obj.getError(obj.Errcode));
            else
                disp(['Loading File "', varargin{1}, '"...']);
            end
            % Hide messages at command window from bin computed
            obj.CMDCODE=1;

            % Using new variable for temp file
            obj.TempInpFile = obj.BinTempfile;
            % Load file only, return
            if nargin==2
                if strcmpi(varargin{2}, 'LOADFILE')
                    obj.libFunctions = libfunctions(obj.LibEPANET);
                    disp(['Input File "', varargin{1}, '" loaded sucessfuly.']);
                    return;
                end
            end
            % Get some link data
            lnkInfo = obj.getLinksInfo;
            getFields_link_info = fields(lnkInfo);
            for i=1:length(getFields_link_info)
                obj.(getFields_link_info{i}) = eval(['lnkInfo.', getFields_link_info{i}]);
            end
            % Get some node data
            ndInfo = obj.getNodesInfo;
            getFields_node_info = fields(ndInfo);
            for i=1:length(getFields_node_info)
                obj.(getFields_node_info{i}) = eval(['ndInfo.', getFields_node_info{i}]);
            end
            if obj.getVersion > 20101
                % Get demand model type and parameters
                demModelInfo = obj.getDemandModel;
                getFields_demModelInfo = fields(demModelInfo);
                for i=1:length(getFields_demModelInfo)
                    obj.(getFields_demModelInfo{i}) = eval(['demModelInfo.', getFields_demModelInfo{i}]);
                end
            end
            %Get all the countable network parameters
            obj.NodeCount = obj.getNodeCount;
            if ~obj.NodeCount
                return;
            end
            obj.NodeTankReservoirCount = obj.getNodeTankReservoirCount;
            obj.LinkCount = obj.getLinkCount;
            obj.PatternCount = obj.getPatternCount;
            obj.CurveCount = obj.getCurveCount;
            obj.CurveIndex = obj.getCurveIndex;
            obj.ControlRulesCount = obj.getControlRulesCount;
            obj.NodeJunctionCount = obj.NodeCount-obj.NodeTankReservoirCount; %obj.getNodeJunctionCount;
            % Get type of the parameters
            obj.LinkType=obj.TYPELINK(obj.LinkTypeIndex+1);
            obj.NodeType=obj.TYPENODE(obj.NodeTypeIndex+1);
            %Get all the countable network parameters
            obj.LinkPipeCount = sum(strcmp(obj.LinkType, 'PIPE')+strcmp(obj.LinkType, 'CVPIPE')); %obj.getLinkPipeCount;
            obj.LinkPumpCount = sum(strcmp(obj.LinkType, 'PUMP')); %obj.getLinkPumpCount;
            obj.NodeReservoirCount = sum(strcmp(obj.NodeType, 'RESERVOIR')); %obj.getNodeReservoirCount;
            obj.NodeTankCount = obj.NodeCount-obj.NodeJunctionCount-obj.NodeReservoirCount; %obj.getNodeTankCount;
            obj.LinkValveCount = obj.LinkCount-obj.LinkPipeCount-obj.LinkPumpCount;%obj.getLinkValveCount;
            %Get all the controls
            obj.Controls = obj.getControls;
            %Get the flow units
            obj.LinkFlowUnits = obj.getFlowUnits;
            %Get all the link data
            obj.LinkNameID = obj.getLinkNameID;
            obj.LinkIndex = 1:obj.LinkCount;   %obj.getLinkIndex;
            obj.LinkPipeIndex = 1:obj.LinkPipeCount; %find(strcmp(obj.LinkType, 'PIPE'));
            obj.LinkPumpIndex = obj.LinkPipeCount+1:obj.LinkPipeCount+obj.LinkPumpCount; %find(strcmp(obj.LinkType, 'PUMP'));
            obj.LinkValveIndex = find(obj.LinkTypeIndex>2);
            obj.LinkPipeNameID = obj.LinkNameID(obj.LinkPipeIndex);
            obj.LinkPumpNameID = obj.LinkNameID(obj.LinkPumpIndex);
            obj.LinkValveNameID = obj.LinkNameID(obj.LinkValveIndex);
            %Get all the node data
            obj.NodeNameID = obj.getNodeNameID;
            tmp(:, 1) = obj.NodeNameID(obj.NodesConnectingLinksIndex(:, 1)');
            tmp(:, 2) = obj.NodeNameID(obj.NodesConnectingLinksIndex(:, 2)');
            obj.NodesConnectingLinksID=tmp;
            obj.NodeIndex = 1:obj.NodeCount; %obj.getNodeIndex;
            obj.NodeReservoirIndex = find(obj.NodeTypeIndex==1); %find(strcmp(obj.NodeType, 'RESERVOIR'));
            obj.NodeTankIndex = find(obj.NodeTypeIndex==2); %find(strcmp(obj.NodeType, 'TANK'));
            obj.NodeJunctionIndex = 1:obj.NodeJunctionCount; %find(strcmp(obj.NodeType, 'JUNCTION'));
            obj.NodeReservoirNameID=obj.NodeNameID(obj.NodeReservoirIndex);
            obj.NodeTankNameID=obj.NodeNameID(obj.NodeTankIndex);
            obj.NodeJunctionNameID=obj.NodeNameID(obj.NodeJunctionIndex);
            obj.NodePatternIndex=obj.getNodePatternIndex;
            obj.NodeBaseDemands = obj.getNodeBaseDemands;
            %Get all tank data
            obj.NodeTankInitialLevel = obj.getNodeTankInitialLevel;
            obj.NodeTankInitialWaterVolume = obj.getNodeTankInitialWaterVolume;
            obj.NodeTankMixingModelCode = obj.getNodeTankMixingModelCode;
            obj.NodeTankMixingModelType = obj.getNodeTankMixingModelType;
            obj.NodeTankMixZoneVolume = obj.getNodeTankMixZoneVolume;
            obj.NodeTankDiameter = obj.getNodeTankDiameter;
            obj.NodeTankMinimumWaterVolume = obj.getNodeTankMinimumWaterVolume;
            obj.NodeTankVolumeCurveIndex = obj.getNodeTankVolumeCurveIndex;
            obj.NodeTankMinimumWaterLevel = obj.getNodeTankMinimumWaterLevel;
            obj.NodeTankMaximumWaterLevel = obj.getNodeTankMaximumWaterLevel;
            obj.NodeTankMinimumFraction = obj.getNodeTankMixingFraction;
            obj.NodeTankBulkReactionCoeff = obj.getNodeTankBulkReactionCoeff;
            %Get all options
            obj.OptionsMaxTrials = obj.getOptionsMaxTrials;
            obj.OptionsAccuracyValue = obj.getOptionsAccuracyValue;
            obj.OptionsQualityTolerance = obj.getOptionsQualityTolerance;
            obj.OptionsEmitterExponent = obj.getOptionsEmitterExponent;
            obj.OptionsPatternDemandMultiplier = obj.getOptionsPatternDemandMultiplier;
            if obj.getVersion > 20101
                obj.OptionsHeadError = obj.getOptionsHeadError;
                obj.OptionsFlowChange = obj.getOptionsFlowChange;
                obj.OptionsHeadLossFormula = obj.getOptionsHeadLossFormula;
            end
            %Get pattern data
            obj.PatternNameID = obj.getPatternNameID;
            obj.PatternIndex = 1:obj.PatternCount; %obj.getPatternIndex;
            obj.PatternLengths = obj.getPatternLengths;
            obj.Pattern = obj.getPattern;
            %Get quality types
            obj.QualityCode = obj.getQualityCode;
            obj.QualityTraceNodeIndex = obj.getQualityTraceNodeIndex;
            %obj.QualityType = obj.getQualityType;
            obj.libFunctions = libfunctions(obj.LibEPANET);
            if sum(strcmp(obj.libFunctions, 'ENgetqualinfo'))
                n = obj.getQualityInfo;
                obj.QualityChemUnits = n.QualityChemUnits;
                obj.QualityChemName= n.QualityChemName;
            end
            %Get time parameters
            obj.TimeSimulationDuration = obj.getTimeSimulationDuration;
            obj.TimeHydraulicStep = obj.getTimeHydraulicStep;
            obj.TimeQualityStep = obj.getTimeQualityStep;
            obj.TimePatternStep = obj.getTimePatternStep;
            obj.TimePatternStart = obj.getTimePatternStart;
            obj.TimeReportingStep = obj.getTimeReportingStep;
            obj.TimeReportingStart = obj.getTimeReportingStart;
            obj.TimeRuleControlStep = obj.getTimeRuleControlStep;
            obj.TimeStatisticsIndex = obj.getTimeStatisticsIndex;
            obj.TimeStatisticsType = obj.getTimeStatisticsType;
            obj.TimeReportingPeriods = obj.getTimeReportingPeriods;
            % Get EPANET version
            obj.Version = obj.getVersion;
            try %EPANET Version 2.1.dll LibEPANET
                obj.TimeStartTime = obj.getTimeStartTime;
                obj.TimeHTime = obj.getTimeHTime;
                obj.TimeHaltFlag = obj.getTimeHaltFlag;
                obj.TimeNextEvent = obj.getTimeNextEvent;
                obj.NodeTankMaximumWaterVolume = obj.getNodeTankMaximumWaterVolume;
                obj.NodeBaseDemands = obj.getNodeBaseDemands;
                obj.NodeDemandCategoriesNumber = obj.getNodeDemandCategoriesNumber;
                obj.PatternAverageValue = obj.getPatternAverageValue;
                n = obj.getStatistic;
                obj.RelativeError = n.RelativeError;
                obj.Iterations = n.Iterations;
                obj.NodeDemandPatternNameID = obj.getNodeDemandPatternNameID;
                obj.NodeDemandPatternIndex = obj.getNodeDemandPatternIndex;
                obj.LinkPumpHeadCurveIndex = obj.getLinkPumpHeadCurveIndex;
                obj.LinkPumpPatternNameID = obj.getLinkPumpPatternNameID;
                obj.LinkPumpPatternIndex = obj.getLinkPumpPatternIndex;
                obj.LinkPumpTypeCode = obj.getLinkPumpTypeCode;
                obj.LinkPumpType = obj.getLinkPumpType;
                obj.LinkPumpPower = obj.getLinkPumpPower;
                obj.CurvesInfo = obj.getCurvesInfo;
            catch
            end
            %Get data from raw file (for information which cannot be
            %accessed by the epanet library)
            value=obj.getNodeCoordinates;
            %Get coordinates
            obj.NodeCoordinates{1} = value{1};
            obj.NodeCoordinates{2} = value{2};
            obj.NodeCoordinates{3} = value{3};
            obj.NodeCoordinates{4} = value{4};

            % US Customary - SI metric
            infoUnits = obj.getUnits;
            getFields_infoUnits = fields(infoUnits);
            for i=1:length(getFields_infoUnits)
                obj.(getFields_infoUnits{i}) = eval(['infoUnits.', getFields_infoUnits{i}]);
            end
            disp(['Input File "', varargin{1}, '" loaded sucessfuly.']);
        end % End of epanet class constructor
        function openAnyInp(obj, varargin)
            % Open as on matlab editor any EPANET input file using built
            % function open. Open current loaded input file (not temporary)
            % Example:
            %   d.openAnyInp;
            %   d.openAnyInp('Net2.inp');
            arg = obj.InputFile;
            if nargin > 1
                arg = varargin{1};
            end
            open(arg);
        end
        function openCurrentInp(obj, varargin)
            % Open EPANET input file who is loaded
            % Example:
            %   d.openCurrentInp;
            open(obj.TempInpFile);
        end
        function Errcode = loadEPANETFile(obj, varargin)
            % Load epanet file when use bin functions.
            % Example:
            %   d.loadEPANETFile(d.TempInpFile);
            obj.solve = 0;
            if nargin==2
                [Errcode] = ENopen(varargin{1}, [varargin{1}(1:end-4), '.txt'], [varargin{1}(1:end-4), '.bin'], obj.LibEPANET);
            else
                [Errcode] = ENopen(varargin{1}, varargin{2}, varargin{3}, obj.LibEPANET);
            end
        end
        function Errcode = runsCompleteSimulation(obj, varargin)
            % Runs a complete hydraulic and water simulation to create
            % binary & report files with name: [NETWORK_temp.txt], [NETWORK_temp.bin]
            % OR you can use argument to runs a complete simulation via ENepanet
            % Example:
            %       d.runsCompleteSimulation;
            %       d.runsCompleteSimulation('results'); % using ENepanet

            if nargin == 1
                obj.solveCompleteHydraulics;
                Errcode = obj.solveCompleteQuality;
                obj.writeReport();
            elseif nargin == 2
                rptfile = [varargin{1}, '.txt'];
                binfile = [varargin{1}, '.bin'];
                obj.Errcode = ENepanet(obj.LibEPANET, obj.BinTempfile, rptfile, binfile);
                Errcode = reloadNetwork(obj);
            end
        end
        function [value] = plot(obj, varargin)
            % Plots network in a new Matlab figure
            % Arguments:
            % 'nodes': yes/no
            % 'links': yes/no
            % 'line' : yes/no
            % 'point': yes/no
            % 'highlightnode': array of node IDs
            % 'highlightlink': array of link IDs
            % 'fontsize': number (px)
            % 'colornode': array of node IDs
            % 'colorlink': array of link id
            % 'axes': axes coordinates
            % 'linksindex': yes
            % 'nodesindex': yes
            % 'legend': show/hide
            % 'extend': yes/no
            % 'legendposition':
                %       'North'              inside plot box near top
                %       'South'              inside bottom
                %       'East'               inside right
                %       'West'               inside left
                %       'NorthEast'          inside top right (default for 2-D plots)
                %       'NorthWest'          inside top left
                %       'SouthEast'          inside bottom right
                %       'SouthWest'          inside bottom left
                %       'NorthOutside'       outside plot box near top
                %       'SouthOutside'       outside bottom
                %       'EastOutside'        outside right
                %       'WestOutside'        outside left
                %       'NorthEastOutside'   outside top right (default for 3-D plots)
                %       'NorthWestOutside'   outside top left
                %       'SouthEastOutside'   outside bottom right
                %       'SouthWestOutside'   outside bottom left
                %       'Best'               least conflict with data in plot
                %       'BestOutside'        least unused space outside plot
            % 'uifigure': app.UIFigure
            % Example:
            %   d.plot('nodes','yes','links','yes','highlightnode',{'10','11'},
            %'highlightlink',{'10'},'fontsize',8);
            %   d.plot('line','no');
            %   d.plot('point','no','linksindex','yes');
            %   d.plot('linksindex','yes','fontsize',8);
            %   d.plot('nodesindex','yes','fontsize',14);
            %   f=figure()
            %   d.plot('axes', h) % d.plot('axes', handles.axes) #gui
            %   d.plot('legend','hide')
            %   d.plot('legendposition','northwest')
            %   d = epanet('Net1.inp');
            %   nodeSet1={'10','11','22'};
            %   nodeSet2={'21','23','31'};
            %   colorNodeSet1={'r','r','r'};
            %   colorNodeSet2={'g','g','g'};
            %   linkSet1={'111','122','121'};
            %   linkSet2={'110','12','113'};
            %   colorLinkSet1={'k','k','k'};
            %   colorLinkSet2={'y','y','y'};
            %   d.plot('nodes','yes','links','yes','highlightnode',[nodeSet1 nodeSet2],'colornode',[colorNodeSet1 colorNodeSet2],...
            %   'highlightlink',[linkSet1 linkSet2],'colorlink',[colorLinkSet1 colorLinkSet2])
            %
            %   % Highlight multiple links with different colors #nodes
            %   linkSet1=d.getLinkNameID(1:4);
            %   colorLinkSet1=repmat({'r'},1,length(linkSet1));
            %   linkSet2=d.getLinkNameID([5,6,7,8]);
            %   colorLinkSet2=repmat({'g'},1,length(linkSet2));
            %   d.plot('highlightlink',[linkSet1 linkSet2],'colorlink',[colorLinkSet1 colorLinkSet2])
            [value] = plotnet(obj, 'bin', 0, varargin{:});
        end
        function value = getControls(obj, varargin)
            % Retrieves the parameters of all control statements.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.getControls              % Retrieves the parameters of all control statements
            %   d.getControls(1).Type      % Retrieves the type of the 1st control
            %   d.getControls(1).LinkID    % Retrieves the ID of the link associated with the 1st control
            %   d.getControls(1).Setting   % Retrieves the setting of the link associated with the 1st control
            %   d.getControls(1).NodeID    % Retrieves the ID of the node associated with the 1st control
            %   d.getControls(1).Value     % Retrieves the value of the node associated with the 1st control
            %   d.getControls(1).Control   % Retrieves the 1st control statement
            %
            % See also setControls, addControls, deleteControls,
            %          getRules, setRules, addRules, deleteRules.
            indices = getControlIndices(obj, varargin);j=1;
            value=struct();
            obj.ControlTypes={};
            for i=indices
                [obj.Errcode, obj.ControlTypesIndex(j), ...
                    obj.ControlLinkIndex(j), obj.ControlSettings(j), ...
                    obj.ControlNodeIndex(j), obj.ControlLevelValues(j)]...
                    = obj.ENgetcontrol(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                obj.ControlTypes(j)=obj.TYPECONTROL(obj.ControlTypesIndex(j)+1);
                value(j).Type = obj.ControlTypes{j};
                %value{i}.TypeIndex = obj.ControlTypesIndex(i);
                value(j).LinkID = obj.getLinkNameID{obj.ControlLinkIndex(j)};
                if all(obj.ControlSettings(j) ~= [0, 1])
                    value(j).Setting = obj.ControlSettings(j);
                else
                    value(j).Setting = obj.TYPESTATUS{obj.ControlSettings(j)+1};
                end
                if obj.ControlNodeIndex(j)
                    value(j).NodeID = obj.getNodeNameID{obj.ControlNodeIndex(j)};
                else
                    value(j).NodeID = NaN;
                end
                value(j).Value = obj.ControlLevelValues(j);
                switch obj.ControlTypes{j}
                    case 'LOWLEVEL'
                        value(j).Control = ['LINK ', value(j).LinkID, ' ', ...
                            value(j).Setting, ' IF NODE ', value(j).NodeID, ...
                            ' BELOW ', num2str(value(j).Value)];
                    case 'HIGHLEVEL'
                        value(j).Control = ['LINK ', value(j).LinkID, ' ', ...
                            value(j).Setting, ' IF NODE ', value(j).NodeID, ...
                            ' ABOVE ', num2str(value(j).Value)];
                    case 'TIMER'
                        value(j).Control = ['LINK ', value(j).LinkID, ' ', ...
                            num2str(value(j).Setting), ' AT TIME ', num2str(value(j).Value)];
                    case 'TIMEOFDAY'
                        value(j).Control = ['LINK ', value(j).LinkID, ' ', ...
                            num2str(value(j).Setting), ' AT CLOCKTIME ', num2str(value(j).Value)];
                end
                j=j+1;
            end
        end
        function value = getRules(obj, ruleIndex)
            % Retrieves the rule - based control statements. (EPANET Version 2.2)
            %
            % % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   rules = d.getRules
            %   rule_first_index = 1;
            %   rule_first = rules(rule_first_index)                   % Retrieves all the statements of the 1st rule - based control
            %   rule_second_index = 2;
            %   rule_second = rules(rule_second_index)                 % Retrieves all the statements of the 2nd rule - based control
            %
            % Example 2:
            %   rule_first = d.getRules(1)                             % Retrieves all the statements of the 1st rule - based control
            %   rule_first_ID = d.getRules(1).Rule_ID                  % Retrieves the ID of the 1st rule - based control
            %   rule_first_premises = d.getRules(1).Premises           % Retrieves all the premises of the 1st rule - based control
            %   rule_first_Then_Actions = d.getRules(1).Then_Actions   % Retrieves all the then actions of the 1st rule - based control
            %   rule_first_Else_Actions = d.getRules(1).Else_Actions   % Retrieves all the else actions of the 1st rule - based control
            %   rule_first_Rule = d.getRules(1).Rule                   % Retrieves the 1st rule - based control
            %
            % See also getRuleInfo, getRuleID, getRuleCount,
            %          setRules, deleteRules, addRules.
            value = struct();
            if nargin==1
                ruleIndex = 1:obj.getRuleCount;
            end
            for i=ruleIndex
                cnt = obj.getRuleInfo.Premises(i);
                premises = cell(cnt, 1);
                for j=1:obj.getRuleInfo.Premises(i)
                    [obj.Errcode, logop, object, objIndex, variable, relop, status, value_premise] =  ENgetpremise(i, j, obj.LibEPANET);
                    if j==1
                        logop = 1;
                    end
                    if object == obj.ToolkitConstants.EN_R_NODE
                        objectNameID = obj.getNodeNameID(objIndex);
                        space=' ';
                    elseif object == obj.ToolkitConstants.EN_R_LINK
                        objectNameID = obj.getLinkNameID(objIndex);
                        space=' ';
                    elseif object == obj.ToolkitConstants.EN_R_SYSTEM
                        objectNameID = ' ';
                        space = '';
                    end
                    if variable >= obj.ToolkitConstants.EN_R_TIME
                        value_premise = datestr((double(value_premise)/86400), 'HH:MM PM');
                    else
                        value_premise = num2str(value_premise);
                    end
                    if status==0
                        ruleStatus = '';
                    else
                        ruleStatus = char([obj.RULESTATUS{status} ' ']);
                        value_premise = '';
                    end
                    premises{j, 1} = [obj.LOGOP{logop}, ' ', obj.RULEOBJECT{object-5 }, space, char(objectNameID), space, obj.RULEVARIABLE{variable+1}, ' ', obj.RULEOPERATOR{relop+1}, ' ', ruleStatus, value_premise];
                    error(obj.getError(obj.Errcode));
                end
                cnt = obj.getRuleInfo.ThenActions(i);
                thenactions = cell(cnt, 1);
                for j=1:obj.getRuleInfo.ThenActions(i)
                    [obj.Errcode, linkIndex, status, setting] = ENgetthenaction(i, j, obj.LibEPANET);
                    if j==1
                        logop = 'THEN';
                    else
                        logop = 'AND';
                    end
                    link_type = char(obj.getLinkType(linkIndex));
                    linkNameID = char(obj.getLinkNameID(linkIndex));
                    if ismember(status, [1,2,3])
                        status = char([' STATUS IS ', char(obj.RULESTATUS{status})]);
                    else
                        status = '';
                    end
                    if setting>=0
                        setting = char([' SETTING IS ', num2str(setting)]);
                    else
                        setting = '';
                    end
                    thenactions{j, 1} = [logop, ' ', link_type, ' ', linkNameID, status, setting];
                    error(obj.getError(obj.Errcode));
                end
                cnt = obj.getRuleInfo.ElseActions(i);
                elseactions = cell(cnt, 1);
                for j=1:obj.getRuleInfo.ElseActions(i)
                    [obj.Errcode, linkIndex, status, setting] = ENgetelseaction(i, j, obj.LibEPANET);
                    if j==1
                        logop = 'ELSE';
                    else
                        logop = 'AND';
                    end
                    link_type = char(obj.getLinkType(linkIndex));
                    linkNameID = char(obj.getLinkNameID(linkIndex));
                    if (status==1) || (status == 2) || (status == 3)
                        status = char([' STATUS IS ', char(obj.RULESTATUS{status})]);
                    else
                        status = '';
                    end
                    if setting>=0
                        setting = char([' SETTING IS ', num2str(setting)]);
                    else
                        setting = '';
                    end
                    elseactions{j, 1} = [logop, ' ', link_type, ' ', linkNameID, status, setting];
                    error(obj.getError(obj.Errcode));
                end
                if nargin==1
                    k = i;
                elseif nargin==2
                    k = 1;
                end
                value(k).Rule_ID = obj.getRuleID{i};
                value(k).Premises=char(premises);
                value(k).Then_Actions=char(thenactions);
                value(k).Else_Actions=char(elseactions);
                value(k).Rule=char(['RULE ' obj.getRuleID{i}; premises; thenactions; elseactions; {['PRIORITY ' num2str(obj.getRuleInfo.Priority(i))]}]);
            end
        end
        function addRules(obj, rule)
            % Adds a new rule-based control to a project. (EPANET Version 2.2)
            %
            % Rule format: Following the format used in an EPANET input file.
            %              'RULE ruleid \n IF object objectid attribute relation attributevalue \n THEN object objectid STATUS/SETTING IS value \n PRIORITY value'
            %
            % See more: 'https://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf' (Page 164)
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.getRuleCount
            %   d.addRules('RULE RULE-1 \n IF TANK 2 LEVEL >= 140 \n THEN PUMP 9 STATUS IS CLOSED \n PRIORITY 1')
            %   d.getRuleCount
            %   d.getRules(1).Rule
            %
            % See also deleteRules, setRules, getRules, getRuleInfo,
            %          setRuleThenAction, setRuleElseAction, setRulePriority.

            %rule_new = split(rule, '\n ');
            rule_new = regexp(rule, '\\n', 'split');
            rule_final = [];
            for i=1:length(rule_new)
                rule_final = [rule_final rule_new{i} char(10)];
            end
            [obj.Errcode] = ENaddrule(rule_final, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRules(obj, ruleIndex, rule)
            % Sets a rule - based control. (EPANET Version 2.2)
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   rule = 'RULE RULE-1 \n IF NODE 2 LEVEL >= 140 \n THEN PIPE 10 STATUS IS CLOSED \n ELSE PIPE 10 STATUS IS OPEN \n PRIORITY 1';
            %   d.addRules(rule)                  % Adds a new rule - based control
            %   d.getRules(1).Rule                % Retrieves the 1st rule - based control
            %   ruleIndex = 1;
            %   rule_new = 'IF NODE 2 LEVEL > 150 \n THEN PIPE 10 STATUS IS OPEN \n ELSE PIPE 11 STATUS IS OPEN \n PRIORITY 2';
            %   d.setRules(ruleIndex, rule_new)   % Sets rule - based control
            %   d.getRules(1).Rule
            %
            % See also setRulePremise, setRuleThenAction, setRuleElseAction,
            %          getRules, addRules, deleteRules.
            rule_new = regexp(rule, '\\n ', 'split');
            i = 1;
            while strcmp(rule_new{i}(1:2), 'IF') || strcmp(rule_new{i}(1:3), 'AND') || strcmp(rule_new{i}(1:2), 'OR')
                obj.setRulePremise(ruleIndex, i, rule_new{i})
                i = i+1;
            end
            j = 1;
            while strcmp(rule_new{i}(1:4), 'THEN') || strcmp(rule_new{i}(1:3), 'AND')
                obj.setRuleThenAction(ruleIndex, j, rule_new{i})
                i = i+1;
                j = j+1;
            end
            j = 1;
            while strcmp(rule_new{i}(1:4), 'ELSE') || strcmp(rule_new{i}(1:3), 'AND')
                obj.setRuleElseAction(ruleIndex, j, rule_new{i})
                i = i+1;
                j = j+1;
            end
            if obj.getRuleInfo.Priority(ruleIndex) ~= 0
                obj.setRulePriority(ruleIndex, str2num(rule_new{i}(end)));
            end
        end
        function setRuleThenAction(obj, ruleIndex, actionIndex, then_action)
            % Sets rule - based control then actions. (EPANET Version 2.2)
            %
            % Input Arguments: Rule Index, Action Index, Then clause.
            %
            % See more: 'https://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf' (Page 164)
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.addRules('RULE RULE-1 \n IF TANK 2 LEVEL >= 140 \n THEN PIPE 10 STATUS IS CLOSED \n ELSE PIPE 10 STATUS IS OPEN \n PRIORITY 1')   % Adds a new rule - based control
            %   rule = d.getRules(1)   % Retrieves the 1st rule - based control
            %   ruleIndex = 1;
            %   actionIndex = 1;
            %   then_action = 'THEN PIPE 11 STATUS IS OPEN';
            %   d.setRuleThenAction(ruleIndex, actionIndex, then_action)
            %   rule = d.getRules(1)
            %
            % See also setRules, setRuleElseAction, setRulePriority,
            %          getRuleInfo, getRules, addRules, deleteRules.
            then_new = regexp(then_action, '\s', 'split');
            if strcmp(then_new{4},'STATUS')
                status = eval(['obj.ToolkitConstants.EN_R_IS_', then_new{6}]);
                setting = -1;
            elseif strcmp(then_new{4},'SETTING')
                status = -1;
                setting = str2num(then_new{6});
            end
            linkIndex = obj.getLinkIndex(then_new{3});
            [obj.Errcode] = ENsetthenaction(ruleIndex, actionIndex, linkIndex, status, setting, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRuleElseAction(obj, ruleIndex, actionIndex, else_action)
            % Sets rule - based control else actions. (EPANET Version 2.2)
            %
            % Input Arguments: Rule Index, Action Index, Link Index, Type, Value.   % Where Type = 'STATUS' or 'SETTING' and Value = the value of STATUS/SETTING
            %
            % See more: 'https://nepis.epa.gov/Adobe/PDF/P1007WWU.pdf' (Page 164)
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.addRules('RULE RULE-1 \n IF TANK 2 LEVEL >= 140 \n THEN PIPE 10 STATUS IS CLOSED \n ELSE PIPE 10 STATUS IS OPEN \n PRIORITY 1')   % Adds a new rule - based control
            %   rule = d.getRules(1)   % Retrieves the 1st rule - based control
            %   ruleIndex = 1;
            %   actionIndex = 1;
            %   else_action = 'ELSE PIPE 11 STATUS IS CLOSED';
            %   d.setRuleElseAction(ruleIndex, actionIndex, else_action)   % Sets the new else - action in the 1st rule - based control, in the 1st else - action.
            %   rule = d.getRules(1)
            %
            % See also setRules, setRuleThenAction, setRulePriority,
            %          getRuleInfo, getRules, addRules, deleteRules.
            else_new = regexp(else_action, '\s', 'split');
            if strcmp(else_new{4},'STATUS')
                status = eval(['obj.ToolkitConstants.EN_R_IS_', else_new{6}]);
                setting = -1;
            elseif strcmp(else_new{4},'SETTING')
                status = -1;
                setting = str2num(else_new{6});
            end
            linkIndex = obj.getLinkIndex(else_new{3});
            [obj.Errcode] = ENsetelseaction(ruleIndex, actionIndex, linkIndex, status, setting, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRulePriority(obj, ruleIndex, priority)
            % Sets rule - based control priority. (EPANET Version 2.2)
            %
            % % The example is based on d=epanet('BWSN_Network_1.inp');
            %
            % Example:
            %   d.getRules(1).Rule                       % Retrieves the 1st rule - based control
            %   ruleIndex = 1;
            %   priority = 2;
            %   d.setRulePriority(ruleIndex, priority)   % Sets the 1st rule - based control priority = 2
            %   d.getRules(1).Rule
            %
            % See also setRules, setRuleThenAction, setRuleElseAction,
            %          getRuleInfo, getRules, addRules, deleteRules.
            [obj.Errcode] = ENsetrulepriority(ruleIndex, priority, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRulePremise(obj, ruleIndex, premiseIndex, premise)
            % Sets the premise of a rule - based control. (EPANET Version 2.2)
            %
            % % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getRules(1).Premises                               % Retrieves the premise of the 1st rule
            %   ruleIndex = 1;
            %   premiseIndex = 1;
            %   premise = 'IF SYSTEM CLOCKTIME >= 8 PM';
            %   d.setRulePremise(ruleIndex, premiseIndex, premise)   % Sets the 1st premise of the 1st rule - based control
            %   d.getRules(1).Premises
            %
            % Example 2:
            %   d.getRules(1).Premises
            %   ruleIndex = 1;
            %   premiseIndex = 1;
            %   premise = 'IF NODE TANK-131 LEVEL > 20';
            %   d.setRulePremise(1, 1, premise)                      % Sets the 1st premise of the 1st rule - based control
            %   d.getRules(1).Premises
            %
            % See also setRulePremiseObejctNameID, setRulePremiseStatus, setRulePremiseValue,
            %          setRules, getRules, addRules, deleteRules.

            %  premise_new = split(premise, ' ');
            premise_new = regexp(premise, '\s', 'split');
            logop_code = ismember(obj.LOGOP, premise_new{1});
            logop = find(logop_code, 1);
            object = eval(['obj.ToolkitConstants.EN_R_', premise_new{2}]);
            if object == obj.ToolkitConstants.EN_R_NODE
                objIndex = obj.getNodeIndex(premise_new{3});
            elseif object == obj.ToolkitConstants.EN_R_LINK
                objIndex = obj.getLinkIndex(premise_new{3});
            elseif object == obj.ToolkitConstants.EN_R_SYSTEM
                objIndex = 0;
            end
            if object == obj.ToolkitConstants.EN_R_SYSTEM
                j = 3; k = 4; m = 5;
            else
                j = 4; k = 5; m = 6;
            end
            variable = eval(['obj.ToolkitConstants.EN_R_', premise_new{j}]);
            opearator_code = ismember(obj.RULEOPERATOR, premise_new{k});
            relop = find(opearator_code, 1) - 1;
            if variable == obj.ToolkitConstants.EN_R_STATUS
                value = -1;
                status = eval(['obj.ToolkitConstants.EN_R_IS_', premise_new{m}]);
            else
                value = str2double(premise_new{m});
                status = 0;
            end
            if object == obj.ToolkitConstants.EN_R_SYSTEM
               if strcmp(premise_new(6), 'AM')
                   value = value*3600;
               elseif strcmp(premise_new(6), 'PM')
                   value = value*3600 + 43200;
               end
            end
            [obj.Errcode] = ENsetpremise(ruleIndex, premiseIndex, logop, object, objIndex, variable, relop, status, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRulePremiseObejctNameID(obj, ruleIndex, premiseIndex, objNameID)
            % Sets the ID of an object in a premise of a rule-based control. (EPANET Version 2.2)
            %
            % % The example is based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getRules(1).Premises
            %   ruleIndex = 1;
            %   premiseIndex = 1;
            %   objNameID = 'TANK-131';
            %   d.setRulePremiseObejctNameID(ruleIndex, premiseIndex, objNameID)   % Sets the node's ID = 'TANK-131' to the 1st premise of the 1st rule - based control
            %   d.getRules(1).Premises
            %
            % See also setRulePremise, setRulePremiseStatus, setRulePremiseValue,
            %          setRules, getRules, addRules, deleteRules.
            [obj.Errcode, ~, object, objIndex, ~, ~, ~, ~] =  ENgetpremise(ruleIndex, premiseIndex, obj.LibEPANET);
            if object == obj.ToolkitConstants.EN_R_NODE
                objIndex = obj.getNodeIndex(objNameID);
            elseif object == obj.ToolkitConstants.EN_R_LINK
                objIndex = obj.getLinkIndex(objNameID);
            end
            [obj.Errcode] = ENsetpremiseindex(ruleIndex, premiseIndex, objIndex, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRulePremiseStatus(obj, ruleIndex, premiseIndex, status)
            % Sets the status being compared to in a premise of a rule-based control. (EPANET Version 2.2)
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.getRules
            %   d.addRules('RULE RULE-1 \n IF LINK 110 STATUS = CLOSED \n THEN PUMP 9 STATUS IS CLOSED \n PRIORITY 1')
            %   d.getRules(1)
            %   ruleIndex = 1;
            %   premiseIndex = 1;
            %   status = 'OPEN';
            %   d.setRulePremiseStatus(ruleIndex, premiseIndex, status)   % Sets the status = 'OPEN' to the 1st premise of the 1st rule - based control
            %   d.getRules(1).Premises
            %
            % See also setRulePremise, setRulePremiseObejctNameID, setRulePremiseValue,
            %          setRules, getRules, addRules, deleteRules.
            if strcmp(status,'OPEN')
                status_code = obj.ToolkitConstants.EN_R_IS_OPEN;
            elseif strcmp(status,'CLOSED')
                status_code = obj.ToolkitConstants.EN_R_IS_CLOSED;
            elseif strcmp(status,'ACTIVE')
                status_code = obj.ToolkitConstants.EN_R_IS_ACTIVE;
            end
            [obj.Errcode] = ENsetpremisestatus(ruleIndex, premiseIndex, status_code, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setRulePremiseValue(obj, ruleIndex, premiseIndex, value)
            % Sets the value being compared to in a premise of a rule-based control. (EPANET Version 2.2)
            %
            % % The example is based on d=epanet('BWSN_Network_1.inp');
            %
            % Example:
            %   d.getRules(1).Premises
            %   ruleIndex = 1;
            %   premiseIndex = 1;
            %   value = 20;
            %   d.setRulePremiseValue(ruleIndex, premiseIndex, value)   % Sets the value = 20 to the 1st premise of the 1st rule - based control
            %   d.getRules(1).Premises
            %
            % See also setRulePremise, setRulePremiseObejctNameID, setRulePremiseStatus,
            %          setRules, getRules, addRules, deleteRules.
            [obj.Errcode] = ENsetpremisevalue(ruleIndex, premiseIndex, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getRuleID(obj, varargin)
            % Retrieves the ID name of a rule-based control given its index. (EPANET Version 2.2)
            %
            % % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getRuleID         % Retrieves the ID name of every rule-based control
            %
            % Example 2:
            %   d.getRuleID(1)      % Retrieves the ID name of the 1st rule-based control
            %
            % Example 3:
            %   d.getRuleID(1:3)    % Retrieves the ID names of the 1st to 3rd rule-based control
            %
            % See also getRules, getRuleInfo, addRules.
            if nargin==1
                index = 1:obj.getRuleCount;
                value = cell(1, length(index));
            elseif nargin==2
                index = varargin{1};
                value = cell(1,1);
            end
            j=1;
            for i=1:length(index)
                [Errcode, value{j}] = ENgetruleID(index(i), obj.LibEPANET);
                error(obj.getError(Errcode));
                j=j+1;
            end
        end
        function value = getRuleInfo(obj, varargin)
            % Retrieves summary information about a rule-based control given it's index. (EPANET Version 2.2)
            %
            % % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getRuleInfo         % Retrieves summary information about every rule-based control
            %
            % Example 2:
            %   d.getRuleInfo(1)      % Retrieves summary information about the 1st rule-based control
            %
            % Example 3:
            %   d.getRuleInfo(1:3)    % Retrieves summary information about the 1st to 3rd rule-based control
            %
            % See also getRuleID, getRules, addRules.
            value=struct();
            if nargin==1
                index = 1:obj.getRuleCount;
            elseif nargin==2
                index = varargin{1};
            end
            value.Index = index;
            for i=1:length(index)
                [Errcode, value.Premises(i), value.ThenActions(i), value.ElseActions(i), value.Priority(i)] = ENgetrule(index(i), obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function Errcode = deleteRules(varargin)
            % Deletes an existing rule-based control given it's index. (EPANET Version 2.2)
            % Returns error code.
            %
            % % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getRuleCount        % Retrieves the number of rules
            %   d.deleteRules;        % Deletes all the rule-based control
            %   d.getRuleCount
            %
            % Example 2:
            %   d.deleteRules(1);     % Deletes the 1st rule-based control
            %   d.getRuleCount
            %
            % Example 3:
            %   d.deleteRules(1:3);   % Deletes the 1st to 3rd rule-based control
            %   d.getRuleCount
            %
            % See also addRules, getRules, setRules, getRuleCount.
            obj = varargin{1};
            if nargin==1
                index = 1:obj.getRuleCount;
            else
                index = varargin{2};
            end
            for i=length(index):-1:1
                [Errcode] = ENdeleterule(index(i), obj.LibEPANET);
                error(obj.getError(Errcode));
            end
        end
        function value = getNodeCount(obj)
            % Retrieves the number of nodes.
            %
            % Example:
            %   d.getNodeCount
            %
            % See also getNodeIndex, getLinkCount.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_NODECOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getNodeTankReservoirCount(obj)
            % Retrieves the number of tanks.
            %
            % Example:
            %   d.getNodeTankReservoirCount
            %
            % See also getNodeTankIndex, getNodeReservoirIndex.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_TANKCOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getLinkCount(obj)
            % Retrieves the number of links.
            %
            % Example:
            %   d.getLinkCount
            %
            % See also getLinkIndex, getNodeCount.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_LINKCOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getPatternCount(obj)
            % Retrieves the number of patterns.
            %
            % Example:
            %   d.getPatternCount
            %
            % See also getPatternIndex, getPattern.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_PATCOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getCurveCount(obj)
            % Retrieves the number of curves.
            %
            % Example:
            %   d.getCurveCount
            %
            % See also getCurveIndex, getCurvesInfo.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_CURVECOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getControlRulesCount(obj)
            % Retrieves the number of controls.
            %
            % Example:
            %   d.getControlRulesCount
            %
            % See also getControls, getRuleCount.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_CONTROLCOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getRuleCount(obj)
            % Retrieves the number of rules. (EPANET Version 2.2)
            %
            % Example:
            %   d.getRuleCount
            %
            % See also getRules, getControlRulesCount.
            [obj.Errcode, value] = obj.ENgetcount(obj.ToolkitConstants.EN_RULECOUNT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getNodeTankCount(obj)
            % Retrieves the number of Tanks.
            %
            % Example:
            %   d.getNodeTankCount
            %
            % See also getNodeReservoirCount, getNodeCount.
            value = sum(strcmp(obj.getNodeType, 'TANK'));
        end
        function value = getNodeReservoirCount(obj)
            % Retrieves the number of Reservoirs.
            %
            % Example:
            %   d.getNodeReservoirCount
            %
            % See also getNodeTankCount, getNodeCount.
            value = sum(strcmp(obj.getNodeType, 'RESERVOIR'));
        end
        function value = getNodeJunctionCount(obj)
            % Retrieves the number of junction nodes.
            %
            % Example:
            %   d.getNodeJunctionCount
            %
            % See also getNodeTankCount, getNodeCount.
            value = obj.getNodeCount - obj.getNodeTankReservoirCount;
        end
        function value = getLinkPipeCount(obj)
            % Retrieves the number of pipes.
            %
            % Example:
            %   d.getLinkPipeCount
            %
            % See also getLinkPumpCount, getLinkCount.
            LinkType1=obj.getLinkType;
            value = sum(strcmp(LinkType1, 'PIPE')+strcmp(LinkType1, 'CVPIPE'));
        end
        function value = getLinkPumpCount(obj)
            % Retrieves the number of pumps.
            %
            % Example:
            %   d.getLinkPumpCount
            %
            % See also getLinkPipeCount, getLinkCount.
            value = sum(strcmp(obj.getLinkType, 'PUMP'));
        end
        function value = getLinkValveCount(obj)
            % Retrieves the number of valves.
            %
            % Example:
            %   d.getLinkValveCount
            %
            % See also getLinkPumpCount, getLinkCount.
            LinkType1=obj.getLinkType;
            pipepump = sum(strcmp(LinkType1, 'PIPE')+strcmp(LinkType1, 'CVPIPE')+strcmp(LinkType1, 'PUMP'));
            value = obj.getLinkCount - pipepump;
        end
        function [errmssg, Errcode] = getError(obj, Errcode)
            % Retrieves the text of the message associated with a particular error or warning code.
            %
            % Example:
            %   error = 250;
            %   d.getError(error)
              [errmssg , Errcode] = obj.ENgeterror(Errcode, obj.LibEPANET);
            warning(errmssg);
            errmssg = '';
        end
        function value = getFlowUnits(obj)
            % Retrieves flow units used to express all flow rates.
            %
            % Example:
            %   d.getFlowUnits
            [obj.Errcode, flowunitsindex] = obj.ENgetflowunits(obj.LibEPANET);
            error(obj.getError(obj.Errcode));
            value=obj.TYPEUNITS{flowunitsindex+1};
        end
        function value = getLinkNameID(obj, varargin)
            % Retrieves the ID label(s) of all links, or the IDs of an index set of links.
            %
            % Example 1:
            %   d.getLinkNameID                % Retrieves the ID's of all links
            %
            % Example 2:
            %   linkIndex = 1;
            %   d.getLinkNameID(linkIndex)     % Retrieves the ID of the link with index = 1
            %
            % Example 3:
            %   linkIndices = 1:3;
            %   d.getLinkNameID(linkIndices)   % Retrieves the IDs of the links with indices = 1, 2, 3
            %
            % See also getNodeNameID, getLinkPipeNameID, getLinkIndex.
            if isempty(varargin)
                cnt = obj.getLinkCount;
                value=cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=obj.ENgetlinkid(i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            else
                k=1;
                value = cell(1, length(varargin{1}));
                if isempty(varargin{1}), varargin{1}=0; end
                for i=varargin{1}
                    [obj.Errcode, value{k}]=obj.ENgetlinkid(i, obj.LibEPANET);
                    if obj.Errcode==204, value=[];  return; end
                    k=k+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function value = getLinkPipeNameID(obj)
            % Retrieves the pipe ID.
            %
            % Example 1:
            %   d.getLinkPipeNameID        % Retrieves the ID's of all pipes
            %
            % Example 2:
            %   d.getLinkPipeNameID{1}     % Retrieves the ID of the 1st pipe
            %
            % Example 3:
            %   d.getLinkPipeNameID{1:3}   % Retrieves the ID of the first 3 pipes
            %
            % See also getLinkNameID, getLinkPumpNameID, getNodeNameID.
            value=obj.getLinkNameID(obj.getLinkPipeIndex);
        end
        function value = getLinkPumpNameID(obj)
            % Retrieves the pump ID.
            %
            % Example 1:
            %   d.getLinkPumpNameID        % Retrieves the ID's of all pumps
            %
            % Example 2:
            %   d.getLinkPumpNameID{1}     % Retrieves the ID of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpNameID{1:2}   % Retrieves the ID of the first 2 pumps
            %
            % See also getLinkNameID, getLinkPipeNameID, getNodeNameID.
            value=obj.getLinkNameID(obj.getLinkPumpIndex);
        end
        function value = getLinkValveNameID(obj)
            % Retrieves the valve ID.
            %
            % Example 1:
            %   d.getLinkValveNameID        % Retrieves the ID's of all valves
            %
            % Example 2:
            %   d.getLinkValveNameID{1}     % Retrieves the ID of the 1st valve
            %
            % Example 3:
            %   d.getLinkValveNameID{1:3}   % Retrieves the ID of the first 3 valves
            %
            % See also getLinkNameID, getLinkPumpNameID, getNodeNameID.
            value=obj.getLinkNameID(obj.getLinkValveIndex);
        end
        function value = getLinkIndex(obj, varargin)
            % Retrieves the indices of all links, or the indices of an ID set of links.
            %
            % Example 1:
            %   d.getLinkIndex           % Retrieves the indices of all links
            %
            % Example 2:
            %   linkID = d.getLinkNameID{1};
            %   d.getLinkIndex(linkID)   % Retrieves the index of the 1st link given it's ID
            %
            % Example 3:
            %   linkID = d.getLinkNameID(1:3);
            %   d.getLinkIndex(linkID)   % Retrieves the indices of the first 3 links given their ID
            %
            % See also getLinkNameID, getLinkPipeIndex, getNodeIndex.
            value = [];
            if isempty(varargin)
                value=1:obj.getLinkCount;
            elseif isempty(varargin{1})
                   return;
            elseif isa(varargin{1}, 'cell')
                k=1;
                value = zeros(1, length(varargin{1}));
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = obj.ENgetlinkindex(varargin{1}{j}, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = obj.ENgetlinkindex(varargin{1}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getLinkPipeIndex(obj)
            % Retrieves the pipe indices.
            %
            % Example:
            %   d.getLinkPipeIndex
            %
            % See also getLinkIndex, getLinkPumpIndex.
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes, 'PIPE'));
        end
        function value = getLinkPumpIndex(obj, varargin)
            % Retrieves the pump indices.
            %
            % Example 1:
            %   d.getLinkPumpIndex        % Retrieves the indices of all pumps
            %
            % Example 2:
            %   d.getLinkPumpIndex(1)     % Retrieves the index of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpIndex(1:2)   % Retrieves the indices of the first 2 pumps
            %
            % See also getLinkIndex, getLinkPipeIndex, getLinkValveIndex.
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes, 'PUMP'));
            if ~isempty(varargin)
                value = value(varargin{1});
            end
        end
        function value = getLinkValveIndex(obj)
            % Retrieves the valve indices.
            %
            % Example:
            %   d.getLinkValveIndex
            %
            % See also getLinkIndex, getLinkPipeIndex, getLinkPumpIndex.
            value = obj.getLinkPipeCount+obj.getLinkPumpCount+1:obj.getLinkCount;
        end
        function value = getNodesConnectingLinksIndex(obj)
            % Retrieves the indexes of the from/to nodes of all links.
            % Duplicate function with getLinkNodesIndex for new version.
            %
            % Example:
            %   d.getNodesConnectingLinksIndex
            %
            % See also getLinkNodesIndex, getNodesConnectingLinksID.
            value = obj.getLinkNodesIndex;
        end
        function value = getLinkNodesIndex(obj, varargin)
            % Retrieves the indexes of the from/to nodes of all links.
            %
            % Example 1:
            %   d.getLinkNodesIndex
            %
            % Example 2:
            %   d.getLinkNodesIndex(2)  % Link index
            %
            % See also getNodesConnectingLinksID.
            if nargin == 2
                cntL = varargin{1};
                indices = cntL;
            else
                cntL = obj.getLinkCount;
                indices = 1:cntL;
                value(cntL, 1:2)=[nan nan];
            end
            j=1;
            for i=indices
                [obj.Errcode, linkFromNode, linkToNode] = obj.ENgetlinknodes(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                value(j, :)= [linkFromNode, linkToNode];
                j = j +1;
            end
        end
        function value = getNodesConnectingLinksID(obj)
            % Retrieves the id of the from/to nodes of all links.
            %
            % Example 1:
            %   d.getNodesConnectingLinksID           % Retrieves the id of the from/to nodes of all links
            %
            % Example 2:
            %   linkIndex = 1;
            %   d.getNodesConnectingLinksID{1, 1:2}   % Retrieves the id of the from/to nodes of the 1st link
            %
            % See also getLinkNodesIndex.
            value={};
            obj.NodesConnectingLinksIndex=obj.getLinkNodesIndex;
            if obj.getLinkCount
                value(:, 1)=obj.getNodeNameID(obj.NodesConnectingLinksIndex(:, 1)');
                value(:, 2)=obj.getNodeNameID(obj.NodesConnectingLinksIndex(:, 2)');
            end
        end
        function value = getLinkType(obj, varargin)
            % Retrieves the link-type code for all links.
            %
            % Example 1:
            %    d.getLinkType      % Retrieves the link-type code for all links
            %
            % Example 2:
            %    d.getLinkType(1)   % Retrieves the link-type code for the first link
            %
            % See also getLinkTypeIndex, getLinksInfo, getLinkDiameter,
            %          getLinkLength, getLinkRoughnessCoeff, getLinkMinorLossCoeff.
            if ~isempty(varargin), varargin=varargin{1}; end
            value=obj.TYPELINK(obj.getLinkTypeIndex(varargin)+1);
        end
        function value = getLinkTypeIndex(obj, varargin)
            % Retrieves the link-type index for all links.
            %
            % Example 1:
            %    d.getLinkTypeIndex      % Retrieves the link-type index for all links
            %
            % Example 2:
            %    d.getLinkTypeIndex(1)   % Retrieves the link-type index for the first link
            %
            % See also getLinkType, getLinksInfo, getLinkDiameter,
            %          getLinkLength, getLinkRoughnessCoeff, getLinkMinorLossCoeff.
            [indices, value] = getLinkIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = obj.ENgetlinktype(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function [value] = getLinksInfo(obj)
            % Retrieves all link info.
            %
            % Example:
            %    d.getLinksInfo
            %
            % See also getLinkType, getLinkTypeIndex, getLinkDiameter,
            %          getLinkLength, getLinkRoughnessCoeff, getLinkMinorLossCoeff.
            value = struct();
            for i=1:obj.getLinkCount
                [~, value.LinkDiameter(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_DIAMETER, obj.LibEPANET);
                [~, value.LinkLength(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_LENGTH, obj.LibEPANET);
                [~, value.LinkRoughnessCoeff(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_ROUGHNESS, obj.LibEPANET);
                [~, value.LinkMinorLossCoeff(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_MINORLOSS, obj.LibEPANET);
                [~, value.LinkInitialStatus(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_INITSTATUS, obj.LibEPANET);
                [~, value.LinkInitialSetting(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_INITSETTING, obj.LibEPANET);
                [~, value.LinkBulkReactionCoeff(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_KBULK, obj.LibEPANET);
                [~, value.LinkWallReactionCoeff(i)] = obj.ENgetlinkvalue(i, obj.ToolkitConstants.EN_KWALL, obj.LibEPANET);
                [~, value.NodesConnectingLinksIndex(i, 1), value.NodesConnectingLinksIndex(i, 2)] = obj.ENgetlinknodes(i, obj.LibEPANET);
                if obj.getVersion > 20101
                    [~, value.LinkTypeIndex(i)] = obj.ENgetlinktype(i, obj.LibEPANET);
                end
            end
        end
        function value = getLinkDiameter(obj, varargin)
            % Retrieves the value of link diameters.
            % Pipe/valve diameter
            %
            % Example 1:
            %    d.getLinkDiameter      % Retrieves the value of all link diameters
            %
            % Example 2:
            %    d.getLinkDiameter(1)   % Retrieves the value of the first link diameter
            %
            % See also getLinkType, getLinksInfo, getLinkLength,
            %          getLinkRoughnessCoeff, getLinkMinorLossCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_DIAMETER, varargin{:});
        end
        function value = getLinkLength(obj, varargin)
            % Retrieves the value of link lengths.
            % Pipe length
            %
            % Example 1:
            %    d.getLinkLength      % Retrieves the value of all link lengths
            %
            % Example 2:
            %    d.getLinkLength(1)   % Retrieves the value of the first link length
            %
            % See also getLinkType, getLinksInfo, getLinkDiameter,
            %          getLinkRoughnessCoeff, getLinkMinorLossCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_LENGTH, varargin{:});
        end
        function value = getLinkRoughnessCoeff(obj, varargin)
            % Retrieves the value of all link roughness coefficients.
            % Pipe roughness coefficient
            %
            % Example 1:
            %    d.getLinkRoughnessCoeff      % Retrieves the value of all link roughness coefficients
            %
            % Example 2:
            %    d.getLinkRoughnessCoeff(1)   % Retrieves the value of the first link roughness coefficient
            %
            % See also getLinkType, getLinksInfo, getLinkDiameter,
            %          getLinkLength, getLinkMinorLossCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_ROUGHNESS, varargin{:});
        end
        function value = getLinkMinorLossCoeff(obj, varargin)
            % Retrieves the value of link minor loss coefficients.
            % Pipe/valve minor loss coefficient
            %
            % Example 1:
            %   d.getLinkMinorLossCoeff      % Retrieves the value of all link minor loss coefficients
            %
            % Example 2:
            %   d.getLinkMinorLossCoeff(1)   % Retrieves the value of the first link minor loss coefficient
            %
            % See also getLinkType, getLinksInfo, getLinkDiameter,
            %          getLinkLength, getLinkRoughnessCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_MINORLOSS, varargin{:});
        end
        function value = getLinkInitialStatus(obj, varargin)
            % Retrieves the value of all link initial status.
            % Initial status (see @ref EN_LinkStatusType)
            %
            % Example 1:
            %   d.getLinkInitialStatus      % Retrieves the value of all link initial status
            %
            % Example 2:
            %   d.getLinkInitialStatus(1)   % Retrieves the value of the first link initial status
            %
            % See also getLinkType, getLinksInfo, getLinkInitialSetting,
            %          getLinkBulkReactionCoeff, getLinkWallReactionCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_INITSTATUS, varargin{:});
        end
        function value = getLinkInitialSetting(obj, varargin)
            % Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves.
            %
            % Example 1:
            %   d.getLinkInitialSetting      % Retrieves the value of all link initial settings
            %
            % Example 2:
            %   d.getLinkInitialSetting(1)   % Retrieves the value of the first link initial setting
            %
            % See also getLinkType, getLinksInfo, getLinkInitialStatus,
            %          getLinkBulkReactionCoeff, getLinkWallReactionCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_INITSETTING, varargin{:});
        end
        function value = getLinkBulkReactionCoeff(obj, varargin)
            % Retrieves the value of all link bulk chemical reaction coefficient.
            %
            % Example 1:
            %   d.getLinkBulkReactionCoeff      % Retrieves the value of all link bulk chemical reaction coefficient
            %
            % Example 2:
            %   d.getLinkBulkReactionCoeff(1)   % Retrieves the value of the first link bulk chemical reaction coefficient
            %
            % See also getLinkType, getLinksInfo, getLinkRoughnessCoeff,
            %          getLinkMinorLossCoeff, getLinkInitialStatus,
            %          getLinkInitialSetting, getLinkWallReactionCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_KBULK, varargin{:});
        end
        function value = getLinkWallReactionCoeff(obj, varargin)
            % Retrieves the value of all pipe wall chemical reaction coefficient.
            %
            % Example 1:
            %   d.getLinkWallReactionCoeff      % Retrieves the value of all pipe wall chemical reaction coefficient
            %
            % Example 2:
            %   d.getLinkWallReactionCoeff(1)   % Retrieves the value of the first pipe wall chemical reaction coefficient
            %
            % See also getLinkType, getLinksInfo, getLinkRoughnessCoeff,
            %          getLinkMinorLossCoeff, getLinkInitialStatus,
            %          getLinkInitialSetting, getLinkBulkReactionCoeff.
            value = get_link_info(obj, obj.ToolkitConstants.EN_KWALL, varargin{:});
        end
        function value = getLinkFlows(obj, varargin)
            % Retrieves the current computed flow rate (read only).
            % Using step-by-step hydraulic analysis
            %
            % Example 1:
            %    d.getLinkFlows      % Retrieves the current computed flow rate for all links
            %
            % Example 2:
            %    d.getLinkFlows(1)   % Retrieves the current computed flow rate for the first link
            %
            % Example 3:
            %    % Hydraulic analysis step-by-step
            %    d.openHydraulicAnalysis;
            %    d.initializeHydraulicAnalysis;
            %    tstep=1;P=[];T_H=[];D=[];H=[];F=[];S=[];
            %    while (tstep>0)
            %       t= d.runHydraulicAnalysis;
            %       P = [P; d.getNodePressure];
            %       D = [D; d.getNodeActualDemand];
            %       H = [H; d.getNodeHydaulicHead];
            %       S = [S; d.getLinkStatus];
            %       F = [F; d.getLinkFlows];
            %       T_H = [T_H; t];
            %       tstep=d.nextHydraulicAnalysisStep;
            %    end
            %    d.closeHydraulicAnalysis;
            %
            % Example 4:
            %    % Hydraulic and Quality analysis step-by-step
            %    d.openHydraulicAnalysis;
            %    d.openQualityAnalysis;
            %    d.initializeHydraulicAnalysis(0);
            %    d.initializeQualityAnalysis(d.ToolkitConstants.EN_NOSAVE);
            %
            %    tstep = 1;
            %    T = []; P = []; F = []; QN = []; QL = [];
            %    while (tstep>0)
            %       t  = d.runHydraulicAnalysis;
            %       qt = d.runQualityAnalysis;
            %
            %       P  = [P; d.getNodePressure];
            %       F  = [F; d.getLinkFlows];
            %
            %       QN = [QN; d.getNodeActualQuality];
            %       QL = [QL; d.getLinkActualQuality];
            %       T  = [T; t];
            %
            %       tstep = d.nextHydraulicAnalysisStep;
            %       qtstep = d.nextQualityAnalysisStep;
            %    end
            %    d.closeQualityAnalysis;
            %    d.closeHydraulicAnalysis;
            %
            % See also getLinkVelocity, getLinkHeadloss, getLinkStatus,
            %          getLinkPumpState, getLinkSettings, getLinkEnergy,
            %          getLinkActualQuality, getLinkPumpEfficiency.

            value = get_link_info(obj, obj.ToolkitConstants.EN_FLOW, varargin{:});
        end
        function value = getLinkVelocity(obj, varargin)
            % Retrieves the current computed flow velocity (read only).
            %
            % Using step-by-step hydraulic analysis
            %
            % Example 1:
            %   d.getLinkVelocity      % Retrieves the current computed flow velocity for all links
            %
            % Example 2:
            %   d.getLinkVelocity(1)   % Retrieves the current computed flow velocity for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkHeadloss, getLinkStatus,
            %          getLinkPumpState, getLinkSettings, getLinkActualQuality.
            value = get_link_info(obj, obj.ToolkitConstants.EN_VELOCITY, varargin{:});
        end
        function value = getLinkHeadloss(obj, varargin)
            % Retrieves the current computed head loss (read only).
            %
            % Using step-by-step hydraulic analysis,
            %
            % Example 1:
            %   d.getLinkHeadloss      % Retrieves the current computed head loss for all links
            %
            % Example 2:
            %   d.getLinkHeadloss(1)   % Retrieves the current computed head loss for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkVelocity, getLinkStatus,
            %          getLinkPumpState, getLinkSettings, getLinkActualQuality.
            value = get_link_info(obj, obj.ToolkitConstants.EN_HEADLOSS, varargin{:});
        end
        function value = getLinkStatus(obj, varargin)
            % Retrieves the current link status (see @ref EN_LinkStatusType) (0 = closed, 1 = open).
            %
            % Using step-by-step hydraulic analysis,
            %
            % Example 1:
            %   d.getLinkStatus      % Retrieves the current link status for all links
            %
            % Example 2:
            %   d.getLinkStatus(1)   % Retrieves the current link status for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkVelocity, getLinkHeadloss,
            %          getLinkPumpState, getLinkSettings.
            value = get_link_info(obj, obj.ToolkitConstants.EN_STATUS, varargin{:});
        end
        function value = getLinkPumpState(obj, varargin)
            % Retrieves the current computed pump state (read only) (see @ref EN_PumpStateType). (EPANET Version 2.2)
            % same as status: open, active, closed
            % Using step-by-step hydraulic analysis,
            %
            % Example 1:
            %   d.getLinkPumpState      % Retrieves the current computed pump state for all links
            %
            % Example 2:
            %   d.getLinkPumpState(1)   % Retrieves the current computed pump state for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkHeadloss, getLinkStatus,
            %          getLinkSettings, getLinkEnergy, getLinkPumpEfficiency.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_STATE, varargin);
        end
        function value = getLinkSettings(obj, varargin)
            % Retrieves the current computed value of all link roughness for pipes
            % or actual speed for pumps or actual setting for valves.
            %
            % Using step-by-step hydraulic analysis,
            %
            % Example 1:
            %   d.getLinkSettings      % Retrieves the current values of settings for all links
            %
            % Example 2:
            %   d.getLinkSettings(1)   % Retrieves the current value of setting for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkVelocity, getLinkHeadloss,
            %          getLinkStatus, getLinkPumpState, getLinkEnergy.
            value = get_link_info(obj, obj.ToolkitConstants.EN_SETTING, varargin{:});
        end
        function value = getLinkEnergy(obj, varargin)
            % Retrieves the current computed pump energy usage (read only).
            %
            % Using step-by-step hydraulic analysis,
            %
            % Example 1:
            %   d.getLinkEnergy      % Retrieves the current computed pump energy usage for all links
            %
            % Example 2:
            %   d.getLinkEnergy(1)   % Retrieves the current computed pump energy usage for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkVelocity, getLinkHeadloss,
            %          getLinkStatus, getLinkPumpState, getLinkPumpEfficiency.
            value = get_link_info(obj, obj.ToolkitConstants.EN_ENERGY, varargin{:});
        end
        function value = getLinkActualQuality(obj, varargin)
            % Retrieves the current computed link quality (read only). (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkActualQuality      % Retrieves the current computed link quality for all links
            %
            % Example 2:
            %   d.getLinkActualQuality(1)   % Retrieves the current computed link quality for the first link
            %
            % Example 3:
            %   check examples/EX14_hydraulic_and_quality_analysis.m
            %
            % See also getLinkFlows, getLinkStatus, getLinkPumpState,
            %          getLinkSettings, getLinkPumpEfficiency.
            value = get_link_info(obj, obj.ToolkitConstants.EN_LINKQUAL, varargin{:});
        end
        function value = getLinkPumpEfficiency(obj, varargin)
            % Retrieves the current computed pump efficiency (read only). (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpEfficiency      % Retrieves the current computed pump efficiency for all links
            %
            % Example 2:
            %   d.getLinkPumpEfficiency(1)   % Retrieves the current computed pump efficiency for the first link
            %
            % For more, you can type `help getLinkFlows` and check examples 3 & 4
            %
            % See also getLinkFlows, getLinkStatus, getLinkPumpState,
            %          getLinkSettings, getLinkEnergy, getLinkActualQuality.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_EFFIC, varargin);
        end
        function value = getLinkPumpPower(obj, varargin)
            % Retrieves the pump constant power rating (read only). (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpPower              % Retrieves the constant power rating of all pumps
            %
            % Example 2:
            %   d.getLinkPumpPower(1)           % Retrieves the constant power rating of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpPower(1:2)         % Retrieves the constant power rating of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpPower(pumpIndex)   % Retrieves the constant power rating of the pumps given their indices
            %
            % See also getLinkPumpHCurve, getLinkPumpECurve,getLinkPumpECost,
            %          getLinkPumpEPat, getLinkPumpPatternIndex, getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_POWER, varargin);
        end
        function value = getLinkPumpHCurve(obj, varargin)
            % Retrieves the pump head v. flow curve index. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpHCurve              % Retrieves the head v. flow curve index of all pumps
            %
            % Example 2:
            %   d.getLinkPumpHCurve(1)           % Retrieves the head v. flow curve index of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpHCurve(1:2)         % Retrieves the head v. flow curve index of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpHCurve(pumpIndex)   % Retrieves the head v. flow curve index of the pumps given their indices
            %
            % See also setLinkPumpHCurve, getLinkPumpECurve,getLinkPumpECost,
            %          getLinkPumpEPat, getLinkPumpPatternIndex, getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_HCURVE, varargin);
        end
        function value = getLinkPumpECurve(obj, varargin)
            % Retrieves the pump efficiency v. flow curve index. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpECurve              % Retrieves the efficiency v. flow curve index of all pumps
            %
            % Example 2:
            %   d.getLinkPumpECurve(1)           % Retrieves the efficiency v. flow curve index of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpECurve(1:2)         % Retrieves the efficiency v. flow curve index of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpECurve(pumpIndex)   % Retrieves the efficiency v. flow curve index of the pumps given their indices
            %
            % See also setLinkPumpECurve, getLinkPumpHCurve, getLinkPumpECost,
            %          getLinkPumpEPat, getLinkPumpPatternIndex, getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_ECURVE, varargin);
        end
        function value = getLinkPumpECost(obj, varargin)
            % Retrieves the pump average energy price. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpECost              % Retrieves the average energy price of all pumps
            %
            % Example 2:
            %   d.getLinkPumpECost(1)           % Retrieves the average energy price of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpECost(1:2)         % Retrieves the average energy price of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpECost(pumpIndex)   % Retrieves the average energy price of the pumps given their indices
            %
            % See also setLinkPumpECost, getLinkPumpPower, getLinkPumpHCurve,
            %          getLinkPumpEPat, getLinkPumpPatternIndex, getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_ECOST, varargin);
        end
        function value = getLinkPumpEPat(obj, varargin)
            % Retrieves the pump energy price time pattern index. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getLinkPumpEPat              % Retrieves the energy price time pattern index of all pumps
            %
            % Example 2:
            %   d.getLinkPumpEPat(1)           % Retrieves the energy price time pattern index of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpEPat(1:2)         % Retrieves the energy price time pattern index of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpEPat(pumpIndex)   % Retrieves the energy price time pattern index of the pumps given their indices
            %
            % See also setLinkPumpEPat, getLinkPumpHCurve, getLinkPumpECurve,
            %          getLinkPumpECost, getLinkPumpPatternIndex, getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_PUMP_EPAT, varargin);
        end
        function value = getLinkPumpPatternIndex(obj, varargin)
            % Retrieves the pump speed time pattern index. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getLinkPumpPatternIndex              % Retrieves the speed time pattern index of all pumps
            %
            % Example 2:
            %   d.getLinkPumpPatternIndex(1)           % Retrieves the speed time pattern index of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpPatternIndex(1:2)         % Retrieves the speed time pattern index of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpPatternIndex(pumpIndex)   % Retrieves the speed time pattern index of the pumps given their indices
            %
            % See also setLinkPumpPatternIndex, getLinkPumpPower, getLinkPumpHCurve,
            %          getLinkPumpECost, getLinkPumpEPat,  getLinkPumpPatternNameID.
            value = get_node_link(obj, 'pump', 'ENgetlinkvalue', obj.ToolkitConstants.EN_LINKPATTERN, varargin);
        end
        function value = getLinkPumpPatternNameID(obj, varargin)
            % Retrieves pump pattern name ID. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getLinkPumpPatternNameID              % Retrieves the pattern name ID of all pumps
            %
            % Example 2:
            %   d.getLinkPumpPatternNameID(1)           % Retrieves the pattern name ID of the 1st pump
            %
            % Example 3:
            %   d.getLinkPumpPatternNameID(1:2)         % Retrieves the pattern name ID of the first 2 pumps
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.getLinkPumpPatternNameID(pumpIndex)   % Retrieves the pattern name ID of the pumps given their indices
            %
            % See also getLinkPumpPower, getLinkPumpHCurve, getLinkPumpECurve,
            %          getLinkPumpECost, getLinkPumpEPat, getLinkPumpPatternIndex.
            if isempty(varargin)
                patindices = obj.getLinkPumpPatternIndex;
            else
                patindices = obj.getLinkPumpPatternIndex(varargin{1});
            end
            value = obj.getPatternNameID(patindices);
        end
        function value = getNodeNameID(obj, varargin)
            % Retrieves the ID label of all nodes or some nodes with a specified index.
            %
            % Example 1:
            %   d.getNodeNameID                  % Retrieves the ID label of all nodes
            %
            % Example 2:
            %   d.getNodeNameID(1)               % Retrieves the ID label of the first node
            %
            % Example 3:
            %   junctionIndex = d.getNodeJunctionIndex;
            %   d.getNodeNameID(junctionIndex)   % Retrieves the ID labels of all junctions give their indices
            %
            % See also getNodeReservoirNameID, getNodeJunctionNameID,
            %          getNodeIndex, getNodeType, getNodesInfo.
            if isempty(varargin)
                cnt = obj.getNodeCount;
                value = cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=obj.ENgetnodeid(i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            else
                if isempty(varargin{1}), varargin{1}=0; end
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value{k}]=obj.ENgetnodeid(i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            end
        end
        function value = getNodeReservoirNameID(obj)
            % Retrieves the reservoir ID label.
            %
            % Example :
            %   d.getNodeReservoirNameID
            %
            % See also getNodeNameID, getNodeJunctionNameID, getNodeIndex,
            %          getNodeReservoirIndex, getNodeType, getNodesInfo.
            %
            value=obj.getNodeNameID(obj.getNodeReservoirIndex);
        end
        function value = getNodeJunctionNameID(obj)
            % Retrieves the junction ID label.
            %
            % Example:
            %   d.getNodeJunctionNameID
            %
            % See also getNodeNameID, getNodeReservoirNameID, getNodeIndex,
            %          getNodeJunctionIndex, getNodeType, getNodesInfo.
            value=obj.getNodeNameID(obj.getNodeJunctionIndex);
        end
        function value = getNodeIndex(obj, varargin)
            % Retrieves the indices of all nodes or some nodes with a specified ID.
            %
            % Example 1:
            %   d.getNodeIndex           % Retrieves the indices of all nodes
            %
            % Example 2:
            %   nameID = d.getNodeNameID(1)
            %   d.getNodeIndex(nameID)   % Retrieves the node index given the ID label of the 1st node
            %
            % See also getNodeNameID, getNodeReservoirIndex, getNodeJunctionIndex,
            %          getNodeType, getNodeTypeIndex, getNodesInfo.
            value = [];
            if isempty(varargin)
                value=1:obj.getNodeCount;
            elseif isempty(varargin{1})
                    return;
            elseif isa(varargin{1}, 'cell')
                k=1;
                value = zeros(1, length(varargin{1}));
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = obj.ENgetnodeindex(varargin{1}{j}, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = obj.ENgetnodeindex(varargin{1}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getNodeReservoirIndex(obj)
            % Retrieves the indices of reservoirs.
            %
            % Example:
            %   d.getNodeReservoirIndex
            %
            % See also getNodeNameID, getNodeIndex, getNodeJunctionIndex,
            %          getNodeType, getNodeTypeIndex, getNodesInfo.
            tmpNodeTypes=obj.getNodeType;
            value = find(strcmp(tmpNodeTypes, 'RESERVOIR'));
        end
        function value = getNodeJunctionIndex(obj)
            % Retrieves the indices of junctions.
            %
            % Example:
            %   d.getNodeJunctionIndex
            %
            % See also getNodeNameID, getNodeIndex, getNodeReservoirIndex,
            %          getNodeType, getNodeTypeIndex, getNodesInfo.
            tmpNodeTypes=obj.getNodeType;
            value = find(strcmp(tmpNodeTypes, 'JUNCTION'));
        end
        function value = getNodeType(obj, varargin)
            % Retrieves the node-type code for all nodes.
            %
            % Example 1:
            %   d.getNodeType      % Retrieves the node-type code for all nodes
            %
            % Example 2:
            %   d.getNodeType(1)   % Retrieves the node-type code for the first node
            %
            % See also getNodeNameID, getNodeIndex,
            %          getNodeTypeIndex, getNodesInfo.
            indices = getNodeIndices(obj, varargin);
            value=obj.TYPENODE(obj.getNodeTypeIndex(indices)+1);
        end
        function value = getNodeTypeIndex(obj, varargin)
            % Retrieves the node-type code index for all nodes.
            %
            % Code meaning:
            %   0 = Junction
            %   1 = Tank
            %   2 = Reservoir
            %
            % Example 1:
            %   d.getNodeTypeIndex      % Retrieves the node-type code index for all nodes
            %
            % Example 2:
            %   d.getNodeTypeIndex(1)   % Retrieves the node-type code index for the first node
            %
            % See also getNodeNameID, getNodeIndex,
            %          getNodeType, getNodesInfo.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = obj.ENgetnodetype(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodesInfo(obj)
            % Retrieves nodes info (elevations, demand patterns, emmitter coeff, initial quality,
            % source quality, source pattern index, source type index, node type index).
            %
            % Example:
            %   d.getNodesInfo
            %
            % See also getNodeElevations, getNodeDemandPatternIndex, getNodeEmitterCoeff,
            %          getNodeInitialQuality, NodeTypeIndex.
            value = struct();
            for i=1:obj.getNodeCount
                [~, value.NodeElevations(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_ELEVATION, obj.LibEPANET);
                [~, value.NodeDemandPatternIndex(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_PATTERN, obj.LibEPANET);
                [~, value.NodeEmitterCoeff(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_EMITTER, obj.LibEPANET);
                [~, value.NodeInitialQuality(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_INITQUAL, obj.LibEPANET);
                [~, value.NodeSourceQuality(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCEQUAL, obj.LibEPANET);
                [~, value.NodeSourcePatternIndex(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCEPAT, obj.LibEPANET);
                [~, value.NodeSourceTypeIndex(i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCETYPE, obj.LibEPANET);
                [~, value.NodeTypeIndex(i)] = obj.ENgetnodetype(i, obj.LibEPANET);
            end
        end
        function value = getNodeElevations(obj, varargin)
            % Retrieves the value of all node elevations.
            %
            % Example 1:
            %   d.getNodeElevations        % Retrieves the value of all node elevations
            %
            % Example 2:
            %   d.getNodeElevations(1)     % Retrieves the value of the first node elevation
            %
            % Example 3:
            %   d.getNodeElevations(5:7)   % Retrieves the value of the 5th to 7th node elevations
            %
            % See also setNodeElevations, getNodesInfo, getNodeNameID,
            %          getNodeType,getNodeEmitterCoeff, getNodeInitialQuality.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_ELEVATION, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getDemandModel(obj, varargin)
            % Retrieves the type of demand model in use and its parameters.
            %
            % Example:
            %   model = d.getDemandModel()
            %
            % See also setDemandModel, getNodeBaseDemands, getNodeDemandCategoriesNumber
            %          getNodeDemandPatternIndex, getNodeDemandPatternNameID.
            [obj.Errcode, value.DemandModelCode, value.DemandModelPmin, value.DemandModelPreq, value.DemandModelPexp] = ENgetdemandmodel(obj.LibEPANET);
            error(obj.getError(obj.Errcode));
            value.DemandModelType = obj.DEMANDMODEL(value.DemandModelCode+1);
        end
        function categoryIndex = addNodeJunctionDemand(obj, varargin)
            % Adds a new demand to a junction given the junction index, base demand, demand time pattern and demand category name. (EPANET Version 2.2)
            % Returns the values of the new demand category index.
            % A blank string can be used for demand time pattern and demand name category to indicate
            % that no time pattern or category name is associated with the demand.
            %
            % Example 1:
            %   % New demand added with the name 'new demand' to the 1st node, with 100 base demand, using the 1st time pattern.
            %   d.addNodeJunctionDemand(1, 100, '1', 'new demand')
            %   d.getNodeJunctionDemandIndex     % Retrieves the indices of all demands for all nodes.
            %   d.getNodeJunctionDemandName{2}   % Retrieves the demand category names of the 2nd demand index for all nodes.
            %
            % Example 2:
            %   % New demands added with the name 'new demand' to the 1st and 2nd node, with 100 base demand, using the 1st time pattern.
            %   d.addNodeJunctionDemand([1, 2], 100, '1', 'new demand')
            %   d.getNodeJunctionDemandIndex     % Retrieves the indices of all demands for all nodes.
            %   d.getNodeJunctionDemandName{2}   % Retrieves the demand category names of the 2nd demand index for all nodes.
            %
            % Example 3:
            %   % New demands added with the name 'new demand' to the 1st and 2nd node, with 100 and 110 base demand respectively, using the 1st time pattern.
            %   d.addNodeJunctionDemand([1, 2], [100, 110], '1', 'new demand')
            %   d.getNodeJunctionDemandIndex     % Retrieves the indices of all demands for all nodes.
            %   d.getNodeJunctionDemandName{2}   % Retrieves the demand category names of the 2nd demand index for all nodes.
            %
            % Example 4:
            %   % New demands added with the name 'new demand' to the 1st and 2nd node, with 100 and 110 base demand respectively, using the 1st time pattern.
            %   d.addNodeJunctionDemand([1, 2], [100, 110], {'1', '1'}, 'new demand')
            %   d.getNodeJunctionDemandIndex     % Retrieves the indices of all demands for all nodes.
            %   d.getNodeJunctionDemandName{2}   % Retrieves the demand category names of the 2nd demand index for all nodes.
            %
            % Example 5:
            %   % New demands added with the names 'new demand1' and 'new demand2' to the 1st and 2nd node, with 100 and 110 base demand respectively, using the 1st and 2nd(if exists) time pattern respectively.
            %   d.addNodeJunctionDemand([1, 2], [100, 110], {'1', '2'}, {'new demand1', 'new demand2'})
            %   d.getNodeJunctionDemandIndex     % Retrieves the indices of all demands for all nodes.
            %   d.getNodeJunctionDemandName{2}   % Retrieves the demand category names of the 2nd demand index for all nodes.
            %
            % See also deleteNodeJunctionDemand, getNodeJunctionDemandIndex, getNodeJunctionDemandName,
            %          setNodeJunctionDemandName, getNodeBaseDemands.
            nodeIndex = varargin{1};
            baseDemand = varargin{2};
            if nargin==3
                demandPattern = '';
                demandName = '';
            elseif nargin==4
                demandPattern = varargin{3};
                demandName = '';
            elseif nargin==5
                demandPattern = varargin{3};
                demandName = varargin{4};
            end
            if isscalar(nodeIndex)
                [obj.Errcode]=ENadddemand(nodeIndex, baseDemand , demandPattern, demandName, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            elseif ~isscalar(nodeIndex)&&  isscalar(baseDemand) && ~iscell(demandPattern) && ~iscell(demandName)
                for i=1:length(nodeIndex)
                    [obj.Errcode]=ENadddemand(nodeIndex(i), baseDemand , demandPattern, demandName, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            elseif ~isscalar(nodeIndex)&& ~isscalar(baseDemand) && ~iscell(demandPattern) && ~iscell(demandName)
                for i=1:length(nodeIndex)
                    [obj.Errcode]=ENadddemand(nodeIndex(i), baseDemand(i) , demandPattern, demandName, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            elseif ~isscalar(nodeIndex) &&  ~isscalar(baseDemand) && iscell(demandPattern) && ~iscell(demandName)
                for i=1:length(nodeIndex)
                    [obj.Errcode]=ENadddemand(nodeIndex(i), baseDemand(i) , demandPattern{i}, demandName, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            elseif ~isscalar(nodeIndex) &&  ~isscalar(baseDemand) && iscell(demandPattern) && iscell(demandName)
                for i=1:length(nodeIndex)
                    [obj.Errcode]=ENadddemand(nodeIndex(i), baseDemand(i) , demandPattern{i}, demandName{i}, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            end
            if ~isscalar(nodeIndex) && ~iscell(demandName)
                demandNameCell=cell(1,length(nodeIndex));
                for j=1:length(nodeIndex)
                    demandNameCell{j}=demandName;
                end
                demandName=demandNameCell;
            end
            categoryIndex = obj.getNodeJunctionDemandIndex(nodeIndex,demandName);
        end
        function Errcode = deleteNodeJunctionDemand(obj, varargin)
            % Deletes a demand from a junction given the junction index and demand index. (EPANET Version 2.2)
            % Returns the remaining(if exist) node demand indices.
            %
            % Example 1:
            %   nodeIndex = 1;
            %   baseDemand = 100;
            %   patternId = '1';
            %   categoryIndex = 1;
            %   d.getNodeJunctionDemandIndex(nodeIndex)                                                   % Retrieves the indices of all demands for the 1st node
            %   d.getNodeJunctionDemandName                                                               % Retrieves the names of all nodes demand category
            %   d.getNodeJunctionDemandName{categoryIndex}(nodeIndex)                                     % Retrieves the name of the 1st demand category of the 1st node
            %   categoryIndex = d.addNodeJunctionDemand(nodeIndex, baseDemand, patternId, 'new demand')   % Adds a new demand to the 1st node and returns the new demand index
            %   d.getNodeJunctionDemandIndex(nodeIndex)                                                   % Retrieves the indices of all demands for the 1st node
            %   d.getNodeJunctionDemandName                                                               % Retrieves the names of all nodes demand category
            %   d.getNodeJunctionDemandName{categoryIndex}(nodeIndex)                                     % Retrieves the name of the 2nd demand category of the 1st node
            %   d.deleteNodeJunctionDemand(1,2)                                                           % Deletes the 2nd demand of the 1st node
            %   d.getNodeJunctionDemandIndex(nodeIndex)                                                   % Retrieves the indices of all demands for the 1st node
            %
            % Example 2:
            %   nodeIndex = 1;
            %   baseDemand = 100;
            %   patternId = '1';
            %   categoryIndex_2 = d.addNodeJunctionDemand(nodeIndex, baseDemand, patternId, 'new demand_2')   % Adds a new demand to the first node and returns the new demand index
            %   categoryIndex_3 = d.addNodeJunctionDemand(nodeIndex, baseDemand, patternId, 'new demand_3')   % Adds a new demand to the first node and returns the new demand index
            %   d.getNodeJunctionDemandName{categoryIndex_2}(nodeIndex)                                       % Retrieves the name of the 2nd demand category of the 1st node
            %   d.deleteNodeJunctionDemand(1)                                                                 % Deletes all the demands of the 1st node
            %   d.getNodeJunctionDemandIndex(nodeIndex)                                                       % Retrieves the indices of all demands for the 1st node
            %
            % Example 3:
            %   nodeIndex = [1, 2, 3];
            %   baseDemand = [100, 110, 150];
            %   patternId = {'1', '1',''};
            %   categoryIndex = d.addNodeJunctionDemand(nodeIndex, baseDemand, ...
            %   patternId, {'new demand_1', 'new demand_2', 'new demand_3'})            % Adds 3 new demands to the first 3 nodes
            %   d.getNodeJunctionDemandName{2}(nodeIndex)
            %   d.deleteNodeJunctionDemand(1:3)                                         % Deletes all the demands of the first 3 nodes
            %   d.getNodeJunctionDemandIndex(nodeIndex)                                 % Retrieves the indices of all demands for the first 3 nodes
            %
            % See also addNodeJunctionDemand, getNodeJunctionDemandIndex, getNodeJunctionDemandName,
            %          setNodeJunctionDemandName, getNodeBaseDemands.
            nodeIndex = varargin{1};
            if nargin==2
                if isscalar(nodeIndex)
                    numDemand=size(obj.getNodeJunctionDemandIndex);
                    for i=1:numDemand(1)
                        [obj.Errcode]=ENdeletedemand(nodeIndex, 1, obj.LibEPANET);
                        error(obj.getError(obj.Errcode));
                    end
                elseif ~isscalar(nodeIndex)
                    for j=1:length(nodeIndex)
                        numDemand=size(obj.getNodeJunctionDemandIndex);
                        for i=1:numDemand(1)
                            [obj.Errcode]=ENdeletedemand(nodeIndex(j), 1, obj.LibEPANET);
                            error(obj.getError(obj.Errcode));
                        end
                    end
                end
            elseif nargin==3
                [obj.Errcode]=ENdeletedemand(nodeIndex, varargin{2}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getNodeJunctionDemandName(obj, varargin)
            % Gets the name of a node's demand category.
            %
            % Example:
            %   model = d.getNodeJunctionDemandName()
            %
            % See also setNodeJunctionDemandName, getNodeBaseDemands,
            %          getNodeDemandCategoriesNumber, getNodeDemandPatternNameID.
            [indices, ~] = getNodeJunctionIndices(obj, varargin);
            numdemands = obj.getNodeDemandCategoriesNumber(indices);
            value = cell(1, max(numdemands));
            cnt = length(indices);
            val = cell(max(numdemands), cnt);j=1;
            for i=indices
                v=1;
                for u=1:numdemands(j)
                    [obj.Errcode, val{v, j}] = ENgetdemandname(i, u, obj.LibEPANET);v=v+1;
                end
                j=j+1;
            end
            for i=1:size(val, 1)
                value{i} = val(i, :);
            end
        end
        function value = getNodeComment(obj, varargin)
            % Retrieves the comment string assigned to the node object.
            %
            % Example 1:
            %   d.getNodeComment        % Retrieves the comment string assigned to all node objects
            %
            % Example 2:
            %   d.getNodeComment(4)     % Retrieves the comment string assigned to the 4th node object
            %
            % Example 3:
            %   d.getNodeComment(1:5)   % Retrieves the comment string assigned to the 1st to 5th node object
            %
            % See also setNodeComment, getNodesInfo,
            %          getNodeNameID, getNodeType.
            if isempty(varargin)
                cnt = obj.getNodeCount;
                value = cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=ENgetcomment(obj.ToolkitConstants.EN_NODE, i, obj.LibEPANET);
                end
            else
                if isempty(varargin{1}), varargin{1}=0; end
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value{k}]=ENgetcomment(obj.ToolkitConstants.EN_NODE, i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            end
        end
        function value = setNodeComment(obj, value, varargin)
            % Sets the comment string assigned to the node object.
            %
            % Example 1:
            %   d.setNodeComment(1, 'This is a node');                    % Sets a comment to the 1st node
            %   d.getNodeComment(1)
            %
            % Example 2:
            %   d.setNodeComment(1:2, {'This is a node', 'Test comm'});   % Sets a comment to the 1st and 2nd node
            %   d.getNodeComment(1:2)
            %
            % See also getNodeComment, getNodesInfo,
            %          setNodeNameID, setNodeCoordinates.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            if length(indices) == 1
                [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_NODE, indices, value, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                for i=indices
                    [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_NODE, i, value{j}, obj.LibEPANET); j=j+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function [Line1, Line2, Line3] = getTitle(obj, varargin)
            % Retrieves the title lines of the project.
            %
            % Example:
            %   [Line1, Line2, Line3] = d.getTitle()    % Retrieves the three title lines of the project
            %
            % See also setTitle.
            [obj.Errcode, Line1, Line2, Line3] = ENgettitle(obj.LibEPANET);
        end
        function value = getNodeBaseDemands(obj, varargin)
            % Retrieves the value of all node base demands.
            %
            % Example 1:
            %   d.getNodeBaseDemands
            %   d.getNodeBaseDemands{1}   % Get categories 1
            %
            % Example 2:
            %   d.getNodeBaseDemands(2)   % Get node base demand with categories for specific node index
            %
            % See also setNodeBaseDemands, getNodeDemandCategoriesNumber,
            %          getNodeDemandPatternIndex, getNodeDemandPatternNameID.
            [indices, ~] = getNodeIndices(obj, varargin);
            numdemands = obj.getNodeDemandCategoriesNumber(indices);
            value = cell(1, max(numdemands));
            cnt = length(indices);
            val = zeros(max(numdemands), cnt);j=1;
            for i=indices
                v=1;
                for u=1:numdemands(j)
                    [obj.Errcode, val(v, j)] = ENgetbasedemand(i, u, obj.LibEPANET);v=v+1;
                end
                j=j+1;
            end
            for i=1:size(val, 1)
                value{i} = val(i, :);
            end
        end
        function value = getNodeDemandCategoriesNumber(obj, varargin)
            % Retrieves the value of all node base demands categorie number. (EPANET Version 2.1)
            %
            % Example 1:
            %	d.getNodeDemandCategoriesNumber      % Retrieves the value of all node base demands categorie number
            %
            % Example 2:
            %	d.getNodeDemandCategoriesNumber(1)   % Retrieves the value of the first node base demand categorie number
            %
            % See also getNodeBaseDemands, getNodeDemandPatternIndex,
            %          getNodeDemandPatternNameID.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = obj.ENgetnumdemands(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeDemandPatternIndex(obj)
            % Retrieves the value of all node base demands pattern index. (EPANET Version 2.1)
            %
            % Example:
            %   d.getNodeDemandPatternIndex
            %   d.getNodeDemandPatternIndex{1}
            %
            % See also getNodeBaseDemands, getNodeDemandCategoriesNumber,
            %          getNodeDemandPatternNameID,
            %          setNodeDemandPatternIndex.
            numdemands = obj.getNodeDemandCategoriesNumber;
            value = cell(1, max(numdemands));
            val = zeros(max(numdemands), obj.getNodeCount);
            for i=obj.getNodeIndex
                v=1;
                for u=1:numdemands(i)
                    [obj.Errcode, val(v, i)] = ENgetdemandpattern(i, u, obj.LibEPANET);
                    v = v+1;
                end
            end
            for i=obj.getNodeReservoirIndex
                if ismember(i, obj.getNodeReservoirIndex)
                    [obj.Errcode, val(v, i)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_PATTERN, obj.LibEPANET);
                end
            end
            for i=1:size(val, 1)
                value{i} = val(i, :);
            end
        end
        function value = getNodeDemandPatternNameID(obj, varargin)
            % Retrieves the value of all node base demands pattern name ID. (EPANET Version 2.1)
            %
            % Example:
            %   d.getNodeDemandPatternNameID
            %   d.getNodeDemandPatternNameID{1}
            %
            % See also getNodeBaseDemands, getNodeDemandCategoriesNumber,
            %          getNodeDemandPatternIndex.
            v = obj.getNodeDemandPatternIndex;
            m = obj.getPatternNameID;
            if ~isempty(varargin)
                numdemands = obj.getNodeDemandCategoriesNumber(varargin{1});
            else
                numdemands = obj.getNodeDemandCategoriesNumber;
            end
            indices = getNodeIndices(obj, varargin);
            value = cell(1, max(numdemands));
            val = cell(max(numdemands), obj.getNodeCount);
            for i=indices
                for u=1:numdemands(i)
                    if v{u}(i)~=0
                        val{u, i}= char(m(v{u}(i)));
                    else
                        val{u, i}= '';
                    end
                end
                if numdemands(i)==0
                    val{1, i}= [];
                end
            end
            for i=1:size(val, 1)
                value{i} = val(i, :);
            end
        end
        function value = getNodeJunctionDemandIndex(obj, varargin)
            % Retrieves the demand index of the junctions. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getNodeJunctionDemandIndex         % Retrieves the demand index of all junctions
            %
            % Example 2:
            %   d.getNodeJunctionDemandIndex(1,'')   % Retrieves the demand index of the 1st junction given it's name (i.e. '')
            %
            % Example 3:
            %   d.getNodeJunctionDemandIndex(1:3)    % Retrieves the demand index of the first 3 junctions
            %
            % Example 4:
            %   % Adds two new demands and retrieves the two new demand indices.
            %   d.addNodeJunctionDemand([1, 2], [100, 110], {'1', '1'}, {'new demand1', 'new demand2'});
            %   d.getNodeJunctionDemandIndex([1,2],{'new demand1','new demand2'})
            %
            % See also getNodeJunctionDemandName, getNodeJunctionIndex, getNodeJunctionNameID,
            %          addNodeJunctionDemand, deleteNodeJunctionDemand, getNodeJunctionCount.
            if nargin==3
                nodeIndex = varargin{1};
                demandName = varargin{2};
                if isscalar(nodeIndex) && ~iscell(demandName)
                    [obj.Errcode, value] = obj.ENgetdemandindex(nodeIndex, demandName, obj.LibEPANET);
                elseif ~isscalar(nodeIndex) && iscell(demandName)
                    value=zeros(1,length(nodeIndex));
                    for i=1:length(nodeIndex)
                        [obj.Errcode, value(i)] = obj.ENgetdemandindex(nodeIndex(i), demandName{i}, obj.LibEPANET);
                    end
                end
            elseif nargin==2
                nodeIndex = varargin{1};
                demandName = obj.getNodeJunctionDemandName;
                if isscalar(nodeIndex)
                    value = zeros(length(demandName), 1);
                    for i=1:length(demandName)
                        demandNameIn = demandName{i};
                        [obj.Errcode, value(i)] = obj.ENgetdemandindex(nodeIndex, demandNameIn{varargin{1}}, obj.LibEPANET);
                    end
                else
                    value = zeros(length(demandName),length(nodeIndex));
                    for i=1:length(demandName)
                        demandNameIn = demandName{i};
                        for j=1:length(nodeIndex)
                            [obj.Errcode, value(i,j)] = obj.ENgetdemandindex(nodeIndex(j), demandNameIn{nodeIndex(j)}, obj.LibEPANET);
                        end
                    end
                end
            elseif nargin==1
				demandName = obj.getNodeJunctionDemandName;
				[indices, ~] = getNodeJunctionIndices(obj, varargin);
                value = zeros(length(demandName),length(indices));
                for i=1:length(demandName)
                    for j=1:length(demandName{i})
                        demandNameIn = demandName{i}{j};
                        [obj.Errcode, value(i,j)] = obj.ENgetdemandindex(j, demandNameIn, obj.LibEPANET);
                    end
                end
            else
                error(obj.getError(250))
            end
        end
        function value = getStatistic(obj)
            % Returns error code. (EPANET Version 2.1)
            %
            % Example:
            %   Input:  none
            %   Output: *iter = # of iterations to reach solution
            %           *relerr = convergence error in solution
            [obj.Errcode, value.Iterations] = obj.ENgetstatistic(obj.ToolkitConstants.EN_ITERATIONS, obj.LibEPANET);
            [obj.Errcode, value.RelativeError] = obj.ENgetstatistic(obj.ToolkitConstants.EN_RELATIVEERROR, obj.LibEPANET);
            [obj.Errcode, value.DeficientNodes] = obj.ENgetstatistic(obj.ToolkitConstants.EN_DEFICIENTNODES, obj.LibEPANET);
            [obj.Errcode, value.DemandReduction] = obj.ENgetstatistic(obj.ToolkitConstants.EN_DEMANDREDUCTION, obj.LibEPANET);
        end
        function value = getNodePatternIndex(obj, varargin)
            % Retrieves the value of all node demand pattern indices.
            %
            % Example 1:
            %   d.getNodePatternIndex      % Retrieves the value of all node demand pattern indices
            %
            % Example 2:
            %   d.getNodePatternIndex(1)   % Retrieves the value of the first node demand pattern index
            %
            % See also getNodeBaseDemands, getNodeDemandCategoriesNumber,
            %          getNodeDemandPatternIndex, getNodeDemandPatternNameID.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_PATTERN, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeEmitterCoeff(obj, varargin)
            % Retrieves the value of all node emmitter coefficients.
            %
            % Example 1:
            %   d.getNodeEmitterCoeff      % Retrieves the value of all node emmitter coefficients
            %
            % Example 2:
            %   d.getNodeEmitterCoeff(1)   % Retrieves the value of the first node emmitter coefficient
            %
            % See also setNodeEmitterCoeff, getNodesInfo, getNodeElevations.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_EMITTER, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeInitialQuality(obj, varargin)
            % Retrieves the value of all node initial quality.
            %
            % Example 1:
            %   d.getNodeInitialQuality      % Retrieves the value of all node initial quality
            %
            % Example 2:
            %   d.getNodeInitialQuality(1)   % Retrieves the value of the first node initial quality
            %
            % See also setNodeInitialQuality, getNodesInfo, getNodeSourceQuality.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_INITQUAL, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeSourceQuality(obj, varargin)
            % Retrieves the value of all node source quality.
            %
            % Example 1:
            %   d.getNodeSourceQuality      % Retrieves the value of all node source quality
            %
            % Example 2:
            %   d.getNodeSourceQuality(1)   % Retrieves the value of the first node source quality
            %
            % See also setNodeSourceQuality, getNodeInitialQuality, getNodeSourcePatternIndex,
            %          getNodeSourceTypeIndex, getNodeSourceType.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCEQUAL, obj.LibEPANET);
                if isnan(value(j)), value(j)=0; end
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end
                j=j+1;
            end
        end
        function value = getNodeSourcePatternIndex(obj, varargin)
            % Retrieves the value of all node source pattern index.
            %
            % Example 1:
            %   d.getNodeSourcePatternIndex      % Retrieves the value of all node source pattern index
            %
            % Example 2:
            %   d.getNodeSourcePatternIndex(1)   % Retrieves the value of the first node source pattern index
            %
            % See also setNodeSourcePatternIndex, getNodeSourceQuality,
            %          getNodeSourceTypeIndex, getNodeSourceType.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCEPAT, obj.LibEPANET);
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end
                j=j+1;
            end
        end
        function value = getNodeSourceTypeIndex(obj, varargin)
            % Retrieves the value of all node source type index.
            %
            % Example 1:
            %   d.getNodeSourceTypeIndex      % Retrieves the value of all node source type index
            %
            % Example 2:
            %   d.getNodeSourceTypeIndex(1)   % Retrieves the value of the first node source type index
            %
            % See also getNodeSourceQuality, getNodeSourcePatternIndex, getNodeSourceType.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCETYPE, obj.LibEPANET);
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end
                j=j+1;
            end
        end
        function value = getNodeSourceType(obj, varargin)
            % Retrieves the value of all node source type.
            %
            % Example 1:
            %   d.getNodeSourceType      % Retrieves the value of all node source type
            %
            % Example 2:
            %   d.getNodeSourceType(1)   % Retrieves the value of the first node source type
            %
            % See also setNodeSourceType, getNodeSourceQuality,
            %          getNodeSourcePatternIndex, getNodeSourceTypeIndex.
            [indices, ~] = getNodeIndices(obj,varargin);j=1;
            value = cell(1, length(indices));
            for i=indices
                [obj.Errcode, temp] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCETYPE, obj.LibEPANET);
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end
                if ~isnan(temp)
                    value{j}=obj.TYPESOURCE(temp+1);
                else
                    value{j}=[];
                end
                j=j+1;
            end
        end
        function value = getNodeTankInitialLevel(obj, varargin)
            % Retrieves the value of all tank initial water levels.
            %
            % Example 1:
            %   d.getNodeTankInitialLevel       % Retrieves the value of all tank initial water levels
            %
            % Example 2:
            %   d.getNodeTankInitialLevel(11)   % Retrieves the value of the eleventh node(tank) water level
            %
            % See also setNodeTankInitialLevel, getNodeTankInitialWaterVolume, getNodeTankVolume,
            %          getNodeTankMaximumWaterLevel, getNodeTankMinimumWaterLevel.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_TANKLEVEL, varargin);
        end
        function value = getNodeActualDemand(obj, varargin)
            % Retrieves the computed value of all node actual demands.
            %
            % Example 1:
            %   d.getNodeActualDemand      % Retrieves the computed value of all node actual demands
            %
            % Example 2:
            %   d.getNodeActualDemand(1)   % Retrieves the computed value of the first node actual demand
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemandSensingNodes, getNodeHydaulicHead, getNodePressure,
            %          getNodeActualQuality, getNodeMassFlowRate, getNodeActualQualitySensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_DEMAND, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeActualDemandSensingNodes(obj, varargin)
            % Retrieves the computed demand values at some sensing nodes.
            %
            % Example:
            %   d.getNodeActualDemandSensingNodes(1)   % Retrieves the computed demand value of the first sensing node
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemand, getNodeHydaulicHead, getNodePressure,
            %          getNodeActualQuality, getNodeMassFlowRate, getNodeActualQualitySensingNodes.
            value=zeros(1, length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.Errcode, value(v)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_DEMAND, obj.LibEPANET);v=v+1;
            end
        end
        function value = getNodeHydaulicHead(obj, varargin)
            % Retrieves the computed values of all node hydraulic heads.
            %
            % Example 1:
            %   d.getNodeHydaulicHead      % Retrieves the computed value of all node hydraulic heads
            %
            % Example 2:
            %   d.getNodeHydaulicHead(1)   % Retrieves the computed value of the first node hydraulic head
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemand, getNodeActualDemandSensingNodes, getNodePressure,
            %          getNodeActualQuality, getNodeMassFlowRate, getNodeActualQualitySensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_HEAD, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodePressure(obj, varargin)
            % Retrieves the computed values of all node pressures.
            %
            % Example 1:
            %    d.getNodePressure      % Retrieves the computed values of all node pressures
            %
            % Example 2:
            %    d.getNodePressure(1)   % Retrieves the computed value of the first node pressure
            %
            % Example 3:
            %    % Hydraulic analysis step-by-step
            %    d.openHydraulicAnalysis;
            %    d.initializeHydraulicAnalysis;
            %    tstep=1;P=[];T_H=[];D=[];H=[];F=[];S=[];
            %    while (tstep>0)
            %       t= d.runHydraulicAnalysis;
            %       P = [P; d.getNodePressure];
            %       D = [D; d.getNodeActualDemand];
            %       H = [H; d.getNodeHydaulicHead];
            %       S = [S; d.getLinkStatus];
            %       F = [F; d.getLinkFlows];
            %       T_H = [T_H; t];
            %       tstep=d.nextHydraulicAnalysisStep;
            %    end
            %    d.closeHydraulicAnalysis;
            %
            % Example 4:
            %    % Hydraulic and Quality analysis step-by-step
            %    d.openHydraulicAnalysis;
            %    d.openQualityAnalysis;
            %    d.initializeHydraulicAnalysis(0);
            %    d.initializeQualityAnalysis(d.ToolkitConstants.EN_NOSAVE);
            %
            %    tstep = 1;
            %    T = []; P = []; F = []; QN = []; QL = [];
            %    while (tstep>0)
            %       t  = d.runHydraulicAnalysis;
            %       qt = d.runQualityAnalysis;
            %
            %       P  = [P; d.getNodePressure];
            %       F  = [F; d.getLinkFlows];
            %
            %       QN = [QN; d.getNodeActualQuality];
            %       QL = [QL; d.getLinkActualQuality];
            %       T  = [T; t];
            %
            %       tstep = d.nextHydraulicAnalysisStep;
            %       qtstep = d.nextQualityAnalysisStep;
            %    end
            %    d.closeQualityAnalysis;
            %    d.closeHydraulicAnalysis;
            %
            % See also getNodeActualDemand, getNodeActualDemandSensingNodes, getNodeHydaulicHead
            %          getNodeActualQuality, getNodeMassFlowRate, getNodeActualQualitySensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_PRESSURE, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeActualQuality(obj, varargin)
            % Retrieves the computed values of the actual quality for all nodes.
            %
            % Example 1:
            %   d.getNodeActualQuality      % Retrieves the computed values of the actual quality for all nodes
            %
            % Example 2:
            %   d.getNodeActualQuality(1)   % Retrieves the computed value of the actual quality for the first node
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemand, getNodeActualDemandSensingNodes, getNodePressure,
            %          getNodeHydaulicHead, getNodeMassFlowRate, getNodeActualQualitySensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_QUALITY, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeMassFlowRate(obj, varargin)
            % Retrieves the computed mass flow rates per minute of chemical sources for all nodes.
            %
            % Example 1:
            %   d.getNodeMassFlowRate      % Retrieves the computed mass flow rates per minute of chemical sources for all nodes
            %
            % Example 2:
            %   d.getNodeMassFlowRate(1)   % Retrieves the computed mass flow rates per minute of chemical sources for the first node
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemand, getNodeActualDemandSensingNodes, getNodePressure,
            %          getNodeHydaulicHead, getNodeActualQuality, getNodeActualQualitySensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_SOURCEMASS, obj.LibEPANET);
                j=j+1;
            end
        end
        function value = getNodeActualQualitySensingNodes(obj, varargin)
            % Retrieves the computed quality values at some sensing nodes
            %
            % Example:
            %   d.getNodeActualQualitySensingNodes(1)   % Retrieves the computed quality value at the first node
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also getNodeActualDemand, getNodeActualDemandSensingNodes, getNodePressure,
            %          getNodeHydaulicHead, getNodeActualQuality, getNodeMassFlowRate.
            value=zeros(1, length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.Errcode, value(v)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_QUALITY, obj.LibEPANET);v=v+1;
            end
        end
        function value = getNodeTankInitialWaterVolume(obj, varargin)
            % Retrieves the tank initial water volume.
            %
            % Example 1:
            %   d.getNodeTankInitialWaterVolume              % Retrieves the initial water volume of all tanks
            %
            % Example 2:
            %   d.getNodeTankInitialWaterVolume(1)           % Retrieves the initial water volume of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankInitialWaterVolume(1:2)         % Retrieves the initial water volume of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankInitialWaterVolume(tankIndex)   % Retrieves the initial water volume of the tanks given their indices
            %
            % See also getNodeTankInitialLevel,  getNodeTankVolume,
            %          getNodeTankMaximumWaterVolume, getNodeTankMinimumWaterVolume.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_INITVOLUME, varargin);
        end
        function value = getNodeTankMixingModelCode(obj, varargin)
            % Retrieves the tank mixing model code.
            %
            % Code meaning:
            %   0 = Complete mix model (MIX1)
            %   1 = 2-compartment model (MIX2)
            %   2 = First in, first out model (FIFO)
            %   3 = Last in, first out model (LIFO)
            %
            % Example 1:
            %   d.getNodeTankMixingModelCode              % Retrieves the mixing model code of all tanks
            %
            % Example 2:
            %   d.getNodeTankMixingModelCode(1)           % Retrieves the mixing model code of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMixingModelCode(1:2)         % Retrieves the mixing model code of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMixingModelCode(tankIndex)   % Retrieves the mixing model code of the tanks given their indices
            %
            % See also setNodeTankMixingModelType, getNodeTankMixingModelType, getNodeTankMixZoneVolume.
            value = obj.get_node_tank_mixining_model(varargin{:});
            value = value{1};
        end
        function value = getNodeTankMixingModelType(obj, varargin)
            % Retrieves the tank mixing model type.
            %
            % Types of models that describe water quality mixing in storage tanks:
            %   MIX1 = Complete mix model
            %   MIX2 = 2-compartment model
            %   FIFO = First in, first out model
            %   LIFO = Last in, first out model
            %
            % Example 1:
            %   d.getNodeTankMixingModelType              % Retrieves the mixing model type of all tanks
            %
            % Example 2:
            %   d.getNodeTankMixingModelType(1)           % Retrieves the mixing model type of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMixingModelType(1:2)         % Retrieves the mixing model type of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMixingModelType(tankIndex)   % Retrieves the mixing model type of the tanks given their indices
            %
            % See also setNodeTankMixingModelType, getNodeTankMixingModelCode, getNodeTankMixZoneVolume.
            value = obj.get_node_tank_mixining_model(varargin{:});
            value = value{2};
        end
        function value = getNodeTankMixZoneVolume(obj, varargin)
            % Retrieves the tank mixing zone volume.
            %
            % Example 1:
            %   d.getNodeTankMixZoneVolume              % Retrieves the mixing zone volume of all tanks
            %
            % Example 2:
            %   d.getNodeTankMixZoneVolume(1)           % Retrieves the mixing zone volume of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMixZoneVolume(1:2)         % Retrieves the mixing zone volume of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMixZoneVolume(tankIndex)   % Retrieves the mixing zone volume of the tanks given their indices
            %
            % See also getNodeTankMixingModelCode, getNodeTankMixingModelType.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MIXZONEVOL, varargin);
        end
        function value = getNodeTankDiameter(obj, varargin)
            % Retrieves the tank diameters.
            %
            % Example 1:
            %   d.getNodeTankDiameter              % Retrieves the diameters of all tanks
            %
            % Example 2:
            %   d.getNodeTankDiameter(1)           % Retrieves the diameter of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankDiameter(1:2)         % Retrieves the diameters of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankDiameter(tankIndex)   % Retrieves the diameters of the tanks given their indices
            %
            % See also setNodeTankDiameter, getNodeTankBulkReactionCoeff, getNodeTankInitialLevel,
            %          getNodeTankMixingModelType, getNodeTankVolume, getNodeTankNameID.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_TANKDIAM, varargin);
        end
        function value = getNodeTankMinimumWaterVolume(obj, varargin)
            % Retrieves the tank minimum water volume.
            %
            % Example 1:
            %   d.getNodeTankMinimumWaterVolume              % Retrieves the minimum water volume of all tanks
            %
            % Example 2:
            %   d.getNodeTankMinimumWaterVolume(1)           % Retrieves the minimum water volume of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMinimumWaterVolume(1:2)         % Retrieves the minimum water volume of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMinimumWaterVolume(tankIndex)   % Retrieves the minimum water volume of the tanks given their indices
            %
            % See also setNodeTankMinimumWaterVolume, getNodeTankMaximumWaterVolume, getNodeTankInitialWaterVolume,
            %          getNodeTankInitialLevel,  getNodeTankVolume, getNodeTankMixZoneVolume.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MINVOLUME, varargin);
        end
        function value = getNodeTankVolumeCurveIndex(obj, varargin)
            % Retrieves the tank volume curve index.
            %
            % Example 1:
            %   d.getNodeTankVolumeCurveIndex              % Retrieves the volume curve index of all tanks
            %
            % Example 2:
            %   d.getNodeTankVolumeCurveIndex(1)           % Retrieves the volume curve index of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankVolumeCurveIndex(1:2)         % Retrieves the volume curve index of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankVolumeCurveIndex(tankIndex)   % Retrieves the volume curve index of the tanks given their indices
            %
            % See also getNodeTankVolume, getNodeTankMaximumWaterVolume, getNodeTankMinimumWaterVolume,
            %          getNodeTankInitialWaterVolume, getNodeTankMixZoneVolume.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_VOLCURVE, varargin);
        end
        function value = getNodeTankMinimumWaterLevel(obj, varargin)
            % Retrieves the tank minimum water level.
            %
            % Example 1:
            %   d.getNodeTankMinimumWaterLevel              % Retrieves the minimum water level of all tanks
            %
            % Example 2:
            %   d.getNodeTankMinimumWaterLevel(1)           % Retrieves the minimum water level of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMinimumWaterLevel(1:2)         % Retrieves the minimum water level of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMinimumWaterLevel(tankIndex)   % Retrieves the minimum water level of the tanks given their indices
            %
            % See also setNodeTankMinimumWaterLevel, getNodeTankMaximumWaterLevel, getNodeTankInitialLevel,
            %          getNodeTankMaximumWaterVolume, getNodeTankMinimumWaterVolume, getNodeTankVolume.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MINLEVEL, varargin);
        end
        function value = getNodeTankMaximumWaterLevel(obj, varargin)
            % Retrieves the tank maximum water level.
            %
            % Example 1:
            %   d.getNodeTankMaximumWaterLevel              % Retrieves the maximum water level of all tanks
            %
            % Example 2:
            %   d.getNodeTankMaximumWaterLevel(1)           % Retrieves the maximum water level of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMaximumWaterLevel(1:2)         % Retrieves the maximum water level of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMaximumWaterLevel(tankIndex)   % Retrieves the maximum water level of the tanks given their indices
            %
            % See also setNodeTankMaximumWaterLevel, getNodeTankMinimumWaterLevel, getNodeTankInitialLevel,
            %          getNodeTankMaximumWaterVolume, getNodeTankMinimumWaterVolume, getNodeTankVolume.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MAXLEVEL, varargin);
        end
        function value = getNodeTankMixingFraction(obj, varargin)
            % Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank.
            %
            % Example 1:
            %    d.getNodeTankMixingFraction             % Retrieves the mixing fraction of all tanks
            %
            % Example 2:
            %   d.getNodeTankMixingFraction(1)           % Retrieves the mixing fraction of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMixingFraction(1:2)         % Retrieves the mixing fraction of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMixingFraction(tankIndex)   % Retrieves the mixing fraction of the tanks given their indices
            %
            % See also setNodeTankMixingFraction, getNodeTankData.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MIXFRACTION, varargin);
        end
        function value = getNodeTankBulkReactionCoeff(obj, varargin)
            % Retrieves the tank bulk rate coefficient.
            %
            % Example 1:
            %   d.getNodeTankBulkReactionCoeff              % Retrieves the bulk rate coefficient of all tanks
            %
            % Example 2:
            %   d.getNodeTankBulkReactionCoeff(1)           % Retrieves the bulk rate coefficient of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankBulkReactionCoeff(1:2)         % Retrieves the bulk rate coefficient of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankBulkReactionCoeff(tankIndex)   % Retrieves the bulk rate coefficient of the tanks given their indices
            %
            % See also setNodeTankBulkReactionCoeff, getNodeTankData.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_TANK_KBULK, varargin);
        end
        function value = getNodeTankVolume(obj, varargin)
            % Retrieves the tank volume. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getNodeTankVolume              % Retrieves the volume of all tanks
            %
            % Example 2:
            %   d.getNodeTankVolume(1)           % Retrieves the volume of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankVolume(1:2)         % Retrieves the volume of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankVolume(tankIndex)   % Retrieves the volume of the tanks given their indices
            %
            % See also getNodeTankData.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_TANKVOLUME, varargin);
        end
        function value = getNodeTankMaximumWaterVolume(obj, varargin)
            % Retrieves the tank maximum water volume. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getNodeTankMaximumWaterVolume              % Retrieves the maximum water volume of all tanks
            %
            % Example 2:
            %   d.getNodeTankMaximumWaterVolume(1)           % Retrieves the maximum water volume of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankMaximumWaterVolume(1:2)         % Retrieves the maximum water volume of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankMaximumWaterVolume(tankIndex)   % Retrieves the maximum water volume of the tanks given their indices
            %
            % See also getNodeTankMinimumWaterVolume, getNodeTankData.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_MAXVOLUME, varargin);
        end
        function value = getNodeTankCanOverFlow(obj, varargin)
            % Retrieves the tank can overflow (= 1) or not (= 0). (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getNodeTankCanOverFlow              % Retrieves the can overflow of all tanks
            %
            % Example 2:
            %   d.getNodeTankCanOverFlow(1)           % Retrieves the can overflow of the 1st tank
            %
            % Example 3:
            %   d.getNodeTankCanOverFlow(1:2)         % Retrieves the can overflow of the first 2 tanks
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.getNodeTankCanOverFlow(tankIndex)   % Retrieves the can overflow of the tanks given their indices
            %
            % See also setNodeTankCanOverFlow, getNodeTankData.
            value = get_node_link(obj, 'tank', 'ENgetnodevalue', obj.ToolkitConstants.EN_CANOVERFLOW, varargin);
        end
        function value = getNodeDemandDeficit(obj, varargin)
            % Retrieves the amount that full demand is reduced under PDA. (EPANET Version 2.2)
            %
            % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d.setDemandModel('PDA', 0, 0.1, 0.5);   % Sets a type of demand model and its parameters
            %   d.getComputedHydraulicTimeSeries        % Computes hydraulic simulation and retrieve all time-series
            %   d.getNodeDemandDeficit                  % Retrieves the amount that full demand is reduced under PDA
            %
            % See also setDemandModel, getComputedHydraulicTimeSeries,
            %          getNodeActualDemand, getNodeActualDemandSensingNodes.
            [indices, value] = getNodeIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i, obj.ToolkitConstants.EN_DEMANDDEFICIT, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getNodeTankIndex(obj)
            % Retrieves the tank indices.
            %
            % Example:
            %   d.getNodeTankIndex
            %
            % See also getNodeTankCount, getNodeTankNameID.
            tmpNodeTypes=obj.getNodeType;
            value = find(strcmp(tmpNodeTypes, 'TANK'));
        end
        function value = getNodeTankNameID(obj)
            % Retrieves the tank IDs.
            %
            % Example 1:
            %   d.getNodeTankNameID      % Retrieves the IDs of all tanks
            %
            % Example 2:
            %   d.getNodeTankNameID{1}   % Retrieves the ID of the 1st tank
            %
            % See also getNodeTankCount, getNodeTankIndex.
            value=obj.getNodeNameID(obj.getNodeTankIndex);
        end
        function tankData = getNodeTankData(obj, varargin)
            % Retrieves a group of properties for a tank. (EPANET Version 2.2)
            %
            % Tank data that is retrieved:
            %
            % 1) Tank index
            % 2) Elevation
            % 3) Initial Level
            % 4) Minimum Water Level
            % 5) Maximum Water Level
            % 6) Diameter
            % 7) Minimum Water Volume
            % 8) Volume Curve Index
            %
            % Example 1:
            %   tankData = d.getNodeTankData;                  % Retrieves all the data of all tanks
            %   disp(tankData)
            %
            % Example 2:
            %   tankIndex = d.getNodeTankIndex;
            %   tankData = d.getNodeTankData(tankIndex);       % Retrieves all the data given the index/indices of tanks.
            %
            % Example 3:
            %   tankElevation = d.getNodeTankData.Elevation;   % Retrieves the elevations of all tanks.
            %   disp(tankElevation)
            %
            % See also setNodeTankData, getNodeElevations, getNodeTankInitialLevel,
            %          getNodeTankMinimumWaterLevel, getNodeTankDiameter.
            tankData = struct();
            if nargin == 1
                tankIndices = obj.getNodeTankIndex;
            elseif nargin == 2
                tankIndices = varargin{1};
            end
            tankData.Index = tankIndices;
            tankData.Elevation = obj.getNodeElevations(tankIndices);
            tankData.Initial_Level = obj.getNodeTankInitialLevel(tankIndices);
            tankData.Minimum_Water_Level = obj.getNodeTankMinimumWaterLevel(tankIndices);
            tankData.Maximum_Water_Level = obj.getNodeTankMaximumWaterLevel(tankIndices);
            tankData.Diameter = obj.getNodeTankDiameter(tankIndices);
            tankData.Minimum_Water_Volume = obj.getNodeTankMinimumWaterVolume(tankIndices);
            tankData.Volume_Curve_Index = obj.getNodeTankVolumeCurveIndex(tankIndices);
        end
        function value = getOptionsMaxTrials(obj)
            % Retrieves the maximum hydraulic trials allowed for hydraulic convergence.
            %
            % Example:
            %   d.getOptionsMaxTrials
            %
            % See also setOptionsMaxTrials, getOptionsExtraTrials, getOptionsAccuracyValue.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_TRIALS, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsAccuracyValue(obj)
            % Retrieves the total normalized flow change for hydraulic convergence.
            %
            % Example:
            %   d.getOptionsAccuracyValue
            %
            % See also setOptionsAccuracyValue, getOptionsExtraTrials, getOptionsMaxTrials.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_ACCURACY, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsQualityTolerance(obj)
            % Retrieves the water quality analysis tolerance.
            %
            % Example:
            %   d.getOptionsQualityTolerance
            %
            % See also setOptionsQualityTolerance, getOptionsSpecificDiffusivity, getOptionsLimitingConcentration.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_TOLERANCE, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsEmitterExponent(obj)
            % Retrieves the power exponent for the emmitters.
            %
            % Example:
            %   d.getOptionsEmitterExponent
            %
            % See also setOptionsEmitterExponent, getOptionsPatternDemandMultiplier, getOptionsAccuracyValue.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_EMITEXPON, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsPatternDemandMultiplier(obj)
            % Retrieves the global pattern demand multiplier.
            %
            % Example:
            %   d.getOptionsPatternDemandMultiplier
            %
            % See also setOptionsPatternDemandMultiplier, getOptionsEmitterExponent, getOptionsAccuracyValue.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_DEMANDMULT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsHeadError(obj)
            % Retrieves the maximum head loss error for hydraulic convergence. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsHeadError
            %
            % See also setOptionsHeadError, getOptionsEmitterExponent, getOptionsAccuracyValue.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_HEADERROR, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsFlowChange(obj)
            % Retrieves the maximum flow change for hydraulic convergence. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsFlowChange
            %
            % See also setOptionsFlowChange, getOptionsHeadError, getOptionsHeadLossFormula.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_FLOWCHANGE, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsHeadLossFormula(obj)
            % Retrieves the headloss formula. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsHeadLossFormula
            %
            % See also setOptionsHeadLossFormula, getOptionsHeadError, getOptionsFlowChange.
            [obj.Errcode, headloss] = obj.ENgetoption(obj.ToolkitConstants.EN_HEADLOSSFORM, obj.LibEPANET);
            value= obj.TYPEHEADLOSS{headloss+1};
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsGlobalEffic(obj)
            % Retrieves the global efficiency for pumps(percent). (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsGlobalEffic
            %
            % See also setOptionsGlobalEffic, getOptionsGlobalPrice, getOptionsGlobalPattern.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_GLOBALEFFIC, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsGlobalPrice(obj)
            % Retrieves the global average energy price per kW-Hour. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsGlobalPrice
            %
            % See also setOptionsGlobalPrice, getOptionsGlobalEffic, getOptionsGlobalPattern.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_GLOBALPRICE, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsGlobalPattern(obj)
            % Retrieves the index of the global energy price pattern. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsGlobalPattern
            %
            % See also setOptionsGlobalPattern, getOptionsGlobalEffic, getOptionsGlobalPrice.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_GLOBALPATTERN, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsDemandCharge(obj)
            % Retrieves the energy charge per maximum KW usage. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsDemandCharge
            %
            % See also setOptionsDemandCharge, getOptionsGlobalPrice, getOptionsGlobalPattern.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_DEMANDCHARGE, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsSpecificGravity(obj)
            % Retrieves the specific gravity. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsSpecificGravity
            %
            % See also setOptionsSpecificGravity, getOptionsSpecificViscosity, getOptionsHeadLossFormula.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_SP_GRAVITY, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsSpecificViscosity(obj)
            % Retrieves the specific viscosity. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsSpecificViscosity
            %
            % See also setOptionsSpecificViscosity, getOptionsSpecificGravity, getOptionsHeadLossFormula.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_SP_VISCOS, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsExtraTrials(obj)
            % Retrieves the extra trials allowed if hydraulics don't converge. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsExtraTrials
            %
            % See also setOptionsExtraTrials, getOptionsMaxTrials, getOptionsMaximumCheck.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_UNBALANCED, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsCheckFrequency(obj)
            % Retrieves the frequency of hydraulic status checks. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsCheckFrequency
            %
            % See also setOptionsCheckFrequency, getOptionsMaxTrials, getOptionsMaximumCheck.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_CHECKFREQ, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsMaximumCheck(obj)
            % Retrieves the maximum trials for status checking. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsMaximumCheck
            %
            % See also setOptionsMaximumCheck, getOptionsMaxTrials, getOptionsCheckFrequency.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_MAXCHECK, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsDampLimit(obj)
            % Retrieves the accuracy level where solution damping begins. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsDampLimit
            %
            % See also setOptionsDampLimit, getOptionsMaxTrials, getOptionsCheckFrequency.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_DAMPLIMIT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsSpecificDiffusivity(obj)
            % Retrieves the specific diffusivity (relative to chlorine at 20 deg C). (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsSpecificDiffusivity
            %
            % See also setOptionsSpecificDiffusivity, getOptionsSpecificViscosity, getOptionsSpecificGravity.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_SP_DIFFUS, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsPipeBulkReactionOrder(obj)
            % Retrieves the bulk water reaction order for pipes. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsPipeBulkReactionOrder
            %
            % See also setOptionsPipeBulkReactionOrder, getOptionsPipeWallReactionOrder, getOptionsTankBulkReactionOrder.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_BULKORDER, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsPipeWallReactionOrder(obj)
            % Retrieves the wall reaction order for pipes (either 0 or 1). (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsPipeWallReactionOrder
            %
            % See also setOptionsPipeWallReactionOrder, getOptionsPipeBulkReactionOrder, getOptionsTankBulkReactionOrder.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_WALLORDER, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsTankBulkReactionOrder(obj)
            % Retrieves the bulk water reaction order for tanks. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsTankBulkReactionOrder
            %
            % See also setOptionsTankBulkReactionOrder, getOptionsPipeBulkReactionOrder, getOptionsPipeWallReactionOrder.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_TANKORDER, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getOptionsLimitingConcentration(obj)
            % Retrieves the limiting concentration for growth reactions. (EPANET Version 2.2)
            %
            % Example:
            %   d.getOptionsLimitingConcentration
            %
            % See also setOptionsLimitingConcentration, getOptionsPipeBulkReactionOrder, getOptionsPipeWallReactionOrder.
            [obj.Errcode, value] = obj.ENgetoption(obj.ToolkitConstants.EN_CONCENLIMIT, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = getPatternComment(obj, varargin)
            % Retrieves the comment string assigned to the pattern object.
            %
            % Example 1:
            %   d.getPatternComment        % Retrieves the comments of all the patterns
            %
            % Example 2:
            %   d.getPatternComment(1)     % Retrieves the comment of the 1st pattern
            %
            % Example 3:
            %   d.getPatternComment(1:2)   % Retrieves the comments of the first 2 patterns
            %
            % See also setPatternComment, getPattern.
            if isempty(varargin)
                cnt = obj.getPatternCount;
                value = cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=ENgetcomment(obj.ToolkitConstants.EN_TIMEPAT, i, obj.LibEPANET);
                end
            else
                if isempty(varargin{1}), varargin{1}=0; end
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value{k}]=ENgetcomment(obj.ToolkitConstants.EN_TIMEPAT, i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            end
        end
        function value = getPatternNameID(obj, varargin)
            % Retrieves the ID label of all or some time patterns indices.
            %
            % Example 1:
            %   d.getPatternNameID        % Retrieves the IDs of all the patterns
            %
            % Example 2:
            %   d.getPatternNameID(1)     % Retrieves the ID of the 1st pattern
            %
            % Example 3:
            %   d.getPatternNameID(1:2)   % Retrieves the IDs of the first 2 patterns
            %
            % See also setPatternNameID, getPattern.
            patCnt = obj.getPatternCount;
            value = {};
            if patCnt
                if isempty(varargin)
                    value = cell(1, patCnt);
                    for i=1:patCnt
                        [obj.Errcode, value{i}]=obj.ENgetpatternid(i, obj.LibEPANET);
                    end
                else
                    k=1;
                    value = cell(1, length(varargin{1}));
                    for i=varargin{1}
                        [obj.Errcode, value{k}]=obj.ENgetpatternid(i, obj.LibEPANET);
                        k=k+1;
                    end
                end
            end
        end
        function setPatternNameID(obj, index, id)
            % Sets the name ID of a time pattern given it's index and the new ID. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getPatternNameID                                     % Retrieves the name IDs of all the time patterns
            %   d.setPatternNameID(1, 'Pattern1')                      % Sets to the 1st time pattern the new name ID 'Pattern1'
            %   d.getPatternNameID
            %
            % Example 2:
            %   d.setPatternNameID([1, 2], {'Pattern1', 'Pattern2'})   % Sets to the 1st and 2nd time pattern the new name IDs 'Pattern1' and 'Pattern2' respectively
            %   d.getPatternNameID
            %
            % See also getPatternNameID, getPatternIndex, getPatternLengths,
            %          setPatternComment, setPattern.
            if ischar(id)
                id={id};
            end
            for i=1:length(index)
                [obj.Errcode] = obj.ENsetpatternid(index(i), id{i}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getCurveNameID(obj, varargin)
            % Retrieves the IDs of curves. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getCurveNameID        % Retrieves the IDs of all the curves
            %
            % Example 2:
            %   d.getCurveNameID(1)     % Retrieves the ID of the 1st curve
            %
            % Example 3:
            %   d.getCurveNameID(1:2)   % Retrieves the IDs of the first 2 curves
            %
            % See also setCurveNameID, getCurvesInfo.
            curCnt = obj.getCurveCount;
            value = {};
            if curCnt
                if isempty(varargin)
                    value = cell(1, curCnt);
                    for i=1:curCnt
                        [obj.Errcode, value{i}]=obj.ENgetcurveid(i, obj.LibEPANET);
                    end
                else
                    k=1;
                    value = cell(1, length(varargin{1}));
                    for i=varargin{1}
                        [obj.Errcode, value{k}]=obj.ENgetcurveid(i, obj.LibEPANET);
                        k=k+1;
                    end
                end
            end
        end
        function setCurveNameID(obj, index, id)
            % Sets the name ID of a curve given it's index and the new ID. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getCurveNameID                                % Retrieves the name IDs of all the curves
            %   d.setCurveNameID(1,'Curve1')                    % Sets to the 1st curve the new name ID 'Curve1'
            %   d.getCurveNameID
            %
            % Example 2:
            %   d.setCurveNameID([1, 2],{'Curve1', 'Curve2'})   % Sets to the 1st and 2nd curve the new name IDs 'Curve1' and 'Curve2' respectively
            %   d.getCurveNameID
            %
            % See also getCurveNameID, getCurveIndex, getCurveLengths,
            %          setCurve, setCurveComment, getCurveComment.
            if ischar(id)
                id ={id};
            end
            for i=1:length(index)
                [obj.Errcode] = obj.ENsetcurveid(index(i), id{i}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getCurveLengths(obj, varargin)
            % Retrieves number of points in a curve. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getCurveLengths        % Retrieves the number of points in all the curves
            %
            % Example 2:
            %   d.getCurveLengths(1)     % Retrieves the number of points in the 1st curve
            %
            % Example 3:
            %   d.getCurveLengths(1:2)   % Retrieves the number of points in the first 2 curves
            %
            % See also getCurvesInfo, setCurve.
            if isempty(varargin)
                curCnt = obj.getCurveCount;
                tmpCurves = 1:curCnt;
                value = zeros(1, curCnt);
                for i=tmpCurves
                    [obj.Errcode, value(i)]=ENgetcurvelen(i, obj.LibEPANET);
                end
            elseif isa(varargin{1}, 'cell')
                k=1;
                lentmpCurves = length(varargin{1});
                value = zeros(1, lentmpCurves);
                for j=1:lentmpCurves
                    [obj.Errcode, value(k)] = ENgetcurvelen(obj.getCurveIndex(varargin{1}{j}), obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = ENgetcurvelen(obj.getCurveIndex(varargin{1}), obj.LibEPANET);
            elseif isa(varargin{1}, 'numeric')
                k=1;
                value = zeros(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value(k)]=ENgetcurvelen(i, obj.LibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getCurveIndex(obj, varargin)
            % Retrieves the index of a curve with specific ID. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getCurveIndex            % Retrieves the indices of all the curves
            %
            % Example 2:
            %   curveID = d.getCurveNameID(1);
            %   d.getCurveIndex(curveID)   % Retrieves the index of the 1st curve given it's ID
            %
            % Example 3:
            %   curveID = d.getCurveNameID(1:2);
            %   d.getCurveIndex(curveID)   % Retrieves the indices of the first 2 curves given their ID
            %
            % See also getCurveNameID, getCurvesInfo.
            if isempty(varargin)
                value=1:obj.getCurveCount;
            elseif isa(varargin{1}, 'cell')
                k=1;
                value = zeros(1, length(varargin{1}));
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = ENgetcurveindex(varargin{1}{j}, obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = ENgetcurveindex(varargin{1}, obj.LibEPANET);
            end
        end
        function value = getCurveTypeIndex(obj, varargin)
            % Retrieves the curve-type index for all curves.
            %
            % Example 1:
            %   d.getCurveTypeIndex        % Retrieves the curve-type index for all curves
            %
            % Example 2:
            %   d.getCurveTypeIndex(1)     % Retrieves the curve-type index for the 1st curve
            %
            % Example 3:
            %   d.getCurveTypeIndex(1:2)   % Retrieves the curve-type index for the first 2 curves
            %
            % See also getCurveType, getCurvesInfo.
            [indices, value] = getCurveIndices(obj, varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] =obj.ENgetcurvetype(i, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                j=j+1;
            end
        end
        function value = getCurveType(obj, varargin)
            % Retrieves the curve-type for all curves.
            %
            % Example 1:
            %   d.getCurveType        % Retrieves the curve-type for all curves
            %
            % Example 2:
            %   d.getCurveType(1)     % Retrieves the curve-type for the 1st curve
            %
            % Example 3:
            %   d.getCurveType(1:2)   % Retrieves the curve-type for the first 2 curves
            %
            % See also getCurveTypeIndex, getCurvesInfo.
            indices = getCurveIndices(obj, varargin);
            value=obj.TYPECURVE(obj.getCurveTypeIndex(indices)+1);
        end
        function value = getCurveComment(obj, varargin)
            % Retrieves the comment string of a curve.
            %
            % Example 1:
            %   d.getCurveComment        % Retrieves the comment string assigned to all the curves
            %
            % Example 2:
            %   d.getCurveComment(1)     % Retrieves the comment string assigned to the 1st curve
            %
            % Example 3:
            %   d.getCurveComment(1:2)   % Retrieves the comment string assigned to the first 2 curves
            %
            % See also getCurveNameID, getCurveType, getCurvesInfo,
            if isempty(varargin)
                cnt = obj.getCurveCount;
                value = cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=ENgetcomment(obj.ToolkitConstants.EN_CURVE, i, obj.LibEPANET);
                end
            else
                if isempty(varargin{1}), varargin{1}=0; end
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value{k}]=ENgetcomment(obj.ToolkitConstants.EN_CURVE, i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            end
        end
        function setCurve(obj, index, curveVector)
            % Sets x, y values for a specific curve. (EPANET Version 2.1)
            %
            % % The example is based on d=epanet('BWSN_Network_1.inp');
            %
            % Example:
            %   curveIndex = 1;
            %   d.getCurvesInfo.CurveXvalue{curveIndex}   % Retrieves the X values of the 1st curve
            %   d.getCurvesInfo.CurveYvalue{curveIndex}   % Retrieves the Y values of the 1st curve
            %   x_y_1 = [0, 730];
            %   x_y_2 = [1000, 500];
            %   x_y_3 = [1350, 260];
            %   values = [x_y_1; x_y_2; x_y_3];           % X and Y values selected.
            %   d.setCurve(curveIndex, values)            % Sets the X and Y values of the 1st curve
            %   d.getCurvesInfo.CurveXvalue{curveIndex}
            %   d.getCurvesInfo.CurveYvalue{curveIndex}
            %
            % See also setCurveValue, getCurvesInfo.
            nfactors=size(curveVector, 1);%x = number of points in curve
            [obj.Errcode] = ENsetcurve(index, curveVector(:, 1), curveVector(:, 2), nfactors, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = setCurveComment(obj, value, varargin)
            % Sets the comment string of a curve.
            %
            % Example 1:
            %   d.getCurveComment                         % Retrieves the comments of all the curves
            %   curveIndex = 1;
            %   comment = 'This is a curve';
            %   d.setCurveComment(curveIndex, comment);   % Sets a comment to the 1st curve
            %   d.getCurveComment(curveIndex)
            %
            % Example 2:
            %   d.getCurveComment
            %   curveIndex = 1:2;
            %   comment = {'This is the 1st curve', 'This is the 2nd curve'};
            %   d.setCurveComment(curveIndex, comment);   % Sets comments to the first 2 curves
            %   d.getCurveComment(curveIndex)
            %
            % See also getCurveComment, getCurveIndex, getCurvesInfo.
            if nargin==3, indices = value; value=varargin{1}; else indices = getCurveIndices(obj, varargin); end
            j=1;
            if length(indices) == 1
                [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_CURVE, indices, value, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                for i=indices
                    [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_CURVE, i, value{j}, obj.LibEPANET); j=j+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function setCurveValue(obj, index, curvePnt, value)
            % Sets x, y point for a specific point number and curve. (EPANET Version 2.1)
            %
            % % The example is based on d=epanet('BWSN_Network_1.inp');
            %
            % Example:
            %   curveIndex = 1;
            %   d.getCurvesInfo.CurveXvalue{curveIndex}               % Retrieves the X values of the 1st curve
            %   d.getCurvesInfo.CurveYvalue{curveIndex}               % Retrieves the Y values of the 1st curve
            %   curvePoint = 1;                                       % Point of the curve selected
            %   x_y_values = [10, 400];                               % X and Y values selected
            %   d.setCurveValue(curveIndex, curvePoint, x_y_values)   % Sets the X and Y values of the 1st point on the 1st curve
            %   d.getCurvesInfo.CurveXvalue{curveIndex}
            %   d.getCurvesInfo.CurveYvalue{curveIndex}
            %
            % See also getCurveValue, setCurve, getCurvesInfo.
            x=value(1); y=value(2);
            [obj.Errcode] = ENsetcurvevalue(index, curvePnt, x, y, obj.LibEPANET);
        end
        function value = getPatternIndex(obj, varargin)
            % Retrieves the index of all or some time patterns given their IDs.
            %
            % Example 1:
            %   d.getPatternIndex              % Retrieves the indices of all time patterns
            %
            % Example 2:
            %   patternIndex = 1;
            %   patternID = d.getPatternNameID(patternIndex);
            %   d.getPatternIndex(patternID)   % Retrieves the index of the 1st time pattern given it's ID
            %
            % Example 3:
            %   patternIndex = 1:2;
            %   patternID = d.getPatternNameID(patternIndex);
            %   d.getPatternIndex(patternID)   % Retrieves the index of the first 2 time patterns given their IDs
            %
            % See also getPatternNameID, getPattern.
            if isempty(varargin)
                value=1:obj.getPatternCount;
            elseif isa(varargin{1}, 'cell')
                k=1;
                value = zeros(1, length(varargin{1}));
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = obj.ENgetpatternindex(varargin{1}{j}, obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = obj.ENgetpatternindex(varargin{1}, obj.LibEPANET);
            end
        end
        function value = getPatternLengths(obj, varargin)
            % Retrieves the number of time periods in all or some time patterns.
            %
            % Example 1:
            %   d.getPatternLengths                 % Retrieves the number of time periods of all time patterns
            %
            % Example 2:
            %   patternIndex = 1;
            %   d.getPatternLengths(patternIndex)   % Retrieves the number of time periods of the 1st time pattern
            %
            % Example 3:
            %   patternIndex = 1:2;
            %   d.getPatternLengths(patternIndex)   % Retrieves the number of time periods of the first 2 time patterns
            %
            % See also getPatternIndex, getPattern.
            if isempty(varargin)
                patCnt = obj.getPatternCount;
                tmpPatterns=1:patCnt;
                value = zeros(1, patCnt);
                for i=tmpPatterns
                    [obj.Errcode, value(i)]=ENgetpatternlen(i, obj.LibEPANET);
                end
            elseif isa(varargin{1}, 'cell')
                k=1;
                lentmppat = length(varargin{1});
                value = zeros(1, lentmppat);
                for j=1:lentmppat
                    [obj.Errcode, value(k)] = ENgetpatternlen(obj.getPatternIndex(varargin{1}{j}), obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = ENgetpatternlen(obj.getPatternIndex(varargin{1}), obj.LibEPANET);
            elseif isa(varargin{1}, 'numeric')
                k=1;
                value = zeros(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value(k)]=ENgetpatternlen(i, obj.LibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getPattern(obj)
            % Retrieves the multiplier factor for all patterns and all times.
            %
            % Example:
            %   d.getPattern
            %
            % See also getPatternLengths, getPatternValue, getPatternAverageValue.
            tmpmaxlen=max(obj.getPatternLengths);
            value=nan(obj.getPatternCount, tmpmaxlen);
            for i=1:obj.getPatternCount
                tmplength=obj.getPatternLengths(i);
                for j=1:tmplength
                    [obj.Errcode, value(i, j)] = ENgetpatternvalue(i, j, obj.LibEPANET);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i, j)=value(i, j-tmplength);
                    end
                end
            end
        end
        function value = getPatternValue(obj, patternIndex, patternStep)
            % Retrieves the multiplier factor for a certain pattern and time.
            %
            % Example:
            %   patternIndex = 1;
            %   patternStep = 5;
            %   d.getPatternValue(patternIndex, patternStep)   % Retrieves the 5th multiplier factor of the 1st time pattern
            %
            % See also getPattern, getPatternLengths, getPatternAverageValue.
            [obj.Errcode, value] = ENgetpatternvalue(patternIndex, patternStep, obj.LibEPANET);
        end
        function value = getQualityType(obj, varargin)
            % Retrieves the type of water quality analysis type.
            %
            % Example:
            %   d.getQualityType
            %
            % See also getQualityInfo, getQualityCode.
            if any(strcmp(obj.getLibFunctions, 'ENgetqualinfo'))
                [obj.Errcode, ~, value] = ENgetqualinfo(obj.LibEPANET);
            else
                obj.saveInputFile(obj.BinTempfile);
                value = obj.getBinQualType;
            end
        end
        function value = getQualityInfo(obj)
            % Retrieves quality analysis information (type, chemical name, units, trace node ID).
            %
            % Information that is retrieved:
            %   1) Water quality analysis code
            %      0 = No quality analysis
            %      1 = Chemical analysis
            %      2 = Water age analysis
            %      3 = Source tracing
            %   2) Name of the chemical being analyzed
            %   3) Units that the chemical is measured in
            %   4) Index of node traced in a source tracing analysis
            %   5) Quality type
            %
            % Example:
            %   d.getQualityInfo                    % Retrieves all the quality info
            %   d.getQualityInfo.QualityCode        % Retrieves the water quality analysis code
            %   d.getQualityInfo.QualityChemName    % Retrieves the name of the chemical being analyzed
            %   d.getQualityInfo.QualityChemUnits   % Retrieves the units that the chemical is measured in
            %   d.getQualityInfo.TraceNode          % Retrieves the index of node traced in a source tracing analysis
            %   d.getQualityInfo.QualityType        % Retrieves the quality type
            %
            % See also getQualityType, getQualityCode.
            [obj.Errcode, value.QualityCode, value.QualityChemName, value.QualityChemUnits, value.TraceNode] = ENgetqualinfo(obj.LibEPANET);
            value.QualityType = obj.TYPEQUALITY(value.QualityCode+1);
        end
        function value = getQualityCode(obj)
            % Retrieves the code of water quality analysis type.
            %
            % Water quality analysis code:
            %   0 = No quality analysis
            %   1 = Chemical analysis
            %   2 = Water age analysis
            %   3 = Source tracing
            %
            % Example:
            %   d.getQualityCode
            %
            % See also getQualityInfo, getQualityType.
            [obj.Errcode, value, obj.QualityTraceNodeIndex] = ENgetqualtype(obj.LibEPANET);
        end
        function value = getQualityTraceNodeIndex(obj)
            % Retrieves the trace node index of water quality analysis type.
            %
            % Example:
            %   d.getQualityTraceNodeIndex
            %
            % See also getQualityInfo, getQualityType.
            [obj.Errcode, obj.QualityCode, value] = ENgetqualtype(obj.LibEPANET);
        end
        function value = getTimeSimulationDuration(obj)
            % Retrieves the value of simulation duration.
            %
            % Example:
            %   d.getTimeSimulationDuration
            %
            % See also getTimePatternStep, getTimeHydraulicStep.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_DURATION, obj.LibEPANET);
        end
        function value = getTimeHydraulicStep(obj)
            % Retrieves the value of the hydraulic time step.
            %
            % Example:
            %   d.getTimeHydraulicStep
            %
            % See also getTimeQualityStep, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_HYDSTEP, obj.LibEPANET);
        end
        function value = getTimeQualityStep(obj)
            % Retrieves the value of the water quality time step.
            %
            % Example:
            %   d.getTimeQualityStep
            %
            % See also getTimeHydraulicStep, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_QUALSTEP, obj.LibEPANET);
        end
        function value = getTimePatternStep(obj)
            % Retrieves the value of the pattern time step.
            %
            % Example:
            %   d.getTimePatternStep
            %
            % See also getTimePatternStart, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_PATTERNSTEP, obj.LibEPANET);
        end
        function value = getTimePatternStart(obj)
            % Retrieves the value of pattern start time.
            %
            % Example:
            %   d.getTimePatternStart
            %
            % See also getTimePatternStep, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_PATTERNSTART, obj.LibEPANET);
        end
        function value = getTimeReportingStep(obj)
            % Retrieves the value of the reporting time step.
            %
            % Example:
            %   d.getTimeReportingStep
            %
            % See also getTimeReportingPeriods, getTimeReportingStart.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_REPORTSTEP, obj.LibEPANET);
        end
        function value = getTimeReportingStart(obj)
            % Retrieves the value of the reporting start time.
            %
            % Example:
            %   d.getTimeReportingStart
            %
            % See also getTimeReportingPeriods, getTimeReportingStep.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_REPORTSTART, obj.LibEPANET);
        end
        function value = getTimeRuleControlStep(obj)
            % Retrieves the time step for evaluating rule-based controls.
            %
            % Example:
            %   d.getTimeRuleControlStep
            %
            % See also getTimeHydraulicStep.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_RULESTEP, obj.LibEPANET);
        end
        function value = getTimeStatisticsType(obj)
            % Retrieves the type of time series post-processing.
            %
            % Types:
            %   1) NONE:    Reports the full time series for all quantities for all nodes and links (default)
            %   2) AVERAGE: Reports a set of time-averaged results
            %   3) MINIMUM: Reports only the minimum values
            %   3) MAXIMUM: Reports only the maximum values
            %   4) RANGE:   Reports the difference between the minimum and maximum values
            %
            % Example:
            %   d.getTimeStatisticsType
            %
            % See also getTimeStatisticsIndex, getTimeSimulationDuration.
            [obj.Errcode, obj.TimeStatisticsIndex] = ENgettimeparam(obj.ToolkitConstants.EN_STATISTIC, obj.LibEPANET);
            value=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
        end
        function value = getTimeStatisticsIndex(obj)
            % Retrieves the index of the type of time series post-processing.
            %
            % Type of time series post-processing:
            %   0 = 'NONE'
            %   1 = 'AVERAGE'
            %   2 = 'MINIMUM'
            %   3 = 'MAXIMUM'
            %   4 = 'RANGE'
            %
            % Example:
            %   d.getTimeStatisticsIndex
            %
            % See also getTimeStatisticsType, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_STATISTIC, obj.LibEPANET);
        end
        function value = getTimeReportingPeriods(obj)
            % Retrieves the number of reporting periods saved to the binary.
            %
            % Example:
            %   d.getTimeReportingPeriods
            %
            % See also getTimeReportingStart, getTimeReportingStep.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_PERIODS, obj.LibEPANET);
        end
        %%%%% EPANET Version 2.1 %%%%%
        function value = getTimeStartTime(obj)
            % Retrieves the simulation starting time of day.
            %
            % Example:
            %   d.getTimeStartTime
            %
            % See also getTimeSimulationDuration, getTimePatternStart.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_STARTTIME, obj.LibEPANET);
        end
        function value = getTimeHTime(obj)
            % Retrieves the elapsed time of current hydraulic solution.
            %
            % Example:
            %   d.getTimeHTime
            %
            % See also getTimeHydraulicStep, getComputedHydraulicTimeSeries.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_HTIME, obj.LibEPANET);
        end
        function value = getTimeQTime(obj)
            % Retrieves the elapsed time of current quality solution.
            %
            % Example:
            %   d.getTimeQTime
            %
            % See also getTimeQualityStep, getComputedQualityTimeSeries.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_QTIME, obj.LibEPANET);
        end
        function value = getTimeHaltFlag(obj)
            % Retrieves the number of halt flag indicating if the simulation was halted.
            %
            % Example:
            %   d.getTimeHaltFlag
            %
            % See also getTimeStartTime, getTimeSimulationDuration.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_HALTFLAG, obj.LibEPANET);
        end
        function value = getTimeNextEvent(obj)
            % Retrieves the shortest time until a tank becomes empty or full.
            %
            % Example:
            %   d.getTimeNextEvent
            %
            % See also getTimeNextEventTank.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_NEXTEVENT, obj.LibEPANET);
        end
        function value = getTimeNextEventTank(obj)
            % Retrieves the index of tank with shortest time to become empty or full.
            %
            % Example:
            %   d.getTimeNextEventTank
            %
            % See also getTimeNextEvent.
            [obj.Errcode, value] = ENgettimeparam(obj.ToolkitConstants.EN_NEXTEVENTTANK, obj.LibEPANET);
        end
        function value = getCurvesInfo(obj)
            % Retrieves all the info of curves. (EPANET Version 2.1)
            %
            % Returns the following informations:
            %   1) Curve Name ID
            %   2) Number of points on curve
            %   3) X values of points
            %   4) Y values of points
            %
            % Example:
            %   d.getCurvesInfo
            %   d.getCurvesInfo.CurveNameID      % Retrieves the IDs of curves
            %   d.getCurvesInfo.CurveNvalue      % Retrieves the number of points on curve
            %   d.getCurvesInfo.CurveXvalue      % Retrieves the X values of points of all curves
            %   d.getCurvesInfo.CurveXvalue{1}   % Retrieves the X values of points of the 1st curve
            %   d.getCurvesInfo.CurveYvalue      % Retrieves the Y values of points of all curves
            %   d.getCurvesInfo.CurveYvalue{1}   % Retrieves the Y values of points of the 1st curve
            %
            % See also setCurve, getCurveType, getCurveLengths,
            %          getCurveValue, getCurveNameID, getCurveComment.
            value = struct();
            for i=1:obj.getCurveCount
                [obj.Errcode, value.CurveNameID{i}, value.CurveNvalue{i}, value.CurveXvalue{i}, value.CurveYvalue{i}] = ENgetcurve(obj, i, obj.LibEPANET);
            end
        end
        function value = getCounts(obj)
            % Retrieves the number of network components.
            % Nodes, Links, Junctions, Reservoirs, Tanks, Pipes, Pumps,
            % Valves, Curves, SimpleControls, RuleBasedControls, Patterns.
            %
            % Example:
            %   d.getCounts                  % Retrieves the number of all network components
            %   d.getCounts.Nodes            % Retrieves the number of nodes
            %   d.getCounts.SimpleControls   % Retrieves the number of simple controls
            %
            % See also getNodeCount, getNodeJunctionCount,
            %          getLinkCount, getControlRulesCount.
            value.Nodes = obj.getNodeCount;
            value.Links = obj.getLinkCount;
            value.Junctions = obj.getNodeJunctionCount;
            value.Reservoirs = obj.getNodeReservoirCount;
            value.Tanks = obj.getNodeTankCount;
            value.Pipes = obj.getLinkPipeCount;
            value.Pumps = obj.getLinkPumpCount;
            value.Valves = obj.getLinkValveCount;
            value.Curves = obj.getCurveCount;
            value.SimpleControls = obj.getControlRulesCount;
            value.RuleBasedControls = obj.getRuleCount;
            value.Patterns = obj.getPatternCount;

        end
        function value = getConnectivityMatrix(obj, varargin)
            conn = obj.getNodesConnectingLinksID;
            nodesID = obj.getNodeNameID;
            cnt = obj.getNodeCount;
            value = zeros(cnt, cnt);
            for i=1:cnt
                mm = strcmp(nodesID(i), conn);
                mS = mm(:, 1)+mm(:, 2);
                chIndex = find(mS);
                connFinal = conn(chIndex, :);
                mmFinal = mm(chIndex, :);
                Ok = connFinal(~mmFinal);
                nodesIndOk = obj.getNodeIndex(Ok);
                value(i, nodesIndOk) = 1;
            end
        end
        function valueIndex = addCurve(obj, varargin)
            % Adds a new curve appended to the end of the existing curves. (EPANET Version 2.1)
            % Returns the new curve's index.
            %
            % Example:
            %   new_curve_ID = 'NewCurve';                        % ID selected without a space in between the letters
            %   x_y_1 = [0, 730];
            %   x_y_2 = [1000, 500];
            %   x_y_3 = [1350, 260];
            %   values = [x_y_1; x_y_2; x_y_3];                   % X and Y values selected
            %   curve_index = d.addCurve(new_curve_ID, values);   % New curve added
            %   d.getCurvesInfo                                   % Retrieves all the info of curves
            %
            % See also getCurvesInfo, getCurveType, setCurve,
            %          setCurveValue, setCurveNameID, setCurveComment.
            valueIndex = 0;
            if (4>nargin && nargin>1)
                [obj.Errcode] = ENaddcurve(varargin{1}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
                valueIndex = getCurveIndex(obj, varargin{1});
                if nargin==3
                    setCurve(obj, valueIndex, varargin{2});
                end
            end
        end
        function value = getCurveValue(obj, varargin)
            % Retrieves the X, Y values of points of curves. (EPANET Version 2.1)
            %
            % Example 1:
            %   d.getCurveValue                           % Retrieves all the X and Y values of all curves
            %
            % Example 2:
            %   curveIndex = 1;
            %   d.getCurveValue(curveIndex)
            %   d.getCurveValue{curveIndex}               % Retrieves all the X and Y values of the 1st curve
            %
            % Example 3:
            %   curveIndex = 1;
            %   pointIndex = 1;
            %   d.getCurveValue(curveIndex, pointIndex)   % Retrieves the X and Y values of the 1st point of the 1st curve
            %
            % See also setCurveValue, setCurve, getCurvesInfo.
            tmplen = obj.getCurveLengths;
            cur_index = obj.getCurveIndex;
            pnt = 0;
            if nargin>2
                pnt = varargin{2};
            end
            if nargin>1
               index = varargin{1};
            else
               index = cur_index;
            end
            j = 1;
            if find(index == cur_index)
                for cur = cur_index(index)
                    for i=1:tmplen(cur)
                        if pnt
                            [obj.Errcode, value(i, 1), value(i, 2)] = ENgetcurvevalue(cur, pnt, obj.LibEPANET);
                            break;
                        else
                            [obj.Errcode, value{j}(i, 1), value{j}(i, 2)] = ENgetcurvevalue(cur, i, obj.LibEPANET);
                        end
                    end
                    j=j+1;
                end
            else
                obj.Errcode = 206;
                error(obj.getError(obj.Errcode));
            end
        end
        function [curveIndex, pumpIndex] = getLinkPumpHeadCurveIndex(obj)
            % Retrieves the index of a head curve for all pumps. (EPANET Version 2.1)
            %
            % Example:
            %   [curveIndex, pumpIndex] = d.getLinkPumpHeadCurveIndex
            %
            % See also getLinkPumpHCurve, getLinkPumpECurve.
            j = 1;
            pumpIndex = obj.getLinkPumpIndex;
            curveIndex = zeros(1, length(pumpIndex));
            for i=pumpIndex
                [obj.Errcode, curveIndex(j)] = ENgetheadcurveindex(i, obj.LibEPANET);
                j=j+1;
            end
        end
        function value = getLinkPumpTypeCode(obj)
            % Retrieves the code of type of a pump. (EPANET Version 2.1)
            %
            % Type of pump codes:
            %   0 = Constant horsepower
            %   1 = Power function
            %   2 = User-defined custom curve
            %
            % Example:
            %   d.getLinkPumpTypeCode
            %
            % See also getLinkPumpType, getLinkPumpPower.
            j = 1;
            pumpCnt = obj.getLinkPumpCount;
            value = ones(1, pumpCnt);
            if pumpCnt
                for i=obj.getLinkPumpIndex
                    [obj.Errcode, value(j)] = ENgetpumptype(i, obj.LibEPANET);
                    j=j+1;
                end
            end
        end
        function value = getLinkPumpType(obj)
            % Retrieves the type of a pump. (EPANET Version 2.1)
            %
            % Example:
            %   d.getLinkPumpType
            %
            % See also getLinkPumpTypeCode, getLinkPumpPower.
            v = obj.getLinkPumpTypeCode;
            value = obj.TYPEPUMP(v+1);
        end
        function value = getPatternAverageValue(obj)
            % Retrieves the average values of all the time patterns. (EPANET Version 2.1)
            %
            % Example:
            %   d.getPatternAverageValue
            %
            % See also getPattern, setPattern,
            %          getPatternValue, getPatternLengths.
            value = zeros(1, obj.getPatternCount);
            for i=obj.getPatternIndex
                [obj.Errcode, value(i)] = ENgetaveragepatternvalue(i, obj.LibEPANET);
            end
        end
        function value = getENfunctionsImpemented(obj)
            % Retrieves the epanet functions that have been developed.
            %
            % Example:
            %   d.getENfunctionsImpemented
            %
            % See also getLibFunctions, getVersion.
            [tline]=regexp(fileread([mfilename, '.m']), '\n', 'split');
            u=1;obj.Errcode;
            value = cell(1, 1);
            for i=1:length(tline)
                if ~isempty(regexp(tline{i}, 'EN\w', 'match'))
                    enfunc = regexp(tline{i}, '\w*EN\w*', 'match');
                    checkL = isstrprop(enfunc, 'upper');
                    if length(checkL{1})>2
                        if strcmp(enfunc{1}(1:2), 'EN') && ~checkL{1}(3) && ~strcmp(enfunc{1}(1:3), 'EN_')
                            value(u) = enfunc(1);
                            u=u+1;
                        end
                    end
                end
            end
            value = unique(value)';
        end
        function value = getLibFunctions(obj)
            % Retrieves the functions of DLL.
            %
            % Example:
            %   d.getLibFunctions
            %
            % See also getENfunctionsImpemented, getVersion.
            value = obj.libFunctions;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = getVersion(obj)
            % Retrieves the current EPANET version of DLL.
            %
            % Example:
            %   d.getVersion
            %
            % See also getENfunctionsImpemented, getLibFunctions.
            [obj.Errcode, value] = ENgetversion(obj.LibEPANET);
        end
        function value = getLinkPumpSwitches(obj)
            % Retrieves the number of pump switches.
            %
            % Example:
            %   d.getLinkPumpSwitches
            value=[];
            s = obj.getComputedTimeSeries;
            for i=1:obj.getLinkPumpCount
                index = find(diff(s.Status(:, obj.getLinkPumpIndex(i))));
                value(i) = length(index);
            end
        end
        function value = getComputedHydraulicTimeSeries(obj, varargin)
            % Computes hydraulic simulation and retrieves all time-series.
            %
            % Data that is computed:
            %   1) Time              8)  Velocity
            %   2) Pressure          9)  HeadLoss
            %   3) Demand            10) Status
            %   4) DemandDeficit     11) Setting
            %   5) Head              12) Energy
            %   6) TankVolume        13) Efficiency
            %   7) Flow
            %
            % Example 1:
            %   d.getComputedHydraulicTimeSeries          % Retrieves all the time-series data
            %
            % Example 2:
            %   d.getComputedHydraulicTimeSeries.Demand   % Retrieves all the time-series demands
            %   d.getComputedHydraulicTimeSeries.Flow     % Retrieves all the time-series flows
            %
            % Example 3:
            %   data = d.getComputedHydraulicTimeSeries('Time', ...
            %    'Pressure', 'Velocity');                  %  Retrieves all the time-series Time, Pressure, Velocity
            %   time = data.Time;
            %   pressure = data.Pressure;
            %   velocity = data.Velocity;
            %
            % See also getComputedQualityTimeSeries, getComputedTimeSeries.
            obj.openHydraulicAnalysis;
            obj.solve = 1;
            obj.initializeHydraulicAnalysis
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep;
            initnodematrix = zeros(totalsteps, obj.getNodeCount);
            initlinkmatrix = zeros(totalsteps, obj.getLinkCount);
            if size(varargin, 2)==0
                varargin = {'time', 'pressure', 'demand', 'demanddeficit', 'head', 'tankvolume', 'flow', 'velocity', 'headloss', 'status', 'setting', 'energy', 'efficiency', 'state'};
            %if ~sum(strcmpi(fields(obj.ToolkitConstants), 'EN_EFFICIENCY'))
            %    varargin{end} = {''};
            %end
            else
                for i = 1:length(varargin)
                    if isnumeric(varargin{i})
                        sensingnodes = i;
                    end
                end
            end
            if find(strcmpi(varargin, 'time'))
                value.Time = zeros(totalsteps, 1);
            end
            if find(strcmpi(varargin, 'pressure'))
                value.Pressure = initnodematrix;
            end
            if find(strcmpi(varargin, 'demand'))
                value.Demand = initnodematrix;
            end
            if obj.getVersion > 20101
                if find(strcmpi(varargin, 'demanddeficit'))
                    value.DemandDeficit = initnodematrix;
                end
            end
            if find(strcmpi(varargin, 'demandSensingNodes'))
                value.DemandSensingNodes = zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices = varargin{sensingnodes};
            end
            if find(strcmpi(varargin, 'head'))
                value.Head = initnodematrix;
            end
            if find(strcmpi(varargin, 'tankvolume'))
                value.TankVolume = initnodematrix;
            end
            if find(strcmpi(varargin, 'flow'))
                value.Flow = initlinkmatrix;
            end
            if find(strcmpi(varargin, 'velocity'))
                value.Velocity = initlinkmatrix;
            end
            if find(strcmpi(varargin, 'headloss'))
                value.HeadLoss = initlinkmatrix;
            end
            if find(strcmpi(varargin, 'status'))
                value.Status = initlinkmatrix;
            end
            if find(strcmpi(varargin, 'setting'))
                value.Setting = initlinkmatrix;
            end
            if find(strcmpi(varargin, 'energy'))
                value.Energy = initlinkmatrix;
            end
            if obj.getVersion > 20101
                if find(strcmpi(varargin, 'efficiency'))
                    value.Efficiency = initlinkmatrix;
                end
            end
            if find(strcmpi(varargin, 'state'))
                value.State = zeros(totalsteps, obj.getLinkPumpCount);
            end
            clear initlinkmatrix initnodematrix;
            k = 1;tstep = 1;
            while (tstep>0)
                t = obj.runHydraulicAnalysis;
                if find(strcmpi(varargin, 'time'))
                    value.Time(k, :) = t;
                end
                if find(strcmpi(varargin, 'pressure'))
                    value.Pressure(k, :) = obj.getNodePressure;
                end
                if find(strcmpi(varargin, 'demand'))
                    value.Demand(k, :) = obj.getNodeActualDemand;
                end
                if obj.getVersion > 20101
                    if find(strcmpi(varargin, 'demanddeficit'))
                        value.DemandDeficit(k, :) = obj.getNodeDemandDeficit;
                    end
                end
                if find(strcmpi(varargin, 'demandSensingNodes'))
                    value.DemandSensingNodes(k, :) = obj.getNodeActualDemandSensingNodes(varargin{sensingnodes});
                end
                if find(strcmpi(varargin, 'head'))
                    value.Head(k, :) = obj.getNodeHydaulicHead;
                end
                if find(strcmpi(varargin, 'tankvolume'))
                    value.TankVolume(k, :) = [zeros(1, obj.getNodeJunctionCount+obj.getNodeReservoirCount) obj.getNodeTankVolume];
                end
                if find(strcmpi(varargin, 'flow'))
                    value.Flow(k, :) = obj.getLinkFlows;
                end
                if find(strcmpi(varargin, 'velocity'))
                    value.Velocity(k, :) = obj.getLinkVelocity;
                end
                if find(strcmpi(varargin, 'headloss'))
                    value.HeadLoss(k, :) = obj.getLinkHeadloss;
                end
                if find(strcmpi(varargin, 'status'))
                    value.Status(k, :) = obj.getLinkStatus;
                    value.StatusStr(k, :) = obj.TYPESTATUS(value.Status(k, :) + 1);
                end
                if find(strcmpi(varargin, 'setting'))
                    value.Setting(k, :) = obj.getLinkSettings;
                end
                if find(strcmpi(varargin, 'energy'))
                    value.Energy(k, :) = obj.getLinkEnergy;
                end
                if obj.getVersion > 20101
                    if find(strcmpi(varargin, 'efficiency'))
                        value.Efficiency(k, :) = [zeros(1, obj.getLinkPipeCount) obj.getLinkPumpEfficiency zeros(1, obj.getLinkValveCount)];
                    end
                end
                if obj.getVersion > 20101
                    if find(strcmpi(varargin, 'state'))
                        value.State(k, :) = obj.getLinkPumpState;
                        value.StateStr(k, :) = obj.TYPEPUMPSTATE(value.State(k, :) + 1);
                    end
                end
                tstep = obj.nextHydraulicAnalysisStep;
                k = k+1;
            end
            obj.closeHydraulicAnalysis;
        end
        function value = getComputedQualityTimeSeries(obj, varargin)
            % Computes Quality simulation and retrieves all or some time-series.
            %
            % Data that is computed:
            %   1) Time
            %   2) NodeQuality
            %   3) LinkQuality
            %   4) MassFlowRate
            %
            % Example 1:
            %   d.getComputedQualityTimeSeries               % Retrieves all the time-series data
            %
            % Example 2:
            %   d.getComputedQualityTimeSeries.NodeQuality   % Retrieves all the time-series node quality
            %   d.getComputedQualityTimeSeries.LinkQuality   % Retrieves all the time-series link quality
            %
            % Example 3:
            %   data = d.getComputedQualityTimeSeries('Time', ...
            %   'NodeQuality', 'LinkQuality');              %  Retrieves all the time-series Time, NodeQuality, LinkQuality
            %   time = data.Time;
            %   node_quality = data.NodeQuality;
            %   link_quality = data.LinkQuality;
            %
            % See also getComputedHydraulicTimeSeries, getComputedTimeSeries.
            if ~obj.solve
                obj.solveCompleteHydraulics;
                obj.solve = 1;
            end
            obj.openQualityAnalysis;
            obj.initializeQualityAnalysis;
            %tleft=obj.nextQualityAnalysisStep;
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep;
            initnodematrix=zeros(totalsteps, obj.getNodeCount);
            if size(varargin, 2)==0
                varargin={'time', 'quality', 'linkquality', 'mass'};
            else
                for i=1:length(varargin)
                    if isnumeric(varargin{i})
                        sensingnodes=i;
                    end
                end
            end
            if find(strcmpi(varargin, 'time'))
                value.Time=zeros(totalsteps, 1);
            end
            if find(strcmpi(varargin, 'quality'))
                value.NodeQuality=initnodematrix;
            end
            if find(strcmpi(varargin, 'linkquality'))
                value.LinkQuality=zeros(totalsteps, obj.getLinkCount);
            end
            if find(strcmpi(varargin, 'qualitySensingNodes'))
                value.QualitySensingNodes=zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices=varargin{sensingnodes};
            end
            if find(strcmpi(varargin, 'demandSensingNodes'))
                value.DemandSensingNodes=zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices=varargin{sensingnodes};
            end
            if find(strcmpi(varargin, 'mass'))
                value.MassFlowRate=initnodematrix;
            end
            if find(strcmpi(varargin, 'demand'))
                value.Demand=initnodematrix;
            end
            clear initnodematrix;
            k=1;t=1;tleft=1;
            while (tleft>0)||(t<obj.getTimeSimulationDuration)
                t=obj.runQualityAnalysis;
                if find(strcmpi(varargin, 'time'))
                    value.Time(k, :)=t;
                end
                if find(strcmpi(varargin, 'quality'))
                    value.NodeQuality(k, :)=obj.getNodeActualQuality;
                end
                if find(strcmpi(varargin, 'linkquality'))
                    value.LinkQuality(k, :)=obj.getLinkActualQuality;
                end
                if find(strcmpi(varargin, 'mass'))
                    value.MassFlowRate(k, :)=obj.getNodeMassFlowRate;
                end
                if find(strcmpi(varargin, 'demand'))
                    value.Demand(k, :)=obj.getNodeActualDemand;
                end
                if find(strcmpi(varargin, 'qualitySensingNodes'))
                    value.QualitySensingNodes(k, :)=obj.getNodeActualQualitySensingNodes(varargin{2});
                end
                if find(strcmpi(varargin, 'demandSensingNodes'))
                    value.DemandSensingNodes(k, :)=obj.getNodeActualDemandSensingNodes(varargin{sensingnodes});
                end
                tleft = obj.nextQualityAnalysisStep;
                k=k+1;
                if t==obj.getTimeSimulationDuration
                    t=obj.getTimeSimulationDuration+1;
                end
            end
            obj.closeQualityAnalysis;
        end
        function value = getComputedTimeSeries(obj)
            obj.saveInputFile(obj.TempInpFile);
            [fid,binfile, ~] = runEPANETexe(obj);
            value = readEpanetBin(fid, binfile, 0);
            value.StatusStr = obj.TYPEBINSTATUS(value.Status + 1);
            % Remove report bin, files @#
            warning off;
            fclose('all');
            files=dir('@#*');
            if ~isempty(files); delete('@#*'); end
            warning on;
        end
        function value = getComputedTimeSeries_ENepanet(obj)
            obj.saveInputFile(obj.TempInpFile);
            uuID = char(java.util.UUID.randomUUID);
            rptfile=[obj.TempInpFile(1:end-4), '.txt'];
            binfile=['@#', uuID, '.bin'];
            obj.Errcode = ENepanet(obj.LibEPANET, obj.TempInpFile, rptfile, binfile);
            fid = fopen(binfile, 'r');
            value = readEpanetBin(fid, binfile, 0);
            value.StatusStr = obj.TYPEBINSTATUS(value.Status + 1);
            fclose('all');
            files=dir('@#*');
            if ~isempty(files); delete('@#*'); end
            obj.Errcode = ENopen(obj.TempInpFile, [obj.TempInpFile(1:end-4), '.txt'], [obj.TempInpFile(1:end-4), '.bin'], obj.LibEPANET);
        end
        function value = getUnits(obj)
            % Retrieves the Units of Measurement.
            %
            % Example 1:
            %   allUnits = d.getUnits           % Retrieves all the units
            %
            % Example 2:
            %   d.getUnits.NodeElevationUnits   % Retrieves elevation units
            %   d.getUnits.LinkVelocityUnits    % Retrieves velocity units
            %
            % More info: https://github.com/OpenWaterAnalytics/EPANET/wiki/Units-of-Measurement
            %
            % See also getFlowUnits.
            if find(strcmp(obj.getFlowUnits, obj.TYPEUNITS))<6
                value.Units_US_Customary=1;
                value.Units_SI_Metric=0;
            else
                value.Units_SI_Metric=1;
                value.Units_US_Customary=0;
            end

            value.LinkFlowUnits = obj.getFlowUnits;
            if value.Units_US_Customary
                value.NodePressureUnits='psi';
                value.PatternDemandsUnits=value.LinkFlowUnits;
                value.LinkPipeDiameterUnits='inches';
                value.NodeTankDiameterUnits='feet';
                value.EnergyEfficiencyUnits='percent';
                value.NodeElevationUnits='feet';
                value.NodeDemandUnits=value.LinkFlowUnits;
                value.NodeEmitterCoefficientUnits='flow units @ 1 psi drop';
                value.EnergyUnits='kwatt-hours';
                value.LinkFrictionFactorUnits='unitless';
                value.NodeHeadUnits='feet';
                value.LinkLengthsUnits='feet';
                value.LinkMinorLossCoeffUnits='unitless';
                value.LinkPumpPowerUnits='horsepower';
                value.QualityReactionCoeffBulkUnits='1/day (1st-order)';
                value.QualityReactionCoeffWallUnits='mass/sq-ft/day (0-order), ft/day (1st-order)';
                value.LinkPipeRoughnessCoeffUnits='millifeet(Darcy-Weisbach), unitless otherwise';
                value.QualitySourceMassInjectionUnits='mass/minute';
                value.LinkVelocityUnits='ft/sec';
                value.NodeTankVolumeUnits='cubic feet';
                value.QualityWaterAgeUnits='hours';
            else % SI Metric
                value.NodePressureUnits='meters';
                value.PatternDemandsUnits=value.LinkFlowUnits;
                value.LinkPipeDiameterUnits='millimeters';
                value.NodeTankDiameterUnits='meters';
                value.EnergyEfficiencyUnits='percent';
                value.NodeElevationUnits='meters';
                value.NodeDemandUnits=value.LinkFlowUnits;
                value.NodeEmitterCoefficientUnits='flow units @ 1 meter drop';
                value.EnergyUnits='kwatt-hours';
                value.LinkFrictionFactorUnits='unitless';
                value.NodeHeadUnits='meters';
                value.LinkLengthsUnits='meters';
                value.LinkMinorLossCoeffUnits='unitless';
                value.LinkPumpPowerUnits='kwatts';
                value.QualityReactionCoeffBulkUnits='1/day (1st-order)';
                value.QualityReactionCoeffWallUnits='mass/sq-m/day(0-order), meters/day (1st-order)';
                value.LinkPipeRoughnessCoeffUnits='mm(Darcy-Weisbach), unitless otherwise';
                value.QualitySourceMassInjectionUnits='mass/minute';
                value.LinkVelocityUnits='meters/sec';
                value.NodeTankVolumeUnits='cubic meters';
                value.QualityWaterAgeUnits='hours';
            end
        end
        function Errcode = solveCompleteHydraulics(obj)
            % Runs a complete hydraulic simulation with results for all time periods written to the binary Hydraulics file.
            %
            % Example:
            %   d.solveCompleteHydraulics
            %
            % See also solveCompleteQuality.
            obj.solve = 1;
            [Errcode] = ENsolveH(obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function Errcode = solveCompleteQuality(obj)
            % Runs a complete water quality simulation with results at uniform reporting intervals written to EPANET's binary Output file.
            %
            % Example:
            %   d.solveCompleteQuality
            %
            % See also solveCompleteHydraulics.
            [Errcode] = ENsolveQ(obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function index = addPattern(obj, varargin)
            % Adds a new time pattern to the network.
            %
            % Example 1:
            %   d.getPatternNameID                                     % Retrieves the ID labels of time patterns
            %   patternID = 'new_pattern';
            %   patternIndex = d.addPattern(patternID);                % Adds a new time pattern given it's ID
            %   d.getPatternNameID
            %
            % Example 2:
            %   patternID = 'new_pattern';
            %   patternMult = [1.56, 1.36, 1.17, 1.13, 1.08, ...
            %   1.04, 1.2, 0.64, 1.08, 0.53, 0.29, 0.9, 1.11, ...
            %   1.06, 1.00, 1.65, 0.55, 0.74, 0.64, 0.46, ...
            %   0.58, 0.64, 0.71, 0.66];
            %   patternIndex = d.addPattern(patternID, patternMult);   % Adds a new time pattern given it's ID and the multiplier
            %   d.getPatternNameID
            %   d.getPattern
            %
            % See also getPattern, setPattern, setPatternNameID
                %          setPatternValue, setPatternComment.
            [obj.Errcode] = obj.ENaddpattern(varargin{1}, obj.LibEPANET);
            index = getPatternIndex(obj, varargin{1});
            if nargin==2
                setPattern(obj, index, ones(1, max(obj.getPatternLengths)));
            elseif nargin==3
                setPattern(obj, index, varargin{2});
            end
        end
        function index = addNodeJunction(obj, juncID, varargin)
            % Adds a new junction.
            % Returns the index of the new junction.
            %
            % The following data can be set(optional):
            %   1) Coordinates
            %   2) Elevation
            %   3) Primary base demand
            %   4) ID name of the demand's time pattern
            %
            % Example 1:
            %   % Adds a new junction with the default coordinates (i.e. [0, 0])
            %   junctionID = 'newJunction_1';
            %   junctionIndex = d.addNodeJunction(junctionID);
            %   d.plot;
            %
            % Example 2:
            %   % Adds a new junction with coordinates [X, Y] = [20, 10].
            %   junctionID = 'newJunction_2';
            %   junctionCoords = [20 10];
            %   junctionIndex = d.addNodeJunction(junctionID, junctionCoords);
            %   d.plot;
            %
            % Example 3:
            %   % Adds a new junction with coordinates [X, Y] = [20, 20] and elevation = 500.
            %   junctionID = 'newJunction_3';
            %   junctionCoords = [20 20];
            %   junctionElevation = 500;
            %   junctionIndex = d.addNodeJunction(junctionID, junctionCoords, junctionElevation);
            %   d.getNodeElevations(junctionIndex)
            %   d.plot;
            %
            % Example 4:
            %   % Adds a new junction with coordinates [X, Y] = [10, 10], elevation = 500 and demand = 50.
            %   junctionID = 'newJunction_4';
            %   junctionCoords = [10 10];
            %   junctionElevation = 500;
            %   demand = 50;
            %   junctionIndex = d.addNodeJunction(junctionID, junctionCoords, junctionElevation, demand);
            %   d.getNodeBaseDemands{1}(junctionIndex)
            %   d.plot;
            %
            % Example 5:
            %   % Adds a new junction with coordinates [X, Y] = [10, 20], elevation = 500, demand = 50 and pattern ID = the 1st time pattern ID(if exists).
            %   junctionID = 'newJunction_5';
            %   junctionCoords = [10 20];
            %   junctionElevation = 500;
            %   demand = 50;
            %   demandPatternID = d.getPatternNameID{1};
            %   junctionIndex = d.addNodeJunction(junctionID, junctionCoords, junctionElevation, demand, demandPatternID);
            %   d.getNodeDemandPatternNameID{1}(junctionIndex)
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addNodeReservoir,
            %          setNodeComment, deleteNode, setNodeBaseDemands.
            xy = [0 0];
            elev = 0;
            dmnd = 0;
            dmndpat = '';
            if nargin >= 3
                xy = varargin{1};
            end
            if nargin >= 4
                elev = varargin{2};
            end
            if nargin >= 5
                dmnd = varargin{3};
            end
            if nargin == 6
                dmndpat = varargin{4};
            end
            index = ENaddnode(obj, juncID, obj.ToolkitConstants.EN_JUNCTION);
            obj.setNodeCoordinates(index, [xy(1),xy(2)]);
            obj.setNodeJunctionData(index, elev, dmnd, dmndpat);
        end
        function index = addNodeReservoir(obj, resID, varargin)
            % Adds a new reservoir.
            % Returns the index of the new reservoir.
            %
            % Example 1:
            %   % Adds a new reservoir with the default coordinates (i.e. [0, 0])
            %   reservoirID = 'newReservoir_1';
            %   reservoirIndex = d.addNodeReservoir(reservoirID);
            %   d.plot
            %
            % Example 2:
            %   % Adds a new reservoir with coordinates [X, Y] = [20, 30].
            %   reservoirID = 'newReservoir_2';
            %   reservoirCoords = [20 30];
            %   reservoirIndex = d.addNodeReservoir(reservoirID, reservoirCoords);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addNodeJunction,
            %          addLinkPipe, deleteNode, setNodeBaseDemands.
            xy = [0 0];
            if nargin == 3
                xy = varargin{1};
            end
            index = ENaddnode(obj, resID, obj.ToolkitConstants.EN_RESERVOIR);
            obj.setNodeCoordinates(index,[xy(1),xy(2)]);
        end
        function index = addNodeTank(obj, tankID, varargin)
            % Adds a new tank.
            % Returns the index of the new tank.
            %
            % Example 1:
            %   % Adds a new tank with the default coordinates (i.e. [0, 0])
            %   tankID = 'newTank_1';
            %   tankIndex = d.addNodeTank(tankID);
            %   d.plot
            %
            % Example 2:
            %   % Adds a new tank with coordinates [X, Y] = [10, 10].
            %   tankID = 'newTank_2';
            %   tankCoords = [10 10];
            %   tankIndex = d.addNodeTank(tankID, tankCoords);
            %   d.plot
            %
            % Example 3:
            %   % Adds a new tank with coordinates [X, Y] = [20, 20] and elevation = 100.
            %   tankID = 'newTank_3';
            %   tankCoords = [20 20];
            %   elevation = 100;
            %   tankIndex = d.addNodeTank(tankID, tankCoords, elevation);
            %   d.plot;
            %
            % Example 4:
            %   % Adds a new tank with coordinates [X, Y] = [20, 30], elevation = 100, initial level = 130, minimum water level = 110,
            %   % maximum water level = 160, diameter = 60, minimum water volume = 200000, volume curve ID = '';.
            %   tankID = 'newTank_4';
            %   tankCoords = [20 30];
            %   elevation = 100;
            %   initialLevel = 130;
            %   minimumWaterLevel = 110;
            %   maximumWaterLevel = 160;
            %   diameter = 60;
            %   minimumWaterVolume = 200000;
            %   volumeCurveID = '';   % Empty for no curve
            %   tankIndex = d.addNodeTank(tankID, tankCoords, elevation, initialLevel, minimumWaterLevel, ...
            %   maximumWaterLevel, diameter, minimumWaterVolume, volumeCurveID);
            %   d.getNodeTankData(tankIndex)
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addNodeJunction,
            %          addLinkPipe, deleteNode, setNodeBaseDemands.
            xy = [0 0];
            elev = 0;
            intlvl = 0;
            minlvl =  0;
            maxlvl = 0;
            diam = 0;
            minvol = 0;
            volcurve = '';
            if nargin >= 3
                xy = varargin{1};
            end
            if nargin >= 4
                elev = varargin{2};
            end
            if nargin >= 5
                intlvl = varargin{3};
            end
            if nargin >= 6
                minlvl = varargin{4};
            end
            if nargin >= 7
                maxlvl = varargin{5};
            end
            if nargin >= 8
                diam = varargin{6};
            end
            if nargin >= 9
                minvol = varargin{7};
            end
            if nargin >= 10
                volcurve = varargin{8};
            end
            index = ENaddnode(obj, tankID, obj.ToolkitConstants.EN_TANK);
            obj.setNodeCoordinates(index,[xy(1),xy(2)]);
            if diam == 0 && strcmp(obj.getNodeType(index), 'TANK')
                 minvol = (pi * (diam/2)^2) *minlvl;
                 if minvol == 0
                     return;
                 end
            end
            obj.setNodeTankData(index, elev, intlvl, minlvl, maxlvl, diam, minvol, volcurve)
        end
        function index = addLinkPipeCV(obj, cvpipeID, fromNode, toNode, varargin)
            % Adds a new control valve pipe.
            % Returns the index of the new control valve pipe.
            %
            % % Properties that can be set(optional):
            % 1) Length
            % 2) Diameter
            % 3) Roughness Coefficient
            % 4) Minor Loss Coefficient
            %
            % % If no properties are given, the default values are:
            %   length = 330 feet (~100.5 m)
            %   diameter = 10 inches (25.4 cm)
            %   roughness coefficient = 130 (Hazen-Williams formula) or
            %                           0.15 mm (Darcy-Weisbach formula) or
            %                           0.01 (Chezy-Manning formula)
            %   minor Loss Coefficient = 0
            %
            % % The examples are based on d = epanet('NET1.inp');
            %
            % Example 1:
            %   % Adds a new control valve pipe given no properties.
            %   cvPipeID = 'newCVPipe_1';
            %   fromNode = '10';
            %   toNode = '21';
            %   d.getLinkPipeCount                     % Retrieves the number of pipes
            %   cvPipeIndex = d.addLinkPipeCV(cvPipeID, fromNode, toNode);
            %   d.getLinkPipeCount
            %   d.plot;                                % Plots the network in a new MATLAB figure
            %
            % Example 2:
            %   % Adds a new control valve pipe given it's length.
            %   cvPipeID = 'newCVPipe_2';
            %   fromNode = '11';
            %   toNode = '22';
            %   length = 600;
            %   d.getLinkPipeCount
            %   cvPipeIndex = d.addLinkPipeCV(cvPipeID, fromNode, toNode, length);
            %   d.getLinkPipeCount
            %   d.getLinkLength(cvPipeIndex)           % Retrieves the new link's length
            %   d.plot;
            %
            % Example 3:
            %   % Adds a new control valve pipe given it's length, diameter, roughness coefficient and minor loss coefficient.
            %   cvPipeID = 'newCVPipe_3';
            %   fromNode = '31';
            %   toNode = '22';
            %   length = 500;
            %   diameter = 15;
            %   roughness = 120;
            %   minorLossCoeff = 0.2;
            %   d.getLinkPipeCount
            %   cvPipeIndex = d.addLinkPipeCV(cvPipeID, fromNode, toNode, length, diameter, roughness, minorLossCoeff);
            %   d.getLinkPipeCount
            %   d.getLinkLength(cvPipeIndex)
            %   d.getLinkDiameter(cvPipeIndex)         % Retrieves the new link's diameter
            %   d.getLinkRoughnessCoeff(cvPipeIndex)   % Retrieves the new link's roughness coefficient
            %   d.getLinkMinorLossCoeff(cvPipeIndex)   % Retrieves the new link's minor loss coefficient
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addNodeJunction, deleteLink, setLinkDiameter.
            index = ENaddlink(obj, cvpipeID, obj.ToolkitConstants.EN_CVPIPE, fromNode, toNode);
            if nargin >= 5
                obj.setLinkLength(index, varargin{1});
            end
            if nargin >= 6
                obj.setLinkDiameter(index, varargin{2});
            end
            if nargin >= 7
                obj.setLinkRoughnessCoeff(index, varargin{3});
            end
            if nargin == 8
                obj.setLinkMinorLossCoeff(index, varargin{4});
            end
        end
        function index = addLinkPipe(obj, pipeID, fromNode, toNode, varargin)
            % Adds a new pipe.
            % Returns the index of the new pipe.
            %
            % % Properties that can be set(optional):
            % 1) Length
            % 2) Diameter
            % 3) Roughness Coefficient
            % 4) Minor Loss Coefficient
            %
            % % If no properties are given, the default values are:
            %   length = 330 feet (~100.5 m)
            %   diameter = 10 inches (25.4 cm)
            %   roughness coefficient = 130 (Hazen-Williams formula) or
            %                           0.15 mm (Darcy-Weisbach formula) or
            %                           0.01 (Chezy-Manning formula)
            %   minor Loss Coefficient = 0
            %
            % % The examples are based on d = epanet('NET1.inp');
            %
            % Example 1:
            %   % Adds a new pipe given no properties.
            %   pipeID = 'newPipe_1';
            %   fromNode = '10';
            %   toNode = '21';
            %   d.getLinkPipeCount                   % Retrieves the number of links
            %   pipeIndex = d.addLinkPipe(pipeID, fromNode, toNode);
            %   d.getLinkPipeCount
            %   d.plot;                              % Plots the network in a new MATLAB figure
            %
            % Example 2:
            %   % Adds a new pipe given it's length.
            %   pipeID = 'newPipe_2';
            %   fromNode = '11';
            %   toNode = '22';
            %   length = 600;
            %   d.getLinkPipeCount
            %   pipeIndex = d.addLinkPipe(pipeID, fromNode, toNode, length);
            %   d.getLinkPipeCount
            %   d.getLinkLength(pipeIndex)           % Retrieves the new link's length
            %   d.plot;
            %
            % Example 3:
            %   % Adds a new pipe given it's length, diameter, roughness coefficient and minor loss coefficient.
            %   pipeID = 'newPipe_3';
            %   fromNode = '31';
            %   toNode = '22';
            %   length = 500;
            %   diameter = 15;
            %   roughness = 120;
            %   minorLossCoeff = 0.2;
            %   d.getLinkPipeCount
            %   pipeIndex = d.addLinkPipe(pipeID, fromNode, toNode, length, diameter, roughness, minorLossCoeff);
            %   d.getLinkPipeCount
            %   d.getLinkLength(pipeIndex)
            %   d.getLinkDiameter(pipeIndex)         % Retrieves the new link's diameter
            %   d.getLinkRoughnessCoeff(pipeIndex)   % Retrieves the new link's roughness coefficient
            %   d.getLinkMinorLossCoeff(pipeIndex)   % Retrieves the new link's minor loss coefficient
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addLinkPipeCV,
            %          addNodeJunction, deleteLink, setLinkDiameter.
            index = ENaddlink(obj, pipeID, obj.ToolkitConstants.EN_PIPE, fromNode, toNode);
            if nargin >= 5
                obj.setLinkLength(index, varargin{1});
            end
            if nargin >= 6
                obj.setLinkDiameter(index, varargin{2});
            end
            if nargin >= 7
                obj.setLinkRoughnessCoeff(index, varargin{3});
            end
            if nargin == 8
                obj.setLinkMinorLossCoeff(index, varargin{4});
            end
        end
        function index = addLinkPump(obj, pumpID, fromNode, toNode, varargin)
            % Adds a new pump.
            % Returns the index of the new pump.
            %
            % % Properties that can be set(optional):
            % 1) Initial Status
            % 2) Initial Speed setting
            % 3) Power
            % 4) Pattern index
            %
            % % If no properties are given, the default values are:
            %   initial status = 1 (OPEN)
            %   initial speed setting = 1
            %   power = 0
            %   pattern index = 0
            %
            % % The examples are based on d=epanet('NET1.inp');
            %
            % Example 1:
            %   % Adds a new pump given no properties.
            %   pumpID = 'newPump_1';
            %   fromNode = '10';
            %   toNode = '21';
            %   d.getLinkPumpCount                     % Retrieves the number of pumps
            %   pumpIndex = d.addLinkPump(pumpID, fromNode, toNode);
            %   d.getLinkPumpCount
            %   d.plot;                                % Plots the network in a new MATLAB figure
            %
            % Example 2:
            %   % Adds a new pump given it's initial status.
            %   pumpID = 'newPump_2';
            %   fromNode = '31';
            %   toNode = '22';
            %   initialStatus = 0;   % (CLOSED)
            %   d.getLinkPumpCount
            %   pumpIndex = d.addLinkPump(pumpID, fromNode, toNode, initialStatus);
            %   d.getLinkPumpCount
            %   d.getLinkInitialStatus(pumpIndex)      % Retrieves the new pump's initial status
            %   d.plot;
            %
            % Example 3:
            %   % Adds a new pump given it's initial status, initial speed setting, power and pattern index.
            %   pumpID = 'newPump_3';
            %   fromNode = '11';
            %   toNode = '22';
            %   initialStatus = 1;   % (OPEN)
            %   initialSetting = 1.2;
            %   power = 10;
            %   patternIndex = 1;
            %   d.getLinkPumpCount
            %   pumpIndex = d.addLinkPump(pumpID, fromNode, toNode, initialStatus, initialSetting, power, patternIndex);
            %   d.getLinkPumpCount
            %   d.getLinkInitialStatus(pumpIndex)
            %   d.getLinkInitialSetting(pumpIndex)     % Retrieves the new pump's initial setting
            %   d.getLinkPumpPower(pumpIndex)          % Retrieves the new pump's power
            %   d.getLinkPumpPatternIndex(pumpIndex)   % Retrieves the new pump's pattern index
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addNodeJunction, deleteLink, setLinkInitialStatus.
            index = ENaddlink(obj, pumpID, obj.ToolkitConstants.EN_PUMP, fromNode, toNode);
            if nargin >= 5
                obj.setLinkInitialStatus(index, varargin{1});
            end
            if nargin >= 6
                obj.setLinkInitialSetting(index, varargin{2});
            end
            if nargin >= 7
                obj.setLinkPumpPower(index, varargin{3});
            end
            if nargin == 8
                obj.setLinkPumpPatternIndex(index, varargin{4});
            end
        end
        function index = addLinkValve(obj, vID, fromNode, toNode, varargin)
            % Adds a new valve.
            % Returns the index of the new valve.
            %
            % % Properties that can be set(optional):
            % 1) Type: PRV (Pressure reducing valve)
            %          PSV (Pressure sustaining valve)
            %          PBV (Pressure breaker valve)
            %          FCV (Flow control valve)
            %          TCV (Throttle control valve)
            %          GPV (General purpose valve)
            % 2) Diameter
            % 3) Initial Setting
            % 4) Minor Loss Coefficient
            %
            % % If no properties are given, the default values are:
            %   type = 'GPV'
            %   diameter = 10 inches (25.4 cm)
            %   initial setting = 0
            %   minor Loss Coefficient = 0
            %
            % % The examples are based on d = epanet('NET1.inp');
            %
            % Example 1:
            %   % Adds a new valve given no properties.
            %   valveID = 'newValve_1';
            %   fromNode = '10';
            %   toNode = '21';
            %   d.getLinkValveCount                  % Retrieves the number of valves
            %   valveIndex = d.addLinkValve(valveID, fromNode, toNode);
            %   d.getLinkValveCount
            %   d.plot;                              % Plots the network in a new MATLAB figure
            %
            % Example 2:
            %   % Adds a new valve given it's type.
            %   valveID = 'newValve_2';
            %   fromNode = '11';
            %   toNode = '22';
            %   type = 'PRV';
            %   d.getLinkValveCount
            %   valveIndex = d.addLinkValve(valveID, fromNode, toNode, type);
            %   d.getLinkValveCount
            %   d.getLinkType(valveIndex)           % Retrieves the new valve's type
            %   d.plot;
            %
            % Example 3:
            %   % Adds a new valve given it's type, diameter, initial setting and minor loss coefficient.
            %   valveID = 'newValve_3';
            %   fromNode = '31';
            %   toNode = '22';
            %   type = 'FCV';
            %   diameter = 15;
            %   initialSetting = 1;
            %   minorLossCoeff = 0.2;
            %   d.getLinkValveCount
            %   valveIndex = d.addLinkValve(valveID, fromNode, toNode, type, diameter, initialSetting, minorLossCoeff);
            %   d.getLinkValveCount
            %   d.getLinkType(valveIndex)
            %   d.getLinkDiameter(valveIndex)         % Retrieves the new valve's diameter
            %   d.getLinkInitialSetting(valveIndex)   % Retrieves the new valve's initial setting
            %   d.getLinkMinorLossCoeff(valveIndex)   % Retrieves the new valve's minor loss coefficient
            %   d.plot;
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addNodeJunction, deleteLink, setLinkDiameter.
            if nargin >= 5
                vtype = eval(['obj.ToolkitConstants.EN_', varargin{1}]);
            else
                vtype = obj.ToolkitConstants.EN_GPV;
            end
            index = ENaddlink(obj, vID, vtype, fromNode, toNode);
            if nargin >= 6
                obj.setLinkDiameter(index, varargin{2});
            end
            if nargin >= 7
                obj.setLinkInitialSetting(index, varargin{3});
            end
            if nargin == 8
                obj.setLinkMinorLossCoeff(index, varargin{4});
            end
        end
        function index = addLinkValvePRV(obj, vID, fromNode, toNode)
            % Adds a new PRV valve.
            % Returns the index of the new PRV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValvePRV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValvePRV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePSV, deleteLink, setLinkTypeValveFCV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_PRV, fromNode, toNode);
        end
        function index = addLinkValvePSV(obj, vID, fromNode, toNode)
            % Adds a new PSV valve.
            % Returns the index of the new PSV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValvePSV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValvePSV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePRV, deleteLink, setLinkTypeValveGPV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_PSV, fromNode, toNode);
        end
        function index = addLinkValvePBV(obj, vID, fromNode, toNode)
            % Adds a new PBV valve.
            % Returns the index of the new PBV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValvePBV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValvePBV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePRV, deleteLink, setLinkTypeValvePRV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_PBV, fromNode, toNode);
        end
        function index = addLinkValveFCV(obj, vID, fromNode, toNode)
            % Adds a new FCV valve.
            % Returns the index of the new FCV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValveFCV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValveFCV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePRV, deleteLink, setLinkTypeValveTCV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_FCV, fromNode, toNode);
        end
        function index = addLinkValveTCV(obj, vID, fromNode, toNode)
            % Adds a new TCV valve.
            % Returns the index of the new TCV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValveTCV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValveTCV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePRV, deleteLink, setLinkTypeValveFCV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_TCV, fromNode, toNode);
        end
        function index = addLinkValveGPV(obj, vID, fromNode, toNode)
            % Adds a new GPV valve.
            % Returns the index of the new GPV valve.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   valveID = 'newValveGPV';
            %   fromNode = '10';
            %   toNode = '21';
            %   valveIndex = d.addLinkValveGPV(valveID, fromNode, toNode);
            %   d.plot
            %
            % See also plot, setLinkNodesIndex, addLinkPipe,
            %          addLinkValvePRV, deleteLink, setLinkTypeValveFCV.
            index = ENaddlink(obj, vID, obj.ToolkitConstants.EN_GPV, fromNode, toNode);
        end
        function value = getLinkVerticesCount(obj, varargin)
            % Retrieves the number of internal vertex points assigned to a link.
            %
            % Example:
            %   d = epanet('anytown.inp');
            %   d.getLinkVerticesCount
            %
            %   d.getBinLinkVerticesCount
            %   link_id = '10';
            %   d.getLinkVerticesCount(link_id)
            %
            % See also getLinkVertices, setLinkVertices, getBinLinkVertices,
            %          getBinLinkVerticesCount.
            if nargin == 2
                if ischar(varargin{1})
                    varargin{1} = obj.getLinkIndex(varargin{1});
                end
            end
            indices = getLinkIndices(obj, varargin);
            j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetvertexcount(i, obj.LibEPANET);
                j = j +1;
                error(obj.getError(obj.Errcode));
            end
        end
        function data = getLinkVertices(obj, varargin)
            % Retrieves the coordinate's of a vertex point assigned to a link.
            %
            % % The example is based on d=epanet('Net1.inp');
            %
            % Example:
            %   d = epanet('Net1.inp');
            %   linkID = '10';
            %   x = [22, 24, 28];
            %   y = [69, 68, 69];
            %   d.setLinkVertices(linkID, x, y)
            %   linkID = '112';
            %   x = [10, 24, 18];
            %   y = [49, 58, 60];
            %   d.setLinkVertices(linkID, x, y)
            %   d.getLinkVertices{1}
            %   d.getLinkVertices{2}
            %
            % See also getLinkVertices, getLinkVerticesCount, getBinLinkVertices,
            %          getBinLinkVerticesCount.
            data = {};
            if nargin == 2
                if ischar(varargin{1})
                    indices = obj.getLinkIndex(varargin{1});
                end
            else
                indices = obj.getLinkIndex;
            end
            j=1;
            for l=indices
                if obj.getLinkVerticesCount(l) == 0
                     data{j} = [];
                end
                for i=1:obj.getLinkVerticesCount(l)
                    [obj.Errcode, data{j}.x(i), data{j}.y(i)]=ENgetvertex(l, i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
                j = j +1;
            end
        end
        function setLinkVertices(obj, linkID, x, y, varargin)
            % Assigns a set of internal vertex points to a link.
            %
            % % The example is based on d=epanet('Net1.inp');
            %
            % Example:
            %   d = epanet('Net1.inp');
            %   linkID = '10';
            %   x = [22, 24, 28];
            %   y = [69, 68, 69];
            %   d.setLinkVertices(linkID, x, y)
            % See also getLinkVertices, getLinkVerticesCount, getBinLinkVertices,
            %          getBinLinkVerticesCount.
            index = obj.getLinkIndex(linkID);
            [obj.Errcode]=ENsetvertices(index, x, y, size(x, 2), obj.LibEPANET);
            error(obj.getError(obj.Errcode));

        end
        function Errcode = deleteNode(obj, idNode, varargin)
            % Deletes nodes. (EPANET Version 2.2)
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Deletes all controls, rules and links related to the object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object deletion if contained in controls, rules and links
            % Default condition is 0.
            %
            % Example 1:
            %   d.getNodeCount                   % Retrieves the ID label of all nodes
            %   idNode = d.getNodeNameID(1);      % Retrieves the ID label of the 1st node
            %   d.deleteNode(idNode)              % Deletes the 1st node given it's ID
            %   d.getNodeCount
            %
            % Example 2:
            %   idNode = d.getNodeNameID(1);
            %   condition = 1;
            %   d.deleteNode(idNode, condition)   % Attempts to delete a node connected to links (error occurs)
            %
            % Example 3:
            %   index = 1;
            %   d.deleteNode(index)               % Deletes the 1st node given it's index
            %   d.getNodeNameID
            %
            % Example 4:
            %   idNodes = d.getNodeNameID(1:2);
            %   d.getNodeCount
            %   d.deleteNode(idNodes)             % Deletes 2 nodes given their IDs
            %   d.getNodeCount
            %
            % See also addNodeJunction, deleteLink, deleteRules,
            %          setNodeCoordinates, setNodeJunctionData.
            condition = 0;
            if nargin == 3
                condition = varargin{1};
            end
            if ischar(idNode)
                idNode = {idNode};
            end
            if iscell(idNode)
                for j = 1:length(idNode)
                    indexNode = obj.getNodeIndex(idNode(j));
                    [Errcode] = ENdeletenode(obj.LibEPANET, indexNode, condition);
                    error(obj.getError(Errcode));
                end
            else
                [Errcode] = ENdeletenode(obj.LibEPANET, idNode, condition);
                error(obj.getError(Errcode));
            end
        end
        function Errcode = deleteLink(obj, idLink, varargin)
            % Deletes a link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Deletes all controls and rules related to the object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object deletion if contained in controls and rules
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkNameID                     % Retrieves the ID label of all links
            %   idLink = d.getLinkNameID(1);        % Retrieves the ID label of the 1st link
            %   d.deleteLink(idLink)                % Deletes the 1st link given it's ID
            %   d.getLinkNameID
            %
            % Example 2:
            %   idLink = d.getLinkPumpNameID{1};
            %   condition = 1;
            %   d.deleteLink(idLink, condition)     % Attempts to delete a link contained in controls (error occurs)
            %
            % Example 3:
            %   indexLink = 1;
            %   d.deleteLink(indexLink)             % Deletes the 1st link given it's index
            %   d.getLinkNameID
            %
            % See also addLinkPipe, deleteNode, deleteRules,
            %          setNodeCoordinates, setLinkPipeData.
            condition = 0;
            if nargin == 3
                condition = varargin{1};
            end
            if ischar(idLink) || iscell(idLink)
                indexLink = obj.getLinkIndex(idLink);
            else
                indexLink = idLink;
            end
            [Errcode] = ENdeletelink(obj.LibEPANET, indexLink, condition);
            error(obj.getError(Errcode));
            %if obj.Bin, obj.Errcode = reloadNetwork(obj); end
        end
        function Errcode = deletePattern(obj, idPat)
            % Deletes a time pattern from a project.
            %
            % Example 1:
            %   idPat = d.getPatternNameID(1);   % Retrieves the ID of the 1st pattern
            %   d.deletePattern(idPat)           % Deletes the 1st pattern given it's ID
            %   d.getPatternNameID
            %
            % Example 2:
            %   index = 1;
            %   d.deletePattern(index)           % Deletes the 1st pattern given it's index
            %   d.getPatternNameID
            %
            % See also addPattern, setPattern, setPatternNameID,
            %          setPatternValue, setPatternComment.
            if ischar(idPat) || iscell(idPat)
                indexPat = obj.getPatternIndex(idPat);
            else
                indexPat = idPat;
            end
            [Errcode] = ENdeletepattern(obj.LibEPANET, indexPat);
            error(obj.getError(Errcode));
        end
        function Errcode = deleteCurve(obj, idCurve)
            % Deletes a data curve from a project.
            %
            % Example 1:
            %   idCurve = d.getCurveNameID(1);   % Retrieves the ID of the 1st curve
            %   d.deleteCurve(idCurve)           % Deletes a curve given it's ID
            %   d.getCurveNameID
            %
            % Example 2:
            %   index = 1;
            %   d.deleteCurve(index)             % Deletes a curve given it's index
            %   d.getCurveNameID
            %
            % See also addCurve, setCurve, setCurveNameID,
            %          setCurveValue, setCurveComment.
            if ischar(idCurve) || iscell(idCurve)
                indexCurve = obj.getCurveIndex(idCurve);
            else
                indexCurve = idCurve;
            end
            [Errcode] = ENdeletecurve(obj.LibEPANET, indexCurve);
            error(obj.getError(Errcode));
        end
        function setControls(obj, index, control, varargin)
            % Sets the parameters of a simple control statement.
            %
            % % The examples are based on d=epanet('Net1.inp');
            %
            % Example 1:
            %   controlIndex = 1;
            %   d.getControls(controlIndex)            % Retrieves the 1st control
            %   control = 'LINK 9 CLOSED IF NODE 2 ABOVE 180';
            %   d.setControls(controlIndex, control)   % Sets a control given it's index and the control statement
            %   d.getControls(1)
            %
            % Example 2:
            %   controls = d.getControls;
            %   d.setControls(controls)                % Sets multiple controls given as structs with fields
            %
            % Example 3:
            %   control_1 = 'LINK 9 OPEN IF NODE 2 BELOW 110';
            %   control_2 = 'LINK 9 CLOSED IF NODE 2 ABOVE 200';
            %   controls = {control_1, control_2};
            %   d.setControls(controls)                % Sets multiple controls given as cell
            %   d.getControls(1)
            %   d.getControls(2)
            % Example 4:
            %  Notes:
            %       index:     control statement index
            %       control:   control type code
            %       lindex:    index of link being controlled
            %       setting:   value of the control setting
            %       nindex:    index of controlling node
            %       level:     value of controlling water level or pressure for
            %                  level controls or of time of control action
            %                  (in seconds) for time-based controls
            %
            % Control type codes consist of the following:
            %  EN_LOWLEVEL      0   Control applied when tank level or node pressure drops below specified level
            %  EN_HILEVEL       1   Control applied when tank level or node pressure rises above specified level
            %  EN_TIMER         2   Control applied at specific time into simulation
            %  EN_TIMEOFDAY     3   Control applied at specific time of day
            %
            %   Code example:
            %       % d.setControls(index, control, lindex, setting, nindex, level)
            %       d.setControls(1, 0, 13, 0, 11, 30)
            %
            % See also getControls, getControlRulesCount,
            %          addControls, deleteControls.
            if isstruct(index)
                for c=1:length(index)
                    setControlFunction(obj, c, index(c).Control)
                end
            else
                if nargin <= 3
                    if length(index)> 1
                        tmpC = index;
                        index = 1:length(index);
                    elseif length(index)== 1
                        if isnumeric(index)
                            tmpC{1} = control;
                        else
                            tmpC{1} = index{1};
                            index = 1;
                        end
                    end
                    j=1;
                    for i=index
                        setControlFunction(obj, i, tmpC{j});
                        j=j+1;
                    end
                else
                    linkIndex = varargin{1};
                    controlSettingValue = varargin{2};
                    nodeIndex = varargin{3};
                    controlLevel = varargin{4};
                    [obj.Errcode] = ENsetcontrol(index, control, linkIndex, controlSettingValue, nodeIndex, controlLevel, obj.LibEPANET);
                end
            end
        end
        function index = addControls(obj, control, varargin)
            % Adds a new simple control. (EPANET Version 2.2)
            %
            % % The examples are based on d=epanet('Net1.inp');
            %
            % Example 1:
            %   % Close Link 12 if the level in Tank 2 exceeds 20 ft.
            %   index = d.addControls('LINK 12 CLOSED IF NODE 2 ABOVE 20');
            %   d.getControls(index)
            %
            % Example 2:
            %   % Open Link 12 if the pressure at Node 11 is under 30 psi
            %   index = d.addControls('LINK 12 OPEN IF NODE 11 BELOW 30');
            %   d.getControls(index)
            %
            % Example 3:
            %   % Pump 9 speed is set to 1.5 at 16 hours or 57600 seconds into the simulation
            %   index = d.addControls('LINK 9 1.5 AT TIME 16:00');
            %   d.getControls(index)
            %   index = d.addControls('LINK 9 1.5 AT TIME 57600'); %in seconds
            %   d.getControls(index)
            %
            % Example 4:
            %   % Link 12 is closed at 10 am and opened at 8 pm throughout the simulation
            %   index_3 = d.addControls('LINK 12 CLOSED AT CLOCKTIME 10:00');
            %   d.getControls(index_3)
            %   index_4 = d.addControls('LINK 12 OPEN AT CLOCKTIME 20:00');
            %   d.getControls(index_4)
            %
            % Example 5:
            %   % Adds multiple controls given as cell
            %   control_1 = 'LINK 9 OPEN IF NODE 2 BELOW 110';
            %   control_2 = 'LINK 9 CLOSED IF NODE 2 ABOVE 200';
            %   controls = {control_1, control_2};
            %   index = d.addControls(controls)
            %   d.getControls(index)
            %
            % Example 6:
            % Notes:
            %    index:	    return index of the new control.
            %	 type:  	the type of control to add (see EN_ControlType).
            %    linkIndex:	the index of a link to control (starting from 1).
            %    setting:	control setting applied to the link.
            %    nodeIndex:	index of the node used to control the link (0 for EN_TIMER and EN_TIMEOFDAY controls).
            %    level:	    action level (tank level, junction pressure, or time in seconds) that triggers the control.
            %
            % Control type codes consist of the following:
            % EN_LOWLEVEL      0   Control applied when tank level or node pressure drops below specified level
            % EN_HILEVEL       1   Control applied when tank level or node pressure rises above specified level
            % EN_TIMER         2   Control applied at specific time into simulation
            % EN_TIMEOFDAY     3   Control applied at specific time of day
            %
            % Code example:
            % % index = d.addControls(type, linkIndex, setting, nodeIndex, level)
            % index = d.addControls(0, 13, 0, 11, 100)
            % d.getControls(index)
            %
            % See also deleteControls, getControls,
            %          setControls, getControlRulesCount.
            if iscell(control)
                for i=1:length(control)
                    index(i) = addControlFunction(obj, control{i});
                end
            else
                if nargin <= 3
                   index = addControlFunction(obj, control);
                else
                    linkIndex = varargin{1};
                    controlSettingValue = varargin{2};
                    nodeIndex = varargin{3};
                    controlLevel = varargin{4};
                    [obj.Errcode, index] = ENaddcontrol(control, linkIndex, controlSettingValue, nodeIndex, controlLevel, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function Errcode = deleteControls(varargin)
            % Deletes an existing simple control. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getControls                                                 % Retrieves the parameters of all control statements
            %   d.deleteControls                                              % Deletes the existing simple controls
            %   d.getControls
            %
            % Example 2:
            %   index = d.addControls('LINK 9 43.2392 AT TIME 4:00:00');      % Adds a new simple control(index = 3)
            %   d.getControls(index)
            %   d.deleteControls(index);                                      % Deletes the 3rd simple control
            %   d.getControls
            %
            % Example 3:
            %   index_3 = d.addControls('LINK 9 43.2392 AT TIME 4:00:00');    % Adds a new simple control(index = 3)
            %   index_4 = d.addControls('LINK 10 43.2392 AT TIME 4:00:00');   % Adds a new simple control(index = 4)
            %   d.getControls(index_3)
            %   d.getControls(index_4)
            %   d.deleteControls([index_3, index_4]);                         % Deletes the 3rd and 4th simple controls
            %   d.getControls
            %
            % See also addControls, setControls,
            %          getControls, getControlRulesCount.
            obj = varargin{1};
            if nargin==1
                index = 1:obj.getControlRulesCount;
            else
                index = varargin{2};
            end
            for i=length(index):-1:1
                [Errcode] = ENdeletecontrol(index(i), obj.LibEPANET);
                error(obj.getError(Errcode));
            end
        end
        function setLinkPipeData(obj, Index, Length, Diameter, RoughnessCoeff, MinorLossCoeff)
            % Sets a group of properties for a pipe. (EPANET Version 2.2)
            %
            % Properties:
            % 1) Pipe Index
            % 2) Length
            % 3) Diameter
            % 4) Roughness Coefficient
            % 5) Minor Loss Coefficient
            %
            % Example 1:
            %   % Sets to the 1st pipe the following properties:
            %   pipeIndex = 1;
            %   length = 1000;
            %   diameter = 20;
            %   RoughnessCoeff = 110;
            %   MinorLossCoeff = 0.2;
            %   d.getLinksInfo   % Retrieves all link info
            %   d.setLinkPipeData(pipeIndex, length, diameter, RoughnessCoeff, MinorLossCoeff)
            %   d.getLinksInfo
            %
            % Example 2:
            %   % Sets to the 1st and 2nd pipe the following properties:
            %   pipeIndex = [1, 2];
            %   length = [1000, 1500];
            %   diameter = [20, 23];
            %   RoughnessCoeff = [110, 115];
            %   MinorLossCoeff = [0.2, 0.3];
            %   d.getLinksInfo   % Retrieves all link info
            %   d.setLinkPipeData(pipeIndex, length, diameter, RoughnessCoeff, MinorLossCoeff)
            %   d.getLinksInfo
            %
            % See also getLinksInfo, setLinkComment, setLinkDiameter,
            %          setLinkLength, setLinkStatus, setNodeTankData.
            for i=1:length(Index)
                [obj.Errcode] = ENsetpipedata(Index(i), Length(i), Diameter(i), RoughnessCoeff(i), MinorLossCoeff(i), obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkNodesIndex(obj, linkIndex, startNode, endNode)
            % Sets the indexes of a link's start- and end-nodes. (EPANET Version 2.2)
            %
            % Example 1:
            %   % Sets to the 1st link the start-node index = 2 and end-node index = 3
            %   d.getLinkNodesIndex   % Retrieves the indexes of the from/to nodes of all links
            %   linkIndex = 1;
            %   startNode = 2;
            %   endNode = 3;
            %   d.setLinkNodesIndex(linkIndex, startNode, endNode)
            %   d.getLinkNodesIndex
            %
            % Example 2:
            %   % Sets to the 1st link the start-node index = 2 and end-node index = 3
            %   % and to 2nd link the start-node index = 4 and end-node index = 5.
            %   linkIndex = [1 ,2];
            %   startNode = [2, 4];
            %   endNode = [3, 5];
            %   d.setLinkNodesIndex(linkIndex, startNode, endNode)
            %   d.getLinkNodesIndex
            %
            % See also getLinkNodesIndex, setLinkDiameter, setLinkLength,
            %          setLinkNameID, setLinkComment.
            for i=1:length(linkIndex)
                [obj.Errcode] = ENsetlinknodes(linkIndex(i), startNode(i), endNode(i), obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodesConnectingLinksID(obj, linkIndex, startNodeID, endNodeID)
            % Sets the IDs of a link's start- and end-nodes. (EPANET Version 2.2)
            %
            % Example 1:
            %   d.getNodesConnectingLinksID   % Retrieves the ids of the from/to nodes of all links
            %   linkIndex = 2;
            %   startNodeID = '11';
            %   endNodeID = '22';
            %   d.setNodesConnectingLinksID(linkIndex, startNodeID, endNodeID)
            %   d.getNodesConnectingLinksID
            %
            % Example 2:
            %   linkIndex = [2, 3];
            %   startNodeID = {'12', '13'};
            %   endNodeID = {'21', '22'};
            %   d.setNodesConnectingLinksID(linkIndex, startNodeID, endNodeID)
            %   d.getNodesConnectingLinksID
            %
            % See also getLinkNodesIndex, getNodesConnectingLinksID, setLinkLength,
            %          setLinkNameID, setLinkComment.
            startNode = obj.getNodeIndex(startNodeID);
            endNode = obj.getNodeIndex(endNodeID);
            obj.setLinkNodesIndex(linkIndex, startNode, endNode)
        end
        function setLinkDiameter(obj, value, varargin)
            % Sets the values of diameters.
            %
            % Example 1:
            %   d.getLinkDiameter                            % Retrieves the diameters of all links
            %   index_pipe = 1;
            %   diameter = 20;
            %   d.setLinkDiameter(index_pipe, diameter);     % Sets the diameter of the 1st pipe
            %   d.getLinkDiameter(index_pipe)
            %
            % Example 2:
            %   index_pipes = [1, 2];
            %   diameters = [20, 25];
            %   d.setLinkDiameter(index_pipes, diameters);   % Sets the diameters of the first 2 pipes
            %   d.getLinkDiameter(index_pipes)
            %
            % Example 3:
            %   diameters = d.getLinkDiameter;
            %   diameters = diameters*1.5;
            %   d.setLinkDiameter(diameters)                 % Sets the diameters of all the links
            %   d.getLinkDiameter
            %
            % See also setLinkPipeData, setLinkLength,
            %          setLinkBulkReactionCoeff, setLinkTypePipe.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_DIAMETER, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function value = getLinkComment(obj, varargin)
            % Retrieves the comment string assigned to the link object.
            %
            % Example 1:
            %   d.getLinkComment              % Retrieves the comments of all links
            %
            % Example 2:
            %   linkIndex = 1;
            %   d.getLinkComment(linkIndex)   % Retrieves the comment of the 1st link
            %
            % Example 3:
            %   linkIndex = 1:5;
            %   d.getLinkComment(linkIndex)   % Retrieves the comments of the first 5 links
            %
            % See also setLinkComment, getLinkNameID, getLinksInfo.
            if isempty(varargin)
                cnt = obj.getLinkCount;
                value = cell(1, cnt);
                for i=1:cnt
                    [obj.Errcode, value{i}]=ENgetcomment(obj.ToolkitConstants.EN_LINK, i, obj.LibEPANET);
                end
            else
                if isempty(varargin{1}), varargin{1}=0; end
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, value{k}]=ENgetcomment(obj.ToolkitConstants.EN_LINK, i, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                    k=k+1;
                end
            end
        end
        function value = setLinkComment(obj, value, varargin)
            % Sets the comment string assigned to the link object.
            %
            % Example 1:
            %   linkIndex = 1;
            %   d.getLinkComment(linkIndex)
            %   comment = 'This is a link';
            %   d.setLinkComment(linkIndex, comment);   % Sets a comment to the 1st link
            %   d.getLinkComment(linkIndex)
            %
            % Example 2:
            %   linkIndex = [1, 2];
            %   d.getLinkComment(linkIndex)
            %   comment = {'This is link 1', 'This is link 2'};
            %   d.setLinkComment(linkIndex, comment);   % Sets comments to the first 2 links
            %   d.getLinkComment(linkIndex)
            %
            % See also getLinkComment, setLinkNameID, setLinkPipeData.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            if length(indices) == 1
                [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_LINK, indices, value, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                for i=indices
                    [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_LINK, i, value{j}, obj.LibEPANET); j=j+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function index = setLinkTypePipe(obj, id, varargin)
            % Sets the link type pipe for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType
            %   linkid = d.getLinkPumpNameID{1};
            %   index = d.setLinkTypePipe(linkid);             % Changes the 1st pump to pipe given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPumpNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypePipe(linkid, condition)   % Changes the 1st pump to pipe given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPumpNameID, setLinkTypePipeCV,
            %          setLinkTypePump, setLinkTypeValveFCV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_PIPE, condition, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function index = setLinkTypePipeCV(obj, id, varargin)
            % Sets the link type cvpipe(pipe with check valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                  % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                  % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypePipeCV(linkid);              % Changes the 1st pipe to cvpipe given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypePipeCV(linkid, condition);   % Changes the 1st pipe to cvpipe given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValveFCV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_CVPIPE, condition, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function index = setLinkTypePump(obj, id, varargin)
            % Sets the link type pump for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypePump(linkid);              % Changes the 1st pipe to pump given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypePump(linkid, condition);   % Changes the 1st pipe to pump given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePipeCV, setLinkTypeValveFCV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_PUMP, condition, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function index = setLinkTypeValveFCV(obj, id, varargin)
            % Sets the link type valve FCV(flow control valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValveFCV(linkid);              % Changes the 1st pipe to valve FCV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValveFCV(linkid, condition);   % Changes the 1st pipe to valve FCV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValveGPV.
            condition = 0; % default
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_FCV, condition, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function index = setLinkTypeValveGPV(obj, id, varargin)
            % Sets the link type valve GPV(general purpose valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValveGPV(linkid);              % Changes the 1st pipe to valve GPV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValveGPV(linkid, condition);   % Changes the 1st pipe to valve GPV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValveFCV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_GPV, condition, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function index = setLinkTypeValvePBV(obj, id, varargin)
            % Sets the link type valve PBV(pressure breaker valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValvePBV(linkid);              % Changes the 1st pipe to valve PBV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValvePBV(linkid, condition);   % Changes the 1st pipe to valve PBV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValvePRV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_PBV, condition, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function index = setLinkTypeValvePRV(obj, id, varargin)
            % Sets the link type valve PRV(pressure reducing valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValvePRV(linkid);              % Changes the 1st pipe to valve PRV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValvePRV(linkid, condition);   % Changes the 1st pipe to valve PRV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValvePSV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_PRV, condition, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function index = setLinkTypeValvePSV(obj, id, varargin)
            % Sets the link type valve PSV(pressure sustaining valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValvePSV(linkid);              % Changes the 1st pipe to valve PSV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValvePSV(linkid, condition);   % Changes the 1st pipe to valve PSV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValvePBV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_PSV, condition, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function index = setLinkTypeValveTCV(obj, id, varargin)
            % Sets the link type valve TCV(throttle control valve) for a specified link.
            %
            % condition = 0 | if is EN_UNCONDITIONAL: Delete all controls that contain object
            % condition = 1 | if is EN_CONDITIONAL: Cancel object type change if contained in controls
            % Default condition is 0.
            %
            % Example 1:
            %   d.getLinkType(1)                                    % Retrieves the type of the 1st link
            %   linkid = d.getLinkPipeNameID{1};                    % Retrieves the ID of the 1t pipe
            %   index = d.setLinkTypeValveTCV(linkid);              % Changes the 1st pipe to valve TCV given it's ID
            %   d.getLinkType(index)
            %
            % Example 2:
            %   linkid = d.getLinkPipeNameID{1};
            %   condition = 1;
            %   index = d.setLinkTypeValveTCV(linkid, condition);   % Changes the 1st pipe to valve TCV given it's ID and a condition (if possible)
            %   d.getLinkType(index)
            %
            % See also getLinkType, getLinkPipeNameID, setLinkTypePipe,
            %          setLinkTypePump, setLinkTypeValveGPV.
            condition = 0; % default
            if nargin == 3
                condition = varargin{1};
            end
            index = obj.check_if_numeric(id);
            [obj.Errcode, index] = ENsetlinktype(index, obj.ToolkitConstants.EN_TCV, condition, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setLinkLength(obj, value, varargin)
            % Sets the values of lengths.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkLength(index_pipe)                 % Retrieves the length of the 1st link
            %   length_pipe = 100;
            %   d.setLinkLength(index_pipe, length_pipe);   % Sets the length of the 1st link
            %   d.getLinkLength(index_pipe)
            %
            % Example 2:
            %   lengths = d.getLinkLength;                  % Retrieves the lengths of all the links
            %   lengths_new = lengths*1.5;
            %   d.setLinkLength(lengths_new);               % Sets the new lengths of all links
            %   d.getLinkLength
            %
            % See also getLinkLength, setLinkDiameter, setLinkMinorLossCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_LENGTH, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function value = setLinkNameID(obj, value, varargin)
            % Sets the ID name for links.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkNameID(index_pipe)            % Retrieves the ID of the 1st link
            %   linkID = 'New_ID';                     % ID selected without a space in between the letters
            %   d.setLinkNameID(index_pipe, linkID);   % Sets the ID name of the 1st link
            %   d.getLinkNameID(index_pipe)
            %(the size of the cell must equal to the number of links)
            % Example 2:
            %   IDs = {'1', '2', '3', '4'};            % Select the IDs of all the links (the size of the cell must equal the number of links)
            %   d.setLinkNameID(IDs);                  % Sets the ID names of all links
            %   d.getLinkNameID
            %
            % See also getLinkNameID, setLinkComment, setLinkDiameter,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            if length(indices) == 1
                [obj.Errcode] = ENsetlinkid(indices, value, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                for i=indices
                    [obj.Errcode] = ENsetlinkid(i, value{j}, obj.LibEPANET); j=j+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function setLinkRoughnessCoeff(obj, value, varargin)
            % Sets the values of roughness coefficient.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkRoughnessCoeff(index_pipe)           % Retrieves the roughness coefficient of the 1st link
            %   coeff = 105;
            %   d.setLinkRoughnessCoeff(index_pipe, coeff);   % Sets the roughness coefficient of the 1st link
            %   d.getLinkRoughnessCoeff(index_pipe)
            %
            % Example 2:
            %   coeffs = d.getLinkRoughnessCoeff;             % Retrieves the roughness coefficients of all the links
            %   coeffs_new = coeffs + 10;
            %   d.setLinkRoughnessCoeff(coeffs_new);          % Sets the roughness coefficient of all links
            %   d.getLinkRoughnessCoeff
            %
            % See also getLinkRoughnessCoeff, setLinkDiameter, setLinkMinorLossCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_ROUGHNESS, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkMinorLossCoeff(obj, value, varargin)
            % Sets the values of minor loss coefficient.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkMinorLossCoeff(index_pipe)           % Retrieves the minor loss coefficient of the 1st link
            %   coeff = 105;
            %   d.setLinkMinorLossCoeff(index_pipe, coeff);   % Sets the minor loss coefficient of the 1st link
            %   d.getLinkMinorLossCoeff(index_pipe)
            %
            % Example 2:
            %   coeffs = d.getLinkMinorLossCoeff;             % Retrieves the minor loss coefficients of all the links
            %   coeffs_new = coeffs + 0.2;
            %   d.setLinkMinorLossCoeff(coeffs_new);          % Sets the minor loss coefficient of all links
            %   d.getLinkMinorLossCoeff
            %
            % See also getLinkMinorLossCoeff, setLinkDiameter, setLinkRoughnessCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_MINORLOSS, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkInitialStatus(obj, value, varargin)
            % Sets the values of initial status.
            %
            % Note: Cannot set status for a check valve
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkInitialStatus(index_pipe)            % Retrieves the initial status of the 1st link
            %   status = 0;
            %   d.setLinkInitialStatus(index_pipe, status);   % Sets the initial status of the 1st link
            %   d.getLinkInitialStatus(index_pipe)
            %
            % Example 2:
            %   statuses = d.getLinkInitialStatus;            % Retrieves the initial status of all links
            %   statuses_new = zeros(1, length(statuses));
            %   d.setLinkInitialStatus(statuses_new);         % Sets the initial status of all links
            %   d.getLinkInitialStatus
            %
            % See also getLinkInitialStatus, setLinkInitialSetting, setLinkDiameter,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_INITSTATUS, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkInitialSetting(obj, value, varargin)
            % Sets the values of initial settings, roughness for pipes or initial speed for pumps or initial setting for valves.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkInitialSetting(index_pipe)             % Retrieves the initial setting of the 1st link
            %   setting = 80;
            %   d.setLinkInitialSetting(index_pipe, setting);   % Sets the initial setting of the 1st link
            %   d.getLinkInitialSetting(index_pipe)
            %
            % Example 2:
            %   settings = d.getLinkInitialSetting;             % Retrieves the initial setting of all links
            %   settings_new = settings + 40;
            %   d.setLinkInitialSetting(settings_new);          % Sets the initial setting of all links
            %   d.getLinkInitialSetting
            %
            % See also getLinkInitialSetting, setLinkInitialStatus, setLinkRoughnessCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_INITSETTING, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkBulkReactionCoeff(obj, value, varargin)
            % Sets the value of bulk chemical reaction coefficient.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkBulkReactionCoeff(index_pipe)           % Retrieves the bulk chemical reaction coefficient of the 1st link
            %   coeff = 0;
            %   d.setLinkBulkReactionCoeff(index_pipe, coeff);   % Sets the bulk chemical reaction coefficient of the 1st link
            %   d.getLinkBulkReactionCoeff(index_pipe)
            %
            % Example 2:
            %   coeffs = d.getLinkBulkReactionCoeff;             % Retrieves the bulk chemical reaction coefficients of all links
            %   coeffs_new = 0*coeffs;
            %   d.setLinkBulkReactionCoeff(coeffs_new);          % Sets the bulk chemical reaction coefficient of all links
            %   d.getLinkBulkReactionCoeff
            %
            % See also getLinkBulkReactionCoeff, setLinkRoughnessCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_KBULK, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkWallReactionCoeff(obj, value, varargin)
            % Sets the value of wall chemical reaction coefficient.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkWallReactionCoeff(index_pipe)           % Retrieves the wall chemical reaction coefficient of the 1st link
            %   coeff = 0;
            %   d.setLinkWallReactionCoeff(index_pipe, coeff);   % Sets the wall chemical reaction coefficient of the 1st link
            %   d.getLinkWallReactionCoeff(index_pipe)
            %
            % Example 2:
            %   coeffs = d.getLinkWallReactionCoeff;             % Retrieves the wall chemical reaction coefficients of all links
            %   coeffs_new = 0*coeffs;
            %   d.setLinkWallReactionCoeff(coeffs_new);          % Sets the wall chemical reaction coefficient of all links
            %   d.getLinkWallReactionCoeff
            %
            % See also getLinkWallReactionCoeff, setLinkBulkReactionCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_KWALL, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkStatus(obj, value, varargin)
            % Sets the values of current status for links.
            %
            % Note: Cannot set status for a check valve
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkStatus(index_pipe)            % Retrieves the current status of the 1st link
            %   status = 1;
            %   d.setLinkStatus(index_pipe, status);   % Sets the current status of the 1st link
            %   d.getLinkStatus(index_pipe)
            %
            % Example 2:
            %   statuses = d.getLinkStatus;            % Retrieves the current status of all links
            %   statuses_new = zeros(1, length(statuses));
            %   d.setLinkStatus(statuses_new);         % Sets the current status of all links
            %   d.getLinkStatus
            %
            % See also getLinkStatus, setLinkInitialStatus, setLinkDiameter,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_STATUS, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkSettings(obj, value, varargin)
            % Sets the values of current settings, roughness for pipes or initial speed for pumps or initial setting for valves.
            %
            % Example 1:
            %   index_pipe = 1;
            %   d.getLinkSettings(index_pipe)             % Retrieves the current setting of the 1st link
            %   setting = 80;
            %   d.setLinkSettings(index_pipe, setting);   % Sets the current setting of the 1st link
            %   d.getLinkSettings(index_pipe)
            %
            % Example 2:
            %   settings = d.getLinkSettings;             % Retrieves the current setting of all links
            %   settings_new = settings + 40;
            %   d.setLinkSettings(settings_new);          % Sets the current setting of all links
            %   d.getLinkSettings
            %
            % See also getLinkSettings, setLinkStatus, setLinkRoughnessCoeff,
            %          setLinkPipeData, addLink, deleteLink.
            if nargin==3, indices = value; value=varargin{1}; else indices = getLinkIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, obj.ToolkitConstants.EN_SETTING, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setLinkPumpPower(obj, value, varargin)
            % Sets the power for pumps. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   d.getLinkPumpPower                       % Retrieves the power of all pumps
            %   d.setLinkPumpPower(10)                   % Sets the pump power = 10 to every pump
            %   d.getLinkPumpPower
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpPower([10, 15])             % Sets the pump power = 10 and 15 to the 2 pumps
            %   d.getLinkPumpPower
            %
            % Example 3:
            %   d.setLinkPumpPower(1, 10)                % Sets the pump power = 10 to the 1st pump
            %   d.getLinkPumpPower
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpPower(pumpIndex, 10)        % Sets the pump power = 10 to the pumps with index 118 and 119
            %   d.getLinkPumpPower
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpPower(pumpIndex,[10, 15])   % Sets the pump power = 10 and 15 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpPower
            %
            % See also getLinkPumpPower, setLinkPumpHCurve, setLinkPumpECurve,
            %          setLinkPumpECost, setLinkPumpEPat.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_PUMP_POWER, value, varargin)
        end
        function setLinkPumpHCurve(obj, value, varargin)
            % Sets the pump head v. flow curve index. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   d.getLinkPumpHCurve                     % Retrieves the pump head v. flow curve index of all pumps
            %   d.setLinkPumpHCurve(1)                  % Sets the pump head v. flow curve index = 1 to every pump
            %   d.getLinkPumpHCurve
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpHCurve([1, 2])             % Sets the pump head v. flow curve index = 1 and 2 to the 2 pumps
            %   d.getLinkPumpHCurve
            %
            % Example 3:
            %   d.setLinkPumpHCurve(1, 2)               % Sets the pump head v. flow curve index = 2 to the 1st pump
            %   d.getLinkPumpHCurve
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpHCurve(pumpIndex, 1)       % Sets the pump head v. flow curve index = 1 to the pumps with index 118 and 119
            %   d.getLinkPumpHCurve
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpHCurve(pumpIndex,[1, 2])   % Sets the pump head v. flow curve index = 1 and 2 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpHCurve
            %
            % See also getLinkPumpHCurve, setLinkPumpPower, setLinkPumpECurve,
            %          setLinkPumpECost, setLinkPumpEPat.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_PUMP_HCURVE, value, varargin)
        end
        function setLinkPumpECurve(obj, value, varargin)
            % Sets the pump efficiency v. flow curve index. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   d.getLinkPumpECurve                     % Retrieves the pump efficiency v. flow curve index of all pumps
            %   d.setLinkPumpECurve(1)                  % Sets the pump efficiency v. flow curve index = 1 to every pump
            %   d.getLinkPumpECurve
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpECurve([1, 2])             % Sets the pump efficiency v. flow curve index = 1 and 2 to the 2 pumps
            %   d.getLinkPumpECurve
            %
            % Example 3:
            %   d.setLinkPumpECurve(1, 2)               % Sets the pump efficiency v. flow curve index = 2 to the 1st pump
            %   d.getLinkPumpECurve
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpECurve(pumpIndex, 1)       % Sets the pump efficiency v. flow curve index = 1 to the pumps with index 118 and 119
            %   d.getLinkPumpECurve
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpECurve(pumpIndex,[1, 2])   % Sets the pump efficiency v. flow curve index = 1 and 2 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpECurve
            %
            % See also getLinkPumpECurve, setLinkPumpPower, setLinkPumpHCurve,
            %          setLinkPumpECost, setLinkPumpEPat.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_PUMP_ECURVE, value, varargin)
        end
        function setLinkPumpECost(obj, value, varargin)
            % Sets the pump average energy price. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   d.getLinkPumpECost                            % Retrieves the pump average energy price of all pumps
            %   d.setLinkPumpECost(0.10)                      % Sets the pump average energy price = 0.10 to every pump
            %   d.getLinkPumpECost
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpECost([0.10, 0.12])              % Sets the pump average energy price = 0.10 and 0.12 to the 2 pumps
            %   d.getLinkPumpECost
            %
            % Example 3:
            %   d.setLinkPumpECost(1, 0.10)                   % Sets the pump average energy price = 0.10 to the 1st pump
            %   d.getLinkPumpECost
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpECost(pumpIndex, 0.10)           % Sets the pump average energy price = 0.10 to the pumps with index 118 and 119
            %   d.getLinkPumpECost
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpECost(pumpIndex, [0.10, 0.12])   % Sets the pump average energy price = 0.10 and 0.12 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpECost
            %
            % See also getLinkPumpECost, setLinkPumpPower, setLinkPumpHCurve,
            %          setLinkPumpECurve, setLinkPumpEPat.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_PUMP_ECOST, value, varargin)
        end
        function setLinkPumpEPat(obj, value, varargin)
            % Sets the pump energy price time pattern index. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   d.getLinkPumpEPat                     % Retrieves the pump energy price time pattern index of all pumps
            %   d.setLinkPumpEPat(1)                  % Sets the pump energy price time pattern index = 1 to every pump
            %   d.getLinkPumpEPat
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpEPat([1, 2])             % Sets the pump energy price time pattern index = 1 and 2 to the 2 pumps
            %   d.getLinkPumpEPat
            %
            % Example 3:
            %   d.setLinkPumpEPat(1, 2)               % Sets the pump energy price time pattern index = 2 to the 1st pump
            %   d.getLinkPumpEPat
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpEPat(pumpIndex, 1)       % Sets the pump energy price time pattern index = 1 to the pumps with index 118 and 119
            %   d.getLinkPumpEPat
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpEPat(pumpIndex,[1, 2])   % Sets the pump energy price time pattern index = 1 and 2 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpEPat
            %
            % See also getLinkPumpEPat, setLinkPumpPower, setLinkPumpHCurve,
            %          setLinkPumpECurve, setLinkPumpECost.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_PUMP_EPAT, value, varargin)
        end
        function setLinkPumpPatternIndex(obj, value, varargin)
            % Sets the pump speed time pattern index. (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('Net3_trace.inp');
            %
            % Example 1:
            %   d.getLinkPumpPatternIndex                     % Retrieves the pump speed time pattern index of all pumps
            %   d.setLinkPumpPatternIndex(1)                  % Sets the speed time pattern index = 1 to every pump
            %   d.getLinkPumpPatternIndex
            %
            % Example 2:
            %   % The input array must have a length equal to the number of pumps
            %   d.setLinkPumpPatternIndex([1, 2])             % Sets the pump speed time pattern index = 1 and 2 to the 2 pumps
            %   d.getLinkPumpPatternIndex
            %
            % Example 3:
            %   d.setLinkPumpPatternIndex(1, 2)               % Sets the pump speed time pattern index = 2 to the 1st pump
            %   d.getLinkPumpPatternIndex
            %
            % Example 4:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpPatternIndex(pumpIndex, 1)       % Sets the pump speed time pattern index = 1 to the pumps with index 118 and 119
            %   d.getLinkPumpPatternIndex
            %
            % Example 5:
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpPatternIndex(pumpIndex, [1, 2])   % Sets the pump speed time pattern index = 1 and 2 to the pumps with index 118 and 119 respectively
            %   d.getLinkPumpPatternIndex
            %
	    % Example 6:
            %   % To remove the pattern index from the pumps you can use input 0
            %   pumpIndex = d.getLinkPumpIndex;
	    %   d.setLinkPumpPatternIndex(pumpIndex, 0)
	    %
            % See also getLinkPumpPatternIndex, setLinkPumpPower, setLinkPumpHCurve,
            %          setLinkPumpECurve, setLinkPumpECost.
            set_Node_Link(obj, 'pump', 'ENsetlinkvalue', obj.ToolkitConstants.EN_LINKPATTERN, value, varargin)
        end
        function value = setLinkPumpHeadCurveIndex(obj, value, varargin)
	    % Example 1:
	    %   To remove curve index from the pumps you can use input 0
            %   pumpIndex = d.getLinkPumpIndex;
            %   d.setLinkPumpHeadCurveIndex(pumpIndex, 0)
	    %
            % See also setLinkPumpPatternIndex, getLinkPumpPower, setLinkPumpHCurve,
            %          setLinkPumpECurve, setLinkPumpECost.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetheadcurveindex(obj.LibEPANET, i, value(j)); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeElevations(obj, value, varargin)
            % Sets the values of elevation for nodes.
            %
            % Example 1:
            %   index_node = 1;
            %   d.getNodeElevations(index_node)          % Retrieves the elevation of the 1st node
            %   elev = 500;
            %   d.setNodeElevations(index_node, elev);   % Sets the elevation of the 1st node
            %   d.getNodeElevations(index_node)
            %
            % Example 2:
            %   elevs = d.getNodeElevations;             % Retrieves the elevations of all the nodes
            %   elevs_new = elevs + 100;
            %   d.setNodeElevations(elevs_new);          % Sets the elevations of all nodes
            %   d.getNodeElevations
            %
            % See also getNodeElevations, setNodeCoordinates, setNodeBaseDemands,
            %          setNodeJunctionData, addNodeJunction, deleteNode.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, obj.ToolkitConstants.EN_ELEVATION, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeBaseDemands(obj, value, varargin)
            % Sets the values of demand for nodes.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   index_node = 1;
            %   d.getNodeBaseDemands{1}(index_node)                       % Retrieves the demand of the 1st node
            %   demand = 5;
            %   d.setNodeBaseDemands(index_node, demand);                 % Sets the demand of the 1st node
            %   d.getNodeBaseDemands{1}(index_node)
            %
            % Example 2:
            %   nodeIndex = 1:5;
            %   d.getNodeBaseDemands{1}(nodeIndex)                        % Retrieves the demands of first 5 nodes
            %   demands = [10, 5, 15, 20, 5];
            %   d.setNodeBaseDemands(nodeIndex, demands)                  % Sets the demands of first 5 nodes
            %   d.getNodeBaseDemands{1}(nodeIndex)
            %
            % Example 3:
            %   demands = d.getNodeBaseDemands{1};                        % Retrieves the demands of all nodes
            %   demands_new = demands + 15;
            %   d.setNodeBaseDemands(demands_new);                        % Sets the demands of all nodes
            %   d.getNodeBaseDemands{1}
            %
            % For the following examples EPANET Version 2.1 or higher is required.
            %
            % If a category is not given, the default is categoryIndex = 1.
            %   d = epanet('BWSN_Network_1.inp');
            %
            % Example 4:
            %   nodeIndex = 121;
            %   categoryIndex = 2;
            %   d.getNodeBaseDemands{categoryIndex}(nodeIndex)            % Retrieves the demand of the 2nd category of the 121th node
            %   demand = 25;
            %   d.setNodeBaseDemands(nodeIndex, categoryIndex, demand)    % Sets the demand of the 2nd category of the 121th node
            %   d.getNodeBaseDemands{categoryIndex}(nodeIndex)
            %
            % Example 5:
            %   nodeIndex = 1:5;
            %   categoryIndex = 1;
            %   d.getNodeBaseDemands{categoryIndex}(nodeIndex)            % Retrieves the demands of the 1st category of the first 5 nodes
            %   demands = [10, 5, 15, 20, 5];
            %   d.setNodeBaseDemands(nodeIndex, categoryIndex, demands)   % Sets the demands of the 1st category of the first 5 nodes
            %   d.getNodeBaseDemands{categoryIndex}(nodeIndex)
            %
            % See also getNodeBaseDemands, setNodeJunctionDemandName,
            %          setNodeDemandPatternIndex, addNodeJunction, deleteNode.
            set_node_demand_pattern(obj, 'ENsetbasedemand', obj.ToolkitConstants.EN_BASEDEMAND, value, varargin)
        end
        function setNodeCoordinates(obj, value, varargin)
            % Sets node coordinates.
            %
            % Example 1:
            %   nodeIndex = 1;
            %   d.getNodeCoordinates{1}(nodeIndex)             % Retrieves the X coordinates of the 1st node
            %   d.getNodeCoordinates{2}(nodeIndex)             % Retrieves the Y coordinates of the 1st node
            %   coords = [0, 0];
            %   d.setNodeCoordinates(nodeIndex, coords)        % Sets the coordinates of the 1st node
            %   d.getNodeCoordinates{1}(nodeIndex)
            %   d.getNodeCoordinates{2}(nodeIndex)
            %
            % Example 2:
            %   x_values = d.getNodeCoordinates{1};
            %   y_values = d.getNodeCoordinates{2};
            %   new_coords = {x_values + 10, y_values + 10};   % Creates a cell array with the new coordinates
            %   d.setNodeCoordinates(new_coords)               % Sets the coordinates of all nodes
            %   d.getNodeCoordinates{1}
            %   d.getNodeCoordinates{2}
            %
            % See also getNodeCoordinates, setNodeElevations, plot,
            %          addNodeJunction, addNodeTank, deleteNode.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            if ~isempty(varargin)
                for i=indices
                    [obj.Errcode] = ENsetcoord(i, value(1), value(2), obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            else
                for i=1:length(value{1})
                    x=value{1}(i);
                    y=value{2}(i);
                    [obj.Errcode] = ENsetcoord(i, x, y, obj.LibEPANET);
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function setNodeDemandPatternIndex(obj, value, varargin)
            % Sets the values of demand time pattern indices.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   nodeIndex = 1;
            %   d.getNodeDemandPatternIndex{1}(nodeIndex)                               % Retrieves the index of the 1st category's time pattern of the 1st node
            %   patternIndex = 2;
            %   d.setNodeDemandPatternIndex(nodeIndex, patternIndex)                    % Sets the demand time pattern index to the 1st node
            %   d.getNodeDemandPatternIndex{1}(nodeIndex)
            %
            % Example 2:
            %   nodeIndex = 1:5;
            %   d.getNodeDemandPatternIndex{1}(nodeIndex)
            %   patternIndices = [1, 3, 2, 4, 2];
            %   d.setNodeDemandPatternIndex(nodeIndex, patternIndices)                  % Sets the demand time pattern index to the first 5 nodes
            %   d.getNodeDemandPatternIndex{1}(nodeIndex)
            %
            % Example 3:
            %   patternIndices = d.getNodeDemandPatternIndex{1};
            %   patternIndices_new = patternIndices + 1;
            %   d.setNodeDemandPatternIndex(patternIndices_new)                         % Sets all primary demand time pattern indices
            %   d.getNodeDemandPatternIndex{1}
            %
            % For the following examples EPANET Version 2.1 or higher is required.
            %
            % If a category is not given, the default is categoryIndex = 1.
            %
            % Example 4:
            %   nodeIndex = 121;
            %   categoryIndex = 2;
            %   d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)                   % Retrieves the index of the 2nd category's time pattern of the 121th node
            %   patternIndex = 4;
            %   d.setNodeDemandPatternIndex(nodeIndex, categoryIndex, patternIndex)     % Sets the demand time pattern index of the 2nd category of the 121th node
            %   d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
            %
            % Example 5:
            %   nodeIndex = 1:5;
            %   categoryIndex = 1;
            %   d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
            %   patternIndices = [1, 3, 2, 4, 2];
            %   d.setNodeDemandPatternIndex(nodeIndex, categoryIndex, patternIndices)   % Sets the demand time pattern index of the 1st category of the first 5 nodes
            %   d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
            %
            % See also getNodeDemandPatternIndex, getNodeDemandCategoriesNumber
            %          setNodeBaseDemands, addPattern, deletePattern.
            set_node_demand_pattern(obj, 'ENsetdemandpattern', obj.ToolkitConstants.EN_PATTERN, value, varargin)
        end
        function setNodeEmitterCoeff(obj, value, varargin)
            % Sets the values of emitter coefficient for nodes.
            %
            % Example 1:
            %   nodeset = d.getNodeEmitterCoeff                  % Retrieves the value of all nodes emmitter coefficients
            %   nodeset(1) = 0.1;                                % First node emitter coefficient = 0.1
            %   d.setNodeEmitterCoeff(nodeset)                   % Sets the value of all nodes emitter coefficient
            %   d.getNodeEmitterCoeff
            %
            % Example 2:
            %   nodeIndex = 1;
            %   d.getNodeEmitterCoeff(nodeIndex)
            %   emitterCoeff = 0;
            %   d.setNodeEmitterCoeff(nodeIndex, emitterCoeff)   % Sets the value of the 1st node emitter coefficient = 0
            %   d.getNodeEmitterCoeff(nodeIndex)
            %
            % See also getNodeEmitterCoeff, setNodeBaseDemands, setNodeJunctionData.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, obj.ToolkitConstants.EN_EMITTER, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeInitialQuality(obj, value, varargin)
            % Sets the values of initial quality for nodes.
            %
            % Example 1:
            %   nodeset = d.getNodeInitialQuality                    % Retrieves the value of all nodes initial qualities
            %   nodeset(1) = 0.5;                                    % First node initial quality = 0.5
            %   d.setNodeInitialQuality(nodeset)                     % Sets the values of all nodes initial quality
            %   d.getNodeInitialQuality
            %
            % Example 2:
            %   nodeIndex = 1;
            %   d.getNodeInitialQuality(nodeIndex)
            %   initialQuality = 1;
            %   d.setNodeInitialQuality(nodeIndex, initialQuality)   % Sets the value of the 1st node initial quality
            %   d.getNodeInitialQuality(nodeIndex)
            %
            % See also getNodeInitialQuality, getNodeActualQuality, setNodeJunctionData.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, obj.ToolkitConstants.EN_INITQUAL, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeNameID(obj, value, varargin)
            % Sets the ID name for nodes.
            %
            % Example 1:
            %   nodeIndex = 1;
            %   d.getNodeNameID(nodeIndex)            % Retrieves the ID of the 1st node
            %   nameID = 'newID';
            %   d.setNodeNameID(nodeIndex, nameID);   % Sets the ID of the 1st node.
            %   d.getNodeNameID(nodeIndex)
            %
            % Example 2:
            %   nameID = d.getNodeNameID;             % Retrieves the IDs of all nodes
            %   nameID{1} = 'newID_1';
            %   nameID{5} = 'newID_5';
            %   d.setNodeNameID(nameID);              % Sets the IDs of all nodes
            %   d.getNodeNameID
            %
            % See also getNodeNameID, setNodeComment, setNodeJunctionData.
            if nargin==3
                indices = value;
                value = {};
                value{indices}=varargin{1};
            else
                indices = find(~strcmp(obj.getNodeNameID, value));
            end
            for i=indices
                [obj.Errcode] = obj.ENsetnodeid(i, value{i}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeTankData(obj, index, elev, intlvl, minlvl, maxlvl, diam, minvol, volcurve)
            % Sets a group of properties for a tank. (EPANET Version 2.2)
            %
            % Properties:
	        % 1) Tank index
            % 2) Elevation
            % 3) Initial water Level
            % 4) Minimum Water Level
            % 5) Maximum Water Level
            % 6) Diameter (0 if a volume curve is supplied)
            % 7) Minimum Water Volume
            % 8) Volume Curve Index ("" for no curve)
            %
            % The examples are based on d=epanet('Net3_trace.inp')
            %
            % Example 1:
            %   % Sets to the 1st tank the following properties:
            %   tankIndex = 1;   % You can also use tankIndex = 95; (i.e. the index of the tank).
            %   elev = 100;
            %   intlvl = 13;
            %   minlvl =  0.2;
            %   maxlvl = 33;
            %   diam = 80;
            %   minvol = 50000;
            %   volcurve = '';   % For no curve
            %   d.setNodeTankData(tankIndex, elev, intlvl, minlvl, maxlvl, diam, minvol, volcurve)
            %   d.getNodeTankData
            %
            % Example 2:
            %   % Sets to the 1st and 2nd tank the following properties:
            %   tankIndex = [1, 2];   % You can also use tankIndex = [95, 96]; (i.e. the indices of the tanks).
            %   elev = [100, 105];
            %   intlvl = [13, 13.5];
            %   minlvl =  [0.2, 0.25];
            %   maxlvl = [30, 35];
            %   diam = [80, 85];
            %   minvol = [50000, 60000];
            %   volcurve = {'', ''};   % For no curves
            %   d.setNodeTankData(tankIndex, elev, intlvl, minlvl, maxlvl, diam, minvol, volcurve)
            %   d.getNodeTankData
            %
            % See also getNodeTankData, setNodeTankInitialLevel,
            %          setNodeTankMinimumWaterLevel, setNodeTankDiameter.
            if ischar(volcurve)
                volcurve={volcurve};
            else
                volcurve={''};
            end

            if ~ismember(index, obj.getNodeTankIndex)
                tankIndices = obj.getNodeTankIndex;
                index = tankIndices(index);
            end
            for i=1:length(index)
                [obj.Errcode] = ENsettankdata(index(i), elev(i), intlvl(i), minlvl(i), maxlvl(i), diam(i), minvol(i), volcurve{i}, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeTankInitialLevel(obj, value, varargin)
            % Sets the values of initial level for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankInitialLevel                       % Retrieves the initial level of all tanks
            %   d.setNodeTankInitialLevel(10)                   % Sets the initial level = 10 to every tank
            %   d.getNodeTankInitialLevel
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankInitialLevel([10, 15])             % Sets the initial level = 10 and 15 to the 2 tanks
            %   d.getNodeTankInitialLevel
            %
            % Example 3:
            %   d.setNodeTankInitialLevel(1, 10)                % Sets the initial level = 10 to the 1st tank
            %   d.getNodeTankInitialLevel
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankInitialLevel(tankIndex, 10)        % Sets the initial level = 10 to the tanks with index 128 and 129
            %   d.getNodeTankInitialLevel
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankInitialLevel(tankIndex,[10, 15])   % Sets the initial level = 10 and 15 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankInitialLevel
            %
            % See also getNodeTankInitialLevel, setNodeTankMinimumWaterLevel, setNodeTankMaximumWaterLevel,
            %          setNodeTankMinimumWaterVolume, setNodeTankMixingFraction, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_TANKLEVEL, value, varargin)
        end
        function setNodeTankMixingModelType(obj, value, varargin)
            % Sets the mixing model type value for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankMixingModelType                                % Retrieves the  mixing model type of all tanks
            %   d.setNodeTankMixingModelType('MIX2')                        % Sets the  mixing model type = 'MIX2' to every tank
            %   d.getNodeTankMixingModelType
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankMixingModelType({'MIX1', 'LIFO'})              % Sets the  mixing model type = 'MIX1' and 'LIFO' to the 2 tanks
            %   d.getNodeTankMixingModelType
            %
            % Example 3:
            %   d.setNodeTankMixingModelType(1, 'FIFO')                     % Sets the  mixing model type = 'FIFO' to the 1st tank
            %   d.getNodeTankMixingModelType
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMixingModelType(tankIndex, 'MIX1')             % Sets the  mixing model type = 'MIX1' to the tanks with index 128 and 129
            %   d.getNodeTankMixingModelType
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMixingModelType(tankIndex, {'MIX2', 'LIFO'})   % Sets the  mixing model type = 'MIX2' and 'LIFO' to the tanks with index 128 and 129 respectively
            %   d.getNodeTankMixingModelType
            %
            % See also getNodeTankMixingModelType, setNodeTankBulkReactionCoeff, setNodeTankMixingFraction,
            %          setNodeTankMinimumWaterVolume, setNodeTankMinimumWaterLevel, setNodeTankData.
            if nargin == 2
                type = value;
            elseif nargin == 3
                type = varargin{1};
            end
            if iscell(type)
                code = zeros(1, length(type));
                for i = 1:length(type)
                    code(i)=strfind(strcmpi(type{i}, obj.TYPEMIXMODEL), 1)-1;
                end
            elseif ischar(type)
                code=strfind(strcmpi(type, obj.TYPEMIXMODEL), 1)-1;
            end
            if nargin == 2
                value = code;
            elseif nargin == 3
                varargin{1} = code;
            end
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_MIXMODEL, value, varargin)
        end
        function setNodeTankDiameter(obj, value, varargin)
            % Sets the diameter value for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankDiameter                         % Retrieves the diameter of all tanks
            %   d.setNodeTankDiameter(120)                    % Sets the diameter = 120 to every tank
            %   d.getNodeTankDiameter
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankDiameter([110, 130])             % Sets the diameter = 110 and 130 to the 2 tanks
            %   d.getNodeTankDiameter
            %
            % Example 3:
            %   d.setNodeTankDiameter(1, 120)                 % Sets the diameter = 120 to the 1st tank
            %   d.getNodeTankDiameter
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankDiameter(tankIndex, 150)         % Sets the diameter = 150 to the tanks with index 128 and 129
            %   d.getNodeTankDiameter
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankDiameter(tankIndex,[100, 120])   % Sets the diameter = 100 and 120 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankDiameter
            %
            % See also getNodeTankDiameter, setNodeTankInitialLevel, setNodeTankMinimumWaterLevel,
            %          setNodeTankBulkReactionCoeff, setNodeTankCanOverFlow, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_TANKDIAM, value, varargin)
        end
        function setNodeTankMinimumWaterLevel(obj, value, varargin)
            % Sets the minimum water level value for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankMinimumWaterLevel                      % Retrieves the minimum water level of all tanks
            %   d.setNodeTankMinimumWaterLevel(5)                   % Sets the minimum water level = 5 to every tank
            %   d.getNodeTankMinimumWaterLevel
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankMinimumWaterLevel([10, 15])            % Sets the minimum water level = 10 and 15 to the 2 tanks
            %   d.getNodeTankMinimumWaterLevel
            %
            % Example 3:
            %   d.setNodeTankMinimumWaterLevel(1, 5)                % Sets the minimum water level = 5 to the 1st tank
            %   d.getNodeTankMinimumWaterLevel
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMinimumWaterLevel(tankIndex, 10)       % Sets the minimum water level = 10 to the tanks with index 128 and 129
            %   d.getNodeTankMinimumWaterLevel
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMinimumWaterLevel(tankIndex,[5, 15])   % Sets the minimum water level = 5 and 15 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankMinimumWaterLevel
            %
            % See also getNodeTankMinimumWaterLevel, setNodeTankInitialLevel, setNodeTankMaximumWaterLevel,
            %          setNodeTankMinimumWaterVolume, setNodeTankMixingFraction, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_MINLEVEL, value, varargin)
        end
        function setNodeTankMinimumWaterVolume(obj, value, varargin)
            % Sets the minimum water volume value for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankMinimumWaterVolume                           % Retrieves the minimum water volume of all tanks
            %   d.setNodeTankMinimumWaterVolume(1000)                     % Sets the minimum water volume = 1000 to every tank
            %   d.getNodeTankMinimumWaterVolume
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankMinimumWaterVolume([1500, 2000])             % Sets the minimum water volume = 1500 and 2000 to the 2 tanks
            %   d.getNodeTankMinimumWaterVolume
            %
            % Example 3:
            %   d.setNodeTankMinimumWaterVolume(1, 1000)                  % Sets the minimum water volume = 1000 to the 1st tank
            %   d.getNodeTankMinimumWaterVolume
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMinimumWaterVolume(tankIndex, 1500)          % Sets the minimum water volume = 1500 to the tanks with index 128 and 129
            %   d.getNodeTankMinimumWaterVolume
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMinimumWaterVolume(tankIndex,[1000, 2000])   % Sets the minimum water volume = 1000 and 2000 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankMinimumWaterVolume
            %
            % See also getNodeTankMinimumWaterVolume, setNodeTankInitialLevel, setNodeTankMinimumWaterLevel,
            %          setNodeTankMaximumWaterLevel, setNodeTankMixingFraction, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_MINVOLUME, value, varargin)
        end
        function setNodeTankMaximumWaterLevel(obj, value, varargin)
            % Sets the maximum water level value for tanks.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankMaximumWaterLevel                       % Retrieves the maximum water level of all tanks
            %   d.setNodeTankMaximumWaterLevel(35)                   % Sets the maximum water level = 35 to every tank
            %   d.getNodeTankMaximumWaterLevel
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankMaximumWaterLevel([30, 40])             % Sets the maximum water level = 30 and 40 to the 2 tanks
            %   d.getNodeTankMaximumWaterLevel
            %
            % Example 3:
            %   d.setNodeTankMaximumWaterLevel(1, 35)                % Sets the maximum water level = 35 to the 1st tank
            %   d.getNodeTankMaximumWaterLevel
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMaximumWaterLevel(tankIndex, 30)        % Sets the maximum water level = 30 to the tanks with index 128 and 129
            %   d.getNodeTankMaximumWaterLevel
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMaximumWaterLevel(tankIndex,[35, 45])   % Sets the maximum water level = 35 and 45 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankMaximumWaterLevel
            %
            % See also getNodeTankMaximumWaterLevel, setNodeTankInitialLevel, setNodeTankMinimumWaterLevel,
            %          setNodeTankMinimumWaterVolume, setNodeTankMixingFraction, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_MAXLEVEL, value, varargin)
        end
       function setNodeTankCanOverFlow(obj, value, varargin)
            % Sets the tank can-overflow (= 1) or not (= 0). (EPANET Version 2.2)
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankCanOverFlow                     % Retrieves the can-overflow of all tanks
            %   d.setNodeTankCanOverFlow(1)                  % Sets the can-overflow = 1 to every tank
            %   d.getNodeTankCanOverFlow
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankCanOverFlow([1, 0])             % Sets the can-overflow = 1 and 0 to the 2 tanks
            %   d.getNodeTankCanOverFlow
            %
            % Example 3:
            %   d.setNodeTankCanOverFlow(1, 0)               % Sets the can-overflow = 0 to the 1st tank
            %   d.getNodeTankCanOverFlow
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankCanOverFlow(tankIndex, 1)       % Sets the can-overflow = 1 to the tanks with index 128 and 129
            %   d.getNodeTankCanOverFlow
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankCanOverFlow(tankIndex,[0, 1])   % Sets the can-overflow = 0 and 1 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankCanOverFlow
            %
            % See also getNodeTankCanOverFlow, setNodeTankBulkReactionCoeff, setNodeTankMinimumWaterLevel,
            %          setNodeTankMinimumWaterVolume, setNodeTankDiameter, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_CANOVERFLOW, value, varargin)
        end
        function setNodeTankMixingFraction(obj, value, varargin)
            % Sets the tank mixing fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankMixingFraction                     % Retrieves the mixing fraction of all tanks
            %   d.setNodeTankMixingFraction(0)                  % Sets the mixing fraction = 0 to every tank
            %   d.getNodeTankMixingFraction
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankMixingFraction([1, 0])             % Sets the mixing fraction = 1 and 0 to the 2 tanks
            %   d.getNodeTankMixingFraction
            %
            % Example 3:
            %   d.setNodeTankMixingFraction(1, 0)               % Sets the mixing fraction = 0 to the 1st tank
            %   d.getNodeTankMixingFraction
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMixingFraction(tankIndex, 1)       % Sets the mixing fraction = 1 to the tanks with index 128 and 129
            %   d.getNodeTankMixingFraction
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankMixingFraction(tankIndex,[1, 0])   % Sets the mixing fraction = 1 and 0 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankMixingFraction
            %
            % See also getNodeTankMixingFraction, setNodeTankMixingModelType, setNodeTankMinimumWaterLevel,
            %          setNodeTankMinimumWaterVolume, setNodeTankDiameter, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_MIXFRACTION, value, varargin)
        end
        function setNodeTankBulkReactionCoeff(obj, value, varargin)
            % Sets the tank bulk reaction coefficient.
            %
            % The examples are based on d=epanet('BWSN_Network_1.inp');
            %
            % Example 1:
            %   d.getNodeTankBulkReactionCoeff                        % Retrieves the  bulk reaction coefficient of all tanks
            %   d.setNodeTankBulkReactionCoeff(-0.5)                  % Sets the bulk reaction coefficient = -0.5 to every tank
            %   d.getNodeTankBulkReactionCoeff
            %
            % Example 2:
            %   % The input array must have a length equal to the number of tanks
            %   d.setNodeTankBulkReactionCoeff([0, -0.5])             % Sets the bulk reaction coefficient = 0 and -0.5 to the 2 tanks
            %   d.getNodeTankBulkReactionCoeff
            %
            % Example 3:
            %   d.setNodeTankBulkReactionCoeff(1, -0.5)               % Sets the bulk reaction coefficient = -0.5 to the 1st tank
            %   d.getNodeTankBulkReactionCoeff
            %
            % Example 4:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankBulkReactionCoeff(tankIndex, 0)          % Sets the bulk reaction coefficient = 0 to the tanks with index 128 and 129
            %   d.getNodeTankBulkReactionCoeff
            %
            % Example 5:
            %   tankIndex = d.getNodeTankIndex;
            %   d.setNodeTankBulkReactionCoeff(tankIndex,[-0.5, 0])   % Sets the bulk reaction coefficient = -0.5 and 0 to the tanks with index 128 and 129 respectively
            %   d.getNodeTankBulkReactionCoeff
            %
            % See also getNodeTankBulkReactionCoeff, setNodeTankInitialLevel, setNodeTankMixingModelType,
            %          setNodeTankCanOverFlow, setNodeTankDiameter, setNodeTankData.
            set_Node_Link(obj, 'tank', 'ENsetnodevalue', obj.ToolkitConstants.EN_TANK_KBULK, value, varargin)
        end
        function setNodeSourceQuality(obj, value, varargin)
            % Sets the values of quality source strength.
            %
            % Example 1:
            %   nodeIndex = 1;
            %   d.getNodeSourceQuality(nodeIndex)                   % Retrieves the quality source strength of the 1st node
            %   sourceStrength = 10;
            %   d.setNodeSourceQuality(nodeIndex, sourceStrength)   % Sets the quality source strength = 10 to the 1st node
            %   d.getNodeSourceQuality(nodeIndex)
            %
            % Example 2:
            %   nodeIndex = 1:3;
            %   d.getNodeSourceQuality(nodeIndex)                   % Retrieves the quality source strength of the first 3 nodes
            %   sourceStrength = [10, 12, 8];
            %   d.setNodeSourceQuality(nodeIndex, sourceStrength)   % Sets the quality source strength = 10, 12 and 8 to the first 3 nodes
            %   d.getNodeSourceQuality(nodeIndex)
            %
            % See also getNodeSourceQuality, setNodeSourcePatternIndex, setNodeSourceType.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, obj.ToolkitConstants.EN_SOURCEQUAL, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeSourcePatternIndex(obj, value, varargin)
            % Sets the values of quality source pattern index.
            %
            % Example 1:
            %   nodeIndex = 1;
            %   d.getNodeSourcePatternIndex(nodeIndex)                       % Retrieves the quality source pattern index of the 1st node
            %   sourcePatternIndex = 1;
            %   d.setNodeSourcePatternIndex(nodeIndex, sourcePatternIndex)   % Sets the quality source pattern index = 1 to the 1st node
            %   d.getNodeSourcePatternIndex(nodeIndex)
            %
            % Example 2:
            %   nodeIndex = 1:3;
            %   d.getNodeSourcePatternIndex(nodeIndex)                       % Retrieves the quality source pattern index of the first 3 nodes
            %   sourcePatternIndex = [1, 1, 1];
            %   d.setNodeSourcePatternIndex(nodeIndex, sourcePatternIndex)   % Sets the quality source pattern index = 1 to the first 3 nodes
            %   d.getNodeSourcePatternIndex(nodeIndex)
            %
            % See also getNodeSourcePatternIndex, setNodeSourceQuality, setNodeSourceType.
            if nargin==3, indices = value; value=varargin{1}; else indices = getNodeIndices(obj, varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, obj.ToolkitConstants.EN_SOURCEPAT, value(j), obj.LibEPANET); j=j+1;
                error(obj.getError(obj.Errcode));
            end
        end
        function setNodeSourceType(obj, index, value)
            % Sets the values of quality source type.
            %
            % Types of external water quality sources that can be set:
            %   1) CONCEN      Sets the concentration of external inflow entering a node
            %   2) MASS        Injects a given mass/minute into a node
            %   3) SETPOINT    Sets the concentration leaving a node to a given value
            %   4) FLOWPACED   Adds a given value to the concentration leaving a node
            %
            % Example:
            %   nodeIndex = 1;
            %   d.getNodeSourceType{nodeIndex}               % Retrieves the quality source type of the 1st node
            %   sourceType = 'MASS';
            %   d.setNodeSourceType(nodeIndex, sourceType)   % Sets the quality source type = 'MASS' to the 1st node
            %   d.getNodeSourceType{nodeIndex}
            %
            % See also getNodeSourceType, setNodeSourceQuality, setNodeSourcePatternIndex.
            value=find(strcmpi(obj.TYPESOURCE, value)==1)-1;
            [obj.Errcode] = ENsetnodevalue(index, obj.ToolkitConstants.EN_SOURCETYPE, value, obj.LibEPANET);
        end
        function setDemandModel(obj, code, pmin, preq, pexp)
            % Sets the type of demand model to use and its parameters. (EPANET Version 2.2)
            %
            % Input Arguments:
            %   type - Type of demand model
            %       'DDA' = Demand driven analysis (in which case the remaining three parameter values are ignored)
            %       'PDA' = Pressure driven analysis
            %   pmin - Pressure below which there is no demand
            %   preq - Pressure required to deliver full demand
            %   pexp - Pressure exponent in demand function
            %
            % Example:
            %   d.getDemandModel                            % Retrieves the demand model
            %   type = 'PDA';
            %   pmin = 0;
            %   preq = 0.1;
            %   pexp = 0.5;
            %   d.setDemandModel(type, pmin, preq, pexp);   % Sets the demand model
            %   d.getDemandModel
            %
            % See also getDemandModel, setNodeBaseDemands, setNodeJunctionDemandName,
            %          addNodeJunctionDemand, deleteNodeJunctionDemand.
            model_type=find(strcmpi(obj.DEMANDMODEL, code)==1)-1;
            if isempty(model_type)
                error('Please give Demand model type: DDA or PDA');
            end
            [obj.Errcode] = ENsetdemandmodel(model_type, pmin, preq, pexp, obj.LibEPANET);
        end
        function setNodeJunctionDemandName(obj, nodeIndex, demandIndex, demandName)
            % Assigns a name to a node's demand category. (EPANET Version 2.2)
            %
            % Example:
            %   nodeIndex = 1;
            %   demandIndex = 1;
            %   d.getNodeJunctionDemandName{demandIndex}{nodeIndex}                % Retrieves the name of the 1st node, 1st demand category
            %   demandName = 'NEW NAME';
            %   d.setNodeJunctionDemandName(nodeIndex, demandIndex, demandName);   % Sets a new name of the 1st node, 1st demand category
            %   d.getNodeJunctionDemandName{demandIndex}{nodeIndex}
            %
            % See also getNodeJunctionDemandName, setNodeBaseDemands, setDemandModel,
            %          addNodeJunctionDemand, deleteNodeJunctionDemand.
            [obj.Errcode] = ENsetdemandname(nodeIndex, demandIndex, demandName, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setNodeJunctionData(obj, index, elev, dmnd, dmndpat)
            % Sets a group of properties for a junction node. (EPANET Version 2.2)
            %
            % Properties that can be set:
            %   1) index = a junction node's index (starting from 1).
            %   2) elev = the value of the junction's elevation.
            %   3) dmnd = the value of the junction's primary base demand.
            %   4) dmndpat = the ID name of the demand's time pattern ("" for no pattern)
            %
            % Example:
            %   junctionIndex = 1;
            %   elev = 35;
            %   dmnd = 100;
            %   dmndpat = 'NEW_PATTERN';
            %   d.addPattern(dmndpat)                                        % Adds a new pattern
            %   d.setNodeJunctionData(junctionIndex, elev, dmnd, dmndpat);   % Sets the elevation, primary base demand and time pattern of the 1st junction
            %   d.getNodeElevations(junctionIndex)                           % Retrieves the elevation of the 1st junction
            %   d.getNodeBaseDemands(junctionIndex)                          % Retrieves the primary base demand of the 1st junction
            %   d.getNodeDemandPatternNameID{1}(junctionIndex)               % Retrieves the demand pattern ID (primary base demand is the first category)
            %
            % See also setNodeTankData, getNodeElevations, getNodeBaseDemands,
            %          getNodeDemandPatternNameID, addPattern, setNodeJunctionDemandName.
            [obj.Errcode] = ENsetjuncdata(obj.LibEPANET, index, elev, dmnd, dmndpat);
            error(obj.getError(obj.Errcode));
        end
        function setTitle(obj, varargin)
            % Sets the title lines of the project. (EPANET Version 2.2)
            %
            % Example:
            %   line_1 = 'This is a title';
            %   line_2 = 'This is a test line 2';
            %   line_3 = 'This is a test line 3';
            %   d.setTitle(line_1, line_2, line_3);
            %   [Line1, Line2, Line3] = d.getTitle
            %
            % See also getTitle, setLinkComment, setNodeComment.
            line2 = '';
            line3 = '';
            if nargin > 1
                line1 = varargin{1};
            end
            if nargin > 2
                line2 = varargin{2};
            end
            if nargin > 3
                line3 = varargin{3};
            end
            [obj.Errcode] = ENsettitle(line1, line2, line3, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsMaxTrials(obj, value)
            % Sets the maximum hydraulic trials allowed for hydraulic convergence.
            %
            % Example:
            %   d.setOptionsMaxTrials(40)
            %   d.getOptionsMaxTrials
            %
            % See also getOptionsMaxTrials, setOptionsExtraTrials, setOptionsAccuracyValue.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_TRIALS, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsAccuracyValue(obj, value)
            % Sets the total normalized flow change for hydraulic convergence.
            %
            % Example:
            %   d.setOptionsAccuracyValue(0.001)
            %   d.getOptionsAccuracyValue
            %
            % See also getOptionsAccuracyValue, setOptionsExtraTrials, setOptionsMaxTrials.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_ACCURACY, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsQualityTolerance(obj, value)
            % Sets the water quality analysis tolerance.
            %
            % Example:
            %   d.setOptionsQualityTolerance(0.01)
            %   d.getOptionsQualityTolerance
            %
            % See also getOptionsQualityTolerance, setOptionsSpecificDiffusivity, setOptionsLimitingConcentration.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_TOLERANCE, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsEmitterExponent(obj, value)
            % Sets the power exponent for the emmitters.
            %
            % Example:
            %   d.setOptionsEmitterExponent(0.5)
            %   d.getOptionsEmitterExponent
            %
            % See also getOptionsEmitterExponent, setOptionsPatternDemandMultiplier, setOptionsAccuracyValue.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_EMITEXPON, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsPatternDemandMultiplier(obj, value)
            % Sets the global pattern demand multiplier.
            %
            % Example:
            %   d.setOptionsPatternDemandMultiplier(1)
            %   d.getOptionsPatternDemandMultiplier
            %
            % See also getOptionsPatternDemandMultiplier, setOptionsEmitterExponent, setOptionsAccuracyValue.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_DEMANDMULT, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsHeadError(obj, value)
            % Sets the maximum head loss error for hydraulic convergence. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsHeadError(0)
            %   d.getOptionsHeadError
            %
            % See also getOptionsHeadError, setOptionsEmitterExponent, setOptionsAccuracyValue.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_HEADERROR, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsFlowChange(obj, value)
            % Sets the maximum flow change for hydraulic convergence. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsFlowChange(0)
            %   d.getOptionsFlowChange
            %
            % See also getOptionsFlowChange, setOptionsHeadError, setOptionsHeadLossFormula.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_FLOWCHANGE, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsHeadLossFormula(obj, value)
            % Sets the headloss formula. (EPANET Version 2.2)
            % 'HW' = 0, 'DW' = 1, 'CM' = 2
            %
            % Example:
            %   d.setOptionsHeadLossFormula('HW')   % Sets the 'HW' headloss formula
            %   d.getOptionsHeadLossFormula
            %
            % See also getOptionsHeadLossFormula, setOptionsHeadError, setOptionsFlowChange.
            if value=='HW'
                codevalue=0;
            elseif value=='DW'
                codevalue=1;
            elseif value=='CM'
                codevalue=2;
            end
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_HEADLOSSFORM, codevalue, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setOptionsGlobalEffic(obj, value)
            % Sets the global efficiency for pumps(percent). (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsGlobalEffic(75)
            %   d.getOptionsGlobalEffic
            %
            % See also getOptionsGlobalEffic, setOptionsGlobalPrice, setOptionsGlobalPattern.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_GLOBALEFFIC, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsGlobalPrice(obj, value)
            % Sets the global average energy price per kW-Hour. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsGlobalPrice(0)
            %   d.getOptionsGlobalPrice
            %
            % See also getOptionsGlobalPrice, setOptionsGlobalEffic, setOptionsGlobalPattern.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_GLOBALPRICE, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsGlobalPattern(obj, value)
            % Sets the global energy price pattern. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsGlobalPattern(1)
            %   d.getOptionsGlobalPattern
            %
            % See also getOptionsGlobalPattern, setOptionsGlobalEffic, setOptionsGlobalPrice.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_GLOBALPATTERN, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsDemandCharge(obj, value)
            % Sets the energy charge per maximum KW usage. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsDemandCharge(0)
            %   d.getOptionsDemandCharge
            %
            % See also getOptionsDemandCharge, setOptionsGlobalPrice, setOptionsGlobalPattern.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_DEMANDCHARGE, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsSpecificGravity(obj, value)
            % Sets the specific gravity. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsSpecificGravity(1)
            %   d.getOptionsSpecificGravity
            %
            % See also getOptionsSpecificGravity, setOptionsSpecificViscosity, setOptionsHeadLossFormula.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_SP_GRAVITY, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsSpecificViscosity(obj, value)
            % Sets the specific viscosity. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsSpecificViscosity(1)
            %   d.getOptionsSpecificViscosity
            %
            % See also getOptionsSpecificViscosity, setOptionsSpecificGravity, setOptionsHeadLossFormula.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_SP_VISCOS, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsExtraTrials(obj, value)
            % Sets the extra trials allowed if hydraulics don't converge. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsExtraTrials(10);
            %   d.getOptionsExtraTrials
            %
            %   % Set UNBALANCED to STOP
            %   d.setOptionsExtraTrials(-1);
            %
            % See also getOptionsExtraTrials, setOptionsMaxTrials, setOptionsMaximumCheck.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_UNBALANCED, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsCheckFrequency(obj, value)
            % Sets the frequency of hydraulic status checks. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsCheckFrequency(2)
            %   d.getOptionsCheckFrequency
            %
            % See also getOptionsCheckFrequency, setOptionsMaxTrials, setOptionsMaximumCheck.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_CHECKFREQ, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsMaximumCheck(obj, value)
            % Sets the maximum trials for status checking. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsMaximumCheck(10)
            %   d.getOptionsMaximumCheck
            %
            % See also getOptionsMaximumCheck, setOptionsMaxTrials, setOptionsCheckFrequency.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_MAXCHECK, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsDampLimit(obj, value)
            % Sets the accuracy level where solution damping begins. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsDampLimit(0)
            %   d.getOptionsDampLimit
            %
            % See also getOptionsDampLimit, setOptionsMaxTrials, setOptionsCheckFrequency.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_DAMPLIMIT, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsSpecificDiffusivity(obj, value)
            % Sets the specific diffusivity (relative to chlorine at 20 deg C). (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsSpecificDiffusivity(1)
            %   d.getOptionsSpecificDiffusivity
            %
            % See also getOptionsSpecificDiffusivity, setOptionsSpecificViscosity, setOptionsSpecificGravity.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_SP_DIFFUS, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsPipeBulkReactionOrder(obj, value)
            % Sets the bulk water reaction order for pipes. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsPipeBulkReactionOrder(1)
            %   d.getOptionsPipeBulkReactionOrder
            %
            % See also getOptionsPipeBulkReactionOrder, setOptionsPipeWallReactionOrder, setOptionsTankBulkReactionOrder.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_BULKORDER, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsPipeWallReactionOrder(obj, value)
            % Sets the wall reaction order for pipes (either 0 or 1). (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsPipeWallReactionOrder(1)
            %   d.getOptionsPipeWallReactionOrder
            %
            % See also getOptionsPipeWallReactionOrder, setOptionsPipeBulkReactionOrder, setOptionsTankBulkReactionOrder.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_WALLORDER, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsTankBulkReactionOrder(obj, value)
            % Sets the bulk water reaction order for tanks. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsTankBulkReactionOrder(1)
            %   d.getOptionsTankBulkReactionOrder
            %
            % See also getOptionsTankBulkReactionOrder, setOptionsPipeBulkReactionOrder, setOptionsPipeWallReactionOrder.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_TANKORDER, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setOptionsLimitingConcentration(obj, value)
            % Sets the limiting concentration for growth reactions. (EPANET Version 2.2)
            %
            % Example:
            %   d.setOptionsLimitingConcentration(0)
            %   d.getOptionsLimitingConcentration
            %
            % See also getOptionsLimitingConcentration, setOptionsPipeBulkReactionOrder, setOptionsPipeWallReactionOrder.
            [obj.Errcode] = ENsetoption(obj.ToolkitConstants.EN_CONCENLIMIT, value, obj.LibEPANET);
            if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
        end
        function setTimeSimulationDuration(obj, value)
            % Sets the simulation duration (in seconds).
            %
            % Example:
            %   simulationDuration = 172800;   % 172800 seconds = 2days
            %   d.setTimeSimulationDuration(simulationDuration)
            %   d.getTimeSimulationDuration
            %
            % See also getTimeSimulationDuration, getTimeStartTime, getTimeHaltFlag.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_DURATION, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setTimeHydraulicStep(obj, value)
            % Sets the hydraulic time step.
            %
            % Example:
            %   Hstep = 1800;
            %   d.setTimeHydraulicStep(Hstep)
            %   d.getTimeHydraulicStep
            %
            % See also getTimeHydraulicStep, setTimeQualityStep, setTimePatternStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_HYDSTEP, value, obj.LibEPANET);
        end
        function setTimeQualityStep(obj, value)
            % Sets the quality time step.
            %
            % Example:
            %   Qstep = 1800;
            %   d.setTimeQualityStep(Qstep)
            %   d.getTimeQualityStep
            %
            % See also getTimeQualityStep, setTimeHydraulicStep, setTimePatternStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_QUALSTEP, value, obj.LibEPANET);
        end
        function setTimePatternStep(obj, value)
            % Sets the time pattern step.
            %
            % Example:
            %   patternStep = 3600;
            %   d.setTimePatternStep(patternStep)
            %   d.getTimePatternStep
            %
            % See also getTimePatternStep, setTimePatternStart, setTimeHydraulicStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_PATTERNSTEP, value, obj.LibEPANET);
        end
        function setTimePatternStart(obj, value)
            % Sets the time when time patterns begin.
            %
            % Example:
            %   patternStart = 0;
            %   d.setTimePatternStart(patternStart)
            %   d.getTimePatternStart
            %
            % See also getTimePatternStart, setTimePatternStep, setTimeHydraulicStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_PATTERNSTART, value, obj.LibEPANET);
        end
        function setTimeReportingStep(obj, value)
            % Sets the reporting time step.
            %
            % Example:
            %   reportingStep = 3600;
            %   d.setTimeReportingStep(reportingStep)
            %   d.getTimeReportingStep
            %
            % See also getTimeReportingStep, setTimeReportingStart, setTimeRuleControlStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_REPORTSTEP, value, obj.LibEPANET);
        end
        function setTimeReportingStart(obj, value)
            % Sets the time when reporting starts.
            %
            % Example:
            %   reportingStart = 0;
            %   d.setTimeReportingStart(reportingStart)
            %   d.getTimeReportingStart
            %
            % See also getTimeReportingStart, setTimeReportingStep, setTimePatternStart.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_REPORTSTART, value, obj.LibEPANET);
        end
        function setTimeRuleControlStep(obj, value)
            % Sets the rule-based control evaluation time step.
            %
            % Example:
            %   ruleControlStep = 360;
            %   d.setTimeRuleControlStep(ruleControlStep)
            %   d.getTimeRuleControlStep
            %
            % See also getTimeRuleControlStep, setTimeReportingStep, setTimePatternStep.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_RULESTEP, value, obj.LibEPANET);
        end
        function setTimeStatisticsType(obj, value)
            % Sets the statistic type.
            %
            % Types that can be set:
            %   1) 'NONE'
            %   2) 'AVERAGE'
            %   3) 'MINIMUM'
            %   4) 'MAXIMUM'
            %   5) 'RANGE'
            %
            % Example:
            %   d.getTimeStatisticsType
            %   statisticsType = 'AVERAGE';
            %   d.setTimeStatisticsType(statisticsType)
            %   d.getTimeStatisticsType
            %
            % See also getTimeStatisticsType, setTimeReportingStart, setTimeReportingStep.
            tmpindex=find(strcmpi(obj.TYPESTATS, value)==1)-1;
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_STATISTIC, tmpindex, obj.LibEPANET);
        end
        function setTimeStartTime(obj, value)
            % Sets the simulation starting time of day.
            %
            % Example:
            %   startTime = 0;
            %   d.setTimeStartTime(startTime)
            %   d.getTimeStartTime
            %
            % See also getTimeStartTime, setTimeReportingStart, setTimePatternStart.
            [obj.Errcode] = ENsettimeparam(obj.ToolkitConstants.EN_STARTTIME, value, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function value = setPatternComment(obj, value, varargin)
            % Sets the comment string assigned to the pattern object.
            %
            % Example 1:
            %   patternIndex = 1;
            %   patternComment = 'This is a PATTERN';
            %   d.setPatternComment(patternIndex, patternComment);   % Sets the comment of the 1st pattern
            %   d.getPatternComment(patternIndex)                    % Retrieves the comment of the 1st pattern
            %
            % Example 2:
            %   patternIndex = 1:2;
            %   patternComment = {'1st PATTERN', '2nd PATTERN'};
            %   d.setPatternComment(patternIndex, patternComment);   % Sets the comments of the first 2 patterns (if exist)
            %   d.getPatternComment(patternIndex)
            %
            % Example 3:
            %   patternComment = {'1st PATTERN', '2nd PATTERN'};
            %   d.setPatternComment(patternComment);                 % Sets the comments of all the patterns (the length of the cell must equal the number of patterns)
            %   d.getPatternComment
            %
            % See also getPatternComment, setPatternNameID, setPattern.
            if nargin==3, indices = value; value=varargin{1}; else indices = getPatternIndices(obj, varargin); end
            j=1;
            if length(indices) == 1
                [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_TIMEPAT, indices, value, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            else
                for i=indices
                    [obj.Errcode] = ENsetcomment(obj.ToolkitConstants.EN_TIMEPAT, i, value{j}, obj.LibEPANET); j=j+1;
                    error(obj.getError(obj.Errcode));
                end
            end
        end
        function setPattern(obj, index, patternVector)
            % Sets all of the multiplier factors for a specific time pattern.
            %
            % Example:
            %   patternID = 'new_pattern';
            %   patternIndex = d.addPattern(patternID)    % Adds a new time pattern
            %   patternMult = [1.56, 1.36, 1.17, 1.13, 1.08, ...
            %       1.04, 1.2, 0.64, 1.08, 0.53, 0.29, 0.9, 1.11, ...
            %       1.06, 1.00, 1.65, 0.55, 0.74, 0.64, 0.46, ...
            %       0.58, 0.64, 0.71, 0.66];
            %   d.setPattern(patternIndex, patternMult)   % Sets the multiplier factors for the new time pattern
            %   d.getPattern                              % Retrieves the multiplier factor for all patterns and all times
            %
            % See also getPattern, setPatternValue, setPatternMatrix,
            %          setPatternNameID, addPattern, deletePattern.
            nfactors=length(patternVector);
            [obj.Errcode] = ENsetpattern(index, patternVector, nfactors, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setPatternMatrix(obj, patternMatrix)
            % Sets all of the multiplier factors for all time patterns.
            %
            % Example:
            %   patternID_1 = 'new_pattern_1';
            %   patternIndex_1 = d.addPattern(patternID_1)    % Adds a new time pattern
            %   patternID_2 = 'new_pattern_2';
            %   patternIndex_2 = d.addPattern(patternID_2)    % Adds a new time pattern
            %   patternMult = d.getPattern
            %   patternMult(patternIndex_1, 2) = 5;           % The 2nd multiplier = 5 of the 1st time pattern
            %   patternMult(patternIndex_2, 3) = 7;           % The 3rd multiplier = 7 of the 2nd time pattern
            %   d.setPatternMatrix(patternMult)               % Sets all of the multiplier factors for all the time patterns given a matrix
            %   d.getPattern                                  % Retrieves the multiplier factor for all patterns and all times
            %
            % See also getPattern, setPattern, setPatternValue,
            %          setPatternNameID, addPattern, deletePattern.
            nfactors=size(patternMatrix, 2);
            for i=1:size(patternMatrix, 1)
                [obj.Errcode] = ENsetpattern(i, patternMatrix(i, :), nfactors, obj.LibEPANET);
                error(obj.getError(obj.Errcode));
            end
        end
        function setPatternValue(obj, index, patternTimeStep, patternFactor)
            % Sets the multiplier factor for a specific period within a time pattern.
            %
            % Example:
            %   patternID = 'new_pattern';
            %   patternIndex = d.addPattern(patternID)                            % Adds a new time pattern
            %   patternTimeStep = 2;
            %   patternFactor = 5;
            %   d.setPatternValue(patternIndex, patternTimeStep, patternFactor)   % Sets the multiplier factor = 5 to the 2nd time period of the new time pattern
            %   d.getPattern                                                      % Retrieves the multiplier factor for all patterns and all times
            %
            % See also getPattern, setPattern, setPatternMatrix,
            %          setPatternNameID, addPattern, deletePattern.
            [obj.Errcode] = ENsetpatternvalue(index, patternTimeStep, patternFactor, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setQualityType(obj, varargin)
            % Sets the type of water quality analysis called for.
            %
            % Example 1:
            %   d.setQualityType('none')                        % Sets no quality analysis.
            %   d.getQualityInfo                                % Retrieves quality analysis information
            %
            % Example 2:
            %   d.setQualityType('age')                         % Sets water age analysis
            %   d.getQualityInfo
            %
            % Example 3:
            %   d.setQualityType('chem', 'Chlorine')            % Sets chemical analysis given the name of the chemical being analyzed
            %   d.getQualityInfo
            %   d.setQualityType('chem', 'Chlorine', 'mg/Kg')   % Sets chemical analysis given the name of the chemical being analyzed and units that the chemical is measured in
            %   d.getQualityInfo
            %
            % Example 4:
            %   d.setQualityType('Chlorine')                    % Sets chemical analysis given the name of the chemical being analyzed
            %   d.getQualityInfo
            %   d.setQualityType('Chlorine', 'mg/Kg')           % Sets chemical analysis given the name of the chemical being analyzed and units that the chemical is measured in
            %   d.getQualityInfo
            %
            % Example 5:
            %   nodeID = d.getNodeNameID{1};
            %   d.setQualityType('trace', nodeID)               % Sets source tracing analysis given the ID label of node traced in a source tracing analysis
            %   d.getQualityInfo
            %
            % See also getQualityInfo, getQualityType, getQualityCode, getQualityTraceNodeIndex.
            qualcode=obj.ToolkitConstants.EN_NONE;chemname='';chemunits='';tracenode='';
            if find(strcmpi(varargin, 'none')==1)
            elseif find(strcmpi(varargin, 'age')==1)
                qualcode=obj.ToolkitConstants.EN_AGE;
            elseif find(strcmpi(varargin, 'chem')==1)
                qualcode=obj.ToolkitConstants.EN_CHEM;
                chemname=varargin{2};
                if nargin<=3
                    chemunits='mg/L';
                else
                    chemunits=varargin{3};
                end
            elseif find(strcmpi(varargin, 'trace')==1)
                qualcode=obj.ToolkitConstants.EN_TRACE;
                tracenode=varargin{2};
            else
                qualcode=obj.ToolkitConstants.EN_CHEM;
                chemname=varargin{1};
                if nargin<3
                    chemunits='mg/L';
                else
                    chemunits=varargin{2};
                end
            end
            [obj.Errcode] = ENsetqualtype(qualcode, chemname, chemunits, tracenode, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setReportFormatReset(obj)
            % Resets a project's report options to their default values.
            %
            % Example:
            %   d.setReportFormatReset
            %
            % See also setReport, setReportStatus.
            [obj.Errcode]=ENresetreport(obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setReportStatus(obj, value)
            % Sets the level of hydraulic status reporting.
            %
            % Possible status that can be set:
            %   1) 'yes'
            %   2) 'no'
            %   3) 'full'
            %
            % Example:
            %   d.setReportStatus('full')
            %
            % See also setReport, setReportFormatReset.
            statuslevel=find(strcmpi(obj.TYPEREPORT, value)==1)-1;
            [obj.Errcode] = ENsetstatusreport(statuslevel, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
        function setReport(obj, value)
            % Issues a report formatting command. Formatting commands are the same as used in the [REPORT] section of the EPANET Input file.
            % More: https://github.com/OpenWaterAnalytics/EPANET/wiki/%5BREPORT%5D
            %
            % Example 1:
            %   d.setReport('FILE TestReport.txt')
            %
            % Example 2:
            %   d.setReport('STATUS YES')
            %
            % See also setReportFormatReset, setReport.
            [obj.Errcode] = ENsetreport(value, obj.LibEPANET);
        end
        function [Errcode]=setFlowUnitsGPM(obj, varargin)
            % Sets flow units to GPM(Gallons Per Minute).
            %
            % Example:
            %   d.setFlowUnitsGPM;   % d.setFlowUnitsGPM('NET1_GPM.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsLPS, setFlowUnitsMGD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_GPM, 1, varargin); % gallons per minute
        end
        function [Errcode]=setFlowUnitsLPS(obj, varargin)
            % Sets flow units to LPS(Liters Per Second).
            %
            % Example:
            %   d.setFlowUnitsLPS;   % d.setFlowUnitsLPS('NET1_LPS.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsGPM, setFlowUnitsMGD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_LPS, 1, varargin); % liters per second
        end
        function [Errcode]=setFlowUnitsMGD(obj, varargin)
            % Sets flow units to MGD(Million Gallons per Day).
            %
            % Example:
            %   d.setFlowUnitsMGD;   % d.setFlowUnitsMGD('NET1_MGD.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsGPM, setFlowUnitsLPS.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_MGD, 1, varargin); % million gallons per day
        end
        function [Errcode]=setFlowUnitsIMGD(obj, varargin)
            % Sets flow units to IMGD(Imperial Million Gallons per Day).
            %
            % Example:
            %   d.setFlowUnitsIMGD;   % d.setFlowUnitsIMGD('NET1_IMGD.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsMGD, setFlowUnitsCFS.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_IMGD, 1, varargin); % imperial mgd
        end
        function [Errcode]=setFlowUnitsCFS(obj, varargin)
            % Sets flow units to CFS(Cubic Feet per Second).
            %
            % Example:
            %   d.setFlowUnitsCFS;   % d.setFlowUnitsCFS('NET1_CFS.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsAFD, setFlowUnitsIMGD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_CFS, 1, varargin); % cubic feet per second
        end
        function [Errcode]=setFlowUnitsAFD(obj, varargin)
            % Sets flow units to AFD(Acre-Feet per Day).
            %
            % Example:
            %   d.setFlowUnitsAFD;   % d.setFlowUnitsAFD('NET1_AFD.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsCFS, setFlowUnitsIMGD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_AFD, 1, varargin); % acre-feet per day
        end
        function [Errcode]=setFlowUnitsLPM(obj, varargin)
            % Sets flow units to LPM(Liters Per Minute).
            %
            % Example:
            %   d.setFlowUnitsLPM;   % d.setFlowUnitsLPM('NET1_LPM.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsAFD, setFlowUnitsMLD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_LPM, 1, varargin); % liters per minute
        end
        function [Errcode]=setFlowUnitsMLD(obj, varargin)
            % Sets flow units to MLD(Million Liters per Day).
            %
            % Example:
            %   d.setFlowUnitsMLD;   % d.setFlowUnitsMLD('NET1_MLD.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsLPM, setFlowUnitsCMH.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_MLD, 1, varargin); % million liters per day
        end
        function [Errcode]=setFlowUnitsCMH(obj, varargin)
            % Sets flow units to CMH(Cubic Meters per Hour).
            %
            % Example:
            %   d.setFlowUnitsCMH;   % d.setFlowUnitsCMH('NET1_CMH.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsMLD, setFlowUnitsCMD.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_CMH, 1, varargin); % cubic meters per hour
        end
        function [Errcode]=setFlowUnitsCMD(obj, varargin)
            % Sets flow units to CMD(Cubic Meters per Day).
            %
            % Example:
            %   d.setFlowUnitsCMD;   % d.setFlowUnitsCMD('NET1_CMD.inp');
            %   d.getFlowUnits
            %
            % See also setFlowUnitsMLD, setFlowUnitsCMH.
            Errcode = obj.setFlowUnits(obj.ToolkitConstants.EN_CMD, 1, varargin); % cubic meters per day
        end
        function closeNetwork(obj)
            % Closes down the Toolkit system.
            %
            % Example:
            %   d.closeNetwork
            %
            % See also loadEPANETFile, closeHydraulicAnalysis, closeQualityAnalysis.
            [obj.Errcode] = obj.ENclose(obj.LibEPANET);
        end
        function closeHydraulicAnalysis(obj)
            % Closes the hydraulic analysis system, freeing all allocated memory.
            %
            % Example:
            %   d.closeHydraulicAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also openHydraulicAnalysis, saveHydraulicFile, closeQualityAnalysis.
            [obj.Errcode] = obj.ENcloseH(obj.LibEPANET);
        end
        function closeQualityAnalysis(obj)
            % Closes the water quality analysis system, freeing all allocated memory.
            %
            % Example:
            %   d.closeQualityAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also openQualityAnalysis, initializeQualityAnalysis, closeHydraulicAnalysis.
            [obj.Errcode] = obj.ENcloseQ(obj.LibEPANET);
        end
        function saveHydraulicFile(obj, hydname)
            % Saves the current contents of the binary hydraulics file to a file.
            %
            % Example:
            %   filename = 'test.hyd'
            %   d.saveHydraulicFile(filename)
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also useHydraulicFile, initializeHydraulicAnalysis.
            [obj.Errcode]=ENsavehydfile(hydname, obj.LibEPANET);
        end
        function useHydraulicFile(obj, hydname)
            % Uses the contents of the specified file as the current binary hydraulics file.
            %
            % Example:
            %   filename = 'test.hyd'
            %   d.useHydraulicFile(filename)
            %
            % See also saveHydraulicFile, initializeHydraulicAnalysis.
            [obj.Errcode]=ENusehydfile(hydname, obj.LibEPANET);
        end
        function initializeEPANET(obj, unitsType, headLossType)
            % Initializes an EPANET project that isn't opened with an input file
            %
            % Example:
            %   d.initializeEPANET(d.ToolkitConstants.EN_GPM, d.ToolkitConstants.EN_HW)
            %
            % See also initializeHydraulicAnalysis.
            [obj.Errcode]=ENinit(obj.LibEPANET, unitsType, headLossType);
        end
        function initializeHydraulicAnalysis(obj, varargin)
            % Initializes storage tank levels, link status and settings, and the simulation clock time prior to running a hydraulic analysis.
            %
            % Codes:
            %   1) NOSAVE        = 0,  % Don't save hydraulics; don't re-initialize flows
            %   2) SAVE          = 1,  % Save hydraulics to file, don't re-initialize flows
            %   3) INITFLOW      = 10, % Don't save hydraulics; re-initialize flows
            %   4) SAVE_AND_INIT = 11  % Save hydraulics; re-initialize flows
            %
            % Example 1:
            %   d.initializeHydraulicAnalysis   % Uses the default code i.e. SAVE = 1
            %
            % Example 2:
            %   code = 0;                       % i.e. Don't save
            %   d.initializeHydraulicAnalysis(code)
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also saveHydraulicFile, initializeQualityAnalysis.
            code=obj.ToolkitConstants.EN_SAVE;
            if ~isempty(varargin)
                code=varargin{1};
            end
            [obj.Errcode] = ENinitH(code, obj.LibEPANET);
        end
        function initializeQualityAnalysis(obj, varargin)
            % Initializes water quality and the simulation clock time prior to running a water quality analysis.
            %
            % Codes:
            %   1) NOSAVE        = 0,  % Don't save the results to the project's binary output file.
            %   2) SAVE          = 1,  % Save the results to the project's binary output file.
            %
            % Example 1:
            %   d.initializeQualityAnalysis   % Uses the default code i.e. SAVE = 1
            %
            % Example 2:
            %   code = 0;                     % i.e. Don't save
            %   d.initializeQualityAnalysis(code)
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also openQualityAnalysis, initializeHydraulicAnalysis.
            code=obj.ToolkitConstants.EN_SAVE;
            if ~isempty(varargin)
            % obj.ToolkitConstants.EN_SAVE_AND_INIT; obj.ToolkitConstants.EN_NOSAVE;
            % obj.ToolkitConstants.EN_INITFLOW;
                code=varargin{1};
            end
            [obj.Errcode] = ENinitQ(code, obj.LibEPANET);
        end
        function tstep = nextHydraulicAnalysisStep(obj)
            % Determines the length of time until the next hydraulic event occurs in an extended period simulation.
            %
            % Example:
            %   d.nextHydraulicAnalysisStep
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also nextQualityAnalysisStep, runHydraulicAnalysis.
            [obj.Errcode, tstep] = ENnextH(obj.LibEPANET);
        end
        function tstep = nextQualityAnalysisStep(obj)
            % Advances the water quality simulation to the start of the next hydraulic time period.
            %
            % Example:
            %   d.nextQualityAnalysisStep
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also nextHydraulicAnalysisStep, runQualityAnalysis.
            [obj.Errcode, tstep] = ENnextQ(obj.LibEPANET);
        end
        function openHydraulicAnalysis(obj)
            % Opens the hydraulics analysis system.
            %
            % Example:
            %   d.openHydraulicAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also openQualityAnalysis, initializeHydraulicAnalysis.
            [obj.Errcode] = ENopenH(obj.LibEPANET);
        end
        function openQualityAnalysis(obj)
            % Opens the water quality analysis system.
            %
            % Example:
            %   d.openQualityAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also openHydraulicAnalysis, initializeQualityAnalysis.
            [obj.Errcode] = ENopenQ(obj.LibEPANET);
        end
        function tstep = runHydraulicAnalysis(obj)
            % Runs a single period hydraulic analysis, retrieving the current simulation clock time t.
            %
            % Example:
            %   tstep = d.runHydraulicAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also runQualityAnalysis, initializeHydraulicAnalysis.
            [obj.Errcode, tstep] = ENrunH(obj.LibEPANET);
        end
        function tstep = runQualityAnalysis(obj)
            % Makes available the hydraulic and water quality results that occur at the start of
            % the next time period of a water quality analysis, where the start of the period is returned in t.
            %
            % Example:
            %   tstep = d.runQualityAnalysis
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also runHydraulicAnalysis, initializeQualityAnalysis.
            [obj.Errcode, tstep] = ENrunQ(obj.LibEPANET);
        end
        function saveHydraulicsOutputReportingFile(obj)
            % Transfers results of a hydraulic simulation from the binary Hydraulics file
            % to the binary Output file, where results are only reported at uniform reporting intervals.
            %
            % Example:
            %   d.saveHydraulicsOutputReportingFile
            %
            % See also saveHydraulicFile, closeHydraulicAnalysis.
            [obj.Errcode] = ENsaveH(obj.LibEPANET);
        end
        function tleft=stepQualityAnalysisTimeLeft(obj)
            % Advances the water quality simulation one water quality time step.
            % The time remaining in the overall simulation is returned in tleft.
            %
            % Example:
            %   tleft = d.stepQualityAnalysisTimeLeft
            %
            % For more, you can type `help getNodePressure` and check examples 3 & 4.
            %
            % See also runQualityAnalysis, closeQualityAnalysis.
            [obj.Errcode, tleft] = ENstepQ(obj.LibEPANET);
        end
        function [Errcode] = saveInputFile(obj, inpname, varargin)
            % Writes all current network input data to a file using the format of an EPANET input file.
            % Returns an error code.
            %
            % Example:
            %   filename = ('test.inp');
            %   d.saveInputFile(filename)
            %
            % See also unload, saveHydraulicFile.
            % while exist('@#', 'file') ~= 2
            if nargin == 1
                inpname = obj.TempInpFile;
            end
            [Errcode] = ENsaveinpfile('@#', obj.LibEPANET);
            copyfile('@#', inpname);% temporary
            delete('@#');
            %else
            %    [addSectionCoordinates, addSectionRules] = obj.getBinCoordRuleSections(obj.BinTempfile);
            %    [Errcode] = ENsaveinpfile(inpname,obj.LibEPANET);
            %    [~, info_file] = obj.readInpFile;
            %    vertSectionIndex = find(~cellfun(@isempty, regexp(info_file,'VERTICES','match')), 1);
            %    len_sec = length(addSectionCoordinates);
            %    if isempty(vertSectionIndex)
            %        fid = fopen(obj.BinTempfile); % Opens the file for read access
            %        texta = char;
            %        i = 1; ok = 0;
            %        while (i < len_sec)
            %            aline = fgetl(fid);
            %            if ~ok
            %                texta = [texta, aline, char(10)];
            %            end
            %            if strcmp(aline, '[COORDINATES]') || ok
            %               ok = 1;
            %               texta = [texta, addSectionCoordinates{i+1}, char(10)];
            %               i = i +1;
            %            end
            %        end
            %        fid = fopen(obj.BinTempfile, 'w');   % Opens file for writing and discard existing contents
            %        fprintf(fid, texta);   % Writes the new text in the .inp file
            %        fclose('all');
            %    end
            %end
        end
        function writeLineInReportFile(obj, line)
            % Writes a line of text to the EPANET report file.
            %
            % Example:
            %   line = 'Status YES';
            %   d.writeLineInReportFile(line)
            %
            % See also writeReport, copyReport.
            [obj.Errcode] = obj.ENwriteline(line, obj.LibEPANET);
        end
        function writeReport(obj)
            % Writes a formatted text report on simulation results to the Report file.
            %
            % Example:
            %   d = epanet('Net1.inp');
            %   d.solveCompleteHydraulics;
            %   d.solveCompleteQuality;
            %   d.setReportFormatReset
            %   d.setReport('FILE TestReport3.txt');
            %   d.setReport('NODES ALL')
            %   d.setReport('LINKS ALL')
            %   d.writeReport
            %   open('TestReport3.txt')
            %
            % See also copyReport, writeLineInReportFile.
            [obj.Errcode]=ENreport(obj.LibEPANET);
        end
        function copyReport(obj, fileName)
            % Copies the current contents of a project's report file to another file. (EPANET Version 2.2)
            %
            % Example:
            %   fileName = 'Report_copy';
            %   d.copyReport(fileName)
            %
            % See also writeReport, writeLineInReportFile, clearReport.
            [obj.Errcode] = ENcopyreport (fileName, obj.LibEPANET);
        end
        function resultindex = getNodeResultIndex(obj, node_index)
            % Retrieves the order in which a node's results
            % were saved to an output file. (EPANET Version 2.2)
            %
            % Example:
            %   node_index = 3;
            %   result_index = d.getNodeResultIndex(node_index)
            %
            % See also getComputedHydraulicTimeSeries, deleteNode, getLinkResultIndex
            [obj.Errcode, resultindex] = ENgetresultindex(obj.LibEPANET, obj.ToolkitConstants.EN_NODE, node_index);
        end
        function resultindex = getLinkResultIndex(obj, link_index)
            % Retrieves the order in which a link's results
            % were saved to an output file. (EPANET Version 2.2)
            %
            % Example:
            %   link_index = 3;
            %   result_index = d.getLinkResultIndex(link_index)
            %
            % See also getComputedHydraulicTimeSeries, deleteNode, getNodeResultIndex
            [obj.Errcode, resultindex] = ENgetresultindex(obj.LibEPANET, obj.ToolkitConstants.EN_LINK, link_index);
        end
        function clearReport(obj)
            % Clears the contents of a project's report file. (EPANET Version 2.2)
            %
            % Example:
            %   d.clearReport
            %
            % See also writeReport, writeLineInReportFile, copyReport.
            [obj.Errcode] = ENclearreport(obj.LibEPANET);
        end
        function unload(obj)
            % Unload library and close the EPANET Toolkit system.
            %
            % Example:
            %   d.unload
            %
            % See also epanet, saveInputFile, closeNetwork.

            %ENclose(obj.LibEPANET);
            ENMatlabCleanup(obj.LibEPANET);
            fclose('all');
            files=dir('@#*');
            try delete([obj.InputFile(1:end-4), '.txt']), catch; end
            if ~isempty(files); delete('@#*'); end
            if exist([obj.BinTempfile(1:end-4), '.bin'], 'file')==2
                delete([obj.BinTempfile(1:end-4), '.bin']);
            end
            delete(obj.BinTempfile);
            if exist([obj.BinTempfile(1:end-4), '.txt'], 'file')==2
                delete([obj.BinTempfile(1:end-4), '.txt']);
            end
            if exist(obj.MSXTempFile, 'file')==2
                delete(obj.MSXTempFile);
            end
            disp('EPANET Class is unloaded')
        end
        function loadMSXFile(obj, msxname, varargin)
            if isempty(varargin)
                MSXMatlabSetup(obj, msxname);
            else
                MSXMatlabSetup(obj, msxname, varargin);
            end
        end
        function value = getMSXEquationsTerms(obj)
            [value, ~, ~] = getEquations(obj.MSXFile);
        end
        function value = getMSXEquationsPipes(obj)
            [~, value, ~] = getEquations(obj.MSXFile);
        end
        function value = getMSXEquationsTanks(obj)
            [~, ~, value] = getEquations(obj.MSXFile);
        end
        function value = getMSXOptions(obj)
            [value] = get_MSX_Options(obj.MSXFile, '', 1);
        end
        function value = getMSXTimeStep(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'timestep', 0);
            value = value.TimeStep;
        end
        function value = getMSXSolver(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'solver', 0);
            value = value.Solver;
        end
        function value = getMSXAreaUnits(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'area_units', 0);
            value = value.AreaUnits;
        end
        function value = getMSXRateUnits(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'rate_units', 0);
            value = value.RateUnits;
        end
        function value = getMSXRtol(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'rtol', 0);
            value = value.Rtol;
        end
        function value = getMSXAtol(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'atol', 0);
            value = value.Atol;
        end
        function value = getMSXCoupling(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'COUPLING', 0);
            value = value.Coupling;
        end
        function value = getMSXCompiler(obj)
            [value] = get_MSX_Options(obj.MSXFile, 'compiler', 0);
            value = value.Compiler;
        end
        function value = getMSXSpeciesCount(obj)
            % Species, Constants, Parameters, Patterns
            [obj.Errcode, value] = MSXgetcount(3, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXConstantsCount(obj)
            [obj.Errcode, value] = MSXgetcount(6, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXParametersCount(obj)
            [obj.Errcode, value] = MSXgetcount(5, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXPatternsCount(obj)
            [obj.Errcode, value] = MSXgetcount(7, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXSpeciesNameID(obj, varargin)
            if isempty(varargin)
                spcnt = obj.getMSXSpeciesCount;
                value = cell(1, spcnt);
                for i=1:spcnt
                    [obj.Errcode, len] = MSXgetIDlen(3, i, obj.MSXLibEPANET);
                    [obj.Errcode, value{i}]=MSXgetID(3, i, len, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            else
                k=1;
                value = cell(1, length(varargin{1}));
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(3, i, obj.MSXLibEPANET);
                    [obj.Errcode, value{k}]=MSXgetID(3, i, len, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    k=k+1;
                end
            end
        end
        function value = getMSXSpeciesType(obj)
            msxSpCnt = obj.getMSXSpeciesCount;
            value = cell(1, msxSpCnt);
            if msxSpCnt
                for i=1:msxSpCnt
                    [obj.Errcode, value{i}, ~, ~, ~] = MSXgetspecies(i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSpeciesUnits(obj)
            msxSpCnt = obj.getMSXSpeciesCount;
            value = cell(1, msxSpCnt);
            if msxSpCnt
                for i=1:msxSpCnt
                    [obj.Errcode, ~, value{i}, ~, ~] = MSXgetspecies(i, obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXSpeciesATOL(obj)
            value = [];
            msxSpCnt = obj.getMSXSpeciesCount;
            if msxSpCnt
                for i=1:msxSpCnt
                    [obj.Errcode, ~, ~, value, ~] = MSXgetspecies(i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSpeciesRTOL(obj)
            value = [];
            msxSpCnt = obj.getMSXSpeciesCount;
            if msxSpCnt
                for i=1:msxSpCnt
                    [obj.Errcode, ~, ~, ~, value] = MSXgetspecies(i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSpeciesIndex(obj, varargin)
            if isempty(varargin)
                value=1:obj.getMSXSpeciesCount;
                if isempty(value), value=[]; end
            elseif isa(varargin{1}, 'cell')
                for j=1:length(varargin{1})
                    [obj.Errcode, value(j)] = MSXgetindex(3, varargin{1}{j}, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET, 3, varargin{1});
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function value = getMSXConstantsNameID(obj)
            value={};
            if obj.getMSXConstantsCount
                for i=1:obj.getMSXConstantsCount
                    [obj.Errcode, len] = MSXgetIDlen(6, i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    [obj.Errcode, value{i}] = MSXgetID(6, i, len, obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXConstantsValue(obj)
            value=[];
            if obj.getMSXConstantsCount
                for i=1:obj.getMSXConstantsCount
                    [obj.Errcode, value(i)] = MSXgetconstant(i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXConstantsIndex(obj, varargin)
            if isempty(varargin)
                value=1:obj.getMSXConstantsCount;
                if isempty(value), value=[]; end
            elseif isa(varargin{1}, 'cell')
                for j=1:length(varargin{1})
                    [obj.Errcode, value(j)] = MSXgetindex(6, varargin{1}{j}, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET, 6, varargin{1});
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function value = getMSXParametersNameID(obj, varargin)
            if isempty(varargin)
                if ~obj.getMSXParametersCount, value={};return; end
                for i=1:obj.getMSXParametersCount
                    [obj.Errcode, len] = MSXgetIDlen(5, i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    [obj.Errcode, value{i}]=MSXgetID(5, i, len, obj.MSXLibEPANET);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(5, i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    [obj.Errcode, value{k}]=MSXgetID(5, i, len, obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXParametersIndex(obj, varargin)
            if isempty(varargin)
                value=1:obj.getMSXParametersCount;
                if isempty(value), value=[]; end
            elseif isa(varargin{1}, 'cell')
                for j=1:length(varargin{1})
                    [obj.Errcode, value(j)] = MSXgetindex(5, varargin{1}{j}, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = MSXgetindex(5, varargin{1}, obj.MSXLibEPANET);
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function value = getMSXParametersTanksValue(obj)
            value={};
            if ~obj.getMSXParametersCount, value=0;return;end
            if ~length(obj.NodeTankIndex), value=0;return;end
            for i=1:length(obj.NodeTankIndex)
                for j=1:obj.MSXParametersCount
                    [obj.Errcode, value{obj.NodeTankIndex(i)}(j)] = MSXgetparameter(0, obj.NodeTankIndex(i), j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXParametersPipesValue(obj)
            if ~obj.getMSXParametersCount, value=[];return; end
            for i=1:obj.getLinkPipeCount
                for j=1:obj.getMSXParametersCount
                    [obj.Errcode, value{i}(j)] = MSXgetparameter(1, i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXPatternsNameID(obj, varargin)
            if isempty(varargin)
                if ~obj.getMSXPatternsCount, value={};return; end
                for i=1:obj.getMSXPatternsCount
                    [obj.Errcode, len] = MSXgetIDlen(7, i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    [obj.Errcode, value{i}]=MSXgetID(7, i, len, obj.MSXLibEPANET);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(7, i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    [obj.Errcode, value{k}]=MSXgetID(7, i, len, obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXPatternsIndex(obj, varargin)
            if isempty(varargin)
                value=1:obj.getMSXPatternsCount;
                if isempty(value), value=[]; end
            elseif isa(varargin{1}, 'cell')
                for j=1:length(varargin{1})
                    [obj.Errcode, value(j)] = MSXgetindex(obj.MSXLibEPANET, 7, varargin{1}{j});
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET, 7, varargin{1});
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function value = getMSXPatternsLengths(obj, varargin)
            value =[];
            if isempty(varargin)
                if obj.getMSXPatternsCount
                    for i=obj.getMSXPatternsIndex
                        [obj.Errcode, value(i)]=MSXgetpatternlen(i, obj.MSXLibEPANET);
                        if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    end
                end
            elseif isa(varargin{1}, 'cell')
                for j=1:length(varargin{1})
                    [obj.Errcode, value(j)] = MSXgetpatternlen(obj.getMSXPatternsIndex(varargin{1}{j}), obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            elseif isa(varargin{1}, 'char')
                [obj.Errcode, value] = MSXgetpatternlen(obj.getMSXPatternsIndex(varargin{1}), obj.MSXLibEPANET);
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            elseif isa(varargin{1}, 'numeric')
                k=1;
                for i=varargin{1}
                    [obj.Errcode, value(k)]=MSXgetpatternlen(i, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    k=k+1;
                end
            end
        end
        function value = getMSXNodeInitqualValue(obj)
            if ~obj.getMSXSpeciesCount, value{1}(1)=0; return; end
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, value{i}(j)] = MSXgetinitqual(0, i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXLinkInitqualValue(obj)
            if ~obj.getMSXSpeciesCount, value{1}(1)=0; return; end
            for i=1:obj.getLinkCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, value{i}(j)] = MSXgetinitqual(1, i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSources(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, obj.MSXSourceType{i}{j}, obj.MSXSourceLevel{i}(j), obj.MSXSourcePatternIndex{i}(j)] = MSXgetsource(i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                    obj.MSXSourceTypeCode{i}(j)=find(strcmp(obj.MSXTYPESOURCE, obj.MSXSourceType{i}{j}))-2;
                end
            end
            SnodeID=obj.getMSXSourceNodeNameID;
           % value={obj.MSXSourceType, obj.MSXSourceTypeCode, obj.MSXSourceLevel, obj.MSXSourcePatternIndex, SnodeID};
            value.MSXSourceType=obj.MSXSourceType;
            value.MSXSourceTypeCode=obj.MSXSourceTypeCode;
            value.MSXSourceLevel=obj.MSXSourceLevel;
            value.MSXSourcePatternIndex=obj.MSXSourcePatternIndex;
            value.MSXSourceNodeNameID=SnodeID;
        end
        function value = getMSXSourceNodeNameID(obj)
            value = obj.getNodeNameID;
        end
        function value = getMSXSourceType(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, value{i}{j}, ~, ~] = MSXgetsource(i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSourceLevel(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, ~, value{i}(j), ~] = MSXgetsource(i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXSourcePatternIndex(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, ~, ~, value{i}(j)] = MSXgetsource(i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function value = getMSXPattern(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getMSXPatternsLengths);
            value=nan(obj.getMSXPatternsCount, tmpmaxlen);
            for i=1:obj.getMSXPatternsCount
                tmplength=obj.getMSXPatternsLengths(i);
                for j=1:tmplength
                    [obj.Errcode, value(i, j)] = MSXgetpatternvalue(i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i, j)=value(i, j-tmplength);
                    end
                end

            end
        end
        function value = getMSXPatternValue(obj, patternIndex, patternStep)
            % Mass flow rate per minute of a chemical source
            [obj.Errcode, value] = MSXgetpatternvalue(patternIndex, patternStep, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXSpeciesConcentration(obj, type, index, species)
            [obj.Errcode, value] = MSXgetqual(type, index, species, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function value = getMSXComputedQualitySpecie(obj, specie)
            % Return the quality nodes for specific specie
            % Example:
            %       MSX_comp = d.getMSXComputedQualitySpecie('CL2')
            %            MSX_comp.NodeQuality % row: time, col: node index
            %            MSX_comp.LinkQuality % row: time, col: link index
            %            MSX_comp.Time

            if obj.getMSXSpeciesCount==0
                value=0;
                return;
            end
            link_indices = 1:obj.getLinkCount;%for all link index
            node_indices = 1:obj.getNodeCount;%for all node index
            specie_name = obj.getMSXSpeciesIndex(specie);

            value.NodeQuality = nan(1, length(node_indices));
            value.LinkQuality = nan(1, length(node_indices));
            % Obtain a hydraulic solution
            obj.solveMSXCompleteHydraulics(obj.MSXLibEPANET);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.initializeMSXQualityAnalysis(0);
            % Retrieve species concentration at node
            k=1; tleft=1;t=0;
            value.Time(k, :)=0;
            time_step = obj.getMSXTimeStep;
            timeSmle=obj.getTimeSimulationDuration;%bug at time
            while(tleft>0 && obj.Errcode==0 && timeSmle~=t)
                [t, tleft]=obj.stepMSXQualityAnalysisTimeLeft;
                if t<time_step || t==time_step
                    if node_indices(end) < link_indices(end)
                        for lnk=link_indices
                            value.LinkQuality(k, lnk)=obj.getMSXLinkInitqualValue{lnk}(specie_name);
                            if lnk < node_indices(end) + 1
                                value.NodeQuality(k, lnk)=obj.getMSXNodeInitqualValue{lnk}(specie_name);
                            end
                        end
                    else
                        for lnk=node_indices
                            value.NodeQuality(k, lnk)=obj.getMSXNodeInitqualValue{lnk}(specie_name);
                            if lnk < link_indices(end) + 1
                                value.LinkQuality(k, lnk)=obj.getMSXLinkInitqualValue{lnk}(specie_name);
                            end
                        end
                    end
                else
                    if node_indices(end) < link_indices(end)
                        for lnk=link_indices
                            value.LinkQuality(k, lnk)=obj.getMSXSpeciesConcentration(1, lnk, specie_name);%link code 1
                            if lnk < node_indices(end) + 1
                                value.NodeQuality(k, lnk)=obj.getMSXSpeciesConcentration(0, lnk, specie_name);%node code0
                            end
                        end
                    else
                        for lnk=node_indices
                            value.NodeQuality(k, lnk)=obj.getMSXSpeciesConcentration(0, lnk, specie_name);%link code 1
                            if lnk < link_indices(end) + 1
                                value.LinkQuality(k, lnk)=obj.getMSXSpeciesConcentration(1, lnk, specie_name);%node code0
                            end
                        end
                    end
                end
                if k>1
                    value.Time(k, :)=t;
                end
                k=k+1;
            end
        end
        function value = getMSXComputedQualityNode(obj, varargin)
            if obj.getMSXSpeciesCount==0
                value=0;
                return;
            end
            if ~isempty(varargin)
                if length(varargin)==1
                    ss=varargin{1};%index node    %future work with argument ID
                    uu=1:obj.getMSXSpeciesCount;
                elseif length(varargin)==2
                    ss=varargin{1};%index node
                    uu=varargin{2};%index species
                end
            else
                ss=1:obj.getNodeCount;%index node
                uu=1:obj.getMSXSpeciesCount;
            end
            value.Quality = cell(1, length(ss));
            % Obtain a hydraulic solution
            obj.solveMSXCompleteHydraulics(obj.MSXLibEPANET);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.initializeMSXQualityAnalysis(0);
            % Retrieve species concentration at node
            k = 1; tleft = 1; t = 0;
            value.Time(k, :)=0;
            time_step = obj.getMSXTimeStep;
            timeSmle=obj.getTimeSimulationDuration;%bug at time
            while(tleft>0 && obj.Errcode==0 && timeSmle~=t)
                [t, tleft]=obj.stepMSXQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<time_step || t==time_step
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXNodeInitqualValue{nl}(j);
                                g=g+1;
                            end
                            i=i+1;
                        end
                    else
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXSpeciesConcentration(0, nl, j);%node code0
                                g=g+1;
                            end
                            i=i+1;
                        end
                    end
                else
                    if t<time_step || t==time_step
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXNodeInitqualValue{(nl)}(j);
                                g=g+1;
                            end
                            i=i+1;
                        end
                    else
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXSpeciesConcentration(0, (nl), j);%node code0
                                g=g+1;
                            end
                            i=i+1;
                        end
                    end
                end
                if k>1
                    value.Time(k, :)=t;
                end
                k=k+1;
            end
        end
        function value = getMSXComputedQualityLink(obj, varargin)
            if obj.getMSXSpeciesCount==0
                value=0;
                return;
            end
            if ~isempty(varargin)
                if length(varargin)==1
                    ss=varargin{1};%index link
                    uu=1:obj.getMSXSpeciesCount;
                elseif length(varargin)==2
                    ss=varargin{1};%index link
                    uu=varargin{2};%index species
                end
            else
                ss=1:obj.getLinkCount;%index link
                uu=1:obj.getMSXSpeciesCount;
            end
            value.Quality = cell(1, length(ss));
            % Obtain a hydraulic solution
            obj.solveMSXCompleteHydraulics(obj.MSXLibEPANET);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.initializeMSXQualityAnalysis(0);
            % Retrieve species concentration at node
            k = 1; tleft = 1;
            time_step = obj.getMSXTimeStep;
            value.Time(k, :)=0;
            while(tleft>0 && obj.Errcode==0)
                [t, tleft]=obj.stepMSXQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<time_step || t==time_step
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXLinkInitqualValue{nl}(j);
                                g=g+1;
                            end
                            i=i+1;
                        end
                    else
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXSpeciesConcentration(1, nl, j);
                                g=g+1;
                            end
                            i=i+1;
                        end
                    end
                else
                    if t<time_step || t==time_step
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXLinkInitqualValue{(nl)}(j);
                                g=g+1;
                            end
                            i=i+1;
                        end
                    else
                        i=1;
                        for nl=ss
                            g=1;
                            for j=uu
                                value.Quality{i}(k, g)=obj.getMSXSpeciesConcentration(1, (nl), j);%link code1
                                g=g+1;
                            end
                            i=i+1;
                        end
                    end
                end
                if k>1
                    value.Time(k, :)=t;
                end
                k=k+1;
            end
        end
        function plotMSXSpeciesNodeConcentration(obj, varargin)
            s=obj.getMSXComputedQualityNode(varargin{1}, varargin{2});
            nodesID=obj.getNodeNameID;
            SpeciesNameID=obj.getMSXSpeciesNameID;nd=1;
            for l=varargin{1}
                nodeID=nodesID(l);
                figure('Name', ['NODE ', char(nodeID)]);
                for i=varargin{2}
                    specie(:, i)=s.Quality{nd}(:, i);
                    time(:, i)=s.Time;
                end
                nd=nd+1;
                plot(time, specie);
                title(['NODE ', char(nodeID)]);
                ylabel('Quantity');
                xlabel('Time(s)');
                legend(SpeciesNameID(varargin{2}));
            end
        end
        function plotMSXSpeciesLinkConcentration(obj, varargin)
            s=obj.getMSXComputedQualityLink(varargin{1}, varargin{2});
            linksID=obj.getLinkNameID;
            SpeciesNameID=obj.getMSXSpeciesNameID;nd=1;
            for l=varargin{1}
                linkID=linksID(l);
                figure('Name', ['LINK ', char(linkID)]);
                for i=varargin{2}
                    specie(:, i)=s.Quality{nd}(:, i);
                    time(:, i)=s.Time;
                end
                nd=nd+1;
                plot(time, specie);
                title(['LINK ', char(linkID)]);
                ylabel('Quantity');
                xlabel('Time(s)');
                legend(SpeciesNameID(varargin{2}));
            end
        end
        function value = getMSXError(obj, Errcode)
            [obj.Errcode, value] = MSXgeterror(Errcode, obj.MSXLibEPANET);
            if Errcode == 519
                error('Please check the MSX file. Maybe node/link ids do not exist in the input file.');
            end
        end
        function solveMSXCompleteHydraulics(obj, varargin)
            [obj.Errcode] = MSXsolveH(obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function solveMSXCompleteQuality(obj, varargin)
            [obj.Errcode] = MSXsolveQ(obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function writeMSXReport(obj, varargin)
            [obj.Errcode]=MSXreport(obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function [status, result] = writeMSXReportExe(obj, varargin)
            if isempty(varargin)
                rptfile=['@#', char(java.util.UUID.randomUUID), '.txt'];
            else
                rptfile=varargin{1};
            end
            [status, result] = runMSX(obj, rptfile);
        end
        %       function value = getMSXComputedResultsBinary(obj)
        %           uuID = char(java.util.UUID.randomUUID);
        %           binfile=['@#', uuID, '.bin'];
        %           obj.solveMSXCompleteHydraulics;
        %           obj.saveHydraulicsOutputReportingFile;
        %           obj.solveMSXCompleteQuality;
        %           obj.saveMSXQualityFile(binfile);
        %           value = readMSXBinaryFile(binfile);
        %       end
        %       function value = getMSXComputedResultsBinaryExe(obj)
        %           uuID = char(java.util.UUID.randomUUID);
        %           rptfile=['@#', uuID, '.txt'];
        %           binfile=['@#', uuID, '.bin'];
        %           runMSX(obj, rptfile, binfile);
        %           value = readMSXBinaryFile(binfile);
        %       end
        function index = addMSXPattern(obj, varargin)
            index=-1;
            if nargin==2
                [obj.Errcode] = MSXaddpattern(varargin{1}, obj.MSXLibEPANET);
                [obj.Errcode, index] = MSXgetindex(obj.MSXLibEPANET, 7, varargin{1});
            elseif nargin==3
                [obj.Errcode] = MSXaddpattern(varargin{1}, obj.MSXLibEPANET);
                [obj.Errcode, index] = MSXgetindex(obj.MSXLibEPANET, 7, varargin{1});
                setMSXPattern(obj, index, varargin{2});
            end
        end
        function setMSXSources(obj, nodeID, speciesID, sourcetype, concentration, patID)
            node = obj.getNodeIndex(nodeID);
            species = obj.getMSXSpeciesIndex(speciesID);
            type = find(strcmpi(obj.MSXTYPESOURCE, sourcetype))-2;
            pat = obj.getMSXPatternsIndex(patID);
            MSXsetsource(node, species, type, concentration, pat, obj.MSXLibEPANET);
        end
        function setMSXConstantsValue(obj, value)
            for i=1:length(value)
                [obj.Errcode] = MSXsetconstant(i, value(i), obj.MSXLibEPANET);
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function setMSXParametersTanksValue(obj, NodeTankIndex, paramindex, value)
            if ~sum(NodeTankIndex==obj.NodeTankIndex)
                fprintf('>> Invalid Tank Index <<\n');obj.NodeTankIndex
                return;
            end
            [obj.Errcode] = MSXsetparameter(0, NodeTankIndex, paramindex, value, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function setMSXParametersPipesValue(obj, pipeIndex, value)
            for i=1:length(value)
                [obj.Errcode] = MSXsetparameter(1, pipeIndex, i, value(i), obj.MSXLibEPANET);
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function setMSXNodeInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1}, 2)
                    [obj.Errcode] = MSXsetinitqual(0, i, j, value{i}(j), obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function setMSXLinkInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1}, 2)
                    [obj.Errcode] = MSXsetinitqual(1, i, j, value{i}(j), obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                end
            end
        end
        function setMSXPattern(obj, pat, patternVector)
            if ischar(pat)
                pat=obj.getMSXPatternsIndex(pat);
            end
            nfactors=length(patternVector);
            [obj.Errcode] = MSXsetpattern(pat, patternVector, nfactors, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function setMSXPatternMatrix(obj, patternMatrix)
            nfactors=size(patternMatrix, 2);
            for i=1:size(patternMatrix, 1)
                [obj.Errcode] = MSXsetpattern(i, patternMatrix(i, :), nfactors, obj.MSXLibEPANET);
                if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
            end
        end
        function setMSXPatternValue(obj, index, patternTimeStep, patternFactor)
            [obj.Errcode] = MSXsetpatternvalue(index, patternTimeStep, patternFactor, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function setMSXTimeStep(obj, timestep)
            setMSXOptions(obj, 'timestep', timestep);
        end
        function setMSXAreaUnitsFT2(obj)
            setMSXOptions(obj, 'areaunits', 'FT2');
        end
        function setMSXAreaUnitsM2(obj)
            setMSXOptions(obj, 'areaunits', 'M2');
        end
        function setMSXAreaUnitsCM2(obj)
            setMSXOptions(obj, 'areaunits', 'CM2');
        end
        function setMSXRateUnitsSEC(obj)
            setMSXOptions(obj, 'rateunits', 'SEC');
        end
        function setMSXRateUnitsMIN(obj)
            setMSXOptions(obj, 'rateunits', 'MIN');
        end
        function setMSXRateUnitsHR(obj)
            setMSXOptions(obj, 'rateunits', 'HR');
        end
        function setMSXRateUnitsDAY(obj)
            setMSXOptions(obj, 'rateunits', 'DAY');
        end
        function setMSXSolverEUL(obj)
            setMSXOptions(obj, 'solver', 'EUL');
        end
        function setMSXSolverRK5(obj)
            setMSXOptions(obj, 'solver', 'RK5');
        end
        function setMSXSolverROS2(obj)
            setMSXOptions(obj, 'solver', 'ROS2');
        end
        function setMSXCouplingFULL(obj)
            setMSXOptions(obj, 'coupling', 'FULL');
        end
        function setMSXCouplingNONE(obj)
            setMSXOptions(obj, 'coupling', 'NONE');
        end
        function setMSXCompilerNONE(obj)
            setMSXOptions(obj, 'compiler', 'NONE');
        end
        function setMSXCompilerVC(obj)
            setMSXOptions(obj, 'compiler', 'VC');
        end
        function setMSXCompilerGC(obj)
            setMSXOptions(obj, 'compiler', 'GC');
        end
        function setMSXAtol(obj, atol)
            setMSXOptions(obj, 'atol', atol);
        end
        function setMSXRtol(obj, rtol)
            setMSXOptions(obj, 'rtol', rtol);
        end
        function saveMSXQualityFile(obj, outfname)
            [obj.Errcode]=MSXsaveoutfile(outfname, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function useMSXHydraulicFile(obj, hydname)
            [obj.Errcode]=MSXusehydfile(hydname, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function initializeMSXQualityAnalysis(obj, flag)
            [obj.Errcode] = MSXinit(flag, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function [t, tleft]= stepMSXQualityAnalysisTimeLeft(obj)
            [obj.Errcode, t, tleft] = MSXstep(obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function saveMSXFile(obj, msxname)
            [obj.Errcode] = MSXsavemsxfile(msxname, obj.MSXLibEPANET);
            if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
        end
        function msx = writeMSXFile(obj, msx)
            % Checkout example: /examples/EX15_write_msx_file.m
            space=5;
            f = writenewTemp(msx.FILENAME);
            fprintf(f,'[TITLE]\n');
            if isfield(msx,'title')
                fprintf(f,msx.TITLE);
            end

            fprintf(f,'\n\n[OPTIONS]\n');
            options = {'AREA_UNITS', 'RATE_UNITS', 'SOLVER', 'COUPLING', 'COMPILER',...
                'TIMESTEP', 'ATOL', 'RTOL'};
            spaces=blanks(space);

            for i=1:length(options)
                if isfield(msx,options{i})
                    fprintf(f,num2str(options{i}));
                    fprintf(f,spaces);
                    fprintf(f,num2str(msx.(options{i})));
                    fprintf(f,'\n');
                end
            end

            FIELDS = {'SPECIES', 'COEFFICIENTS', 'TERMS', 'PIPES', ...
                'TANKS', 'SOURCES', 'QUALITY', 'GLOBAL', 'PARAMETERS', 'PATTERNS'};
            for sect=FIELDS
                section = char(sect);
                if ~strcmpi(section, 'GLOBAL')
                    fprintf(f, ['\n[', section, ']\n']);
                end
                sp_field = upper(section);
                if isfield(msx, upper(section))
                    species_all = eval(['msx.', sp_field]);
                    for i=1:length(species_all)
                        fprintf(f, char(species_all{i}));
                        fprintf(f, '\n');
                    end
                end
            end

            fprintf(f, '[REPORT]\n');
            fprintf(f, 'NODES ALL\n');
            fprintf(f, 'LINKS ALL\n');
            fclose(f);
        end
        function unloadMSX(obj)
            MSXclose(obj);
            MSXMatlabCleanup(obj);
            fclose('all');
            disp('MSX unloaded');
        end
        function ToolkitConstants = getToolkitConstants(obj)
            if obj.getVersion <= 20101
                fparam = '.h'; % error('This version of EMT support version 2.2 of EPANET.')
            else
                fparam = '_enums.h';
            end
            if isunix
                file = [obj.LibEPANETpath, 'epanet2', fparam];
            else
                file = [obj.LibEPANETpath, obj.LibEPANET, fparam];
            end
            if isdeployed
            %file=[file(1:end-1), 'txt'];%epanet2.h-->epanet2.txt
                  file = 'epanet2_enums.txt';%epanet2_enums.h-->epanet2_enums.txt
            end
            fid = fopen(file);
            tline = fgetl(fid);
            i=1; constants={};
            while ischar(tline)
               if ~isempty(regexp(tline, 'typedef enum', 'match'))
                   tline = fgetl(fid);
                   while isempty(regexp(tline, '}', 'match')) || isempty(tline)
                       while isempty(tline)
                           tline = fgetl(fid);
                       end
                       find_comm = strfind(tline, '//');
                       if ~isempty(find_comm)
                           tline = tline(1:find_comm-1);
                       end
                       n = regexp(tline, {'\w*EN_\w*', '\d*'}, 'match');
                       if length(n{1})>1
                           constants(i) = n{1}(1);
                       else
                           constants(i) = n{1};
                       end
                       codes(i) = str2double(n{2}{end});
                       ToolkitConstants.(constants{i}) = codes(i);
                       i=i+1;
                       tline = fgetl(fid);
                   end
               else
                   n = regexp(tline, {'\w*EN_\w*', '\d*'}, 'match');
                   if sum(cellfun(@isempty, n)), tline = fgetl(fid); continue; end
                   if length(constants)>118, break; end %temporary
                   constants(i) = n{1};
                   codes(i) = str2double(n{2}{1});
                   ToolkitConstants.(constants{i}) = codes(i);
                   i=i+1;
                   tline = fgetl(fid);
               end
            end
            fclose(fid);
        end
        function obj = BinUpdateClass(obj)
            sect=0;i=1;t=1;q=1;
            typecode=0;x=1;b=1;d=1;
            if obj.Bin
                obj.saveInputFile(obj.BinTempfile);
            end
            obj.BinUnits_SI_Metric=0; obj.BinUnits_US_Customary=0;
            [~, info] = obj.readInpFile;
            for hc=1:length(info)
                tline = info{hc};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                    % [JUNCTIONS] section
                    if strcmpi(tok(1:5), '[JUNC')
                        sect=1;k=1;
                        obj.BinNodeJunctionNameID={};
                        obj.BinNodeJunctionIndex=[];
                        obj.BinNodeJunctionElevation=[];
                        obj.BinNodeType={};
                        obj.BinNodeJunctionsBaseDemands=[];
                        obj.BinNodeJunctionsBaseDemandsID={};
                        obj.BinNodeJunDemandPatternNameID={};
                        continue;
                        % [RESERVOIRS] section
                    elseif strcmpi(tok(1:5), '[RESE')
                        sect=2;r=1;
                        obj.BinNodeReservoirNameID={};
                        obj.BinNodeReservoirIndex=[];
                        obj.BinNodeReservoirElevation=[];
                        obj.BinNodeResDemandPatternNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5), '[TANK')
                        sect=3;p=1;
                        obj.BinNodeTankNameID={};
                        obj.BinNodeTankIndex=[];
                        obj.BinNodeTankElevation=[];
                        obj.BinNodeTankInitialLevel=[];
                        obj.BinNodeTankMinimumWaterLevel=[];
                        obj.BinNodeTankMaximumWaterLevel=[];
                        obj.BinNodeTankDiameter=[];
                        obj.BinNodeTankMinimumWaterVolume=[];
                        continue;
                        % [PIPES] section
                    elseif strcmpi(tok(1:5), '[PIPE')
                        sect=4;
                        obj.BinLinkPipeNameID={};
                        obj.BinLinkPipeIndex=[];
                        obj.BinLinkFromNode={};
                        obj.BinLinkToNode={};
                        obj.BinLinkPipeLengths=[];
                        obj.BinLinkPipeDiameters=[];
                        obj.BinLinkPipeRoughness=[];
                        obj.BinLinkPipeMinorLoss=[];
                        obj.BinLinkType={};

                        obj.BinNodeJunctionCount = length(obj.BinNodeJunctionNameID);
                        obj.BinNodeReservoirCount = length(obj.BinNodeReservoirNameID);
                        obj.BinNodeTankCount = length(obj.BinNodeTankNameID);
                        obj.BinNodeTankReservoirCount = obj.BinNodeReservoirCount + obj.BinNodeTankCount;
                        obj.BinNodeNameID=[obj.BinNodeJunctionNameID obj.BinNodeReservoirNameID obj.BinNodeTankNameID];
                        obj.BinNodeCount=obj.BinNodeJunctionCount+obj.BinNodeTankReservoirCount;
                        obj.BinNodeInitialQuality=zeros(1, obj.BinNodeCount);
                        continue;
                        % [PUMPS] section
                    elseif strcmpi(tok(1:5), '[PUMP')
                        sect=5;
                        obj.BinLinkPumpNameID={};
                        obj.BinLinkPumpIndex=[];
                        obj.BinLinkPumpPatterns={};
                        obj.BinLinkPumpPower=[];
                        obj.BinLinkPumpNameIDPower={};
                        continue;
                        % [VALVES] section
                    elseif strcmpi(tok(1:5), '[VALV')
                        sect=6;
                        obj.BinLinkValveNameID={};
                        obj.BinLinkValveIndex=[];
                        obj.BinLinkValveDiameters=[];
                        obj.BinLinkValveType={};
                        obj.BinLinkValveSetting=[];
                        obj.BinLinkValveMinorLoss=[];
                        continue;
                        % [CURVES] section
                    elseif strcmpi(tok(1:5), '[CURV')
                        sect=7;
                        obj.BinCurveAllLines={};
                        obj.BinCurveXvalue=[];
                        obj.BinCurveYvalue=[];
                        obj.BinCurveAllLines={};
                        continue;
                        % [DEMANDS] section
                    elseif strcmpi(tok(1:5), '[DEMA') %&& max(obj.BinNodeJunctionsBaseDemands)==0
                        sect=8;d=1;pd=1;
                        continue;
                        % [STATUS] section
                    elseif strcmpi(tok(1:5), '[STAT')
                        sect=9;d=1;
                        obj.BinLinkInitialStatus={};
                        obj.BinLinkInitialStatusNameID={};
                        obj.BinCountStatuslines=0;
                        continue;
                        % [PATTERNS]
                    elseif strcmpi(tok(1:5), '[PATT')
                        sect=10;
                        obj.BinPatternLengths=[];
                        obj.BinPatternNameID={};
                        obj.BinPatternValue={};
                        obj.BinCountPatternlines=0;d=1;h=1;
                        continue;
                        % [CONTROLS] section
                    elseif strcmpi(tok(1:5), '[CONT')
                        sect=11;d=1;
                        obj.BinControlsInfo={};
                        obj.BinControlLinksID={};
                        obj.BinControlNodesID={};
                        obj.BinPatternCount=length(obj.BinPatternNameID);
                        tmpmaxlen=max(obj.BinPatternLengths);
                        obj.BinPatternMatrix=nan(obj.BinPatternCount, tmpmaxlen);
                        for i=1:obj.BinPatternCount
                            tmplength=obj.BinPatternLengths(i);
                            for j=1:tmplength
                                obj.BinPatternMatrix(i,j)=double(obj.BinPatternValue{i}(j));
                            end
                            if tmplength<tmpmaxlen
                                for j=(tmplength+1):tmpmaxlen
                                    obj.BinPatternMatrix(i, j)=obj.BinPatternMatrix(i, j-tmplength);
                                end
                            end
                        end
                        continue;
                        % [RULES] section
                    elseif strcmpi(tok(1:5), '[RULE')
                        sect=20;d=1;obj.BinRulesCount=0;
                        obj.BinRulesControlsInfo={};
                        obj.BinRulesControlLinksID={};
                        obj.BinRulesControlNodesID={};
                        continue;
                        % [QUALITY] section
                    elseif strcmpi(tok(1:5), '[QUAL')
                        sect=12;d=1;
                        obj.BinCountInitialQualitylines=0;
                        continue;
                        % [SOURCES] section
                    elseif strcmpi(tok(1:5), '[SOUR')
                        sect=13;
                        obj.BinNodeSourcePatternIndex = nan(1, obj.BinNodeCount);
                        obj.BinNodeSourceQuality = nan(1, obj.BinNodeCount);
                        obj.BinNodeSourceTypeIndex = nan(1, obj.BinNodeCount);
                        obj.BinNodeSourceType = cell(1, obj.BinNodeCount);
                        obj.BinNodeSourcePatternNameID = cell(1, obj.BinNodeCount);
                        continue;
                        % [MIXING] section
                    elseif strcmpi(tok(1:5), '[MIXI')
                        sect=14;d=1;
                        obj.BinNodeTankMixModel={};
                        obj.BinNodeTankMixID={};
                        obj.BinNodeTankMinimumFraction=[];
                        continue;
                        % [REACTIONS] section
                    elseif strcmpi(tok(1:5), '[REAC')
                        sect=15;d=1;
                        obj.BinLinkGlobalBulkReactionCoeff=[];
                        obj.BinLinkGlobalWallReactionCoeff=[];
                        obj.BinLinkBulkReactionCoeff=[];
                        obj.BinLinkWallReactionCoeff=[];
                        obj.BinCountReactionlines=0;
                        continue;
                        % [TIMES] section
                    elseif strcmpi(tok(1:5), '[TIME')
                        sect=16;d=1;
                        obj.BinTimeSimulationDuration=[];
                        obj.BinTimeHydraulicStep=[];
                        obj.BinTimeQualityStep=[];
                        obj.BinTimePatternStep=[];
                        obj.BinTimePatternStart=[];
                        obj.BinTimeReportingStep=[];
                        obj.BinTimeReportingStart=[];
                        obj.BinTimeStatisticsIndex=[];
                        obj.BinTimeStatistics='';
                        continue;
                        % [OPTIONS] section
                    elseif strcmpi(tok(1:5), '[OPTI')
                        sect=17;
                        vx = NaN(obj.BinNodeCount, 1);
                        vy = NaN(obj.BinNodeCount, 1);
                        vertx = cell(obj.BinLinkCount, 1);
                        verty = cell(obj.BinLinkCount, 1);
                        indexV=1;
                        obj.BinLinkFlowUnits={};
                        obj.BinOptionsHeadloss={};
                        obj.BinNodePressureUnits={};
                        obj.BinOptionsSpecificGravity=[];
                        obj.BinOptionsViscosity=[];
                        obj.BinOptionsMaxTrials=[];
                        obj.BinOptionsAccuracyValue=[];
                        obj.BinOptionsUnbalanced={};
                        obj.BinOptionsPattern=[];
                        obj.BinOptionsPatternDemandMultiplier=[];
                        obj.BinOptionsEmitterExponent=[];
                        obj.BinQualityType={};
                        obj.BinQualityCode=[];
                        obj.BinQualityTraceNodeIndex=[];
                        obj.BinQualityTraceNodeID={};
                        obj.BinQualityUnits={};
                        obj.BinOptionsDiffusivity=[];
                        obj.BinOptionsQualityTolerance=[];

                        obj.BinLinkPipeCount = length(obj.BinLinkPipeNameID);
                        obj.BinLinkPumpCount = length(obj.BinLinkPumpNameID);
                        obj.BinLinkValveCount = length(obj.BinLinkValveNameID);

                        obj.BinLinkNameID=[obj.BinLinkPipeNameID obj.BinLinkPumpNameID obj.BinLinkValveNameID];
                        obj.BinLinkCount=obj.BinLinkPipeCount+obj.BinLinkPumpCount+obj.BinLinkValveCount;
                        continue;
                        % [COORDINATES] section
                    elseif strcmpi(tok(1:5), '[COOR')
                        sect=18;
                        continue;
                        % [VERTICES] section
                    elseif strcmpi(tok(1:5), '[VERT')
                        sect=19;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                atline = checktlines(tline);
                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    obj.BinNodeJunctionNameID{k}=atline{1};
                    obj.BinNodeJunctionIndex(k)=k;
                    obj.BinNodeJunctionElevation(k)=str2double(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            obj.BinNodeJunctionsBaseDemands(k)=str2double(atline{3});
                            if length(atline)>3
                                if ~sum(atline{4}==';')
                                    obj.BinNodeJunDemandPatternNameID{k}=atline{4};
                                else
                                    obj.BinNodeJunDemandPatternNameID{k}='';
                                end
                            else
                                obj.BinNodeJunDemandPatternNameID{k}='';
                            end
                        end
                    end
                    k=k+1;
                elseif sect==2
                    obj.BinNodeReservoirNameID{r}=atline{1};
                    obj.BinNodeReservoirIndex(r)=k;
                    obj.BinNodeReservoirElevation(r)=str2double(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            obj.BinNodeResDemandPatternNameID{r}=atline{3};
                        else
                            obj.BinNodeResDemandPatternNameID{r}='';
                        end
                    else
                        obj.BinNodeResDemandPatternNameID{r}='';
                    end
                    k=k+1;
                    r=r+1;
                elseif sect==3
                    obj.BinNodeTankNameID{p}=atline{1};
                    obj.BinNodeTankIndex(p)=k;
                    obj.BinNodeTankElevation(p)=str2double(atline{2});
                    obj.BinNodeTankInitialLevel(p)=single(str2double(atline{3}));
                    obj.BinNodeTankMinimumWaterLevel(p)=str2double(atline{4});
                    obj.BinNodeTankMaximumWaterLevel(p)=single(str2double(atline{5}));
                    obj.BinNodeTankDiameter(p)=str2double(atline{6});
                    obj.BinNodeTankMinimumWaterVolume(p)=single(str2double(atline{7}));
                    k=k+1;
                    p=p+1;
                    % Links
                elseif sect==4
                    obj.BinLinkPipeNameID{t}=atline{1};
                    obj.BinLinkPipeIndex(t)=t;
                    obj.BinLinkFromNode{t}=atline{2};
                    obj.BinLinkToNode{t}=atline{3};
                    obj.BinLinkPipeLengths(t)=str2double(atline{4});
                    obj.BinLinkPipeDiameters(t)=str2double(atline{5});
                    obj.BinLinkPipeRoughness(t)=str2double(atline{6});
                    obj.BinLinkPipeMinorLoss(t)=str2double(atline{7});
                    BinCNameID={};
                    if length(atline)>7
                        obj.BinLinkPipeStatus{t}=atline{8};
                    else
                        obj.BinLinkPipeStatus{t}='Open';
                    end
                    t=t+1;
                elseif sect==5
                    obj.BinLinkPumpNameID{q}=atline{1};
                    obj.BinLinkPumpIndex(q)=t;
                    obj.BinLinkFromNode{t}=atline{2};
                    obj.BinLinkToNode{t}=atline{3};
                    if strcmp(regexp(tline, '\w*HEAD*\w', 'match'), 'HEAD')
                        obj.BinLinkPumpCurveNameID{q}=atline{5};
                    elseif strcmp(regexp(tline, '\w*POWER*\w', 'match'), 'POWER')
                        obj.BinLinkPumpPower(q)=str2double(atline{5});
                        obj.BinLinkPumpNameIDPower{q}=atline{1};
                    end
                    if length(atline)>6
                        obj.BinLinkPumpPatterns{q}=atline{7};
                    end
                    t=t+1;
                    q=q+1;
                elseif sect==6
                    obj.BinLinkValveNameID{i}=atline{1};
                    obj.BinLinkValveIndex(i)=t;
                    obj.BinLinkFromNode{t}=atline{2};
                    obj.BinLinkToNode{t}=atline{3};
                    obj.BinLinkValveDiameters(i)=str2double(atline{4});
                    obj.BinLinkValveType{i}=atline{5};
                    obj.BinLinkValveSetting(i)=str2double(atline{6});
                    if length(atline)>6
                        obj.BinLinkValveMinorLoss(i)=str2double(atline{7});
                    end
                    t=t+1;
                    i=i+1;
                    % Curves
                elseif sect==7
                    ee=regexp(tline, '\w*EFFICIENCY*\w', 'match');
                    nn=regexp(tline, '\w*VOLUME*\w', 'match');
                    kk=regexp(tline, '\w*HEADLOSS*\w', 'match');

                    if strcmp(ee, 'EFFICIENCY'), typecode=1;   % EFFICIENCY
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif strcmp(nn, 'VOLUME'), typecode=2;   % VOLUME
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif strcmp(kk, 'HEADLOSS'), typecode=3; % HEADLOSS
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif (~length(strcmp(nn, 'VOLUME')) || ~length(strcmp(ee, 'EFFICIENCY')) || ~length(strcmp(kk, 'HEADLOSS'))) &&  (tok(1)==';'), typecode=0; % HEADLOSS
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    else
                        obj.BinCurveTypes(x)=typecode;
                    end
                    a = textscan(tline, '%s %f %f');l=1;
                    aa = regexp(tline, '\s*', 'split');
                    if isempty(aa{l})
                        l=l+1;
                    end
                    BinCNameID{x}=aa{l};
                    obj.BinCurveXvalue(x)=a{2};
                    obj.BinCurveYvalue(x)=a{3};
                    obj.BinCurveAllLines{b}=tline;
                    x=x+1;b=b+1;
                    % Demands
                elseif sect==8
                    indd=find(strcmpi(obj.BinNodeNameID, atline{1}));
                    if ~isempty(obj.BinNodeJunctionsBaseDemandsID)
                        if strcmp(obj.BinNodeJunctionsBaseDemandsID{end}, atline{1})
                            pd=pd+1;obj.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2}));
                        else
                            obj.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2})); pd=1;
                        end
                    else
                        obj.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2})); pd=1;
                    end
                    if length(atline)>2
                        if ~isempty(obj.BinNodeJunctionsBaseDemandsID)
                            if strcmp(obj.BinNodeJunctionsBaseDemandsID{end}, atline{1})
                                obj.BinNodeJunDemandPatternNameID{indd}=atline{3};
                            else
                                obj.BinNodeJunDemandPatternNameID{indd}=atline{3};
                            end
                        else
                            obj.BinNodeJunDemandPatternNameID{indd}=atline{3};
                        end
                    else
                        obj.BinNodeJunDemandPatternNameID{indd}='';
                    end
                    obj.BinNodeJunctionsBaseDemandsID{d}=atline{1};
                    d=d+1;
                    % Status
                elseif sect==9
                    obj.BinLinkInitialStatus{d}=atline{2};
                    obj.BinLinkInitialStatusNameID{d}=atline{1};
                    obj.BinCountStatuslines=d;
                    d=d+1;
                    % Patterns
                elseif sect==10
                    obj.BinPatternNameID{d}=atline{1};
                    dd=length(obj.BinPatternNameID);
                    obj.BinPatternNameID=unique(obj.BinPatternNameID);
                    d=length(obj.BinPatternNameID);
                    if dd>1 && dd~=d
                        obj.BinPatternLengths(d)=length(atline(2:end))+obj.BinPatternLengths(d);
                        obj.BinPatternValue{d}=single([obj.BinPatternValue{d} str2num(char(atline(2:end)))']);
                    else
                        obj.BinPatternLengths(d)=length(atline(2:end));
                        obj.BinPatternValue{d}=single(str2num(char(atline(2:end)))');
                    end
                    obj.BinCountPatternlines=h;
                    d=d+1;h=h+1;
                    % Controls
                elseif sect==11
                    obj.BinControlsInfo{d}=atline;
                    obj.BinControlLinksID{d}=atline{2};
                    t = regexp(tline, '\w*TIME\w*', 'match');
                    if length(t)==0
                        obj.BinControlNodesID{d}=atline{6};
                    end
                    d=d+1;
                    % Rules
                elseif sect==20
                    if strcmpi(atline{1}, {'RULE'})
                        obj.BinRulesCount=obj.BinRulesCount+1;d=1;
                    end
                    obj.BinRulesControlsInfo{obj.BinRulesCount}{d}=atline;

                    if sum(strcmpi(atline{2}, {'LINK', 'PIPE', 'PUMP', 'VALVE'}))
                        obj.BinRulesControlLinksID{obj.BinRulesCount}{d}=atline{3};
                    elseif sum(strcmpi(atline{2}, {'NODE', 'JUNCTION', 'RESERVOIR', 'TANK'}))
                        obj.BinRulesControlNodesID{obj.BinRulesCount}{d}=atline{3};
                    end
                    d=d+1;
                    % Quality
                elseif sect==12
                    h=find(strcmpi(obj.BinNodeNameID, atline{1}));
                    obj.BinNodeInitialQuality(h)=str2double(atline{2});
                    obj.BinCountInitialQualitylines=d;
                    d=d+1;
                    % Sources
                elseif sect==13
                    indexPat=find(strcmpi(obj.BinPatternNameID, atline{1}));
                    indexNode=find(strcmpi(obj.BinNodeNameID, atline{1}));
                    if length(atline)>3
                        obj.BinNodeSourcePatternIndex(indexPat)= find(strcmpi(obj.BinPatternNameID, atline{4}));
                        obj.BinNodeSourcePatternNameID{indexNode}=atline{4};
                    end
                    obj.BinNodeSourceQuality(indexNode)=str2double(atline{3});
                    obj.BinNodeSourceTypeIndex(indexNode)=find((strcmpi(obj.TYPESOURCE, atline{2})-1)>-1)-1;
                    obj.BinNodeSourceType{indexNode}=obj.TYPESOURCE{obj.BinNodeSourceTypeIndex(indexNode)+1};
                    % Mixing
                elseif sect==14
                    obj.BinNodeTankMixID{d}=atline{1};
                    obj.BinNodeTankMixModel{d}=atline{2};
                    obj.BinNodeTankMinimumFraction(d)=str2double(atline{3});
                    d=d+1;
                    % Reactions
                elseif sect==15
                    if strcmpi(atline{1}, 'GLOBAL') && strcmpi(atline{2}, 'BULK')
                        obj.BinLinkGlobalBulkReactionCoeff=str2double(atline{3});
                    elseif strcmpi(atline{1}, 'GLOBAL') && strcmpi(atline{2}, 'WALL')
                        obj.BinLinkGlobalWallReactionCoeff=str2double(atline{3});
                        obj.BinLinkWallReactionCoeff=obj.BinLinkGlobalWallReactionCoeff*ones(1, obj.BinLinkCount);
                        obj.BinLinkBulkReactionCoeff=obj.BinLinkGlobalBulkReactionCoeff*ones(1, obj.BinLinkCount);
                        obj.BinLinkWallReactionCoeff(obj.BinLinkPumpIndex)=0;
                        obj.BinLinkWallReactionCoeff(obj.BinLinkValveIndex)=0;
                        obj.BinLinkBulkReactionCoeff(obj.BinLinkPumpIndex)=0;
                        obj.BinLinkBulkReactionCoeff(obj.BinLinkValveIndex)=0;
                    end
                    if strcmpi(atline{1}, 'BULK')
                        LinkIndex = find(strcmpi(obj.BinLinkNameID, atline{2}));
                        obj.BinLinkBulkReactionCoeff(LinkIndex)=str2double(atline{3});
                    elseif strcmpi(atline{1}, 'WALL')
                        LinkIndex = find(strcmpi(obj.BinLinkNameID, atline{2}));
                        obj.BinLinkWallReactionCoeff(LinkIndex)=str2double(atline{3});
                    end
                    obj.BinCountReactionlines=d;
                    d=d+1;
                    % Times
                elseif sect==16
                    r=atline{2};
                    if length(atline)>2
                        if find(~strcmpi(atline{end}, {'HOURS', 'MIN', 'SECONDS', 'MINUTES', 'DAYS'})==0)
                            r=atline{end-1};
                        else
                            r=atline{end};
                        end
                    end
                    obj = getTimes(obj, r, atline, obj);
                    % Options
                elseif sect==17
                    obj = getOptionsValues(obj, obj, atline);
                    % Coordinates
                elseif sect==18
                    A = textscan(tline, '%s %f %f');
                    % get the node index
                    a=strcmp(A{1}, obj.BinNodeNameID);
                    index=strfind(a, 1);
                    if isempty(index), return; end
                    vx(index) = A{2};
                    vy(index) = A{3};
                    % Vertices
                elseif sect==19
                    A = textscan(tline, '%s %f %f');
                    if isempty(A), return;  end
                    vertx(indexV) = A(2);
                    verty(indexV) = A(3);
                    indexV = indexV + 1;
                end
            end
            if ~sum(obj.BinLinkBulkReactionCoeff)
                if isempty(obj.BinLinkGlobalBulkReactionCoeff), obj.BinLinkGlobalBulkReactionCoeff=0;end
                obj.BinLinkBulkReactionCoeff=obj.BinLinkGlobalBulkReactionCoeff*ones(1, obj.BinLinkCount);
                obj.BinLinkBulkReactionCoeff(obj.BinLinkPumpIndex)=0;
                obj.BinLinkBulkReactionCoeff(obj.BinLinkValveIndex)=0;
            end
            if ~sum(obj.BinLinkWallReactionCoeff)
                if isempty(obj.BinLinkGlobalWallReactionCoeff), obj.BinLinkGlobalWallReactionCoeff=0;end
                obj.BinLinkWallReactionCoeff=obj.BinLinkGlobalWallReactionCoeff*ones(1, obj.BinLinkCount);
                obj.BinLinkWallReactionCoeff(obj.BinLinkPumpIndex)=0;
                obj.BinLinkWallReactionCoeff(obj.BinLinkValveIndex)=0;
            end
            obj.BinNodeCoordinates{1} = vx;
            obj.BinNodeCoordinates{2} = vy;
            obj.BinNodeCoordinates{3} = vertx;
            obj.BinNodeCoordinates{4} = verty;
            obj.BinNodeType(obj.BinNodeJunctionIndex)=obj.TYPENODE(1);
            obj.BinNodeType(obj.BinNodeReservoirIndex)=obj.TYPENODE(2);
            obj.BinNodeType(obj.BinNodeTankIndex)=obj.TYPENODE(3);

            obj.BinLinkType(obj.BinLinkPipeIndex)=obj.TYPELINK(2);
            obj.BinLinkType(obj.BinLinkPumpIndex)=obj.TYPELINK(3);
            obj.BinLinkType(obj.BinLinkValveIndex)=obj.BinLinkValveType;

            obj.BinLinkSettings = [obj.BinLinkPipeRoughness zeros(1, obj.BinLinkPumpCount) obj.BinLinkValveSetting]';
            obj.BinNodeElevations = single([obj.BinNodeJunctionElevation obj.BinNodeReservoirElevation obj.BinNodeTankElevation]);
            obj.BinLinkDiameters = single([obj.BinLinkPipeDiameters zeros(1, obj.BinLinkPumpCount) obj.BinLinkValveDiameters]);
            obj.BinLinkLengths = single([obj.BinLinkPipeLengths zeros(1, obj.BinLinkPumpCount) zeros(1, obj.BinLinkValveCount)]);
            obj.BinLinkRoughnessCoeff = [obj.BinLinkPipeRoughness zeros(1, obj.BinLinkPumpCount) zeros(1, obj.BinLinkValveCount)];
            % obj.BinNodeJunctionsBaseDemands(length(obj.BinNodeJunctionsBaseDemands):obj.BinNodeJunctionCount)=0;
            obj.BinNodeBaseDemands = single([obj.BinNodeJunctionsBaseDemands zeros(1, obj.BinNodeReservoirCount) zeros(1, obj.BinNodeTankCount)]);
            % obj.BinNodeJunDemandPatternNameID=[obj.BinNodeJunDemandPatternNameID obj.BinNodeResDemandPatternNameID];
            % for i=obj.BinNodeTankIndex
            %    obj.BinNodeDemandPatternNameID{i}='';
            %      end

            b={};
            for i=1:obj.BinLinkCount
                ind=find((strcmp(obj.BinLinkInitialStatusNameID, obj.BinLinkNameID{i}))==1);
                if isempty(ind), ind=0; end
                if ind~=0
                    if sum(obj.BinLinkPumpIndex==i)
                        bb{i}=obj.BinLinkInitialStatus{ind};
                        r=obj.BinLinkNameID(i);
                        b{i}=r{1};
                    else
                        b{i}=obj.BinLinkNameID{i};
                        bb{i}=obj.BinLinkInitialStatus{ind};
                    end
                else
                    b{i}=obj.BinLinkNameID{i};
                    bb{i}='Open';
                end
            end
            obj.BinLinkInitialStatusNameID=b;
            obj.BinLinkInitialStatus=bb;
            if ~isempty(obj.BinLinkPumpIndex)
                obj.BinLinkPumpStatus=obj.BinLinkInitialStatus(obj.BinLinkPumpIndex);
                obj.BinLinkPumpStatusNameID=obj.BinLinkInitialStatusNameID(obj.BinLinkPumpIndex);
            end
            if ~isempty(obj.BinLinkPumpIndex)
                obj.BinLinkValveStatus=obj.BinLinkInitialStatus(obj.BinLinkValveIndex);
                obj.BinLinkValveStatusNameID=obj.BinLinkInitialStatusNameID(obj.BinLinkValveIndex);
            end

            if ~isempty(BinCNameID)
                j=1;
                for i=1:length(BinCNameID)
                    if ~isempty(BinCNameID{i})
                        nn(j)=BinCNameID(i);j=j+1;
                    end
                end
                obj.BinCurveNameID=unique(nn);
                obj.BinCurveCount=length(obj.BinCurveNameID);
            end
            obj.BinControlRulesCount=length(obj.BinControlsInfo);
            obj.BinUnits=getBinUnits(obj);
            % US Customary - SI metric
            if find(strcmp(obj.BinLinkFlowUnits, obj.TYPEUNITS))<6
                obj.BinUnits_US_Customary=1;
            else
                obj.BinUnits_SI_Metric=1;
            end
        end
        function [Errcode]=setBinNodeInitialQuality(obj, varargin)
            parameter=varargin{1};
            zz=obj.BinNodeCount-obj.BinCountInitialQualitylines+1;
            sections={'[QUALITY]', '[SOURCES]'};
            [Errcode]=setBinParam2(obj, parameter, sections, zz);
        end
        function [Errcode]=setBinLinkReactionCoeff(obj, varargin)
            wall=[];Errcode=0;
            bulk=[];
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'wall'
                        wall=varargin{2*i};
                    case 'bulk'
                        bulk=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
           links=obj.getBinLinkNameID;
           BinLinkCount=length(links.BinLinkNameID);
           [tlines]=regexp( fileread(obj.BinTempfile), '\n', 'split');
           fid = writenewTemp(obj.BinTempfile);start=0;
           for i=1:length(tlines)
               tt=regexp(tlines{i}, '\s*', 'split');
               tok = strtok(tlines{i});m=1;
               % Skip blank Clines and comments
               if isempty(tok), continue;
               elseif isempty(tt{m})
                   m=m+1;
               end
               if strcmpi(tt{m}, 'GLOBAL') && strcmpi(tt{m+1}, 'WALL')
                   start=i;
               end
               if strcmp(tt{m}, '[MIXING]')
                   stop1=i;
               end
               if strcmp(tt{m}, '[ENERGY]')
                   stop2=i;
               end
           end
           stop=max([stop1 stop2]);
           for kk=start+1:stop-1
              tlines{kk}='';
           end
           for u=1:start
               nnlines{u}=tlines{u};
           end
           for u=start+1:start+1+BinLinkCount*2
               nnlines{u}=[];
           end
           for k=start+1:length(tlines)
               nnlines{u}=tlines{k};
               u=u+1;
           end
           tlines=nnlines;clear nnlines;
           for i=start:stop
               % Get first token in the line
               tok = strtok(tlines{i});
               if isempty(tok), continue; end
               if strcmp(tok(1), ';')
               else
                   if ~isempty(wall)
                       for e=1:BinLinkCount
                           clear atlines;
                           atlines{1} = 'WALL';
                           atlines{2} = links.BinLinkNameID{e};
                           atlines{3} = num2str(wall(e));
                           newlines=[];
                           for pp=1:length(atlines)
                               newlines = [newlines, atlines{pp}, blanks(12)];
                           end
                           tlines{i+e}=newlines;
                       end
                   end
                   if ~isempty(bulk)
                       for e=1:obj.BinLinkCount
                           clear atlines;
                           atlines{1} = 'BULK';
                           atlines{2} = links.BinLinkNameID{e};
                           atlines{3} = num2str(bulk(e));
                           newlines=[];
                           for pp=1:length(atlines)
                               newlines = [newlines, atlines{pp}, blanks(12)];
                           end
                           tlines{i+BinLinkCount+e}=newlines;
                       end
                   end
               end
               break;
           end
           fprintf(fid, '%s\n', tlines{:});
           fclose(fid);
           if obj.Bin==1
               Errcode=reloadNetwork(obj);
           end
        end
        function [Errcode]=setBinQualType(obj, chemname, chemunits, varargin)
            sections={'[OPTIONS]', '[REPORT]'};
            indexParameter=1;
            if nargin<3
                chemunits='mg/L';
            end
            parameter=['Quality', blanks(5), chemname, blanks(5), chemunits];
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinQualityChem(obj, varargin)
            sections={'[OPTIONS]', '[REPORT]'};
            indexParameter=1;
            parameter='Quality            	chem   mg/L';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinQualityNone(obj, varargin)
            sections={'[OPTIONS]', '[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	None';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinQualityAge(obj, varargin)
            sections={'[OPTIONS]', '[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	Age';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinQualityTrace(obj, varargin)
            sections={'[OPTIONS]', '[COORDINATES]'};
            indexParameter=1;
            if ~sum(strcmp(varargin{1}, obj.BinNodeNameID))
                warning('Invalid property found.');
                return
            end
            parameter=['Quality            	Trace ', varargin{1}];
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeSimulationDuration(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=1;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeHydraulicStep(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeQualityStep(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimePatternStep(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimePatternStart(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=5;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeReportingStep(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=6;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeReportingStart(obj, varargin)
            parameter=varargin{1};
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=7;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeStatisticsNone(obj, varargin)
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=8;
            parameter='None';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeStatisticsAverage(obj, varargin)
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=8;
            parameter='AVERAGE';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeStatisticsMinimum(obj, varargin)
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=8;
            parameter='MINIMUM';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeStatisticsMaximum(obj, varargin)
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=8;
            parameter='MAXIMUM';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinTimeStatisticsRange(obj, varargin)
            sections={'[TIMES]', '[REPORT]'};
            indexParameter=8;
            parameter='RANGE';
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankElevation(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankInitialLevel(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankMinimumWaterLevel(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankMaximumWaterLevel(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=5;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankDiameter(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=6;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeTankMinimumWaterVolume(obj, varargin)
            parameter=varargin{1};
            sections={'[TANKS]', '[PIPES]'};
            indexParameter=7;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [value] = getBinLimitingPotential(obj)
            [~, value] = limitingPotential(obj, 'get');
        end
        function [Errcode]=setBinLimitingPotential(obj, newlimiting)
            Errcode = limitingPotential(obj, '', newlimiting);
        end
        function [Errcode]=setBinLinkGlobalWallReactionCoeff(obj, varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]', '[MIXING]'};
            [Errcode]=setBinParam(obj, 3, parameter, sections);
        end
        function [Errcode]=setBinLinkGlobalBulkReactionCoeff(obj, varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]', '[MIXING]'};
            [Errcode]=setBinParam(obj, 1, parameter, sections);
        end
        function BinClose(obj)
            fclose all;
            if exist([obj.BinTempfile(1:end-4), '.bin'])==2
                delete([obj.BinTempfile(1:end-4), '.bin'])
            end
            if ~libisloaded(obj.LibEPANET)
                delete(obj.BinTempfile);
            end
            if exist([obj.BinTempfile(1:end-4), '.msx'])==2
                delete([obj.BinTempfile(1:end-4), '.msx'])
            end
        end
        function [Errcode]=setBinLinkValvesParameters(obj, varargin)
            % Initiality
            Diameter=[];
            Type={};
            Setting=[];
            MinorLoss=[];
            Status={};
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'diameter'
                        Diameter=varargin{2*i};
                    case 'type'
                        Type=varargin{2*i};
                        for u=1:length(Type)
                            aa(u)=sum(strcmp(Type{u}, obj.TYPELINK));
                        end
                        if ~sum(aa)
                            warning('Invalid property found.');Errcode=-1;
                            return;
                        end
                    case 'setting'
                        Setting=varargin{2*i};
                    case 'minorloss'
                        MinorLoss=varargin{2*i};
                    case 'status'
                        if sum(strcmpi(varargin{2*i}, 'closed')+strcmpi(varargin{2*i}, 'open')+strcmpi(varargin{2*i}, 'none')+strcmpi(varargin{2*i}, 'nonestatus'))==obj.LinkValveCount
                            Status=varargin{2*i};
                        else
                            warning('Invalid argument found.');Errcode=-1;
                            return;
                        end
                        zz=abs(obj.BinLinkPumpCount+obj.BinLinkValveCount-obj.BinCountStatuslines);
                        sections={'[STATUS]', '[PATTERNS]', 'valve'};
                        setBinParam2(obj, Status, sections, zz);
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
            [tlines]=regexp( fileread(obj.BinTempfile), '\n', 'split');
            fid = writenewTemp(obj.BinTempfile);
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1}, '[VALVES]')
                    start=i;
                end
                if strcmp(tt{1}, '[DEMANDS]')
                    stop=i;
                end
            end
           if ~isempty(Status)
               for kk=start+1:stop-1
                  tlines{kk}='';
               end
               for u=1:start
                   nnlines{u}=tlines{u};
               end
               for u=start+1:start+1+zz
                   nnlines{u}=[];
               end
               for k=start+1:length(tlines)
                   nnlines{u}=tlines{k};
                   u=u+1;
               end
           end
            ll=1;clear atlines;
            if start
               for i=start:stop
                   % Get first token in the line
                   tok = strtok(tlines{i});
                   if isempty(tok), tok='1'; end
                   if strcmp(tok(1), ';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       clear atlines;

                       if ll<obj.BinLinkValveCount+1
                           atlines{1}=obj.BinLinkNameID{obj.BinLinkValveIndex(ll)};
                           atlines{2}=obj.BinLinkFromNode{obj.BinLinkValveIndex(ll)};
                           atlines{3}=obj.BinLinkToNode{obj.BinLinkValveIndex(ll)};
                           if ~isempty(Diameter)%Diameters
                               atlines{4} = num2str(Diameter(ll));
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(Type)%Type
                               atlines{5} = num2str(Type{ll});
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(Setting)%Setting
                               atlines{6} = num2str(Setting(ll));
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(MinorLoss)%MinorLoss
                               atlines{7} = num2str(MinorLoss(ll));
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       ll=ll+1;
                   end
               end
            end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid);
            if obj.Bin==1
                Errcode=reloadNetwork(obj);
            end
        end
        function [Errcode]=setBinNodeResDemandPatternNameID(obj, varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]', '[TANKS]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeReservoirElevation(obj, varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]', '[TANKS]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeJunctionElevation(obj, varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]', '[RESERVOIRS]', '[DEMANDS]', '[STATUS]', '[EMITTERS]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeJunctionsBaseDemands(obj, varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]', '[RESERVOIRS]', '[DEMANDS]', '[STATUS]', '[EMITTERS]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinNodeJunDemandPatternNameID(obj, varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]', '[RESERVOIRS]', '[DEMANDS]', '[STATUS]', '[EMITTERS]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinPattern(obj, varargin)
            idpattern=varargin{1};
            %ind=obj.getPatternIndex(idpattern);
            values=varargin{2}; %obj.getPatternValue{ind};
            sections={'[PATTERNS]', '[CURVES]'};
            v=obj.getBinPatternsInfo;
            patterns=v.BinPatternNameID;
            if sum(strcmp(idpattern, patterns))
                Errcode=setBinParam(obj, idpattern, values, sections);
            else
                warning('Invalid argument found.');Errcode=-1;
                return
            end
        end
        function [Errcode]=addBinPattern(obj, varargin)
            newidpattern=varargin{1};
            values=varargin{2};
            sections={'[PATTERNS]', '[CURVES]'};
            v=obj.getBinPatternsInfo;
            patterns=v.BinPatternNameID;zz=0;
            if ~sum(strcmp(newidpattern, patterns))
                for i=1:length(patterns)
                    if mod(length(v.BinPatternValue{i}), 6)==0
                        zz=zz+length(v.BinPatternValue{i})/6;
                    else
                        zz=zz+1;
                        if mod(length(values), 6)
                            zz=zz+1;
                        end
                    end
                end
                Errcode=setBinParam2(obj, values, sections, zz, newidpattern);
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
            end
        end
        function [Errcode]=setBinNodeSourceQuality(obj, varargin)
            sections={'[SOURCES]', '[MIXING]'};
            values=varargin{1};
            [Errcode]=setBinParam(obj, 11, values, sections);
        end
        function saveBinInpFile(obj, varargin)
            if ~isempty(varargin)
                copyfile(obj.BinTempfile, varargin{1}); return;
            end
            [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
            f = writenewTemp(obj.BinTempfile);
            % /*Write [TITLE] section */
               for i=1:length(tlines)
                   tok = strtok(tlines{i});
                   if sum(tlines{i}=='[') && ~strcmp(tok, '[TITLE]')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [JUNCTIONS] section */
               fprintf(f, '\n[JUNCTIONS]');
               sps=blanks(10);
               for u=obj.BinNodeJunctionIndex
                   fprintf(f, '\n%s%s%f', obj.BinNodeNameID{u}, sps, obj.BinNodeElevations(u));
               end
            % % /*Write [RESERVOIRS] section */
               fprintf(f, '\n[RESERVOIRS]');
               for u=obj.BinNodeReservoirIndex
                   fprintf(f, '\n%s%s%f', obj.BinNodeNameID{u}, sps, obj.BinNodeElevations(u));
               end
            % % /*Write [TANKS] section */
               fprintf(f, '\n[TANKS]');b=1;
               for u=obj.BinNodeTankIndex
               %  InitLevel   	MinLevel    	MaxLevel    	Diameter    	MinVol      	VolCurve
                   fprintf(f, '\n%s%s%f%s%f%s%f%s%f%s%f%s%f', obj.BinNodeNameID{u}, sps, obj.BinNodeElevations(u), sps, obj.BinNodeTankInitialLevel(b), sps, obj.BinNodeTankMinimumWaterLevel(b), ...
                       sps, obj.BinNodeTankMaximumWaterLevel(b), sps, obj.BinNodeTankDiameter(b), sps, obj.BinNodeTankMinimumWaterVolume(b));
                   b=b+1;
               end
            % % /*Write [PIPES] section */
               fprintf(f, '\n[PIPES]');
               for u=obj.BinLinkPipeIndex
               % ;ID              	Node1           	Node2           	Length      	Diameter    	Roughness   	MinorLoss   	Status
                   fprintf(f, '\n%s%s%s%s%s%s%f%s%f%s%f%s%f%s%s', obj.BinLinkNameID{u}, sps, obj.BinLinkFromNode{u}, sps, obj.BinLinkToNode{u}, sps, obj.BinLinkLengths(u), sps, obj.BinLinkDiameters(u), ...
                       sps, obj.BinLinkPipeRoughness(u), sps, obj.BinLinkPipeMinorLoss(u), sps, obj.BinLinkInitialStatus{u});
               end
            % % /*Write [PUMPS] section */
               fprintf(f, '\n[PUMPS]');
               par={';ID', ';Node', ';Junction', ';Demand', ';Type', ';Tank', ';Link'};
               for pp=1:length(par)
                   if find(strcmp(strtok(tlines), par{pp}))
                       for i=find(strcmp(strtok(tlines), par{pp}))
                          tlines{i}='';
                       end
                   end
               end
               for i=find(strcmp(strtok(tlines), '[PUMPS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [VALVES] section */
               fprintf(f, '\n[VALVES]');
               for i=find(strcmp(strtok(tlines), '[VALVES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [DEMANDS] section */
               fprintf(f, '\n[DEMANDS]');
               for u=1:obj.BinNodeJunctionCount
                   if isempty(obj.BinNodeJunDemandPatternNameID{u})
                       fprintf(f, '\n%s%s%f%s%s', obj.BinNodeNameID{u}, sps, obj.BinNodeBaseDemands(u), sps, obj.BinPatternNameID{1});
                   else
                       fprintf(f, '\n%s%s%f%s%s', obj.BinNodeNameID{u}, sps, obj.BinNodeBaseDemands(u), sps, obj.BinNodeJunDemandPatternNameID{u});
                   end
               end
               % % /*Write [EMITTERS] section */
               fprintf(f, '\n[EMITTERS]');
               for i=find(strcmp(strtok(tlines), '[EMITTERS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [STATUS] section */
                fprintf(f, '\n[STATUS]');
                for i=find(strcmp(strtok(tlines), '[STATUS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [PATTERNS] section */
                fprintf(f, '\n[PATTERNS]');
                for i=find(strcmp(strtok(tlines), '[PATTERNS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [CURVES] section */
                fprintf(f, '\n[CURVES]');
                for i=find(strcmp(strtok(tlines), '[CURVES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [CONTROLS] section */
                fprintf(f, '\n[CONTROLS]');
                for i=find(strcmp(strtok(tlines), '[CONTROLS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                      break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [QUALITY] section */
                fprintf(f, '\n[QUALITY]');
                for i=find(strcmp(strtok(tlines), '[QUALITY]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
                end
            % % /*Write [SOURCES] section */
               fprintf(f, '\n[SOURCES]');
               for i=find(strcmp(strtok(tlines), '[SOURCES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [MIXING] section */
               fprintf(f, '\n[MIXING]');
               for i=find(strcmp(strtok(tlines), '[MIXING]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   for u=1:obj.BinNodeTankCount
                       if isempty(obj.BinNodeTankMixModel)
                           obj.BinNodeTankMixModel{u}='MIXED';
                           obj.BinNodeTankMinimumFraction(u)=0;
                           obj.BinNodeTankMixID{u}=obj.BinNodeTankNameID{u};
                           fprintf(f, '\n%s%s%s%s%f', obj.BinNodeTankNameID{u}, sps, obj.BinNodeTankMixModel{u}, sps, obj.BinNodeTankMinimumFraction);
                       end
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [REACTIONS] section */
               fprintf(f, '\n[REACTIONS]');
               ff=find(strcmp(strtok(tlines), '[REACTIONS]'));
               for u=1:length(ff)
                   for i=ff(u)+1:length(tlines)
                       if sum(tlines{i}=='[')
                           break;
                       end
                       fprintf(f, '\n%s', tlines{i});
                   end
               end
            % % /*Write [ENERGY] section */
               fprintf(f, '\n[ENERGY]');
               for i=find(strcmp(strtok(tlines), '[ENERGY]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [TIMES] section */
               fprintf(f, '\n[TIMES]');
               for i=find(strcmp(strtok(tlines), '[TIMES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [OPTIONS] section */
            fprintf(f, '\n[OPTIONS]');
               for i=find(strcmp(strtok(tlines), '[OPTIONS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [REPORT] section */
            fprintf(f, '\n[REPORT]');
               for i=find(strcmp(strtok(tlines), '[REPORT]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [TAGS] section */
            fprintf(f, '\n[TAGS]');
               for i=find(strcmp(strtok(tlines), '[TAGS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [RULES] section */
            fprintf(f, '\n[RULES]');
               for i=find(strcmp(strtok(tlines), '[RULES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [COORDINATES] section */
            fprintf(f, '\n[COORDINATES]');
               for i=find(strcmp(strtok(tlines), '[COORDINATES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [VERTICES] section */
            fprintf(f, '\n[VERTICES]');
               for i=find(strcmp(strtok(tlines), '[VERTICES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [LABELS] section */
            fprintf(f, '\n[LABELS]');
               for i=find(strcmp(strtok(tlines), '[LABELS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
            % % /*Write [BACKDROP] section */
            fprintf(f, '\n[BACKDROP]');
               for i=find(strcmp(strtok(tlines), '[BACKDROP]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f, '\n%s', tlines{i});
               end
               fprintf(f, '\n[END]');
        end
        function [node_index, link_index] = addBinNodeJunction(obj, nodeID, varargin)
            % Adds a new junction with a pipe/pump/valve to the network.
            %
            % Example 1: Add junction with pipe
            %   d=epanet('NET1.inp');
            %   nodeID = 'new';
            %   coordinates = [20, 50];
            %   elevation = 500;
            %   demand = 100;
            %   patternID = '1';
            %   patternCategoryID = '1';
            %   quality = 0.5;
            %   linkID = 'new_pipe';
            %   from = nodeID;
            %   to = '22';
            %   length = 6000;     % Feet
            %   diameter = 20;     % Inches
            %   roughness = 100;   % Darcy-Weisbach formula(millifeet)
            %   minorLoss = 0;
            %   status = 'Open';
            %
            %   [node_index, link_index] = d.addBinNodeJunction(nodeID, coordinates, elevation, demand, patternID, patternCategoryID, quality,...
            %     {'PIPE', linkID, from, to, length, diameter, roughness, minorLoss, status});
            %
            %
            % See also addBinNodeReservoir, addBinNodeTank, addBinPipe,
            %          addBinPump, getBinNodeIndex, getBinNodesInfo.
            if obj.Bin, obj.Errcode = obj.saveInputFile(obj.BinTempfile); end
            cnt = length(nodeID);
            coords = zeros(cnt, 2);
            elev = zeros(cnt, 1);
            demand = zeros(cnt, 1);
            patternID = repmat({''}, 1, cnt);
            category = repmat({''}, 1, cnt);
            quality = zeros(cnt, 1);
            if nargin >= 3
                coords = varargin{1};
            end
            if nargin >= 4
                elev = varargin{2};
            end
            if nargin > 5
                demand = varargin{3};
            end
            if nargin >= 6
                patternID = varargin{4};
                if ~iscell(patternID)
                    patternID = {patternID};
                end
            end
            if nargin >= 7
                category = varargin{5};
                if ~iscell(category)
                    category = {category};
                end
            end
            if nargin >= 8
                quality = varargin{6};
            end
            node_index = addBinNode(obj, 1, nodeID, coords, elev, demand, patternID, category, quality);
            if nargin == 9
                if strcmp(varargin{7}{1}, 'PIPE')
                    link_index = addBinLinkPipe(obj, varargin{7}{2:end});
                elseif strcmp(varargin{7}{1}, 'PUMP')
                    link_index = addBinLinkPump(obj, varargin{7}{2:end});
                else
                    link_index = addBinLinkValve(obj, varargin{7}{2:end});
                end
            end
        end
        function [node_index, link_index] = addBinNodeReservoir(obj, nodeID, varargin)
            % Adds a new reservoir to the network.
            %
            % Example 1:
            %   d=epanet('NET1.inp');
            %   nodeID = 'new';
            %   coordinates = [20, 50];
            %   head = 800;
            %   patternID = '1';
            %   quality = 1;
            %
            % % If properties are not given, the default values are zero.
            %
            % % Adds a new reservoir with the default coordinates (i.e. [0, 0])
            %   node_index = d.addBinNodeReservoir(nodeID)
            %   d.plot;
            %
            % % Adds a new reservoir with coordinates [X, Y] = [20, 50].
            %   node_index = d.addBinNodeReservoir(nodeID, coordinates)
            %   d.plot;
            %   x_value = d.getBinNodeCoordinates{1}(node_index)
            %   y_value = d.getBinNodeCoordinates{2}(node_index)
            %
            % % Adds a new reservoir with coordinates [X, Y] = [20, 50] and head = 500.
            %   node_index = d.addBinNodeReservoir(nodeID, coordinates, head)
            %   d.plot;
            %
            % % Adds a new reservoir with coordinates [X, Y] = [20, 50], head = 800 and pattern ID = '1'.
            %   node_index = d.addBinNodeReservoir(nodeID, coordinates, head, patternID)
            %   d.plot;
            %
            % % Adds a new reservoir with coordinates [X, Y] = [20, 50], head = 800, pattern ID = '1' and quality = 1.
            %   node_index = d.addBinNodeReservoir(nodeID, coordinates, head, patternID, quality)
            %   d.plot;
            %
            % Example 2:
            %   d=epanet('NET1.inp');
            %   nodeID = {'new_1', 'new_2'};
            %   coordinates = [20, 50; 20, 40];
            %   head = [800, 900];
            %   patternID = {'1', '1'};
            %   quality = [1, 1];
            %
            % % Adds 2 new reservoirs with coordinates [X, Y] = [20, 50] and [20, 40], head = 800 and 900, pattern ID = '1' and '1', and quality = 1 and 1.
            %   node_index = d.addBinNodeReservoir(nodeID, coordinates, head, patternID, quality)
            %   d.plot;
            %
            % See also addBinNodeJunction, addBinNodeTank, addBinPipe,
            %          addBinPump, getBinNodeIndex, getBinNodesInfo.
            coords = zeros(length(nodeID), 2);
            head = zeros(length(nodeID), 1);
            patternID = repmat({''}, 1, length(nodeID));
            quality = zeros(length(nodeID), 1);
            if nargin >= 3
                coords = varargin{1};
            end
            if nargin >= 4
                head = varargin{2};
            end
            if nargin >= 5
                patternID = varargin{3};
                if ~iscell(patternID)
                    patternID = {patternID};
                end
            end
            if nargin >= 6
                quality = varargin{4};
            end
            node_index = addBinNode(obj, 2, nodeID, coords, head, patternID, quality);
            if nargin == 7
                if strcmp(varargin{5}{1}, 'PIPE')
                    link_index = addBinLinkPipe(obj, varargin{5}{2:end});
                elseif strcmp(varargin{5}{1}, 'PUMP')
                    link_index = addBinLinkPump(obj, varargin{5}{2:end});
                else
                    link_index = addBinLinkValve(obj, varargin{5}{2:end});
                end
             end
        end
        function [node_index, link_index] = addBinNodeTank(obj, nodeID, varargin)
            % Adds a new tank to the network.
            %
            % Example 1:
            %   d=epanet('NET1.inp');
            %   nodeID = 'new';
            %   coordinates = [20, 50];
            %   elev = 800;
            %   diameter = 50;
            %   initial_level = 120;
            %   min_level = 100;
            %   max_level = 150;
            %   min_volume = 0;
            %   volume_curve = '';
            %   quality = 1;
            %
            % % If properties are not given, the default values are diameter = 50, initlevel = 10, maxlevel = 20 and the remaining are set to zero.
            %
            % % Adds a new tank with the default coordinates (i.e. [0, 0])
            %   node_index = d.addBinNodeTank(nodeID)
            %   d.plot;
            %
            % % Adds a new tank with coordinates [X, Y] = [20, 50].
            %   node_index = d.addBinNodeTank(nodeID, coordinates)
            %   d.plot;
            %   x_value = d.getBinNodeCoordinates{1}(node_index)
            %   y_value = d.getBinNodeCoordinates{2}(node_index)
            %
            % % Adds a new tank with coordinates [X, Y] = [20, 50], elevation = 800 and diameter = 50.
            %   node_index = d.addBinNodeTank(nodeID, coordinates, elev, diameter)
            %   d.plot;
            %   d.getBinNodesInfo.BinNodeTankElevation
            %   d.getBinNodesInfo.BinNodeTankDiameter
            %
            % % Adds a new tank with coordinates [X, Y] = [20, 50], elevation = 800, diameter = 50, initial_level = 120, minimum level = 100 and maximum level = 150.
            %   node_index = d.addBinNodeTank(nodeID, coordinates, elev, diameter, initial_level, min_level, max_level)
            %   d.plot;
            %   d.getBinNodesInfo.BinNodeTankInitialLevel
            %   d.getBinNodesInfo.BinNodeTankMinimumWaterLevel
            %   d.getBinNodesInfo.BinNodeTankMaximumWaterLevel
            %
            % % Adds a new tank with coordinates [X, Y] = [20, 50], elevation = 800, diameter = 50, initial_level = 120,
            % %  minimum level = 100, maximum level = 150, minimum volume = 0 and volume curve = '1'.
            %   node_index = d.addBinNodeTank(nodeID, coordinates, elev, diameter, initial_level, min_level, max_level, min_volume, volume_curve)
            %   d.plot;
            %   d.getBinNodesInfo.BinNodeTankMinimumWaterVolume
            %
            % % Adds a new tank with coordinates [X, Y] = [20, 50],  elevation = 800, diameter = 50, initial_level = 120,
            % %  minimum level = 100, maximum level = 150, minimum volume = 0, volume curve = '1' and quality = 1.
            %   node_index = d.addBinNodeTank(nodeID, coordinates, elev, diameter, initial_level, min_level, max_level, min_volume, volume_curve, quality)
            %   d.plot;
            %
            % Example 2:
            %   d=epanet('NET1.inp');
            %   nodeID = {'new_1', 'new_2'};
            %   coordinates = [20, 50; 20, 40];
            %   elev = [800, 900];
            %   diameter = [50, 60];
            %   initial_level = [120, 110];
            %   min_level = [100, 90];
            %   max_level = [150, 160];
            %   min_volume = [0, 0];
            %   volume_curve = {'', ''};
            %   quality = [1, 0.5];
            %
            % % Adds 2 new tanks with coordinates [X, Y] = [20, 50] and [20, 40],  elevation = 800 and 900, diameter = 50 and 60, initial_level = 120 and 110,
            % %  minimum level = 100 and 90, maximum level = 150 and 160, minimum volume = 0 and 0, volume curve = none and quality = 1 and 0.5.
            %   node_index = d.addBinNodeTank(nodeID, coordinates, elev, diameter, initial_level, min_level, max_level, min_volume, volume_curve, quality)
            %   d.plot;
            %
            % See also addBinNodeJunction, addBinNodeReservoir, addBinPipe,
            %          addBinPump, getBinNodeIndex, getBinNodesInfo.
            coords = zeros(length(nodeID), 2);
            elev = zeros(length(nodeID), 1);
            diameter = zeros(length(nodeID), 1);
            diameter(:) = 50;
            initlevel = zeros(length(nodeID), 1);
            initlevel(:) = 10;
            minlevel = zeros(length(nodeID), 1);
            minlevel(:) = 0;
            maxlevel = zeros(length(nodeID), 1);
            maxlevel(:) = 20;
            minvol = zeros(length(nodeID), 1);
            volcurve = repmat({''}, 1, length(nodeID));
            quality = zeros(length(nodeID), 1);
            if nargin >= 3
                coords = varargin{1};
            end
            if nargin >= 4
                elev = varargin{2};
            end
            if nargin >= 5
                diameter = varargin{3};
            end
            if nargin >= 6
                initlevel = varargin{4};
            end
            if nargin >= 7
                minlevel = varargin{5};
            end
            if nargin >= 8
                maxlevel = varargin{6};
            end
            if nargin >= 9
                minvol = varargin{7};
            end
            if nargin >= 10
                volcurve = varargin{8};
                if ~iscell(volcurve)
                    volcurve = {volcurve};
                end
                for i = 1:length(volcurve)
                    if ~isempty(volcurve{i})
                        if ~sum(strcmp(volcurve{i}, obj.getBinCurvesInfo.BinCurveNameID))
                            warning(['Curve ', volcurve{i}, ' does not exist.'])
                            node_index=-1;
                            return;
                        end
                    end
                end
            end
            if nargin >= 11
                quality = varargin{9};
            end
            node_index = addBinNode(obj, 3, nodeID, coords, elev, diameter, initlevel, minlevel, maxlevel, minvol, volcurve, quality);
            if nargin == 12
                if strcmp(varargin{10}{1}, 'PIPE')
                    link_index = addBinLinkPipe(obj, varargin{10}{2:end});
                elseif strcmp(varargin{10}{1}, 'PUMP')
                    link_index = addBinLinkPump(obj, varargin{10}{2:end});
                else
                    link_index = addBinLinkValve(obj, varargin{10}{2:end});
                end
            end
        end
        function [link_index] = addBinLinkPipe(obj, linkID, from, to, varargin)
            % Adds a new pipe to the nework.
            %
            % The default values are:
            %   length = 330;      % Feet
            %   diameter = 10;     % Inches
            %   roughness = 0.5;   % Darcy-Weisbach formula(millifeet)
            %   minorLoss = 0;
            %   status = 'Open';
            %
            % Example 1:
            %   d=epanet('NET1.inp');
            %   linkID = 'new_pipe';
            %   from = '11';
            %   to = '22';
            %   length = 6000;     % Feet
            %   diameter = 20;     % Inches
            %   roughness = 100;   % Darcy-Weisbach formula(millifeet)
            %   minorLoss = 0;
            %   status = 'Open';
            %
            % % Adds a new pipe from junction '11' to '22' with the default values.
            %   link_index = d.addBinLinkPipe(linkID, from, to)
            %   d.plot;
            %
            % % Adds a new pipe from junction '11' to '22' with length = 6000.
            %   link_index = d.addBinLinkPipe(linkID, from, to, length)
            %   d.plot;
            %
            % % Adds a new pipe from junction '11' to '22' with length = 6000 and diameter = 20.
            %   link_index = d.addBinLinkPipe(linkID, from, to, length, diameter)
            %   d.plot;
            %
            % % Adds a new pipe from junction '11' to '22' with length = 6000, diameter = 20, roughness = 100, minorLoss = 0 and status = Open.
            %   link_index = d.addBinLinkPipe(linkID, from, to, length, diameter, roughness, minorLoss, status)
            %   d.plot;
            %
            % Example 2:
            %   d=epanet('NET1.inp');
            %   linkID = {'new_pipe_1', 'new_pipe_2'};
            %   from = {'11', '13'};
            %   to = {'22', '22'};
            %   length = [6000, 7000];    % Feet
            %   diameter = [20, 23];      % Inches
            %   roughness = [100, 105];   % Darcy-Weisbach formula(millifeet)
            %   minorLoss = [0, 0.2];
            %   status = {'Open', 'Closed'};
            %
            % % Adds 2 new pipes from junctions '11'and '13' to '22' with length = 6000 and 7000, diameter = 20 and 23,
            % % roughness = 100 and 105, minorLoss = 0 and 0.2, and status = Open and Closed.
            %   link_index = d.addBinLinkPipe(linkID, from, to, length, diameter, roughness, minorLoss, status)
            %   d.plot;
            %
            % See also addBinLinkPump, addBinLinkValve, getBinLinksInfo,
            %          addBinNodeJunction, addBinNodeReservoir, addBinNodeTank.
            lengthp = zeros(length(linkID), 1);
            lengthp(:) = 330;
            diameter = zeros(length(linkID), 1);
            diameter(:) = 10;
            roughness = zeros(length(linkID), 1);
            roughness(:) = 0.5;
            minorloss = zeros(length(linkID), 1);
            minorloss(:) = 0;
            status = repmat({'Open'}, 1, length(linkID));
            if nargin >= 5
                lengthp = varargin{1};
            end
            if nargin >= 6
                diameter = varargin{2};
            end
            if nargin >= 7
                roughness = varargin{3};
            end
            if nargin >= 8
                minorloss = varargin{4};
            end
            if nargin >= 9
                status = varargin{5};
            end
            link_index = addBinLink(obj, 'PIPE', linkID, from, to, lengthp, diameter, roughness, minorloss, status);
        end
        function [link_index] = addBinLinkPump(obj, linkID, from, to, varargin)
            % Adds a new pump to the network.
            %
            % If no propertie is given, the default value is 'SPEED 1.0'.
            %
            % Example 1:
            %   d=epanet('NET1.inp');
            %   linkID = 'new_pump';
            %   from = '11';
            %   to = '22';
            %   propertie = 'HEAD 1';
            %
            % % Adds a new pump from junction '11' to '22' with the default value of propertie(i.e. 'SPEED 1.0').
            %   linkIndex = d.addBinLinkPump(linkID, from, to)
            %   d.plot;
            %
            % % Adds a new pump from junction '11' to '22' with propertie = 'HEAD 1'.
            %   linkIndex = d.addBinLinkPump(linkID, from, to, propertie)
            %   d.plot;
            %
            % Example 2:
            %   d=epanet('NET1.inp');
            %   linkID = {'new_pump_1', 'new_pump_2'};
            %   from = {'11', '13'};
            %   to = {'22', '22'};
            %   properties = {'HEAD 1', 'SPEED 1.2'};
            %
            % % Adds 2 new pumps from junctions '11' and '13' to '22' with properties 'HEAD 1' and 'SPEED 1.2'.
            %   linkIndex = d.addBinLinkPump(linkID, from, to, properties)
            %   d.plot;
            %
            % See also addBinLinkPipe, addBinLinkValve, getBinLinksInfo,
            %          addBinNodeJunction, addBinNodeReservoir, addBinNodeTank.
            propertie = repmat({'SPEED 1.0'}, 1, length(linkID));
            if nargin == 5
                propertie = varargin{1};
            end
            link_index = addBinLink(obj, 'PUMP', linkID, from, to, propertie);
        end
        function [link_index] = addBinLinkValve(obj, linkID, from, to, varargin)
            % Adds a new valve to the network.
            %
            % % If no properties are given, the default values are:
            %   type = 'GPV'
            %   diameter = 10 inches (25.4 cm)
            %   initial setting = 0
            %   minor Loss Coefficient = 0
            %
            % Example 1:
            %   d=epanet('NET1.inp');
            %   linkID = 'new_valve';
            %   from = '11';
            %   to = '22';
            %   type = 'PRV';
            %   diameter = 6;
            %   init_setting = 70;
            %   minorLoss = 0;
            %
            % % Adds a new valve form junction '11' to '22' with the default values.
            %   linkIndex = d.addBinLinkValve(linkID, from, to)
            %   d.plot;
            %
            % % Adds a new valve form junction '11' to '22' with type = 'PRV'.
            %   linkIndex = d.addBinLinkValve(linkID, from, to, type)
            %   d.plot;
            %
            % % Adds a new valve form junction '11' to '22' with type = 'PRV' and diameter = 6.
            %   linkIndex = d.addBinLinkValve(linkID, from, to, type, diameter)
            %   d.plot;
            %
            % % Adds a new valve form junction '11' to '22' with type = 'PRV', diameter = 6, initial setting = 70 and minor loss coefficient = 0.
            %   linkIndex = d.addBinLinkValve(linkID, from, to, type, diameter, init_setting, minorLoss)
            %   d.plot;
            %
            % Example 2:
            %   d=epanet('NET1.inp');
            %   linkID = {'new_valve_1', 'new_valve_2'};
            %   from = {'11', '12'};
            %   to = {'22', '23'};
            %   type = {'PRV', 'PRV'};
            %   diameter = [6, 10];
            %   init_setting = [70, 55];
            %   minorLoss = [0, 0];
            %
            % % Adds 2 new valves from junctions '11' and '12' to '22' and '23' with types = 'PRV',
            % % diameters = 6 and 10, initial settings = 70 and 55 and minor loss coefficients = 0.
            %   linkIndex = d.addBinLinkValve(linkID, from, to, type, diameter, init_setting, minorLoss)
            %   d.plot;
            %
            % See also addBinLinkPipe, addBinLinkPump, getBinLinksInfo,
            %          addBinNodeJunction, addBinNodeReservoir, addBinNodeTank.
            type = repmat({'GPV'}, 1, length(linkID));
            diameter = zeros(length(linkID), 1);
            diameter(:) = 10;
            init_set = zeros(length(linkID), 1);
            init_set(:) = 0;
            minorloss = zeros(length(linkID), 1);
            minorloss(:) = 0;
            if nargin >= 5
                type = varargin{1};
            end
            if nargin >= 6
                diameter = varargin{2};
            end
            if nargin >= 7
                init_set = varargin{3};
            end
            if nargin >= 8
                minorloss = varargin{4};
            end
            link_index = addBinLink(obj, 'VALVE', linkID, from, to, type, diameter, init_set, minorloss);
        end
        function [Errcode]=addBinCVPipe(obj, newLink, fromNode, toNode, newLength, newDiameter, newRoughness)
            [Errcode]=addLink(obj, 0, newLink, fromNode, toNode, newLength, newDiameter, newRoughness);
        end
        function [Errcode]=addBinPipe(obj, newLink, fromNode, toNode, newLength, newDiameter, newRoughness)
            [Errcode]=addLink(obj, 1, newLink, fromNode, toNode, newLength, newDiameter, newRoughness);
        end
        function [Errcode]=addBinPump(obj, newPumpID, fromNode, toNode, varargin)
            Errcode=-1;
            if length(varargin)==4
                newCurveIDofPump=varargin{1};
                newCurveXvalue=varargin{2};
                newCurveYvalue=varargin{3};
                newCurveType=varargin{4};
                if strcmpi(newCurveType, 'PUMP')
                    addBinCurvePump(obj, newCurveIDofPump, newCurveXvalue, newCurveYvalue);%Flow-Head
                elseif strcmpi(newCurveType, 'EFFICIENCY')
                    addBinCurveEfficiency(obj, newCurveIDofPump, newCurveXvalue, newCurveYvalue);%Flow-Efficiency
                elseif strcmpi(newCurveType, 'VOLUME')
                    addBinCurveVolume(obj, newCurveIDofPump, newCurveXvalue, newCurveYvalue);%Heigh-Volume
                elseif strcmpi(newCurveType, 'HEADLOSS')
                    addBinCurveHeadloss(obj, newCurveIDofPump, newCurveXvalue, newCurveYvalue);%Flow-Headloss
                end
                [Errcode]=addLink(obj, 2, newPumpID, fromNode, toNode, newCurveIDofPump, newCurveXvalue, newCurveYvalue, newCurveType);
            elseif length(varargin)==1
                power=varargin{1};
                [Errcode]=addLink(obj, 2, newPumpID, fromNode, toNode, power);
            end
        end
        function [Errcode]=addBinValvePRV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 3, newLink, fromNode, toNode, diameter, setting); % Pressure Reducing Valve
        end
        function [Errcode]=addBinValvePSV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 4, newLink, fromNode, toNode, diameter, setting); % Pressure Sustaining Valve
        end
        function [Errcode]=addBinValvePBV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 5, newLink, fromNode, toNode, diameter, setting); % Pressure Breaker Valve
        end
        function [Errcode]=addBinValveFCV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 6, newLink, fromNode, toNode, diameter, setting); % Flow Control Valve
        end
        function [Errcode]=addBinValveTCV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 7, newLink, fromNode, toNode, diameter, setting); % Throttle Control Valve
        end
        function [Errcode]=addBinValveGPV(obj, newLink, fromNode, toNode, diameter, setting)
            [Errcode]=addLink(obj, 8, newLink, fromNode, toNode, diameter, setting); % General Purpose Valve
        end
        function [Errcode]=addBinCurvePump(obj, newCurveID, varargin)
            CurveX=varargin{1};
            CurveY=varargin{2};
            [Errcode]=addCurve(obj, newCurveID, CurveX, CurveY, 0);  %ID Flow-OptionsHeadloss
        end
        function [Errcode]=addBinCurveEfficiency(obj, newCurveID, CurveX, CurveY)
            [Errcode]=addCurve(obj, newCurveID, CurveX, CurveY, 1);  %ID Flow-Efficiency
        end
        function [Errcode]=addBinCurveVolume(obj, newCurveID, CurveX, CurveY)
            [Errcode]=addCurve(obj, newCurveID, CurveX, CurveY, 2);  %ID Heigh-Volume
        end
        function [Errcode]=addBinCurveHeadloss(obj, newCurveID, CurveX, CurveY)
            [Errcode]=addCurve(obj, newCurveID, CurveX, CurveY, 3);  %ID Flow-OptionsHeadloss
        end
        function [Errcode]=addBinControl(obj, x, status, y_t_c, param, z, varargin)
            if nargin==6
                [Errcode]=addNewControl(obj, x, status, y_t_c, param, z);
            elseif nargin==5
                [Errcode]=addNewControl(obj, x, status, y_t_c, param);
            elseif nargin==4
                [Errcode]=addNewControl(obj, x, status, y_t_c);
            else
                [Errcode]=addNewControl(obj, x); % add many controls
            end
        end
        function [Errcode]=removeBinNodeID(obj, NodeID)
            [Errcode]=rmNode(obj, NodeID);
        end
        function [Errcode]=removeBinCurveID(obj, CurveID)
            [Errcode]=rmCurveID(obj, CurveID);
        end
        function [Errcode]=removeBinLinkID(obj, LinkID)
            [Errcode]=rmLink(obj, LinkID);
        end
        function [Errcode]=removeBinControlLinkID(obj, ID)
            [Errcode]=rmControl(obj, 1, ID);
        end
        function [Errcode]=removeBinRulesControlLinkID(obj, ID)
            [Errcode]=rmRulesControl(obj, 1, ID);
        end
        function [Errcode]=removeBinRulesControlNodeID(obj, ID)
            [Errcode]=rmRulesControl(obj, 0, ID);
        end
        function [Errcode]=removeBinControlNodeID(obj, ID)
            [Errcode]=rmControl(obj, 0, ID);
        end
        function [Errcode]=setBinFlowUnitsGPM(obj, varargin)
            Errcode = obj.setFlowUnits('GPM', 0, varargin); % gallons per minute
        end
        function [Errcode]=setBinFlowUnitsLPS(obj, varargin)
            Errcode = obj.setFlowUnits('LPS', 0, varargin); %liters per second
        end
        function [Errcode]=setBinFlowUnitsMGD(obj, varargin)
            Errcode = obj.setFlowUnits('MGD', 0, varargin); %million gallons per day
        end
        function [Errcode]=setBinFlowUnitsIMGD(obj, varargin)
            Errcode = obj.setFlowUnits('IMGD', 0, varargin); %Imperial mgd
        end
        function [Errcode]=setBinFlowUnitsCFS(obj, varargin)
            Errcode = obj.setFlowUnits('CFS', 0, varargin); %cubic feet per second
        end
        function [Errcode]=setBinFlowUnitsAFD(obj, varargin)
            Errcode = obj.setFlowUnits('AFD', 0, varargin); %acre-feet per day
        end
        function [Errcode]=setBinFlowUnitsLPM(obj, varargin)
            Errcode = obj.setFlowUnits('LPM', 0, varargin); %liters per minute
        end
        function [Errcode]=setBinFlowUnitsMLD(obj, varargin)
            Errcode = obj.setFlowUnits('MLD', 0, varargin); %million liters per day
        end
        function [Errcode]=setBinFlowUnitsCMH(obj, varargin)
            Errcode = obj.setFlowUnits('CMH', 0, varargin); %cubic meters per hour
        end
        function [Errcode]=setBinFlowUnitsCMD(obj, varargin)
            Errcode = obj.setFlowUnits('CMD', 0, varargin); %cubic meters per day
        end
        function [Errcode]=setBinHeadlossHW(obj)
            [Errcode]=Options(obj, '', 'H-W');  %Hazen-Wiliams
        end
        function [Errcode]=setBinHeadlossDW(obj)
            [Errcode]=Options(obj, '', 'D-W');  %Darcy-Weisbach
        end
        function [Errcode]=setBinHeadlossCM(obj)
            [Errcode]=Options(obj, '', 'C-M');  %Chezy-Manning
        end
        function [Errcode]=setBinLinkPipeLengths(obj, varargin)
            parameter=varargin{1};
            indexParameter=4;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinLinkPipeDiameters(obj, varargin)
            parameter=varargin{1};
            indexParameter=5;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinLinkPipeRoughness(obj, varargin)
            parameter=varargin{1};
            indexParameter=6;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinLinkPipeMinorLoss(obj, varargin)
            parameter=varargin{1};
            indexParameter=7;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinLinkPipeStatus(obj, varargin)
           indexParameter=8;
           if sum(strcmpi(varargin{1}, 'closed')+strcmpi(varargin{1}, 'open')+strcmpi(varargin{1}, 'cv'))==obj.BinLinkPipeCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
           end
           sections={'[PIPES]', '[PUMPS]'};
           [Errcode]=setBinParam(obj, indexParameter, parameter, sections);
        end
        function [Errcode]=setBinLinkPumpStatus(obj, varargin)
            if sum(strcmpi(varargin{1}, 'closed')+strcmpi(varargin{1}, 'open'))==obj.BinLinkPumpCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
            end
            zz=abs(obj.BinLinkPumpCount+1+obj.BinLinkValveCount-obj.BinCountStatuslines);
            sections={'[STATUS]', '[PATTERNS]', 'pump'};
            [Errcode]=setBinParam2(obj, parameter, sections, zz);
        end
        function [Errcode]=setBinLinkPipesParameters(obj, varargin)
            % Initiality
            Lengths=[];Errcode=0;
            Diameters=[];
            Roughness=[];
            Minorloss=[];
            Status={};
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'length'
                        Lengths=varargin{2*i};
                    case 'diameter' % Highlight Node
                        Diameters=varargin{2*i};
                    case 'roughness' % Highlight Link
                        Roughness=varargin{2*i};
                    case 'minorloss' % font size
                        Minorloss=varargin{2*i};
                    case 'status' % color
                        if sum(strcmpi(varargin{2*i}, 'closed')+strcmpi(varargin{2*i}, 'open'))==obj.BinLinkPipeCount
                            Status=varargin{2*i};
                        else
                            warning('Invalid argument found.');Errcode=-1;
                            return;
                        end
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
            [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
            fid = writenewTemp(obj.BinTempfile);
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1}, '[PIPES]')
                    pipes=i;
                end
                if strcmp(tt{1}, '[PUMPS]')
                    pumps=i;
                end
            end
            ll=1;clear atlines;
           for i=pipes:pumps
               % Get first token in the line
               tok = strtok(tlines{i});
               % if isempty(tok), continue; end
               % Skip blank Clines and comments
               if isempty(tok), continue; end
               if strcmp(tok(1), ';')
               elseif sum(tlines{i}=='[')
               % skip
               else
                   clear atlines;
                   atlines = checktlines(tlines{i});

                   if ~isempty(Diameters)%Diameters
                       atlines{5} = num2str(Diameters(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Lengths)%Lengths
                       atlines{4} = num2str(Lengths(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Roughness)%Roughness
                       atlines{6} = num2str(Roughness(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Minorloss)%Minorloss
                       atlines{7} = num2str(Minorloss(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Status)%Status
                       atlines{8} = Status{ll};
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid);
            if obj.Bin==1
                Errcode=reloadNetwork(obj);
            end
        end
        function [Errcode]=setBinNodeJunctionsParameters(obj, varargin)
            % Initiality
            Elevations=[];
            BaseDemands=[];
            patterns=[];
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'elevation'
                        Elevations=varargin{2*i};
                    case 'basedemand'
                        BaseDemands=varargin{2*i};
                    case 'demandpattern'
                        patterns=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
            [tlines]=regexp(fileread([obj.BinTempfile]), '\n', 'split');
            fid = writenewTemp(obj.BinTempfile);
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1}, '[JUNCTIONS]')
                    junctions=i;
                end
                if strcmp(tt{1}, '[RESERVOIRS]')
                    reservoirs=i;
                end
                if strcmp(tt{1}, '[DEMANDS]')
                    demands=i;
                end
                if strcmp(tt{1}, '[STATUS]')
                    status=i;
                end
            end
            ll=1;
            if junctions
               for i=junctions:reservoirs
                   % Get first token in the line
                   tok = strtok(tlines{i});
                   % Skip blank Clines and comments
                   if isempty(tok), continue; end
                   if strcmp(tok(1), ';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       clear atlines;
                       atlines = checktlines(tlines{i});

                       if ~isempty(Elevations)%Elevations
                           atlines{2} = num2str(Elevations(ll));
                           newlines=[];
                           for pp=1:length(atlines)
                               newlines = [newlines, atlines{pp}, blanks(10)];
                           end
                           tlines{i}=newlines;
                       end
                       if ~isempty(BaseDemands) && length(atlines)>2
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{3} = num2str(BaseDemands(ll));
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       if ~isempty(patterns) && length(atlines)>3
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{4} = num2str(patterns{ll});
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       ll=ll+1;
                   end
               end
            end
            if demands
                ll=1;
               for i=demands:status
                   % Get first token in the line
                   tok = strtok(tlines{i});
                   % Skip blank Clines and comments
                   if isempty(tok), continue; end
                   if strcmp(tok(1), ';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       clear atlines;
                       atlines = checktlines(tlines{i});

                       if ~isempty(BaseDemands) && length(atlines)>2%BaseDemands
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{2} = num2str(BaseDemands(ll));
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       if ~isempty(patterns) && length(atlines)>2
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{3} = num2str(patterns{ll});
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp}, blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       ll=ll+1;
                   end
               end
           end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid);
            if obj.Bin==1
                Errcode=reloadNetwork(obj);
            end
        end
        function [Errcode]=setBinNodeTankParameters(obj, varargin)
            % Initiality
            elevations=[];Errcode=0;
            InitLevel=[];
            MinLevel=[];
            MaxLevel=[];
            Diameter=[];
            MinVol=[];
            MixModel={};
            MixFraction=[];
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'elevation'
                        elevations=varargin{2*i};
                    case 'initlevel'
                        InitLevel=varargin{2*i};
                    case 'minlevel'
                        MinLevel=varargin{2*i};
                    case 'maxlevel'
                        MaxLevel=varargin{2*i};
                    case 'diameter'
                        Diameter=varargin{2*i};
                    case 'minvol'
                        MinVol=varargin{2*i};
                    case 'mixmodel'
                        MixModel=varargin{2*i};
                        %v=obj.getBinNodesInfo;
                        mixm=upper(MixModel(find(cellfun('isempty', (varargin{2*i}))==0)));
                        for u=1:length(mixm)
                            if ~sum(strcmp(mixm(u), {'MIX1', 'FIFO', 'LIFO', 'MIXED', '2COMP'}))
                              warning('Invalid argument found.');Errcode=-1;
                              return;
                            end
                        end
                    case 'mixfraction'
                        MixFraction=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
            [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
            fid = writenewTemp(obj.BinTempfile);
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1}, '[TANKS]')
                    tanks=i;
                end
                if strcmp(tt{1}, '[PIPES]')
                    pipes=i;
                end
                if strcmp(tt{1}, '[MIXING]')
                    start1=i;
                end
                if strcmp(tt{1}, '[REACTIONS]')
                    stop1=i;
                end
            end
            ll=1;
           for i=tanks:pipes
               % Get first token in the line
               tok = strtok(tlines{i});
               if isempty(tok), continue; end
               % Skip blank Clines and comments
               if strcmp(tok(1), ';')
               elseif sum(tlines{i}=='[')
               elseif isempty(tok)
               % skip
               else
                   clear atlines;
                   atlines = checktlines(tlines{i});

                   if ~isempty(elevations)%elevations
                       atlines{2} = num2str(elevations(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(InitLevel)
                       atlines{3} = num2str(InitLevel(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MinLevel)
                       atlines{4} = num2str(MinLevel(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MaxLevel)
                       atlines{5} = num2str(MaxLevel(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Diameter)
                       atlines{6} = num2str(Diameter(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MinVol)
                       atlines{7} = num2str(MinVol(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
           if ~isempty(MixModel)
               ll=find(cellfun('isempty', (MixModel))==0);
               ll=ll(1);
           else
               ll=1;
           end
          for i=start1+1:stop1-1
               % Get first token in the line
               tok = strtok(tlines{i});
               if isempty(tok), continue; end
               if isempty(tlines{i})
               elseif strcmp(tok(1), ';')
               else
                   clear atlines;
                   atlines = checktlines(tlines{i});

                   if ~isempty(MixModel) && ll<=obj.BinNodeTankIndex
                       atlines{2} = num2str(MixModel{ll});
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, '              	'];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MixFraction) && ll<=obj.BinNodeTankIndex
                       atlines{3} = num2str(MixFraction(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, '              	'];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
          end
          fprintf(fid, '%s\n', tlines{:});
          fclose(fid);
           if obj.Bin==1
               Errcode=reloadNetwork(obj);
           end
        end
        function [Errcode]=setBinNodeReservoirParameters(obj, varargin)
            % Initiality
            elevations=[];Errcode=0;
            patterns={};
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'elevation'
                        elevations=varargin{2*i};
                    case 'pattern'
                        patterns=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end
            [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
            fid = writenewTemp(obj.BinTempfile);
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1}, '[RESERVOIRS]')
                    reservoirs=i;
                end
                if strcmp(tt{1}, '[TANKS]')
                    tanks=i;
                end
            end
           ll=1;clear atlines;
           for i=reservoirs:tanks
               % Get first token in the line
               tok = strtok(tlines{i});
               % Skip blank Clines and comments
               if isempty(tok), continue; end
               if strcmp(tok(1), ';')
               elseif sum(tlines{i}=='[')
               % skip
               else
                   clear atlines;
                   atlines = checktlines(tlines{i});

                   if ~isempty(elevations)%elevations
                       atlines{2} = num2str(elevations(ll));
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(patterns)
                       atlines{3} = num2str(patterns{ll});
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid);
            if obj.Bin==1
                Errcode=reloadNetwork(obj);
            end
        end
        function [value] = getBinSections(obj)
            % Open epanet input file
            [info]=regexp( fileread(obj.BinTempfile), '\n', 'split');
            sect=0;
            value=struct;
            p=zeros(1, 24);
            for i=1:length(info)
                tline = info{i};
                if ~ischar(tline),   break,   end
                tok = strtok(tline);
                if isempty(tok), continue, end
                if (tok(1) == '[')
                    if strcmpi(tok(1:5), '[VALV'), sect=24;
                    elseif strcmpi(tok(1:5), '[PIPE'), sect=23;
                    elseif strcmpi(tok(1:5), '[TANK'), sect=22;
                    elseif strcmpi(tok(1:5), '[VERT'), sect=21;
                    elseif strcmpi(tok(1:5), '[RESE'), sect=20;
                    elseif strcmpi(tok(1:5), '[JUNC'), sect=19;
                    elseif strcmpi(tok(1:5), '[PUMP'), sect=17;
                    elseif strcmpi(tok(1:5), '[OPTI'), sect=16;
                    elseif strcmpi(tok(1:5), '[REPO'), sect=15;
                    elseif strcmpi(tok(1:5), '[TIME'), sect=14;
                    elseif strcmpi(tok(1:5), '[REAC'), p(13)=0; sect=13;
                    elseif strcmpi(tok(1:5), '[ENER'), sect=12;
                    elseif strcmpi(tok(1:5), '[DEMA'), sect=11;
                    elseif strcmpi(tok(1:5), '[STAT'), sect=10;
                    elseif strcmpi(tok(1:5), '[EMIT'), sect=9;
                    elseif strcmpi(tok(1:5), '[CONT'), sect=8;
                    elseif strcmpi(tok(1:5), '[PATT'), sect=7;
                    elseif strcmpi(tok(1:5), '[CURV'), sect=6;
                    elseif strcmpi(tok(1:5), '[QUAL'), sect=5;
                    elseif strcmpi(tok(1:5), '[SOUR'), sect=4;
                    elseif strcmpi(tok(1:5), '[MIXI'), sect=3;
                    elseif strcmpi(tok(1:5), '[COOR'), sect=1;
                    elseif strcmpi(tok(1:5), '[RULE'), sect=2;
                    elseif (tok(1) == '['), sect=0; continue;
                    end
                end
                if sect==0, continue;
                elseif sect==1
                    p(1)=p(1)+1;
                    value.Coordinates{p(1)}=tline;
                elseif sect==2
                    p(2)=p(2)+1;
                    value.Rules{p(2)}=tline;
                elseif sect==3
                    p(3)=p(3)+1;
                    value.Mixing{p(3)}=tline;
                elseif sect==4
                    p(4)=p(4)+1;
                    value.Sources{p(4)}=tline;
                elseif sect==5
                    p(5)=p(5)+1;
                    value.Quality{p(5)}=tline;
                elseif sect==6
                    p(6)=p(6)+1;
                    value.Curves{p(6)}=tline;
                elseif sect==7
                    p(7)=p(7)+1;
                    value.Patterns{p(7)}=tline;
                elseif sect==8
                    p(8)=p(8)+1;
                    value.Controls{p(8)}=tline;
                elseif sect==9
                    p(9)=p(9)+1;
                    value.Emitters{p(9)}=tline;
                elseif sect==10
                    p(10)=p(10)+1;
                    value.Status{p(10)}=tline;
                elseif sect==11
                    p(11)=p(11)+1;
                    value.Demands{p(11)}=tline;
                elseif sect==12
                    p(12)=p(12)+1;
                    value.Energy{p(12)}=tline;
                elseif sect==13
                    if (tok(1) == '[')
                        p(18)=p(18)+1;p(13)=p(13)+1;
                        value.OptReactions{p(18)}=tline;
                        value.Reactions{p(13)}=tline; continue;
                    end
                    if sum(strcmpi(tok, {'order', 'global', 'limiting', 'roughness'}))
                        value.OptReactions{p(13)}=tline;p(13)=p(13)+1;
                    elseif sum(strcmpi(tok, {'BULK', 'WALL', 'TANK'}))
                        value.Reactions{p(18)}=tline;p(18)=p(18)+1;
                    end
                elseif sect==14
                    p(14)=p(14)+1;
                    value.Times{p(14)}=tline;
                elseif sect==15
                    p(15)=p(15)+1;
                    value.Report{p(15)}=tline;
                elseif sect==16
                    p(16)=p(16)+1;
                    value.Options{p(16)}=tline;
                elseif sect==17
                    p(17)=p(17)+1;
                    value.Pumps{p(17)}=tline;
                elseif sect==19
                    p(19)=p(19)+1;
                    value.Junctions{p(19)}=tline;
                elseif sect==20
                    p(20)=p(20)+1;
                    value.Reservoirs{p(20)}=tline;
                elseif sect==21
                    p(21)=p(21)+1;
                    value.Vertices{p(21)}=tline;
                elseif sect==22
                    p(22)=p(22)+1;
                    value.Tanks{p(22)}=tline;
                elseif sect==23
                    p(23)=p(23)+1;
                    value.Pipes{p(23)}=tline;
                elseif sect==24
                    p(24)=p(24)+1;
                    value.Valves{p(24)}=tline;
                end
            end
        end
        function value = getBinNodeSourceInfo(obj, varargin)
            sections = {'[SOURCES]' '[MIXING]'};
            value = getBinParam(obj, sections);
        end
        function value = getBinPatternsInfo(obj, varargin)
            sections = {'[PATTERNS]' '[CURVES]'};
            value = getBinParam(obj, sections);
        end
        function value = getBinNodeIndex(obj, varargin)
            v=obj.getBinNodeNameID;
            if isempty(varargin)
                value=1:length(v.BinNodeNameID);
            else
                value=find(strcmpi(v.BinNodeNameID, varargin{1}));
            end
        end
        function value = getBinLinkIndex(obj, varargin)
            if isempty(varargin)
                value = 1:length(obj.getBinLinkNameID.BinLinkNameID);
            else
                value = find(strcmp(obj.getBinLinkNameID.BinLinkNameID, varargin{1}));
            end
        end
        function value = getBinPatternIndex(obj, varargin)
            value = find(strcmpi(obj.getBinPatternsInfo.BinPatternNameID, varargin{1}));
        end
        function value = getBinCurvesInfo(obj)
            [value.BinCurveNameID, value.BinCurveXvalue, value.BinCurveYvalue, value.BinCurveAllLines, value.BinCurveTypes, value.BinCurveCount, value.BinCTypes]=CurveInfo(obj);
        end
        function [value] = Binplot(obj, varargin)
            %Plots network in a new Matlab figure
            %Arguments:
            % 'nodes': yes/no
            % 'links': yes/no
            % 'line' : yes/no
            % 'point': yes/no
            % 'highlightnode': array of node IDs
            % 'highlightlink': array of link IDs
            % 'fontsize': number (px)
            % 'colornode': array of node IDs
            % 'colorlink': array of link id
            % 'axes': axes coordinates
            % 'linksindex': yes
            % 'nodesindex': yes
            % Example:
            % d.Binplot('nodes', 'yes', 'links', 'yes', 'highlightnode', {'10', '11'}, 'highlightlink', {'10'}, 'fontsize', 8);
            % d.Binplot('line', 'no');
            % d.Binplot('point', 'no', 'linksindex', 'yes');
            % d.Binplot('linksindex', 'yes', 'fontsize', 8);
            % d.Binplot('nodesindex', 'yes', 'fontsize', 14);
            [value] = plotnet(obj, 'bin', 1, varargin{:});
        end
        function value = getBinNumberReportingPeriods(obj, varargin)
            value = getBinComputedTimeSeries(obj, 27);
        end
        function value = getBinSimulationDuration(obj, varargin)
            value = getBinComputedTimeSeries(obj, 28);
        end
        function value = getBinElevationEachNode(obj, varargin)
            value = getBinComputedTimeSeries(obj, 1);
        end
        function value = getBinLengthEachLink(obj, varargin)
            value = getBinComputedTimeSeries(obj, 2);
        end
        function value = getBinDiameterEachLink(obj, varargin)
            value = getBinComputedTimeSeries(obj, 3);
        end
        function value = getBinComputedPumpIndexListLinks(obj, varargin)
            value = getBinComputedTimeSeries(obj, 4);
        end
        function value = getBinComputedPumpUtilization(obj, varargin)
            value = getBinComputedTimeSeries(obj, 5);
        end
        function value = getBinComputedAverageEfficiency(obj, varargin)
            value = getBinComputedTimeSeries(obj, 6);
        end
        function value = getBinComputedAverageKwattsOrMillionGallons(obj, varargin)
            value = getBinComputedTimeSeries(obj, 7);
        end
        function value = getBinComputedAverageKwatts(obj, varargin)
            value = getBinComputedTimeSeries(obj, 8);
        end
        function value = getBinComputedPeakKwatts(obj, varargin)
            value = getBinComputedTimeSeries(obj, 9);
        end
        function value = getBinComputedAverageCostPerDay(obj, varargin)
            value = getBinComputedTimeSeries(obj, 10);
        end
        function value = getBinComputedNodeDemand(obj, varargin)
            value = getBinComputedTimeSeries(obj, 11);
        end
        function value = getBinComputedNodeHead(obj, varargin)
            value = getBinComputedTimeSeries(obj, 12);
        end
        function value = getBinComputedNodePressure(obj, varargin)
            value = getBinComputedTimeSeries(obj, 13);
        end
        function value = getBinComputedNodeQuality(obj, varargin)
            value = getBinComputedTimeSeries(obj, 14);
        end
        function value = getBinComputedLinkFlow(obj, varargin)
            value = getBinComputedTimeSeries(obj, 15);
        end
        function value = getBinComputedLinkVelocity(obj, varargin)
            value = getBinComputedTimeSeries(obj, 16);
        end
        function value = getBinComputedLinkHeadloss(obj, varargin)
            value = getBinComputedTimeSeries(obj, 17);
        end
        function value = getBinComputedLinkQuality(obj, varargin)
            value = getBinComputedTimeSeries(obj, 18);
        end
        function value = getBinComputedLinkStatus(obj, varargin)
            value = getBinComputedTimeSeries(obj, 19);
        end
        function value = getBinComputedLinkSetting(obj, varargin)
            value = getBinComputedTimeSeries(obj, 20);
        end
        function value = getBinComputedLinkReactionRate(obj, varargin)
            value = getBinComputedTimeSeries(obj, 21);
        end
        function value = getBinComputedLinkFrictionFactor(obj, varargin)
            value = getBinComputedTimeSeries(obj, 22);
        end
        function value = getBinComputedAverageBulkReactionRate(obj, varargin)
            value = getBinComputedTimeSeries(obj, 23);
        end
        function value = getBinComputedAverageWallReactionRate(obj, varargin)
            value = getBinComputedTimeSeries(obj, 24);
        end
        function value = getBinComputedAverageTankReactionRate(obj, varargin)
            value = getBinComputedTimeSeries(obj, 25);
        end
        function value = getBinComputedAverageSourceInflow(obj, varargin)
            value = getBinComputedTimeSeries(obj, 26);
        end
        function value = getBinComputedAllParameters(obj, varargin)
            [fid, binfile, rptfile] = runEPANETexe(obj);
            value = readEpanetBin(fid, binfile);

            % Remove report bin, files @#
            warning off;
            fclose('all');
            files=dir('@#*');
            if ~isempty(files); delete('@#*'); end
            warning on;
        end
        function [info, tline, allines] = readInpFile(obj, varargin)
            [info, tline, allines] = readAllFile(obj.BinTempfile);
        end
        function value = getBinNodesInfo(obj)
            valueL = obj.getBinLinksInfo;
            % Open epanet input file
            [~, info] = obj.readInpFile;
            sect=0;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                    % [JUNCTIONS] section
                    if strcmpi(tok(1:5), '[JUNC')
                        sect=1;
                        value.BinNodeJunctionNameID={};
                        value.BinNodeType={};
                        value.BinNodeJunctionIndex=[];
                        value.BinNodeJunctionElevation=[];
                        value.BinNodeJunctionsBaseDemands=[];
                        value.BinNodeJunDemandPatternNameID={};
                        value.BinNodeJunctionsBaseDemandsID={};k=1;
                        continue;
                    elseif strcmpi(tok(1:5), '[RESE')
                        sect=2;r=1;
                        value.BinNodeReservoirNameID={};
                        value.BinNodeReservoirIndex=[];
                        value.BinNodeReservoirElevation=[];
                        value.BinNodeResDemandPatternNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5), '[TANK')
                        sect=3;p=1;
                        value.BinNodeTankNameID={};
                        value.BinNodeTankIndex=[];
                        value.BinNodeTankElevation=[];
                        value.BinNodeTankInitialLevel=[];
                        value.BinNodeTankMinimumWaterLevel=[];
                        value.BinNodeTankMaximumWaterLevel=[];
                        value.BinNodeTankDiameter=[];
                        value.BinNodeTankMinimumWaterVolume=[];
                        continue;
                    elseif strcmpi(tok(1:5), '[PIPE')
                        value.BinNodeJunctionCount = length(value.BinNodeJunctionNameID);
                        value.BinNodeReservoirCount = length(value.BinNodeReservoirNameID);
                        value.BinNodeTankCount = length(value.BinNodeTankNameID);
                        value.BinNodeTankReservoirCount = value.BinNodeReservoirCount + value.BinNodeTankCount;
                        value.BinNodeCount = value.BinNodeJunctionCount+value.BinNodeTankReservoirCount;
                        value.BinNodeNameID=[value.BinNodeJunctionNameID value.BinNodeReservoirNameID value.BinNodeTankNameID];
                        vx = NaN(value.BinNodeCount, 1);
                        vy = NaN(value.BinNodeCount, 1);
                        vertx = cell(valueL.BinLinkCount, 1);
                        verty = cell(valueL.BinLinkCount, 1);
                        nvert = zeros(valueL.BinLinkCount, 1);
                        value.BinNodeInitialQuality=zeros(1, value.BinNodeCount);
                        sect=4;
                        continue;
                    elseif strcmpi(tok(1:5), '[DEMA') %&& max(value.BinNodeJunctionsBaseDemands)==0
                        sect=8;d=1;
                        %value.BinNodeJunDemandPatternNameID={};
                        continue;
                        % [QUALITY] section
                    elseif strcmpi(tok(1:5), '[QUAL')
                        sect=12;d=1;
                        value.BinCountInitialQualitylines=0;
                        continue;
                        % [MIXING] section
                    elseif strcmpi(tok(1:5), '[MIXI')
                        sect=14;d=1;
                        value.BinNodeTankMixModel={};
                        value.BinNodeTankMixID={};
                        value.BinNodeTankMinimumFraction=[];
                        continue;
                        % [COORDINATES] section
                    elseif strcmpi(tok(1:5), '[COOR')
                        sect=17;
                        continue;
                        % [VERTICES] section
                    elseif strcmpi(tok(1:5), '[VERT')
                        sect=18;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                atline = checktlines(tline);

                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    value.BinNodeJunctionNameID{k}=atline{1};
                    value.BinNodeJunctionIndex(k)=k;
                    value.BinNodeJunctionElevation(k)=str2double(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            value.BinNodeJunctionsBaseDemands(k)=single(str2double(atline{3}));
                            if length(atline)>3
                                if ~sum(atline{4}==';')
                                    value.BinNodeJunDemandPatternNameID{k}=atline{4};
                                else
                                    value.BinNodeJunDemandPatternNameID{k}='';
                                end
                            else
                                value.BinNodeJunDemandPatternNameID{k}='';
                            end
                        end
                    end
                    k=k+1;
                elseif sect==2
                    value.BinNodeReservoirNameID{r}=atline{1};
                    value.BinNodeReservoirIndex(r)=k;
                    value.BinNodeReservoirElevation(r)=str2double(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            value.BinNodeResDemandPatternNameID{r}=atline{3};
                        else
                            value.BinNodeResDemandPatternNameID{r}='';
                        end
                    else
                        value.BinNodeResDemandPatternNameID{r}='';
                    end
                    k=k+1;
                    r=r+1;
                elseif sect==3
                    value.BinNodeTankNameID{p}=atline{1};
                    value.BinNodeTankIndex(p)=k;
                    value.BinNodeTankElevation(p)=str2double(atline{2});
                    value.BinNodeTankInitialLevel(p)=single(str2double(atline{3}));
                    value.BinNodeTankMinimumWaterLevel(p)=str2double(atline{4});
                    value.BinNodeTankMaximumWaterLevel(p)=single(str2double(atline{5}));
                    value.BinNodeTankDiameter(p)=str2double(atline{6});
                    value.BinNodeTankMinimumWaterVolume(p)=single(str2double(atline{7}));
                    k=k+1;
                    p=p+1;
                    % Demands
                elseif sect==8
                    indd=find(strcmpi(value.BinNodeNameID, atline{1}));
                    if ~isempty(value.BinNodeJunctionsBaseDemandsID)
                        if strcmp(value.BinNodeJunctionsBaseDemandsID{end}, atline{1})
                            value.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2}));
                        else
                            value.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2}));
                        end
                    else
                        value.BinNodeJunctionsBaseDemands(indd)=single(str2double(atline{2}));
                    end
                    value.BinNodeJunctionsBaseDemandsID{d}=atline{1};
                    if length(atline)>2
                        value.BinNodeJunDemandPatternNameID{indd}=atline{3};
                    else
                        value.BinNodeJunDemandPatternNameID{indd}='';
                    end
                    d=d+1;
                    % Quality
                elseif sect==12
                    Hh=find(strcmpi(value.BinNodeNameID, atline{1}));
                    value.BinNodeInitialQuality(Hh)=str2double(atline{2});
                    value.BinCountInitialQualitylines=d;
                    d=d+1;
                elseif sect==14
                    value.BinNodeTankMixID{d}=atline{1};
                    value.BinNodeTankMixModel{d}=atline{2};
                    value.BinNodeTankMinimumFraction(d)=str2double(atline{3});
                    d=d+1;
                elseif sect==17
                    A = textscan(tline, '%s %f %f');
                    % get the node index
                    a=strcmp(A{1}, value.BinNodeNameID);
                    index=strfind(a, 1);
                    if isempty(index), continue; end
                    vx(index) = A{2};
                    vy(index) = A{3};
                    % Vertices
                elseif sect==18
                    A = textscan(tline, '%s %f %f');
                    index =  find(strcmp(valueL.BinLinkNameID, A{1}));
                    if isempty(index), continue; end
                    nvert(index) = nvert(index) + 1;
                    vertx{index}(nvert(index)) = A{2};
                    verty{index}(nvert(index)) = A{3};
                end
            end
            value.BinNodeType(value.BinNodeJunctionIndex)=obj.TYPENODE(1);
            value.BinNodeType(value.BinNodeReservoirIndex)=obj.TYPENODE(2);
            value.BinNodeType(value.BinNodeTankIndex)=obj.TYPENODE(3);
            value.BinNodeCoordinates{1} = vx;
            value.BinNodeCoordinates{2} = vy;
            value.BinNodeCoordinates{3} = vertx;
            value.BinNodeCoordinates{4} = verty;
            value.BinNodeElevations = single([value.BinNodeJunctionElevation value.BinNodeReservoirElevation value.BinNodeTankElevation]);
            value.BinNodeBaseDemands = single([value.BinNodeJunctionsBaseDemands zeros(1, value.BinNodeReservoirCount) zeros(1, value.BinNodeTankCount)]);
            value.BinNodeIndex=obj.getBinNodeIndex;
        end
        function [valueQualityType] = getBinQualType(obj)
            % Open epanet input file
            [info]=regexp( fileread(obj.BinTempfile), '\n', 'split');
            sect=0;valueQualityType={};
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                        % [OPTIONS] section
                    if strcmpi(tok(1:5), '[OPTI')
                        sect=3;
                    end
                end
                if sect==0
                    continue;
                elseif sect==3
                    if strcmp(regexp(tline, '\w*QUALITY*\w', 'match'), 'QUALITY')
                        tl=regexp(tline, '\s*', 'split');mm=2;
                        if isempty(tl{1})
                            mm=3;
                        end
                        valueQualityType=tl(mm);
                    end
                end
            end
        end
        function [valueCoord, valueRule] = getBinCoordRuleSections(obj, file)
            % Open epanet input file
            [info]=regexp( fileread(file), '\n', 'split');
            sect=0;d=1;valueCoord={};valueRule={};dRule=1;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                        % [COORDINATES] section
                    if strcmpi(tok(1:5), '[COOR')
                        sect=1;
                        valueCoord{d}=tline;
                        continue;
                        % [RULES] section
                    elseif strcmpi(tok(1:5), '[RULE')
                        sect=2;
                        valueRule{dRule}=tline;
                        continue;
                    elseif strcmpi(tok(1:5), '[ENER')
                        sect=0;
                        continue;
                    end
                end
                if sect==0
                    continue;
                elseif sect==1
                    d=d+1;
                    valueCoord{d}=tline;
                elseif sect==2
                    dRule=dRule+1;
                    valueRule{dRule}=tline;
                end
            end
        end
        function value = getBinNodeCoordinates(obj)
            value = BinNodeCoords(obj, 0);
        end
        function value = getNodeCoordinates(obj, varargin)
            cnt = obj.getLinkCount;
            if ~cnt
                warning('Not enough network coordinates.');
                return;
            end
            vert_function = 0;
            coord_function = 0;
            if sum(strcmp(obj.libFunctions, 'ENgetvertexcount'))
                vertices = obj.getLinkVertices;
                for i=1:cnt
                    if ~isempty(vertices{i})
                        vertx{i} = vertices{i}.x;
                        verty{i} = vertices{i}.y;
                    else
                        vertx{i} = [];
                        verty{i} = [];
                    end
                end
                vert_function = 1;
            end
            if sum(strcmp(obj.libFunctions, 'ENgetcoord'))
                try
                    indices = getNodeIndices(obj, varargin);j=1;
                    for i=indices
                        [obj.Errcode, vx(j), vy(j)]=ENgetcoord(i, obj.LibEPANET); j=j+1;
                        error(obj.getError(obj.Errcode));
                    end
                catch
                end
                coord_function = 1;
            end
            n1_value = [];
            n2_value = [];
            if coord_function == 0
                n1_value = BinNodeCoords(obj, 0);
            elseif vert_function == 0
                n2_value = BinNodeCoords(obj, 1);
            end

            if isempty(varargin)
                if isempty(n1_value)
                    value{1} = vx;
                    value{2} = vy;
                else
                    value{3} = n1_value{3};
                    value{4} = n1_value{4};
                end
                if isempty(n2_value)
                    value{3} = vertx;
                    value{4} = verty;
                else
                    value{3} = n2_value{3};
                    value{4} = n2_value{4};
                end
            else
                ind_var = 1:length(varargin{1});
                value = [vx(ind_var) vy(ind_var)];
            end
        end
        function value = getBinNodeNameID(obj)
            % Open epanet input file
            % Open epanet input file
            [~, info] = obj.readInpFile;sect=0;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                    % [JUNCTIONS] section
                    if strcmpi(tok(1:5), '[JUNC')
                        sect=1;k=1;
                        value.BinNodeJunctionNameID={};
                        continue;
                    elseif strcmpi(tok(1:5), '[RESE')
                        sect=2;r=1;
                        value.BinNodeReservoirNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5), '[TANK')
                        sect=3;p=1;
                        value.BinNodeTankNameID={};
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        value.BinNodeNameID=[value.BinNodeJunctionNameID value.BinNodeReservoirNameID value.BinNodeTankNameID];
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                atline = checktlines(tline);
                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    value.BinNodeJunctionNameID{k}=atline{1};
                    k=k+1;
                elseif sect==2
                    value.BinNodeReservoirNameID{r}=atline{1};
                    r=r+1;
                elseif sect==3
                    value.BinNodeTankNameID{p}=atline{1};
                    p=p+1;
                end
            end
        end
        function addBinLinkVertices(obj, linkID, x, y)
            % Adds interior vertex points to network links.
            %
            % % The examples are based on d=epanet('NET1.inp');
            %
            % Example 1:
            %   d = epanet('NET1.inp');
            %   linkID='10';
            %   x = 28;                                % One X coordinate selected for a vertex.
            %   y = 68;                                % One Y coordinate selected for a vertex.
            %   d.addBinLinkVertices(linkID, x, y)     % Adds one vertex to the link with ID label = '10'
            %   d.getBinLinkVertices(linkID)           % Retrieves the link's vertex
            %
            % Example 2:
            %   d = epanet('NET1.inp');
            %   linkID_1 = '10';
            %   x = [20, 25];                          % Two X coordinates selected for a vertex.
            %   y = [66, 67];                          % Two Y coordinates selected for a vertex.
            %   d.addBinLinkVertices(linkID_1, x, y)   % Adds two vertices to the link with ID label = '10'
            %
            %   linkID_2='11';
            %   x = [33, 38, 43, 45, 48];
            %   y = [74, 76, 76, 73, 74];
            %   d.addBinLinkVertices(linkID_2, x, y)   % Adds multiple vertices to the link with ID label = '11'
            %
            %   d.getBinLinkVertices(linkID_1)         % Retrieves the link's vertices.
            %   d.getBinLinkVertices(linkID_2)
            %
            % See also deleteBinLinkVertices, setBinLinkVertices, getBinLinkVertices,
            %          getBinLinkVerticesCount, addLinkPipe, addNodeJunction.
            if obj.Bin, obj.Errcode = obj.saveInputFile(obj.BinTempfile); end
            fid = fopen(obj.BinTempfile); % Opens the file for read access
            %
            % Creates the string that will be set under the [VERTICES] section
            %
            str = linkID;
            for i=1:length(x)
                if i>1
                    str = [str, linkID];
                end
                str = [str, blanks(10), num2str(x(i)), blanks(10) num2str(y(i)), char(10)];
            end
            %
            % Creates the entire text that will replace the .inp file
            %
            texta = char;
            while ~feof(fid)
                aline = fgetl(fid);
                texta = [texta, aline, char(10)];
                if strcmp(aline, '[VERTICES]')
                   fline = fgetl(fid);
                   while ~isempty(strfind(fline, '['))   %contains
                       texta = [texta, fline, char(10)];
                       fline = fgetl(fid);
                   end
                   texta = [texta, str];
                   break
                end
            end
            while isempty(strfind(fline, '[END]'))   %contains
                texta = [texta, fline, char(10)];
                fline = fgetl(fid);
            end
            texta = [texta, '[END]'];
            fid = fopen(obj.BinTempfile, 'w');   % Opens file for writing and discard existing contents
            fprintf(fid, texta);   % Writes the new text in the .inp file
            fclose('all');
            if obj.Bin, obj.Errcode = reloadNetwork(obj); end
        end
        function Errcode = deleteBinLinkVertices(obj, varargin)
            % Deletes interior vertex points of network links.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d = epanet('NET1.inp');
            %   linkID = '10';
            %   x = [20, 25];                                  % Two X coordinates selected for a vertex
            %   y = [66, 67];                                  % Two Y coordinates selected for a vertex
            %   d.addBinLinkVertices(linkID, x, y)             % Adds two vertices to the link with ID label = '10'
            %
            %   d.getBinLinkVertices(linkID)                   % Retrieves all vertices of a link given it's ID label
            %
            %   d.deleteBinLinkVertices                        % Deletes all vertices from all links
            %   d.deleteBinLinkVertices(linkID)                % Deletes all vertices from a link given it's ID label
            %   vertexIndex = 1;
            %   d.deleteBinLinkVertices(linkID, vertexIndex)   % Deletes a certain vertex from a link given the link ID label and the vertex index
            %
            %   d.getBinLinkVertices(linkID)
            %
            % See also addBinLinkVertices, getBinLinkVertices, getBinLinkVerticesCount,
            %          setBinLinkVertices, addLinkPipe, addNodeJunction.
            cnt = obj.getBinLinkVerticesCount;
            Errcode = 0;
            if cnt == 0
                Errcode = 'No vertices found in the network';
                error(Errcode)
            end
            fid = fopen(obj.BinTempfile); % Opens the file for read access
            %
            % Creates the string that will be set under the [VERTICES] section
            % Creates the entire text that will replace the .inp file
            %
            texta = char;
            while ~feof(fid)
                aline = fgetl(fid);
                texta = [texta, aline, char(10)];
                if strcmp(aline, '[VERTICES]')
                   texta = [texta, fgetl(fid), char(10)];
                   % If no input arguments given, deletes all vertices
                   if nargin == 1
                       for i = 1:cnt
                           fgetl(fid);
                       end
                   % If one input is given(i.e. link ID), deletes all vertices of this link
                   elseif nargin == 2
                       linkID = varargin{1};
                       for i = 1:cnt
                           aline = fgetl(fid);
                           if isempty(strfind(aline, linkID))
                                texta = [texta, aline, char(10)];
                           end
                       end
                   % If 2 inputs are given(i.e. link ID & integer), deletes the certain vertex of that link
                   elseif nargin == 3
                       linkID = varargin{1};
                       j = 1;
                       for i = 1:cnt
                           aline = fgetl(fid);
                           if isempty(strfind(aline, linkID))
                                texta = [texta, aline, char(10)];
                           else
                               if j ~= varargin{2}
                                   texta = [texta, aline, char(10)];
                               end
                               j = j + 1;
                           end
                       end
                   end
                end
            end
            fid = fopen(inpfile, 'w');   % Opens file for writing and discard existing contents
            fprintf(fid, texta);   % Writes the new text in the .inp file
            fclose('all');
        end
        function value = getBinLinkVerticesCount(obj, varargin)
            % Retrieves the number of vertices.
            %
            % Example:
            %   d = epanet('NET1.inp');
            %   linkID_1 = '10';
            %   x = [20, 25];                          % Two X coordinates selected for a vertex.
            %   y = [66, 67];                          % Two Y coordinates selected for a vertex.
            %   d.addBinLinkVertices(linkID_1, x, y)   % Adds two vertices to the link with ID label = '10'
            %
            %   linkID_2='11';
            %   x = [33, 38, 43, 45, 48];
            %   y = [74, 76, 76, 73, 74];
            %   d.addBinLinkVertices(linkID_2, x, y)   % Adds multiple vertices to the link with ID label = '11'
            %
            %   d.getBinLinkVerticesCount
            %
            %   d.getBinLinkVerticesCount(linkID_2)
            %
            % See also getBinLinkVertices, getLinkCount, getNodeCount.
            fid = fopen(obj.BinTempfile); % Opens the file for read access
            value = 0;
            while ~feof(fid)
                aline = fgetl(fid);
                if strcmp(aline, '[VERTICES]')
                    while true
                        aline = fgetl(fid);
                        if ~isempty(strfind(aline, '['))
                            if strcmpi(aline, '[END]')
                                break;
                            end
                            if nargin == 1
                                value = value - 2;
                            end
                            break
                        end
                        if nargin == 2
                            if ~isempty(strfind(aline, varargin{1}))
                                value = value + 1;
                            end
                        else
                            value = value + 1;
                        end
                    end
                end
            end
            fclose('all');
        end
        function data = getBinLinkVertices(obj, varargin)
            % Retrieves the link vertices.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d = epanet('NET1.inp');
            %   linkID = '10';
            %   x = [20, 25];                                % Two X coordinates selected for a vertex
            %   y = [66, 67];                                % Two Y coordinates selected for a vertex
            %   d.addBinLinkVertices(linkID, x, y)           % Adds two vertices to the link with ID label = '10'
            %
            %   d.getBinLinkVertices                         % Retrieves all vertices of all links and stores them in cells
            %   d.getBinLinkVertices(linkID)                 % Retrieves all vertices of a link given it's ID label
            %
            % See also setBinLinkVertices, addBinLinkVertices, deleteBinLinkVertices,
            %          getBinLinkVerticesCount, getNodeCoordinates.

            % reload the network
            fid = fopen(obj.BinTempfile); % Opens the file for read access
            while ~feof(fid)
                aline = fgetl(fid);
                if strcmp(aline, '[VERTICES]')
                    j = 1;
                    while true
                        aline = fgetl(fid);
                        if ~isempty(strfind(aline, '['))
                            break;
                        end
                        bline = strsplit(aline);
                        bline = bline(~cellfun(@isempty, bline));
                        if isempty(bline)
                            continue
                        end
                        if strcmp(bline{1}(1), ';')
                            continue
                        end
                        data{j} = bline;
                        ids{j} = bline{1};
                        j = j +1;
                    end
                end
            end
            if nargin == 2
                indices = find(strcmp(ids, varargin{1}));
                data = data(indices);
            end
           fclose('all');
        end
        function setBinLinkVertices(obj, linkID, x, y, varargin)
            % Sets interior vertex points of network links.
            %
            % % The example is based on d=epanet('NET1.inp');
            %
            % Example:
            %   d = epanet('NET1.inp');
            %   linkID_1 = '10';
            %   x = [20, 25];                          % Two X coordinates selected for a vertex.
            %   y = [66, 67];                          % Two Y coordinates selected for a vertex.
            %   d.addBinLinkVertices(linkID_1, x, y)   % Adds two vertices to the link with ID label = '10'
            %
            %   linkID_2 = '11';
            %   x = [33, 38, 43, 45, 48];
            %   y = [74, 76, 76, 73, 74];
            %   d.addBinLinkVertices(linkID_2, x, y)   % Adds multiple vertices to the link with ID label = '11'
            %
            %   d.getBinLinkVertices(linkID_1)         % Retrieves the link's vertices.
            %   d.getBinLinkVertices(linkID_2)
            %
            %   % Deletes all vertices and adds the new vertices.
            %   x = [22, 24, 28];
            %   y = [69, 68, 69];
            %   d.setBinLinkVertices(linkID_1, x, y)
            %   d.getBinLinkVertices(linkID_1)
            %
            %   % Replaces a certain vertex given it's index with a new vertex.
            %   x  = 39;
            %   y = 75;
            %   vertexIndex = 2;
            %   d.setBinLinkVertices(linkID_2, x, y, vertexIndex)
            %   d.getBinLinkVertices(linkID_2)
            %
            % See also addBinLinkVertices, deleteBinLinkVertices, getBinLinkVertices,
            %          getBinLinkVerticesCount, addLinkPipe, addNodeJunction.
            if obj.Bin, obj.Errcode = obj.saveInputFile(obj.BinTempfile, 1); end
            if nargin == 4
                obj.deleteBinLinkVertices(linkID);
                obj.addBinLinkVertices(linkID, x, y);
            end
            if nargin == 5
                cnt = obj.getBinLinkVerticesCount;
                filepath = regexp(obj.TempInpFile, '\\', 'split');   % Finds the .inp file
                inpfile = filepath{end};
                fid = fopen(inpfile); % Opens the file for read access
                texta = char;
                while ~feof(fid)
                    aline = fgetl(fid);
                    texta = [texta, aline, char(10)];
                    if strcmp(aline, '[VERTICES]')
                       texta = [texta, fgetl(fid), char(10)];
                       j = 1;
                       for i = 1:cnt
                           bline = fgetl(fid);
                           if ~isempty(strfind(bline, linkID)) && j == varargin{1}
                                texta = [texta, linkID, blanks(10), num2str(x), blanks(10), num2str(y), char(10)];
                           else
                                texta = [texta, bline, char(10)];
                           end
                           if ~isempty(strfind(bline, linkID))
                                j = j + 1;
                           end
                       end
                    end
                end
                fid = fopen(inpfile, 'w');   % Opens file for writing and discard existing contents
                fprintf(fid, texta);   % Writes the new text in the .inp file
                fclose('all');
                if obj.Bin, obj.Errcode = reloadNetwork(obj); end
            end
        end
        function value = getBinLinkNameID(obj)
            sect=0;i=1;t=1;q=1;
            % Open epanet input file
            [~, info] = obj.readInpFile;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                       % [PIPES] section
                    if strcmpi(tok(1:5), '[PIPE')
                        sect=1;
                        value.BinLinkPipeNameID={};
                        continue;
                        % [PUMPS] section
                    elseif strcmpi(tok(1:5), '[PUMP')
                        sect=2;
                        value.BinLinkPumpNameID={};
                        continue;
                        % [VALVES] section
                    elseif strcmpi(tok(1:5), '[VALV')
                        sect=3;
                        value.BinLinkValveNameID={};
                        continue;
                    elseif strcmpi(tok(1:5), '[REAC')
                        sect=4;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                atline = checktlines(tline);
                if sect==0
                    continue;
                    % Links
                elseif sect==1
                    value.BinLinkPipeNameID{t}=atline{1};
                    t=t+1;
                elseif sect==2
                    value.BinLinkPumpNameID{q}=atline{1};
                    q=q+1;
                elseif sect==3
                    value.BinLinkValveNameID{i}=atline{1};
                    i=i+1;
                end
            end
            value.BinLinkNameID = [value.BinLinkPipeNameID value.BinLinkPumpNameID value.BinLinkValveNameID];
        end
        function value = getBinLinksInfo(obj)
            sect=0;i=1;t=1;q=1;d=1;value=[];
            % Open epanet input file
            [~, info] = obj.readInpFile;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                [value, cont, sect, i, t, q, d] = getLV(tok, value, sect, tline, i, t, q, d);
                if cont==0
                    break;
                end
            end
            if ~sum(value.BinLinkBulkReactionCoeff)
                value.BinLinkBulkReactionCoeff=value.BinLinkGlobalBulkReactionCoeff*ones(1, value.BinLinkCount);
                value.BinLinkBulkReactionCoeff(value.BinLinkPumpIndex)=0;
                value.BinLinkBulkReactionCoeff(value.BinLinkValveIndex)=0;
            end
            if ~sum(value.BinLinkWallReactionCoeff)
                value.BinLinkWallReactionCoeff=value.BinLinkGlobalWallReactionCoeff*ones(1, value.BinLinkCount);
                value.BinLinkWallReactionCoeff(value.BinLinkPumpIndex)=0;
                value.BinLinkWallReactionCoeff(value.BinLinkValveIndex)=0;
            end
            value.BinLinkSettings = [value.BinLinkPipeRoughness zeros(1, value.BinLinkPumpCount) value.BinLinkValveSetting]';
            value.BinLinkDiameters = single([value.BinLinkPipeDiameters zeros(1, value.BinLinkPumpCount) value.BinLinkValveDiameters]);
            value.BinLinkLengths = single([value.BinLinkPipeLengths zeros(1, value.BinLinkPumpCount) zeros(1, value.BinLinkValveCount)]);
            value.BinLinkRoughnessCoeff = [value.BinLinkPipeRoughness zeros(1, value.BinLinkPumpCount) zeros(1, value.BinLinkValveCount)];
            value.BinLinkType(value.BinLinkPipeIndex)=obj.TYPELINK(2);
            value.BinLinkType(value.BinLinkPumpIndex)=obj.TYPELINK(3);
            value.BinLinkType(value.BinLinkValveIndex)=value.BinLinkValveType;

            b={};
            for i=1:value.BinLinkCount
                ind=find((strcmp(value.BinLinkInitialStatusNameID, value.BinLinkNameID{i}))==1);
                if isempty(ind), ind=0; end
                if ind~=0
                    if sum(value.BinLinkPumpIndex==i)
                        bb{i}=value.BinLinkInitialStatus{ind};
                        r=value.BinLinkNameID(i);
                        b{i}=r{1};
                    else
                        b{i}=value.BinLinkNameID{i};
                        bb{i}=value.BinLinkInitialStatus{ind};
                    end
                else
                    b{i}=value.BinLinkNameID{i};
                    bb{i}='Open';
                end
            end
            value.BinLinkInitialStatusNameID=b;
            value.BinLinkInitialStatus=bb;
            value.BinLinkPumpStatus=value.BinLinkInitialStatus(value.BinLinkPumpIndex);
            value.BinLinkPumpStatusNameID=value.BinLinkInitialStatusNameID(value.BinLinkPumpIndex);
            value.BinLinkValveStatus=value.BinLinkInitialStatus(value.BinLinkValveIndex);
            value.BinLinkValveStatusNameID=value.BinLinkInitialStatusNameID(value.BinLinkValveIndex);
            value.BinLinkIndex=obj.getBinLinkIndex;
        end
        function value = getBinControlsInfo(obj)
            % Open epanet input file
            [~, info] = obj.readInpFile;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                        % [CONTROLS] section
                    if strcmpi(tok(1:5), '[CONT')
                        sect=1;d=1;
                        value.BinControlsInfo={};
                        value.BinControlLinksID={};
                        value.BinControlNodesID={};
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end

                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    clear atline;
                    atline = checktlines(tline);
                    value.BinControlsInfo{d}=atline;
                    if length(atline)>1
                        value.BinControlLinksID{d}=atline{2};
                        t = regexp(tline, '\w*TIME\w*', 'match');
                        if length(t)==0
                            value.BinControlNodesID{d}=atline{6};
                        end
                    end
                    d=d+1;
                end
            end
            value.BinControlRulesCount=length(value.BinControlsInfo);
        end
        function value = getBinRulesControlsInfo(obj)
            % Open epanet input file
            value={};
            [~, info] = obj.readInpFile;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                        % [RULES] section
                    if strcmpi(tok(1:5), '[RULE')
                        sect=1;d=1;
                        value.BinRulesControlsInfo={};
                        value.BinRulesControlLinksID={};
                        value.BinRulesControlNodesID={};
                        value.BinRulesCount=0;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end

                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    clear atline;
                    atline = checktlines(tline);
                    if strcmpi(atline{1}, {'RULE'})
                        d=1;value.BinRulesCount=value.BinRulesCount+1;
                    end
                    value.BinRulesControlsInfo{value.BinRulesCount}{d}=atline;
                    if sum(strcmpi(atline{2}, {'LINK', 'PIPE', 'PUMP', 'VALVE'}))
                        value.BinRulesControlLinksID{value.BinRulesCount}{d}=atline{3};
                    elseif sum(strcmpi(atline{2}, {'NODE', 'JUNCTION', 'RESERVOIR', 'TANK'}))
                        value.BinRulesControlNodesID{value.BinRulesCount}{d}=atline{3};
                    end
                    d=d+1;
                end
            end
        end
        function value = getBinOptionsInfo(obj)
            % Open epanet input file
            [~, info] = obj.readInpFile;
            value.BinUnits_SI_Metric=0; value.BinUnits_US_Customary=0;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                    % [OPTIONS] section
                    if strcmpi(tok(1:5), '[OPTI')
                        sect=1;
                        value.BinLinkFlowUnits={};
                        value.BinOptionsHeadloss={};
                        value.BinNodePressureUnits={};
                        value.BinOptionsSpecificGravity=[];
                        value.BinOptionsViscosity=[];
                        value.BinOptionsMaxTrials=[];
                        value.BinOptionsAccuracyValue=[];
                        value.BinOptionsUnbalanced={};
                        value.BinOptionsPattern=[];
                        value.BinOptionsPatternDemandMultiplier=[];
                        value.BinOptionsEmitterExponent=[];
                        value.BinQualityType={};
                        value.BinQualityCode=[];
                        value.BinQualityTraceNodeIndex=[];
                        value.BinQualityTraceNodeID={};
                        value.BinQualityUnits={};
                        value.BinOptionsDiffusivity=[];
                        value.BinOptionsQualityTolerance=[];
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                if sect==0
                    continue;
                    % Options
                elseif sect==1
                    clear atline;
                    atline = checktlines(tline);
                    value = getOptionsValues(obj, value, atline);
                end
            end
            % US Customary - SI metric
            if find(strcmp(value.BinLinkFlowUnits, obj.TYPEUNITS))<6
                value.BinUnits_US_Customary=1;
            else
                value.BinUnits_SI_Metric=1;
            end
        end
        function value = getBinUnits(obj)
            value.BinLinkFlowUnits=obj.getBinOptionsInfo.BinLinkFlowUnits;
            if obj.getBinOptionsInfo.BinUnits_US_Customary==1;
                value.BinQualityUnits='mg/L or ug/L';
                value.BinNodePressureUnits='psi';
                value.BinPatternDemandsUnits=value.BinLinkFlowUnits;
                value.BinLinkPipeDiameterUnits='inches';
                value.BinNodeTankDiameterUnits='feet';
                value.BinEnergyEfficiencyUnits='percent';
                value.BinNodeElevationUnits='feet';
                value.BinNodeDemandUnits=value.BinLinkFlowUnits;
                value.BinNodeEmitterCoefficientUnits='flow units @ 1 psi drop';
                value.BinEnergyUnits='kwatt-hours';
                value.BinLinkFrictionFactorUnits='unitless';
                value.BinNodeHeadUnits='feet';
                value.BinLinkLengthsUnits='feet';
                value.BinLinkMinorLossCoeffUnits='unitless';
                value.BinLinkPumpPowerUnits='horsepower';
                value.BinQualityReactionCoeffBulkUnits='1/day (1st-order)';
                value.BinQualityReactionCoeffWallUnits='mass/sq-ft/day (0-order), ft/day (1st-order)';
                value.BinLinkPipeRoughnessCoeffUnits='millifeet(Darcy-Weisbach), unitless otherwise';
                value.BinQualitySourceMassInjectionUnits='mass/minute';
                value.BinLinkVelocityUnits='ft/sec';
                value.BinNodeTankVolumeUnits='cubic feet';
                value.BinQualityWaterAgeUnits='hours';
            else % SI Metric
                value.BinQualityUnits='mg/L or ug/L';
                value.BinNodePressureUnits='meters';
                value.BinPatternDemandsUnits=value.BinLinkFlowUnits;
                value.BinLinkPipeDiameterUnits='millimeters';
                value.BinNodeTankDiameterUnits='meters';
                value.BinEnergyEfficiencyUnits='percent';
                value.BinNodeElevationUnits='meters';
                value.BinNodeDemandUnits=value.BinLinkFlowUnits;
                value.BinNodeEmitterCoefficientUnits='flow units @ 1 meter drop';
                value.BinEnergyUnits='kwatt-hours';
                value.BinLinkFrictionFactorUnits='unitless';
                value.BinNodeHeadUnits='meters';
                value.BinLinkLengthsUnits='meters';
                value.BinLinkMinorLossCoeffUnits='unitless';
                value.BinLinkPumpPowerUnits='kwatts';
                value.BinQualityReactionCoeffBulkUnits='1/day (1st-order)';
                value.BinQualityReactionCoeffWallUnits='mass/sq-m/day(0-order), meters/day (1st-order)';
                value.BinLinkPipeRoughnessCoeffUnits='mm(Darcy-Weisbach), unitless otherwise';
                value.BinQualitySourceMassInjectionUnits='mass/minute';
                value.BinLinkVelocityUnits='meters/sec';
                value.BinNodeTankVolumeUnits='cubic meters';
                value.BinQualityWaterAgeUnits='hours';
            end
        end
        function value = getBinTimesInfo(obj)
            % Open epanet input file
            [~, info] = obj.readInpFile;
            sect=0;
            for h=1:length(info)
                tline = info{h};
                if ~ischar(tline),   break,   end
                % Get first token in the line
                tok = strtok(tline);
                % Skip blank Clines and comments
                if isempty(tok), continue, end
                if (tok(1) == ';'), continue, end
                if (tok(1) == '[')
                        % [TIMES] section
                    if strcmpi(tok(1:5), '[TIME')
                        sect=1;
                        value.BinTimeSimulationDuration=[];
                        value.BinTimeHydraulicStep=[];
                        value.BinTimeQualityStep=[];
                        value.BinTimePatternStep=[];
                        value.BinTimePatternStart=[];
                        value.BinTimeReportingStep=[];
                        value.BinTimeReportingStart=[];
                        value.BinTimeStatisticsIndex=[];
                        value.BinTimeStatistics='';
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4), '[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end

                if sect==0
                    continue;
                    % Times
                elseif sect==1
                    clear atline;
                    atline = checktlines(tline);
                    r=atline{2};
                    if length(atline)>2
                        if find(~strcmpi(atline{end}, {'HOURS', 'MIN', 'SECONDS', 'MINUTES', 'DAYS'})==0)
                            r=atline{end-1};
                        else
                            r=atline{end};
                        end
                    end
                    value = getTimes(obj, r, atline, value);
                end
            end
        end
        function value = getCMDCODE(obj)
            value = obj.CMDCODE;
        end
        function setCMDCODE(obj, code)
            value = code;
            if code~=0 && code~=1,  value = obj.CMDCODE; end
            obj.CMDCODE = value;
        end
    end
end
function [Errcode, value] = ENgetnodevalue(index, paramcode, LibEPANET)
  value=single(0);
  index=int32(index);
  paramcode=int32(paramcode);
  [Errcode, value]=calllib(LibEPANET, 'ENgetnodevalue', index, paramcode, value);
  if Errcode==240, value=NaN; end
  value = double(value);
end
function [Errcode, value] = ENgetbasedemand(index, numdemands, LibEPANET)
%epanet20100
[Errcode, value]=calllib(LibEPANET, 'ENgetbasedemand', index, numdemands, 0);
end


function [Errcode, len] = ENgetpatternlen(index, LibEPANET)
[Errcode, len]=calllib(LibEPANET, 'ENgetpatternlen', index, 0);
end
function [Errcode, value] = ENgetpatternvalue(index, period, LibEPANET)
[Errcode, value]=calllib(LibEPANET, 'ENgetpatternvalue', index, period, 0);
end
function [Errcode, qualcode, tracenode] = ENgetqualtype(LibEPANET)
[Errcode, qualcode, tracenode]=calllib(LibEPANET, 'ENgetqualtype', 0, 0);
end
function [Errcode, timevalue] = ENgettimeparam(paramcode, LibEPANET)
[Errcode, timevalue]=calllib(LibEPANET, 'ENgettimeparam', paramcode, 0);
end
function [Errcode, LibEPANET] = ENgetversion(LibEPANET)
[Errcode, LibEPANET]=calllib(LibEPANET, 'ENgetversion', 0);
if Errcode
      obj.ENgeterror(Errcode, LibEPANET);
end
end
function [Errcode] = ENinit(LibEPANET, unitsType, headLossType)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENinit', '', '', unitsType, headLossType);
end
function [Errcode] = ENinitH(flag, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENinitH', flag);
end
function [Errcode] = ENinitQ(saveflag, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENinitQ', saveflag);
end
function ENMatlabCleanup(LibEPANET)
% Load library
if libisloaded(LibEPANET)
    unloadlibrary(LibEPANET);
else
    errstring =['Library ', LibEPANET, '.dll was not loaded.'];
    disp(errstring);
end
end
function ENLoadLibrary(LibEPANETpath, LibEPANET, varargin)
if ~libisloaded(LibEPANET)
    warning('off', 'MATLAB:loadlibrary:TypeNotFound');
    if ~isdeployed
        if isunix
            loadlibrary(LibEPANET, [LibEPANETpath, LibEPANET, '.h']);
        else
            loadlibrary([LibEPANETpath, LibEPANET], [LibEPANETpath, LibEPANET, '.h']);
        end
    else
        loadlibrary('epanet2', @mxepanet); %loadlibrary('epanet2', 'epanet2.h', 'mfilename', 'mxepanet.m');
    end
    warning('on', 'MATLAB:loadlibrary:TypeNotFound');
    if ~isempty(varargin), return; end
end
if libisloaded(LibEPANET)
    [~, version]=calllib(LibEPANET, 'ENgetversion', 0);
    LibEPANETString = ['EPANET version {', num2str(version), '} loaded'];
    fprintf(LibEPANETString);
else
    warning('There was an error loading the EPANET library (DLL).')
end
end
function [Errcode, tstep] = ENnextH(LibEPANET)
[Errcode, tstep]=calllib(LibEPANET, 'ENnextH', int32(0));
end
function [Errcode, tstep] = ENnextQ(LibEPANET)
[Errcode, tstep]=calllib(LibEPANET, 'ENnextQ', int32(0));
tstep = double(tstep);
end
function [Errcode] = ENopen(inpname, repname, binname, LibEPANET) %DE
    Errcode=calllib(LibEPANET, 'ENopen', inpname, repname, binname);
    if Errcode && Errcode~=200
         [~, errmsg] = calllib(LibEPANET, 'ENgeterror', Errcode, char(32*ones(1, 79)), 79);
       disp(errmsg);
    end
end
function [Errcode] = ENepanet(LibEPANET, tempfile, rptfile, binfile)
[Errcode] = calllib(LibEPANET, 'ENepanet', tempfile, rptfile, binfile, lib.pointer);
end
function [Errcode] = ENopenH(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENopenH');
end
function [Errcode] = ENopenQ(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENopenQ');
end
function [Errcode] = ENreport(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENreport');
end
function [Errcode] = ENcopyreport(filename, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENcopyreport', filename);
end
function [Errcode] = ENclearreport(LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENclearreport');
end
function [Errcode, value] = ENgetresultindex(LibEPANET, objecttype, index)
% EPANET Version 2.2
[Errcode, value]=calllib(LibEPANET, 'ENgetresultindex', objecttype, index, int32(0));
end
function [Errcode] = ENresetreport(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENresetreport');
end
function [Errcode, t] = ENrunH(LibEPANET)
[Errcode, t]=calllib(LibEPANET, 'ENrunH', int32(0));
t = double(t);
end
function [Errcode, t] = ENrunQ(LibEPANET)
t=int32(0);
[Errcode, t]=calllib(LibEPANET, 'ENrunQ', t);
end
function [Errcode] = ENsaveH(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsaveH');
end
function [Errcode] = ENsavehydfile(fname, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsavehydfile', fname);
end
function [Errcode] = ENsaveinpfile(inpname, LibEPANET)
Errcode=calllib(LibEPANET, 'ENsaveinpfile', inpname);
end
function [Errcode] = ENsetcontrol(cindex, ctype, lindex, setting, nindex, level, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetcontrol', cindex, ctype, lindex, setting, nindex, level);
end
function [Errcode] = ENaddrule(rule, LibEPANET)
% EPANET Version 2.2
[Errcode, ~]=calllib(LibEPANET, 'ENaddrule', rule);
end
function [Errcode] = ENdeleterule(index, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENdeleterule', index);
end
function [Errcode, logop, object, objIndex, variable, relop, status, value] = ENgetpremise(ruleIndex, premiseIndex, LibEPANET)
% EPANET Version 2.2
[Errcode, logop, object, objIndex, variable, relop, status, value]=calllib(LibEPANET, 'ENgetpremise', ruleIndex, premiseIndex,0,0,0,0,0,0,0);
end
function [Errcode, linkIndex, status, setting] = ENgetthenaction(ruleIndex, actionIndex, LibEPANET)
% EPANET Version 2.2
[Errcode, linkIndex, status, setting]=calllib(LibEPANET, 'ENgetthenaction', ruleIndex, actionIndex,0,0,0);
end
function [Errcode] = ENsetthenaction(ruleIndex, actionIndex, linkIndex, status, setting, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetthenaction', ruleIndex, actionIndex, linkIndex, status, setting);
end
function [Errcode, linkIndex, status, setting] = ENgetelseaction(ruleIndex, actionIndex, LibEPANET)
% EPANET Version 2.2
[Errcode, linkIndex, status, setting]=calllib(LibEPANET, 'ENgetelseaction', ruleIndex, actionIndex,0,0,0);
end
function [Errcode] = ENsetelseaction(ruleIndex, actionIndex, linkIndex, status, setting, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetelseaction', ruleIndex, actionIndex, linkIndex, status, setting);
end
function [Errcode] = ENsetrulepriority(ruleIndex, priority, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetrulepriority', ruleIndex, priority);
end
function [Errcode] = ENsetpremise(ruleIndex, premiseIndex, logop, object, objIndex, variable, relop, status, value, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetpremise', ruleIndex, premiseIndex, logop, object, objIndex, variable, relop, status, value);
end
function [Errcode] = ENsetpremiseindex(ruleIndex, premiseIndex, objIndex, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetpremiseindex', ruleIndex, premiseIndex, objIndex);
end
function [Errcode] = ENsetpremisestatus(ruleIndex, premiseIndex, status, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetpremisestatus', ruleIndex, premiseIndex, status);
end
function [Errcode] = ENsetpremisevalue(ruleIndex, premiseIndex, value, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetpremisevalue', ruleIndex, premiseIndex, value);
end
function [Errcode, id] = ENgetruleID(index, LibEPANET)
% EPANET Version 2.2
[Errcode, id]=calllib(LibEPANET, 'ENgetruleID', index,'');
end
function [Errcode, nPremises, nThenActions, nElseActions, priority] = ENgetrule(index, LibEPANET)
% EPANET Version 2.2
[Errcode, nPremises, nThenActions, nElseActions, priority]=calllib(LibEPANET, 'ENgetrule', index, 0, 0, 0, 0);
end
function [Errcode, index] = ENsetlinknodes(index, startnode, endnode, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetlinknodes', index, startnode, endnode);
end
function [Errcode] = ENsetlinkvalue(index, paramcode, value, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetlinkvalue', index, paramcode, value);
end
function [Errcode, index] = ENsetpipedata(index, length, diam, rough, mloss, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetpipedata', index, length, diam, rough, mloss);
end
function [Errcode, index] = ENsetlinktype(id, paramcode, actionCode, LibEPANET)
% EPANET Version 2.2
[Errcode, index]=calllib(LibEPANET, 'ENsetlinktype', id, paramcode, actionCode);
end
function [Errcode] = ENsetnodevalue(index, paramcode, value, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetnodevalue', index, paramcode, value);
end
function [Errcode] = ENsettankdata(index, elev, initlvl, minlvl, maxlvl, diam, minvol, volcurve, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsettankdata', index, elev, initlvl, minlvl, maxlvl, diam, minvol, volcurve );
end
function [Errcode] = ENsetoption(optioncode, value, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetoption', optioncode, value);
end
function [Errcode] = ENsetpattern(index, factors, nfactors, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetpattern', index, factors, nfactors);
end
function [Errcode] = ENsetpatternvalue(index, period, value, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetpatternvalue', index, period, value);
end
function [Errcode] = ENsetqualtype(qualcode, chemname, chemunits, tracenode, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetqualtype', qualcode, chemname, chemunits, tracenode);
end
function [Errcode] = ENsetreport(command, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetreport', command);
end
function [Errcode] = ENsetstatusreport(statuslevel, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsetstatusreport', statuslevel);
end
function [Errcode] = ENsettimeparam(paramcode, timevalue, LibEPANET)
paramcode=int32(paramcode);
timevalue=int32(timevalue);
[Errcode]=calllib(LibEPANET, 'ENsettimeparam', paramcode, timevalue);
end
function [Errcode] = ENsolveH(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsolveH');
end
function [Errcode] = ENsolveQ(LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENsolveQ');
end
function [Errcode, tleft] = ENstepQ(LibEPANET)
tleft=int32(0);
[Errcode, tleft]=calllib(LibEPANET, 'ENstepQ', tleft);
tleft=double(tleft);
end
function [Errcode] = ENusehydfile(hydfname, LibEPANET)
[Errcode]=calllib(LibEPANET, 'ENusehydfile', hydfname);
end
function [Errcode] = ENsetcurve(index, x, y, nfactors, LibEPANET)
% EPANET Version 2.1
[Errcode]=calllib(LibEPANET, 'ENsetcurve', index, x, y, nfactors);
end
function [Errcode, x, y] = ENgetcurvevalue(index, period, LibEPANET)
% EPANET Version 2.1
[Errcode, x, y]=calllib(LibEPANET, 'ENgetcurvevalue', index, period, 0, 0);
end
function [Errcode, x, y] = ENsetcurvevalue(index, pnt, x, y, LibEPANET)
% EPANET Version 2.1
% index  = curve index
% pnt    = curve's point number
% x      = curve x value
% y      = curve y value
% sets x, y point for a specific point and curve
[Errcode]=calllib(LibEPANET, 'ENsetcurvevalue', index, pnt, x, y);
end
function [Errcode, index] = ENgetcurveindex(id, LibEPANET)
% EPANET Version 2.1
[Errcode, ~, index]=calllib(LibEPANET, 'ENgetcurveindex', id, 0);
end
function [Errcode] = ENaddcurve(cid, LibEPANET)
% EPANET Version 2.1
[Errcode]=calllib(LibEPANET, 'ENaddcurve', cid);
end
function [Errcode, ids, nvalue, xvalue, yvalue] = ENgetcurve(obj, value, LibEPANET)
[Errcode, ids, nvalue, xvalue, yvalue]=calllib(LibEPANET, 'ENgetcurve', value, char(32*ones(1, 31)), 0, zeros(1, obj.getCurveLengths(value))', zeros(1, obj.getCurveLengths(value))');
end
function [Errcode, len] = ENgetcurvelen(index, LibEPANET)
% EPANET Version 2.1
[Errcode, len]=calllib(LibEPANET, 'ENgetcurvelen', index, 0);
end
function [Errcode, value] = ENgetheadcurveindex(pumpindex, LibEPANET)
% EPANET Version 2.1
[Errcode, value]=calllib(LibEPANET, 'ENgetheadcurveindex', pumpindex, 0);
end
function [Errcode, value] = ENgetpumptype(pumpindex, LibEPANET)
% EPANET Version 2.1
[Errcode, value]=calllib(LibEPANET, 'ENgetpumptype', pumpindex, 0);
end
function [Errcode, value] = ENgetaveragepatternvalue(index, LibEPANET)
% return  average pattern value
% EPANET Version 2.1
[Errcode, value]=calllib(LibEPANET, 'ENgetaveragepatternvalue', index, 0);
end
function [Errcode, x, y] = ENgetcoord(index, LibEPANET)
% EPANET Version 2.1
[Errcode, x, y]=calllib(LibEPANET, 'ENgetcoord', index, 0, 0);
end
function [Errcode] = ENsetcoord(index, x, y, LibEPANET)
% EPANET Version 2.1
[Errcode]=calllib(LibEPANET, 'ENsetcoord', index, x, y);
end
function [Errcode, x, y] = ENgetvertex(index, vertex, LibEPANET)
% EPANET Version 2.2
[Errcode, x, y]=calllib(LibEPANET, 'ENgetvertex', index, vertex, 0, 0);
end
function [Errcode] = ENsetvertices(index, vertex, x, y, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetvertices', index, vertex, x, y);
end
function [Errcode, count] = ENgetvertexcount(index, LibEPANET)
[Errcode, count]=calllib(LibEPANET, 'ENgetvertexcount', index, 0);
end
function [Errcode] = ENadddemand(nodeIndex, baseDemand, demandPattern, demandName, LibEPANET)
% EPANET Version 2.2
[Errcode, ~, ~]=calllib(LibEPANET, 'ENadddemand', nodeIndex, baseDemand , demandPattern, demandName);
end
function [Errcode] = ENdeletedemand(nodeIndex, demandIndex, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENdeletedemand', nodeIndex, demandIndex);
end
function [Errcode] = ENsetbasedemand(index, demandIdx, value, LibEPANET)
% EPANET Version 2.1
[Errcode]=calllib(LibEPANET, 'ENsetbasedemand', index, demandIdx, value);
end
function [Errcode] = ENsetdemandpattern(index, demandIdx, patInd, LibEPANET)
% New version
[Errcode]=calllib(LibEPANET, 'ENsetdemandpattern', index, demandIdx, patInd);
end
function [Errcode, qualcode, chemname, chemunits, tracenode] = ENgetqualinfo(LibEPANET)
chm=char(32*ones(1, 31));
[Errcode, qualcode, chemname, chemunits, tracenode]=calllib(LibEPANET, 'ENgetqualinfo', 0, chm, chm, 0);
end
function [index, Errcode] = ENaddnode(obj, nodeid, nodetype)
% dev-net-builder
[Errcode, ~, index]=calllib(obj.LibEPANET, 'ENaddnode', nodeid, nodetype, 0);
error(obj.getError(Errcode));
end
function [index, Errcode] = ENaddlink(obj, linkid, linktype, fromnode, tonode)
% dev-net-builder
[Errcode, ~, ~, ~, index]=calllib(obj.LibEPANET, 'ENaddlink', linkid, linktype, fromnode, tonode, 0);
error(obj.getError(Errcode));
end
function [Errcode] = ENdeletenode(LibEPANET, indexNode, condition)
% dev-net-builder
[Errcode]=calllib(LibEPANET, 'ENdeletenode', indexNode, condition);
end
function [Errcode] = ENdeletelink(LibEPANET, indexLink, condition)
% dev-net-builder
[Errcode]=calllib(LibEPANET, 'ENdeletelink', indexLink, condition);
end
function [Errcode] = ENsetheadcurveindex(LibEPANET, pumpindex, curveindex)
% dev-net-builder
[Errcode]=calllib(LibEPANET, 'ENsetheadcurveindex', pumpindex, curveindex);
end
function [Errcode, cindex] = ENaddcontrol(ctype, lindex, setting, nindex, level, LibEPANET)
[Errcode, cindex]=calllib(LibEPANET, 'ENaddcontrol', ctype, lindex, setting, nindex, level, 0);
end
function [Errcode] = ENdeletecontrol(index, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENdeletecontrol', index);
end
function [Errcode, type, pmin, preq, pexp] = ENgetdemandmodel(LibEPANET)
% EPANET Version 2.2
[Errcode, type, pmin, preq, pexp]=calllib(LibEPANET, 'ENgetdemandmodel', 0, 0, 0, 0);
end
function [Errcode] = ENsetdemandmodel(type, pmin, preq, pexp, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetdemandmodel', type, pmin, preq, pexp);
end
function [Errcode] = ENsetdemandname(node_index, demand_index, demand_name, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetdemandname', node_index, demand_index, demand_name);
end
function [Errcode, demand_name] = ENgetdemandname(node_index, demand_index, LibEPANET)
% EPANET Version 2.2
demand_name = char(32*ones(1, 31));
[Errcode, demand_name]=calllib(LibEPANET, 'ENgetdemandname', node_index, demand_index, demand_name);
end
function [Errcode, line1, line2, line3] = ENgettitle(LibEPANET)
% EPANET Version 2.2
c = char(32*ones(1, 79));
[Errcode, line1, line2, line3]=calllib(LibEPANET, 'ENgettitle', c, c, c);
end
function [Errcode] = ENsettitle(line1, line2, line3, LibEPANET)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsettitle', line1, line2, line3);
end
function [Errcode] = ENsetcomment(object, index, comment, LibEPANET)
% EPANET Version 2.2
% Object a type of object (either EN_NODE, EN_LINK, EN_TIMEPAT or EN_CURVE)
[Errcode]=calllib(LibEPANET, 'ENsetcomment', object, index, comment);
end
function [Errcode, comment] = ENgetcomment(object, index, LibEPANET)
% EPANET Version 2.2
comment = char(32*ones(1, 79));
[Errcode, comment]=calllib(LibEPANET, 'ENgetcomment', object, index, comment);
end
function [Errcode] = ENdeletepattern(LibEPANET, indexPat)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENdeletepattern', indexPat);
end
function [Errcode] = ENdeletecurve(LibEPANET, indexCurve)
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENdeletecurve', indexCurve);
end
function [Errcode] = ENsetjuncdata(LibEPANET, index, elev, dmnd, dmndpat)
%   @brief Sets a group of properties for a junction node.
%   @param index a junction node's index (starting from 1).
%   @param elev the value of the junction's elevation.
%   @param dmnd the value of the junction's primary base demand.
%   @param dmndpat the ID name of the demand's time pattern ("" for no pattern)
%   @return an error code.
% EPANET Version 2.2
[Errcode]=calllib(LibEPANET, 'ENsetjuncdata', index, elev, dmnd, dmndpat);
end

% function [Errcode, nPremises, nTrueActions, nFalseActions, priority] = EN_getrule(cindex, LibEPANET)
%     [Errcode, nPremises, nTrueActions, nFalseActions, priority]=calllib(LibEPANET, 'ENgetrule', cindex, 0, 0, 0, 0);
%     if Errcode
  %         obj.ENgeterror(Errcode, LibEPANET);
%     end
% end
function [Errcode] = ENsetflowunits(LibEPANET, code)
[Errcode]=calllib(LibEPANET, 'ENsetflowunits', code);
end
function [obj] = MSXMatlabSetup(obj, msxname, varargin)
arch = computer('arch');
pwdepanet = fileparts(which(mfilename));
if strcmp(arch, 'win64')
    obj.MSXLibEPANETPath = [pwdepanet, '\64bit\'];
elseif strcmp(arch, 'win32')
    obj.MSXLibEPANETPath = [pwdepanet, '\32bit\'];
end
if isunix
    obj.MSXLibEPANETPath = [pwdepanet, '/glnx/'];
end
if ~isempty(varargin)
    if varargin{1}{1}~=1
        if nargin==3
            obj.MSXLibEPANETPath=char(varargin{1});
            obj.MSXLibEPANETPath=[fileparts(obj.MSXLibEPANETPath), '\'];
            if isempty(varargin{1})
                obj.MSXLibEPANETPath='';
            end
        end
    end
end
obj.MSXLibEPANET='epanetmsx'; % Get DLL LibEPANET (e.g. epanet20012x86 for 32-bit)
if ~libisloaded(obj.MSXLibEPANET)
    loadlibrary([obj.MSXLibEPANETPath, obj.MSXLibEPANET], [obj.MSXLibEPANETPath, [obj.MSXLibEPANET, '.h']]);
end

obj.MSXFile = which(char(msxname));
%Save the temporary msx file
mm=0;
if ~isempty(varargin)
    if varargin{1}{1}==1
        mm=1; %for set (write) msx functions
    end
end
if mm==1
    if ~iscell(varargin{1})
        obj.MSXTempFile=obj.MSXFile;
    end
else
    obj.MSXTempFile=[obj.MSXFile(1:end-4), '_temp.msx'];
    copyfile(obj.MSXFile, obj.MSXTempFile);
end
%Open the file
[obj.Errcode] = MSXopen(obj);
obj.MSXEquationsTerms = obj.getMSXEquationsTerms;
obj.MSXEquationsPipes = obj.getMSXEquationsPipes;
obj.MSXEquationsTanks = obj.getMSXEquationsTanks;

obj.MSXSpeciesCount = obj.getMSXSpeciesCount;
obj.MSXConstantsCount = obj.getMSXConstantsCount;
obj.MSXParametersCount = obj.getMSXParametersCount;
obj.MSXPatternsCount = obj.getMSXPatternsCount;
obj.MSXSpeciesIndex = obj.getMSXSpeciesIndex;
obj.MSXSpeciesNameID = obj.getMSXSpeciesNameID;
obj.MSXSpeciesType = obj.getMSXSpeciesType;
obj.MSXSpeciesUnits = obj.getMSXSpeciesUnits;
obj.MSXSpeciesATOL = obj.getMSXSpeciesATOL;
obj.MSXSpeciesRTOL = obj.getMSXSpeciesRTOL;
obj.MSXConstantsNameID = obj.getMSXConstantsNameID;
obj.MSXConstantsValue  = obj.getMSXConstantsValue;
obj.MSXConstantsIndex = obj.getMSXConstantsIndex;
obj.MSXParametersNameID = obj.getMSXParametersNameID;
obj.MSXParametersIndex = obj.getMSXParametersIndex;
obj.MSXParametersTanksValue = obj.getMSXParametersTanksValue;
obj.MSXParametersPipesValue = obj.getMSXParametersPipesValue;
obj.MSXPatternsNameID = obj.getMSXPatternsNameID;
obj.MSXPatternsIndex = obj.getMSXPatternsIndex;
obj.MSXPatternsLengths = obj.getMSXPatternsLengths;
obj.MSXNodeInitqualValue = obj.getMSXNodeInitqualValue;
obj.MSXLinkInitqualValue = obj.getMSXLinkInitqualValue;
obj.MSXSources = obj.getMSXSources;
obj.MSXSourceType = obj.getMSXSourceType;
obj.MSXSourceLevel = obj.getMSXSourceLevel;
obj.MSXSourcePatternIndex = obj.getMSXSourcePatternIndex;
obj.MSXSourceNodeNameID = obj.getMSXSourceNodeNameID;
obj.MSXPattern = obj.getMSXPattern;

end
function [Errcode] = MSXopen(obj, varargin)
[Errcode] = calllib(obj.MSXLibEPANET, 'MSXopen', obj.MSXTempFile);
if Errcode
    MSXerror(Errcode, obj.MSXLibEPANET);
end
if (Errcode == 520)
    disp('current MSX project will be closed and the new project will be opened');
    [Errcode] = MSXclose(obj);
    if Errcode
        MSXerror(Errcode, obj.MSXLibEPANET);
    else
        [Errcode] = calllib(obj.MSXLibEPANET, 'MSXopen', obj.MSXTempFile);
        if Errcode
            MSXerror(Errcode, obj.MSXLibEPANET);
        end
    end
end
end
function [Errcode] = MSXclose(obj)
[Errcode] = calllib(obj.MSXLibEPANET, 'MSXclose');
if Errcode
    MSXerror(Errcode, obj.MSXLibEPANET);
end
end
function [e] = MSXerror(Errcode, MSXLibEPANET)
len=80;
errstring=char(32*ones(1, len+1));
[e, errstring] = calllib(MSXLibEPANET, 'MSXgeterror', Errcode, errstring, len);
disp(errstring);
end
function [Errcode, count] = MSXgetcount(code, MSXLibEPANET)
count=0;
[Errcode, count] = calllib(MSXLibEPANET, 'MSXgetcount', code, count);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, index] = MSXgetindex(MSXLibEPANET, varargin)
index =0;
if ~isnumeric(varargin{1})
    varargin{1}=varargin{2};
    varargin{2}=varargin{3};
end
[Errcode, ~, index]=calllib(MSXLibEPANET, 'MSXgetindex', varargin{1}, varargin{2}, index);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, id] = MSXgetID(type, index, len, MSXLibEPANET)
id=char(32*ones(1, len+1));
[Errcode, id]=calllib(MSXLibEPANET, 'MSXgetID', type, index, id, len);
id=id(1:len);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, len] = MSXgetIDlen(type, index, MSXLibEPANET)
len=0;
[Errcode, len]=calllib(MSXLibEPANET, 'MSXgetIDlen', type, index, len);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, type, units, atol, rtol] = MSXgetspecies(index, MSXLibEPANET)
type=0; rtol=0; atol=0;
units=char(32*ones(1, 16));
[Errcode, type, units, atol, rtol]=calllib(MSXLibEPANET, 'MSXgetspecies', index, type, units, atol, rtol);
switch type
    case 0
        type='BULK';   % for a bulk water species
    case 1
        type='WALL';   % for a pipe wall surface species
end
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetconstant(index, MSXLibEPANET)
value=0;
[Errcode, value]=calllib(MSXLibEPANET, 'MSXgetconstant', index, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetparameter(type, index, param, MSXLibEPANET)
value=0;
[Errcode, value]=calllib(MSXLibEPANET, 'MSXgetparameter', type, index, param, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, patlen] = MSXgetpatternlen(patindex, MSXLibEPANET)
patlen=0;
[Errcode, patlen]=calllib(MSXLibEPANET, 'MSXgetpatternlen', patindex, patlen);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetpatternvalue(patindex, period, MSXLibEPANET)
value=0;
[Errcode, value]=calllib(MSXLibEPANET, 'MSXgetpatternvalue', patindex, period, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetinitqual(obj, index, species, MSXLibEPANET)
value=0;
[Errcode, value]=calllib(MSXLibEPANET, 'MSXgetinitqual', obj, index, species, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, type, level, pat] = MSXgetsource(node, species, MSXLibEPANET)
type=0;
level=0;
pat=0;
[Errcode, type, level, pat]=calllib(MSXLibEPANET, 'MSXgetsource', node, species, type, level, pat);
switch type     % type codes
    case -1
        type='NOSOURCE';  % for no source
    case 0
        type='CONCEN';    % for a concentration source
    case 1
        type='MASS';      % for a mass booster source
    case 2
        type='SETPOINT';  % for a setpoint source
    case 3
        type='FLOWPACED'; % for a flow paced source
end
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function MSXMatlabCleanup(obj)
% Unload library
if libisloaded(obj.MSXLibEPANET)
    unloadlibrary(obj.MSXLibEPANET);
else
    errstring =['Library ', obj.MSXLibEPANET, '.dll was not loaded'];
    disp(errstring);
end
end
function [Errcode] = MSXsaveoutfile(outfname, MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET, 'MSXsaveoutfile', outfname);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsavemsxfile(msxname, MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET, 'MSXsavemsxfile', msxname);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetconstant(index, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetconstant', index, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetparameter(type, index, param, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetparameter', type, index, param, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetinitqual(type, index, species, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetinitqual', type, index, species, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetpattern(index, factors, nfactors, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetpattern', index, factors, nfactors);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetpatternvalue(pat, period, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetpatternvalue', pat, period, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsolveQ(MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsolveQ');
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsolveH(MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsolveH');
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXaddpattern(patid, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXaddpattern', patid);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXusehydfile(hydfname, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXusehydfile', hydfname);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode, t, tleft] = MSXstep(MSXLibEPANET)
t=int32(0);
tleft=int32(0);
[Errcode, t, tleft]=calllib(MSXLibEPANET, 'MSXstep', t, tleft);
t = double(t);
tleft = double(tleft);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXinit(flag, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXinit', flag);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXreport(MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET, 'MSXreport');
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [e, errmsg] = MSXgeterror(Errcode, MSXLibEPANET)
errmsg = char(32*ones(1, 80));
[e, errmsg] = calllib(MSXLibEPANET, 'MSXgeterror', Errcode, errmsg, 80);
if e
    MSXerror(e);
end
end
function [Errcode, value] = MSXgetqual(type, index, species, MSXLibEPANET)
value=0;
[Errcode, value]=calllib(MSXLibEPANET, 'MSXgetqual', type, index, species, value);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Errcode] = MSXsetsource(node, species, type, level, pat, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET, 'MSXsetsource', node, species, type, level, pat);
if Errcode
    MSXerror(Errcode, MSXLibEPANET);
end
end
function [Terms, Pipes, Tanks] = getEquations(msxname)
% Open epanet input file
[fid, message]=fopen(msxname, 'rt');
if fid < 0
    disp(message)
    return
end

Terms={};
Pipes={};
Tanks={};

sect=0; i=1; t=1; k=1;
% Read each line from msx file.
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end

    % Get first token in the line
    tok = strtok(tline);

    % Skip blank Clines and comments
    if isempty(tok), continue, end
    if (tok(1) == ';'), continue, end

    if (tok(1) == '[')
        % [TERMS] section
        if strcmpi(tok(1:5), '[TERM')
            sect = 1;
            continue;
            % [PIPES] section
        elseif strcmpi(tok(1:5), '[PIPE')
            sect = 2;
            continue;
            % [TANKS]
        elseif strcmpi(tok(1:5), '[TANK')
            sect = 3;
            continue;
            % [END]
        elseif strcmpi(tok(1:4), '[END')
            break;
        else
            sect = 0;
            continue;
        end
    end

    if sect == 0
        continue;

        % Terms
    elseif sect == 1
        Terms{i} = tline;
        i=i+1;
        % Pipes
    elseif sect == 2
        Pipes{t} = tline;
        t=t+1;
        % Tanks
    elseif sect == 3
        Tanks{k} = tline;
        k=k+1;
    end
    %
end
fclose(fid);
end
function value = get_MSX_Options(msxname, param, getall)

if isempty(msxname)
    warning('Please load MSX File.');
    return;
end

% Open epanet input file
[fid, message] = fopen(msxname, 'rt');
if fid < 0
    disp(message)
    return
end
% DEFAULT OPTIONS
value.AreaUnits='FT2';
value.RateUnits='HR';
value.Solver='EUL';
value.TimeStep=300;
value.Atol=0.01;
value.Rtol=0.001;
value.Coupling='NONE';
value.Compiler='NONE';
sect=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    % Get first token in the line
    tok = strtok(tline);
    % Skip blank Clines and comments
    if isempty(tok), continue, end
    if (tok(1) == ';'), continue, end
    if (tok(1) == '[')
        % [OPTIONS] section
        if strcmpi(tok(1:5), '[OPTI')
            sect=1;
            continue;
            % [END]
        elseif strcmpi(tok(1:4), '[REP')
            break;
        else
            sect = 0;
            continue;
        end
    end

    if sect == 0
        continue;

        % Options
    elseif sect == 1
        atline = checktlines(tline);
        res = atline{2};
        if strcmpi(atline{1}, param) || (isempty(param))
            switch atline{1}
                case 'TIMESTEP'
                   value.TimeStep = str2double(res);
                case 'AREA_UNITS'
                    value.AreaUnits=res;
                case 'RATE_UNITS'
                    value.RateUnits=res;
                case 'SOLVER'
                    value.Solver=res;
                case 'RTOL'
                    value.Rtol=str2double(res);
                case 'ATOL'
                    value.Atol=str2double(res);
                case 'COUPLING'
                    value.Coupling=res;
                case 'COMPILER'
                    value.Compiler=res;
            end
            if (getall == 0 && ~isempty(param))
                fclose(fid);
                return
            end
        end
    end
end
fclose(fid);
end
function [axesid] = plotnet(obj, varargin)
% Initiality
highlightnode=0;
highlightlink=0;
highlightnodeindex=[];
highlightlinkindex=[];
legendIndices=[];
l=zeros(1, 6);
Node=char('no');
Link=char('no');
NodeInd=0;
LinkInd=0;
fontsize=10;
selectColorNode={''};
selectColorLink={''};
axesid=0;
lline='yes';
npoint='yes';
extend='no';
legendposition = 'northeast';
slegend = 'show';
for i=1:(nargin/2)
    argument =lower(varargin{2*(i-1)+1});
    switch argument
        case 'nodes' % Nodes
            if ~strcmpi(varargin{2*i}, 'yes') && ~strcmpi(varargin{2*i}, 'no')
                warning('Invalid argument.');
                return
            end
            Node=varargin{2*i};
        case 'links' % Nodes
            if ~strcmpi(varargin{2*i}, 'yes') && ~strcmpi(varargin{2*i}, 'no')
                warning('Invalid argument.');
                return
            end
            Link=varargin{2*i};
        case 'nodesindex' % Nodes
            if ~strcmpi(varargin{2*i}, 'yes')
                warning('Invalid argument.');
                return
            end
            NodeInd=varargin{2*i};
        case 'linksindex' % Links
            if ~strcmpi(varargin{2*i}, 'yes')
                warning('Invalid argument.');
                return
            end
            LinkInd=varargin{2*i};
        case 'highlightnode' % Highlight Node
            highlightnode=varargin{2*i};
        case 'highlightlink' % Highlight Link
            highlightlink=varargin{2*i};
        case 'fontsize' % font size
            fontsize=varargin{2*i};
        case 'colornode' % color
            selectColorNode=varargin{2*i};
        case 'colorlink' % color
            selectColorLink=varargin{2*i};
        case 'point' % color
            if ~strcmpi(varargin{2*i}, 'yes') && ~strcmpi(varargin{2*i}, 'no')
                warning('Invalid argument.');
                return
            end
            npoint=varargin{2*i};
        case 'line' % remove line
            if ~strcmpi(varargin{2*i}, 'yes') && ~strcmpi(varargin{2*i}, 'no')
                warning('Invalid argument.');
                return
            end
            lline=varargin{2*i};
        case 'axes' % axes id
            try
                axesid=axes('Parent', varargin{2*i});
            catch
                axesid=varargin{2*i};
            end
        case 'uifigure' % figure
            fig=varargin{2*i};
        case 'bin'
            bin=varargin{2*i};
        case 'extend' % extend option
            extend=varargin{2*i};
        case 'legendposition' % extend option
            legendposition=varargin{2*i};
        case 'legend'
            slegend=varargin{2*i};
        otherwise
            error('Invalid property founobj.');
    end
end

if axesid==0
   drawnow;
   fig=figure;
   axesid=axes('Parent', fig);
end

if cellfun('isempty', selectColorNode)==1
    init={'r'};
    for i=1:length(highlightnode)
        selectColorNode=[init selectColorNode];
    end
end
if cellfun('isempty', selectColorLink)==1
    init={'r'};
    for i=1:length(highlightlink)
        selectColorLink=[init selectColorLink];
    end
end

% get info BIN function
if bin==1
    b=obj.getBinLinksInfo;
    v.linknameid=b.BinLinkNameID;
    v.linkcount=b.BinLinkCount;
    v.linkfromnode=b.BinLinkFromNode;
    v.linktonode=b.BinLinkToNode;
    v.linkindex=b.BinLinkIndex;
    v.pumpindex=b.BinLinkPumpIndex;
    v.valveindex=b.BinLinkValveIndex;
    v.nodesconnlinks = [v.linkfromnode;v.linktonode]';
    b=obj.getBinNodesInfo;
    v.nodenameid=b.BinNodeNameID;
    v.nodecoords=obj.getBinNodeCoordinates;
    v.nodecount=b.BinNodeCount;
    v.nodeindex=b.BinNodeIndex;
    v.resindex=b.BinNodeReservoirIndex;
    v.tankindex=b.BinNodeTankIndex;
elseif bin==0
    % get info EN functions
    v.nodenameid=obj.getNodeNameID;
    v.linknameid=obj.getLinkNameID;
    if isempty(v.nodenameid) || isempty(v.linknameid)
        warning('Not enough network nodes/links.');
        return;
    end
    v.nodesconnlinks=obj.getNodesConnectingLinksID;
    if sum(strcmp(obj.libFunctions, 'ENgetcoord'))
        v.nodecoords=obj.getNodeCoordinates;
    else
        v.nodecoords=obj.getBinNodeCoordinates;
    end
    v.pumpindex=obj.getLinkPumpIndex;
    v.valveindex=obj.getLinkValveIndex;
    v.resindex=obj.getNodeReservoirIndex;
    v.tankindex=obj.getNodeTankIndex;
    v.linkcount=obj.getLinkCount;
    v.nodecount=obj.getNodeCount;
    v.linkindex=obj.getLinkIndex;
    v.nodeindex=obj.getNodeIndex;
end

if isnan(v.nodecoords{1}(2))
   warning('Do not exist coordinates.'); close(g);
   return
end
% Get node names and x, y coordiantes
if isa(highlightnode, 'cell')
    for i=1:length(highlightnode)
        n = strcmp(v.nodenameid, highlightnode{i});
        if sum(n)==0
            warning('Undefined node with id "%s" in function call therefore the index is zero.', char(highlightnode{i}));
        else
            highlightnodeindex(i) = strfind(n, 1);
        end
    end
end

if isa(highlightlink, 'cell')
    for i=1:length(highlightlink)
        n = strcmp(v.linknameid, highlightlink{i});
        if sum(n)==0
            warning('Undefined link with id "%s" in function call therefore the index is zero.', char(highlightlink{i}));
        else
            highlightlinkindex(i) = strfind(n, 1);
        end
    end
end

if (strcmpi(lline, 'yes'))
    hold(axesid, 'on')
    for i=1:v.linkcount
        FromNode=strfind(strcmp(v.nodesconnlinks(i, 1), v.nodenameid), 1);
        ToNode=strfind(strcmp(v.nodesconnlinks(i, 2), v.nodenameid), 1);

        if FromNode
            x1 = double(v.nodecoords{1}(FromNode));
            y1 = double(v.nodecoords{2}(FromNode));
        end
        if ToNode
            x2 = double(v.nodecoords{1}(ToNode));
            y2 = double(v.nodecoords{2}(ToNode));
        end

        hh=strfind(highlightlinkindex, i);

        if length(hh) && ~isempty(selectColorLink)
            line([x1 v.nodecoords{3}{i} x2], [y1 v.nodecoords{4}{i} y2], 'LineWidth', 1, 'Color', [.5 .5 .5], 'Parent', axesid);
        end
        if ~length(hh)
            h(:, 1)=line([x1 v.nodecoords{3}{i} x2], [y1 v.nodecoords{4}{i} y2], 'LineWidth', 1, 'Parent', axesid);
            if ~l(1), legendIndices = [legendIndices 1]; l(1)=1; end
        end

        % Plot Pumps
        if sum(strfind(v.pumpindex, i))
            colornode = 'm';
            if length(hh) && isempty(selectColorLink)
                colornode = 'r';
            end
            h(:, 2)=plot((x1+x2)/2, (y1+y2)/2, 'mv', 'LineWidth', 2, 'MarkerEdgeColor', 'm', ...
                'MarkerFaceColor', 'm', ...
                'MarkerSize', 5, 'Parent', axesid);
            if ~l(2), legendIndices = [legendIndices 2]; l(2)=1; end
            plot((x1+x2)/2, (y1+y2)/2, 'mv', 'LineWidth', 2, 'MarkerEdgeColor', colornode, ...
                'MarkerFaceColor', colornode, ...
                'MarkerSize', 5, 'Parent', axesid);
        end

        % Plot Valves
        if sum(strfind(v.valveindex, i))
            colornode = 'k';
            if length(hh) && isempty(selectColorLink)
                colornode = 'r';
            end
            h(:, 3)=plot((x1+x2)/2, (y1+y2)/2, 'k*', 'LineWidth', 2, 'MarkerEdgeColor', colornode, ...
                'MarkerFaceColor', colornode, 'MarkerSize', 7, 'Parent', axesid);
            if ~l(3), legendIndices = [legendIndices 3]; l(3)=1; end
        end

        if length(hh) && isempty(selectColorLink)
            line([x1, x2], [y1, y2], 'LineWidth', 1, 'Color', 'r', 'Parent', axesid);
            text((x1+x2)/2, (y1+y2)/2, v.linknameid(i), 'Fontsize', fontsize, 'Parent', axesid);
        elseif length(hh) && ~isempty(selectColorLink)
            try tt=length(selectColorLink{hh}); catch; tt=2; end
           if tt>1
                if length(selectColorLink(hh))==1
                    nm{1}=selectColorLink(hh);
                else
                    nm=selectColorLink(hh);
                end
                if iscell(nm{1})
                    line([x1 v.nodecoords{3}{i} x2], [y1 v.nodecoords{4}{i} y2], 'LineWidth', 1, 'Color', nm{1}{1}, 'Parent', axesid);
                else
                    line([x1 v.nodecoords{3}{i} x2], [y1 v.nodecoords{4}{i} y2], 'LineWidth', 1, 'Color', nm{1}, 'Parent', axesid);
                end
            else
                line([x1 v.nodecoords{3}{i} x2], [y1 v.nodecoords{4}{i} y2], 'LineWidth', 1, 'Color', char(selectColorLink(hh)), 'Parent', axesid);
            end
        end
        % Show Link id
        if (strcmpi(Link, 'yes')) %&& ~length(hh))
            text((x1+x2)/2, (y1+y2)/2, v.linknameid(i), 'Fontsize', fontsize, 'Parent', axesid);
        end
        % Show Link Index
        if (strcmpi(LinkInd, 'yes')) %&& ~length(hh))
            text((x1+x2)/2, (y1+y2)/2, num2str(v.linkindex(i)), 'Fontsize', fontsize, 'Parent', axesid);
        end
    end
end

if (strcmpi(npoint, 'yes'))
    % Coordinates for node FROM
    hold(axesid, 'on')
    for i=1:v.nodecount
        [x] = double(v.nodecoords{1}(i));
        [y] = double(v.nodecoords{2}(i));

        hh=strfind(highlightnodeindex, i);
        if ~length(hh)
            h(:, 4)=plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'b', ...
            'MarkerFaceColor', 'b', ...
            'MarkerSize', 5, 'Parent', axesid);
            if ~l(4), legendIndices = [legendIndices 4]; l(4)=1; end
        end

        % Plot Reservoirs
        if sum(strfind(v.resindex, i))
            colornode = 'g';
            if length(hh) && isempty(selectColorNode)
                colornode = 'r';
            end
            h(:, 5)=plot(x, y, 's', 'LineWidth', 2, 'MarkerEdgeColor', 'g', ...
                'MarkerFaceColor', 'g', ...
                'MarkerSize', 13, 'Parent', axesid);
            if ~l(5), legendIndices = [legendIndices 5]; l(5)=1; end
            plot(x, y, 's', 'LineWidth', 2, 'MarkerEdgeColor', colornode, ...
                'MarkerFaceColor', colornode, ...
                'MarkerSize', 13, 'Parent', axesid);
        end
        % Plot Tanks
        if sum(strfind(v.tankindex, i))
            colornode='c';
            if length(hh) && isempty(selectColorNode)
                colornode='r';
            elseif length(hh) && ~isempty(selectColorNode)
                colornode= 'c';
            end
            h(:, 6)=plot(x, y, 'p', 'LineWidth', 2, 'MarkerEdgeColor', 'c', ...
                'MarkerFaceColor', 'c', ...
                'MarkerSize', 16, 'Parent', axesid);
            if ~l(6), legendIndices = [legendIndices 6]; l(6)=1; end

            plot(x, y, 'p', 'LineWidth', 2, 'MarkerEdgeColor', colornode, ...
                'MarkerFaceColor', colornode, ...
                'MarkerSize', 16, 'Parent', axesid);
        end

        if length(hh) && isempty(selectColorNode)
            plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'r', ...
                'MarkerFaceColor', 'r', ...
                'MarkerSize', 5, 'Parent', axesid);
            text(x, y, v.nodenameid(i), 'Fontsize', fontsize, 'Parent', axesid)%'BackgroundColor', [.7 .9 .7], 'Margin', margin/4);
        elseif length(hh) && ~isempty(selectColorNode)
            try tt=length(selectColorNode{hh}); catch, tt=2; end
           if tt>1
                if length(selectColorNode(hh))==1
                    nm{1}=selectColorNode(hh);
                    nmplot=nm{1}{1};
                else
                    nm=selectColorNode(hh);
                    nmplot=nm{1};
                end
                if iscell(nm{1})
                    plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, 'MarkerFaceColor', nmplot, 'MarkerSize', 5, 'Parent', axesid);
                else
                    plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, 'MarkerFaceColor', nmplot, 'MarkerSize', 5, 'Parent', axesid);
                end
                if sum(find(i==v.resindex))
                   plot(x, y, 's', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, ...
                   'MarkerFaceColor', nmplot, ...
                   'MarkerSize', 13, 'Parent', axesid);
                end
                if sum(find(i==v.tankindex))
                   plot(x, y, 'p', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, ...
                   'MarkerFaceColor', nmplot, ...
                   'MarkerSize', 16, 'Parent', axesid);
                end
           else
                nmplot=char(selectColorNode(hh));
                plot(x, y, 'o', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, 'MarkerFaceColor', nmplot, ...
                    'MarkerSize', 5, 'Parent', axesid);
                if sum(find(i==v.resindex))
                   plot(x, y, 's', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, ...
                   'MarkerFaceColor', nmplot, ...
                   'MarkerSize', 13, 'Parent', axesid);
                end
                if sum(find(i==v.tankindex))
                   plot(x, y, 'p', 'LineWidth', 2, 'MarkerEdgeColor', nmplot, ...
                   'MarkerFaceColor', nmplot, ...
                   'MarkerSize', 16, 'Parent', axesid);
                end
           end
        end
        % Show Node id
        if (strcmpi(Node, 'yes')) %&& ~length(hh))
            text(x, y, v.nodenameid(i), 'Fontsize', fontsize, 'Parent', axesid);%'BackgroundColor', [.7 .9 .7], 'Margin', margin/4);
        end
        % Show Node index
        if (strcmpi(NodeInd, 'yes')) %&& ~length(hh))
            text(x, y, num2str(v.nodeindex(i)), 'Fontsize', fontsize, 'Parent', axesid);%'BackgroundColor', [.7 .9 .7], 'Margin', margin/4);
        end
    end
end

% Legend Plots
if strcmpi(slegend, 'show')
    if isempty(highlightnodeindex) || isempty(highlightnodeindex)
        legendString={'Pipes', 'Pumps', 'Valves', ...
            'Junctions', 'Reservoirs', 'Tanks'};
        legendIndices=sort(legendIndices, 'descend');
        if exist('h','var')
            try
                legend(h(legendIndices), legendString(legendIndices), 'Location', legendposition, 'AutoUpdate', 'off');
            catch
                legend(h(legendIndices), legendString(legendIndices), 'Location', legendposition);
            end
        end
    end
elseif strcmpi(slegend, 'hide')
    %skip
else
    error('Invalid property founobj(legend: "hide", "show")')
end

% Axis OFF and se Background
[xmax, ~]=max(v.nodecoords{1});
[xmin, ~]=min(v.nodecoords{1});
[ymax, ~]=max(v.nodecoords{2});
[ymin, ~]=min(v.nodecoords{2});

if ~isnan(ymax)
    if ymax==ymin
        xlim(axesid, [xmin-((xmax-xmin)*.1), xmax+((xmax-xmin)*.1)]);
        ylim(axesid, [ymin-.1, ymax+.1]);
    elseif xmax==xmin
        xlim(axesid, [xmin-.1, xmax+.1]);
        ylim(axesid, [ymin-(ymax-ymin)*.1, ymax+(ymax-ymin)*.1]);
    else
        xlim(axesid, [xmin-((xmax-xmin)*.1), xmax+((xmax-xmin)*.1)]);
        ylim(axesid, [ymin-(ymax-ymin)*.1, ymax+(ymax-ymin)*.1]);
    end
else
    warning('Undefined coordinates.');
end
axis(axesid, 'off');
try
    whitebg(fig, 'w');
catch
end
if strcmpi(extend, 'yes')
    set(axesid, 'position', [0 0 1 1], 'units', 'normalized');
end
end
function [info_file, tline, allines] = readAllFile(inpname)
    fid = fopen(inpname, 'rt');%or msxname
    allines = textscan(fid, '%s', 'delimiter', '\n');
    [tline]=regexp( fileread(inpname), '\n', 'split');
    for i=1:length(tline)
        str=regexp( tline{i}, '\s', 'split');
        info_file{i} = str(~cellfun('isempty', str));
    end
    fclose(fid);
end
function [Errcode, value] = limitingPotential(obj, param, varargin)
    [tlines]=regexp( fileread(obj.BinTempfile), '\n', 'split');
    Errcode=0;value=0;
    if strcmp(param, 'get')
        for i=1:length(tlines)
           tmp{i}=regexp(tlines{i}, '\s*', 'split');
           atlines=tmp{i};
           atlines(strcmp('', atlines)) = [];
           newlines{i}=tlines{i};
           if ~isempty(atlines)
               if strcmpi(atlines{1}, 'limiting')
                   value = str2double(atlines{3});return;
               end
           end
        end
    else
        fid = writenewTemp(obj.BinTempfile);
        for i=1:length(tlines)
           tmp{i}=regexp(tlines{i}, '\s*', 'split');
           atlines=tmp{i};
           atlines(strcmp('', atlines)) = [];
           newlines{i}=tlines{i};
           getLimit = obj.getBinLimitingPotential;
           if length(atlines)==3 && isempty(getLimit)
               if strcmpi(atlines{1}, 'global') && strcmpi(atlines{2}, 'wall')
                  index=i;
                  newlines{i}=tlines{i};
                  newlines{i+1}=['Limiting', blanks(3), 'Potential', blanks(3), num2str(varargin{1})];
                  break;
              end
           end
        end
        if isempty(getLimit)
            for i=index+2:length(tlines)+1
                newlines{i}=tlines{i-1};
            end
            fprintf(fid, '%s\n', newlines{:});
            if obj.Bin==1
                Errcode=reloadNetwork(obj);
            end
        else
            fprintf(fid, '%s\n', tlines{:});
        end
        fclose(fid);
    end
end
function [Errcode]=setBinParam(obj, indexParameter, parameter, sections, varargin)
    ok=0;Errcode=0;
    if ~isempty(parameter) && (strcmpi(sections{1}, '[SOURCES]')) && indexParameter==11
        indices=find(parameter.BinNodeSourceQuality>-1);
        sources=obj.getBinNodeNameID.BinNodeNameID(indices);ok=1;
    end
    if strcmp(sections{1}, '[PATTERNS]')
        if ischar(indexParameter)
            param=indexParameter;
            indexParameter=21;
            pat=obj.getBinPatternsInfo;
        end
    end
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
    cntIDpat=0;start=0;stop1=0;stop11=0;itsOkQual=0;stop2=0;stop22=0;
    for i=1:length(tlines)
        tt=regexp(tlines{i}, '\s*', 'split');
        if strcmp(tt{1}, sections{1})
            start=i;
        end
        if ~strcmp(sections{1}, '[REACTIONS]')
            if strcmp(tt{1}, sections{2})
                stop=i;
            end
        else
            if strcmp(tt{1}, '[MIXING]')
                stop1=i;
            end
            if strcmp(tt{1}, '[ENERGY]')
                stop11=i;
            end
            stop=max([stop1 stop11]);
        end
        cnts=obj.BinLinkPipeCount;
        if length(sections)>2
            if strcmp(tt{1}, sections{3})
                start2=i;
            end
            if strcmp(tt{1}, sections{5})
                stop2=i;
            elseif strcmp(tt{1}, sections{4})
                stop22=i;
            end
            stop_2=max([stop2 stop22]);
            cnts=obj.BinNodeJunctionCount;
        end
    end
    if strcmpi(sections{1}, '[RESERVOIRS]')
        cnts=obj.BinNodeReservoirCount;
    elseif strcmpi(sections{1}, '[TANKS]')
        cnts=obj.BinNodeTankCount;
    end
    fid = writenewTemp(obj.BinTempfile);
    ll=1;
    for i=start:stop
       % Get first token in the line
       tok = strtok(tlines{i});
       if isempty(tok), continue; end
       % Skip blank Clines and comments
       if strcmp(tok(1), ';') && ok==0
       elseif sum(tlines{i}=='[') && ok==0
       elseif isempty(tok) && ok==0
       % skip
       else
           clear atlines;
           atlines = checktlines(tlines{i});
           if ll<cnts+1
               if (~isempty(parameter) && indexParameter~=3 && length(sections)<3 || indexParameter==2) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[REACTIONS]')) && (~strcmpi(sections{1}, '[TIMES]')) && (~strcmpi(sections{1}, '[OPTIONS]')) && (~strcmpi(sections{1}, '[PATTERNS]')) && ~strcmpi(sections{1}, '[TANKS]')
                   if indexParameter ~= 8
                       atlines{indexParameter} = num2str(parameter(ll));
                   else
                       atlines{indexParameter} = num2str(parameter{ll});
                   end
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i}=newlines;
               end
               if (~isempty(parameter) && length(atlines)>2 && indexParameter ~= 8) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[PATTERNS]')) && (~strcmpi(sections{1}, '[REACTIONS]')) && (~strcmpi(sections{1}, '[TIMES]')) && (~strcmpi(sections{1}, '[OPTIONS]')) && ~strcmpi(sections{1}, '[TANKS]')
                   if length(atlines)<indexParameter, atlines{indexParameter}={''}; end
                   if ~isempty(atlines{indexParameter}) && ~sum(atlines{3}==';')
                       if indexParameter==4 && length(sections)>2
                          atlines{indexParameter} = num2str(parameter{ll});
                       else
                          atlines{indexParameter} = num2str(parameter(ll));
                       end
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
               end
               if ~isempty(parameter) && (strcmpi(sections{1}, '[RESERVOIRS]') || strcmpi(sections{1}, '[TANKS]')) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[TIMES]')) && (~strcmpi(sections{1}, '[OPTIONS]')) && (~strcmpi(sections{1}, '[PATTERNS]'))
                   if (indexParameter==2 && strcmpi(sections{1}, '[RESERVOIRS]')) || strcmpi(sections{1}, '[TANKS]')
                      atlines{indexParameter} = num2str(parameter(ll));
                   elseif strcmpi(sections{1}, '[RESERVOIRS]')
                      atlines{indexParameter} = num2str(parameter{ll});
                   end
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i}=newlines;
               end
           end
           if ~isempty(parameter) && (strcmpi(sections{1}, '[REACTIONS]')) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[TIMES]')) && (~strcmpi(sections{1}, '[OPTIONS]')) && (~strcmpi(sections{1}, '[PATTERNS]'))
               if strcmpi(atlines{1}, 'global')
                  if strcmpi(atlines{2}, 'bulk') && indexParameter==1
                     atlines{3} = num2str(parameter);
                  end
                  if strcmpi(atlines{2}, 'wall') && indexParameter==3
                     atlines{3} = num2str(parameter);
                  end
               end
               newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp}, blanks(10)];
               end
               tlines{i}=newlines;
           end
           mins=1;
           if ~isempty(parameter) && (strcmpi(sections{1}, '[TIMES]')) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[OPTIONS]')) && (~strcmpi(sections{1}, '[PATTERNS]'))
               if strcmpi(atlines{1}, 'DURATION') && indexParameter==1
                    [mm, mins]=sec2hrs(parameter);
                    atlines{2} = num2str(mm);
               elseif strcmpi(atlines{1}, 'HYDRAULIC') && indexParameter==2
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'QUALITY') && indexParameter==3 && ~strcmpi(atlines{2}, 'TRACE')
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'PATTERN') && strcmpi(atlines{2}, 'TIMESTEP') && indexParameter==4
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'PATTERN') && strcmpi(atlines{2}, 'START') && indexParameter==5
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'REPORT') && strcmpi(atlines{2}, 'TIMESTEP') && indexParameter==6
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'REPORT') && strcmpi(atlines{2}, 'START') && indexParameter==7
                    [mm, mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmpi(atlines{1}, 'STATISTIC') && indexParameter==8
                   atlines{2} = parameter;
               end
               if mins==0 && length(atlines)>3
                   atlines{4}='';
               end
                   newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp}, blanks(10)];
               end
               tlines{i}=newlines;
           end
           if ~isempty(parameter) && (strcmpi(sections{1}, '[OPTIONS]')) && (~strcmpi(sections{1}, '[SOURCES]')) && (~strcmpi(sections{1}, '[PATTERNS]'))
               if strcmpi(atlines{1}, 'QUALITY') && indexParameter==1 && itsOkQual==0
                   clear atlines;
                   atlines{1}=parameter;itsOkQual=1;
               end
               newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp}, blanks(10)];
               end
               tlines{i}=newlines;
           end
           if ~isempty(parameter) && (strcmpi(sections{1}, '[PATTERNS]')) && (~strcmpi(sections{1}, '[SOURCES]'))
               idpat=param;
               if strcmp(atlines{1}, idpat) && cntIDpat==0
                   atlines=[idpat blanks(12)];
                   cntIDpat=1;
                   newlines=atlines;
                   zz=0;lll=0;
                   lengthparam=length(parameter);
                   mlen=pat.BinPatternValue(find(strcmp(idpat, pat.BinPatternNameID)));
                   if lengthparam<length(mlen{1})
                       for j=(lengthparam+1):length(mlen{1})
                           parameter(1, j)=parameter(1, j-lengthparam);
                       end
                   end
                   for pp=1:size(parameter, 1)
                       if mod(lengthparam, 6)==0
                           zz=zz+lengthparam/6;
                       else
                           zz=zz+1;
                           if mod(lengthparam, 6)
                               zz=zz+1;
                           end
                       end
                       m=1;
                       for k=lll+1:zz
                           if mod(lengthparam, 6) && k==zz
                               newlines = [idpat blanks(15) num2str(parameter(m:end))];
                           else
                               newlines = [idpat blanks(15) num2str(parameter(m:m+6-1))];
                           end
                           m=m+6;
                           tlines{i+k-1}=newlines;
                       end
                       lll=zz;
                   end
               elseif cntIDpat==1
                   if strcmp(atlines{1}, idpat)
                       newlines='';
                   else
                   newlines=[atlines{1} blanks(12)];
                   cntIDpat=1;
                   for pp=2:length(atlines)
                       newlines = [newlines, num2str(atlines{pp}), blanks(12)];
                   end
                   end
               else
                   newlines=tlines{i};
               end
               if isempty(parameter) && (~strcmpi(sections{1}, '[PATTERNS]')) && (strcmpi(sections{1}, '[SOURCES]'))
                  tlines{i}=newlines;
               end
           end
           if ~isempty(parameter) && (strcmpi(sections{1}, '[SOURCES]')) && indexParameter==11
               for kk=start+1:stop-1
                  tlines{kk}='';
               end
               for u=1:start
                   nnlines{u}=tlines{u};
               end
               for u=start+1:start+1+length(sources)
                   nnlines{u}=[];
               end
               for k=start+1:length(tlines)
                   nnlines{u}=tlines{k};
                   u=u+1;
               end
               tlines=nnlines;clear nnlines;
               for u=1:length(sources)
                   if sum(strcmp(sources, atlines{1})) || ok==1
                       atlines{1}=sources{u};
                       atlines{2}=parameter.BinNodeSourceType{indices(u)};
                       if length(atlines)>2 || ok==1
                           atlines{3}=num2str(parameter.BinNodeSourceQuality(indices(u)));
                       end
                       if length(atlines)>3 || ok==1
                           atlines{4}=parameter.BinNodeSourcePatternNameID{indices(u)};
                       end
                       newlines=[atlines{1} blanks(10)];
                       for pp=2:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i+u}=newlines;
                   end
               end
               indexParameter=0;break;
           end
           ll=ll+1;
        end
    end
    if length(sections)>2
       ll=1;newInd=indexParameter-1;
       for i=start2:stop_2
           % Get first token in the line
           tok = strtok(tlines{i});
           % Skip blank Clines and comments
           if isempty(tok), continue; end
           if strcmp(tok(1), ';')
           elseif sum(tlines{i}=='[')
           elseif isempty(tok)
           % skip
           else
               clear atlines;
               atlines = checktlines(tlines{i});
               if ~isempty(parameter) && length(atlines)>1 && indexParameter~=2%BaseDemands
                   if ~length(strfind(cell2mat(atlines), ';')) %~isempty(atlines{3}) &&
                       if length(parameter)<ll, continue; end
                       if newInd==3
                          atlines{newInd} = num2str(parameter{ll});
                       else
                          atlines{newInd} = num2str(parameter(ll));
                       end
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp}, blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
               end
               ll=ll+1;
           end
        end
    end
    clear parameter;
    fprintf(fid, '%s\n', tlines{:});
    fclose(fid);
    if obj.Bin, Errcode=reloadNetwork(obj); end
end
function [mm, mins] = sec2hrs(parameter)
   mm='';hrs=0;mins=0;
   if parameter >= 3600
       hrs=floor(parameter/3600);
       mm=[num2str(hrs), ':'];
   end
   if parameter >= 60
       mins=((parameter - 3600*hrs)/60);
   end
   if hrs
       mm=[mm sprintf('%d', (parameter-3600*hrs-60*mins))];
   elseif hrs==0 && mins==0
       mm=parameter;
       if mm<10
           mm=['00:00:0' num2str(mm)];
       else
           mm=['00:00:' num2str(mm)];
       end
   else
       mm=[sprintf('%.20f', mins) '       min'];
   end
end
function [Errcode]=setBinParam2(obj, parameter, sections, zz, varargin)
    Errcode=0;
    if strcmp(sections{1}, '[STATUS]')
        value =obj.getBinLinksInfo;
        if strcmp(sections{3}, 'pump')
            nameID=value.BinLinkPumpStatusNameID;
            cntlv=obj.BinLinkPumpCount;
        elseif strcmp(sections{3}, 'valve')
            nameID=value.BinLinkValveStatusNameID;
            cntlv=obj.BinLinkValveCount;
            if strcmpi(parameter, 'NONE'), Errcode=-1;return;end
        end
    elseif strcmp(sections{1}, '[PATTERNS]')
        value=obj.getBinPatternsInfo;
        if ~isempty(value.BinPatternValue)
            paramAll=[value.BinPatternValue parameter];
        else
            paramAll{1}=[value.BinPatternValue parameter];
        end
        patternsid=[value.BinPatternNameID varargin];
    end
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
    fid = writenewTemp(obj.BinTempfile);
    for i=1:length(tlines)
        tt=regexp(tlines{i}, '\s*', 'split');
        tok = strtok(tlines{i});m=1;
        % Skip blank Clines and comments
        if isempty(tok), continue;
        elseif isempty(tt{m})
            m=m+1;
        end
        if strcmp(tt{m}, sections{1})
            start=i;
        end
        if strcmp(tt{m}, sections{2})
            stop=i;
        end
    end
   for kk=start+1:stop-1
      tlines{kk}='';
   end
   for u=1:start
       nnlines{u}=tlines{u};
   end
   for u=start+1:start+2+zz
       nnlines{u}=[];
   end
   for k=start+1:length(tlines)
       nnlines{u}=tlines{k};
       u=u+1;
   end
   tlines=nnlines;clear nnlines;
   for i=start+1:stop+zz
       % Get first token in the line
       tok = strtok(tlines{i});
       if isempty(tok), tok='1'; end
       if strcmp(tok(1), ';')
       else
           clear atlines;
           if ~isempty(parameter) && strcmp(sections{1}, '[STATUS]')
               for ee=1:cntlv
                   atlines{1} = nameID{ee};
                   atlines{2} = num2str(parameter{ee});
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i+ee}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1}, '[QUALITY]')
               for e=1:obj.BinNodeCount
                   atlines{1} = obj.BinNodeNameID{e};
                   atlines{2} = num2str(parameter(e));
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i+e}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1}, '[PATTERNS]')
               zz=0;ll=0;
               for e=1:length(patternsid)
                   if e<length(patternsid)
                       if mod(length(value.BinPatternValue{e}), 6)==0
                           zz=zz+length(value.BinPatternValue{e})/6;
                       else
                           zz=zz+1;
                       end
                   else
                       zz=zz+1;
                       if mod(length(paramAll{e}), 6)
                           zz=zz+1;
                       end
                   end
                   m=1;
                   for k=ll+1:zz
                       if mod(length(paramAll{e}), 6) && k==zz
                           newlines = [patternsid{e} blanks(15) num2str(paramAll{e}(m:end))];
                       else
                           newlines = [patternsid{e} blanks(15) num2str(paramAll{e}(m:m+6-1))];
                       end
                       m=m+6;
                       tlines{i+k}=newlines;
                   end
                   ll=zz;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1}, 'Global') && varargin{1}==1
               for e=1:obj.BinLinkCount
                   atlines{1} = 'WALL';
                   atlines{2} = obj.BinLinkNameID{e};
                   atlines{3} = num2str(parameter(e));
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i+e}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1}, 'Global') && varargin{1}==2
               for e=1:obj.BinLinkCount
                   atlines{1} = 'BULK';
                   atlines{2} = obj.BinLinkNameID{e};
                   atlines{3} = num2str(parameter(e));
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp}, blanks(10)];
                   end
                   tlines{i+e}=newlines;
               end
           end
       end
       break;
   end
   fprintf(fid, '%s\n', tlines{:});
   fclose(fid);
    if obj.Bin==1
        Errcode=reloadNetwork(obj);
    end
end
function value = getBinParam(obj, sections, varargin)
    warning off;
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
    if strcmp(sections{1}, '[SOURCES]')
        value.BinNodeSourcePatternIndex = nan(1, obj.BinNodeCount);
        value.BinNodeSourceQuality = nan(1, obj.BinNodeCount);
        value.BinNodeSourceTypeIndex = nan(1, obj.BinNodeCount);
        value.BinNodeSourceType = cell(1, obj.BinNodeCount);
        value.BinNodeSourcePatternNameID = cell(1, obj.BinNodeCount);
    else
        value=[];
    end
    v{1}='';
    for i=1:length(tlines)
        tt=regexp(tlines{i}, '\s*', 'split');
        tok = strtok(tlines{i});m=1;
        % Skip blank Clines and comments
        if isempty(tok), continue;
        elseif isempty(tt{m})
            m=m+1;
        end
        if strcmp(tt{m}, sections{1})
            start=i;
        end
        if strcmp(tt{m}, sections{2})
            stop=i;
        end
    end
    d=1;
   for i=start+1:stop-1
       % Get first token in the line
       tok = strtok(tlines{i});
       if isempty(tok), continue; end
       if strcmp(tok(1), ';')
       else
           clear atline;
           atline = checktlines(tlines{i});

            if strcmp(sections{1}, '[STATUS]')
                if sum(strcmp(who, 'atline'))
                    value.BinLinkInitialStatus{d}=atline{2};
                    value.BinLinkInitialStatusNameID{d}=atline{1};
                    d=d+1;
                end
           end
           if strcmp(sections{1}, '[PATTERNS]')
                if sum(strcmp(who, 'atline'))
                    value.BinPatternNameID{d}=atline{1};
                    value.BinPatternNameID=unique(value.BinPatternNameID);
                    w{d}=atline{1};
                    dd=length(w);
                    w=unique(w);
                    d=length(w);
                    if dd>1 && dd~=d
                        v{d}=[v{d} str2num(char(atline(2:end)))'];
                    else
                        v{d}=str2num(char(atline(2:end)))';
                    end
                    value.BinPatternValue=v; % single
                    d=d+1;
                    value.BinCountPatternlines=d;
                end
           end
           if strcmp(sections{1}, '[SOURCES]')
               if sum(strcmp(who, 'atline'))
                   if length(atline)>2
                       indexPat=getBinNodeIndex(obj, atline{1});
                       indexNode=getBinNodeIndex(obj, atline{1});
                       if length(atline)>3
                           value.BinNodeSourcePatternIndex(indexPat)=getBinPatternIndex(obj, atline{4});
                           value.BinNodeSourcePatternNameID{indexNode}=atline{4};
                       end
                       value.BinNodeSourceQuality(indexNode)=str2double(atline{3});
                       value.BinNodeSourceTypeIndex(indexNode)=find((strcmpi(obj.TYPESOURCE, atline{2})-1)>-1)-1;
                       value.BinNodeSourceType{indexNode}=obj.TYPESOURCE{value.BinNodeSourceTypeIndex(indexNode)+1};
                   end
               end
           end
       end
   end
   if strcmp(sections{1}, '[STATUS]')
        n=obj.getBinLinksInfo;
        value.BinLinkInitialStatusNameID=n.BinLinkInitialStatusNameID;
        value.BinLinkInitialStatus=n.BinLinkInitialStatus;
        value.BinLinkPumpStatus=n.BinLinkPumpStatus;
        value.BinLinkPumpStatusNameID=n.BinLinkPumpStatusNameID;
        value.BinLinkValveStatus=n.BinLinkValveStatus;
        value.BinLinkValveStatusNameID=n.BinLinkValveStatusNameID;
   end
   warning on;
end
function [Errcode]=addCurve(obj, newCurveID, varargin)
v=obj.getBinCurvesInfo;Errcode=0;
CurveX=varargin{1};
CurveY=varargin{2};
typecode=varargin{3};
% PUMP 0 EFFICIENCY 1 VOLUME 2 HEADLOSS 3
for i=1:length(CurveX)
    if i+1<length(CurveX)+1
        if CurveX(i)>=CurveX(i+1)
            if strfind([0 1 3], typecode)
                warning('Flow values are not in ascending order.');
                Errcode=-1;
                return;
            elseif typecode==2
                warning('Heigh values are not in ascending order.');
                Errcode=-1;
                return;
            end
        end
    end
end

% Check if new ID already exists
if ismember(newCurveID, v.BinCurveNameID)
    warning('Curve "%s" already exists.', newCurveID);Errcode=-1; return;
end
sect=0;
% Open and read inpname
% Read all file and save in variable info
[~, info] = obj.readInpFile;
% write
fid2 = writenewTemp(obj.BinTempfile);
sps=blanks(18);
nn=0;yy=0;
for t = 1:length(info)
    a = regexp(info{t}, '\s*', 'split');
    if isempty(a)
        % skip
    elseif isempty(info{t})
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u}, '[CURVES]')
                fprintf(fid2, '[CURVES]');
                sect=1; break;
            end
            if (sum(info{t}=='[') && nn==0)
                if yy==0
                    if sect==0
                        fprintf(fid2, '[CURVES]\n;ID                X-Value            Y-Value\n');
                    end
                    if typecode==0
                        fprintf(fid2, ';PUMP: PUMP:%sX-Value%sY-Value\n', sps, sps); yy=1;
                    elseif typecode==1
                        fprintf(fid2, ';PUMP: EFFICIENCY:\n'); yy=1;
                    elseif typecode==2
                        fprintf(fid2, ';PUMP: VOLUME:\n'); yy=1;
                    elseif typecode==3
                        fprintf(fid2, ';PUMP: HEADLOSS:\n'); yy=1;
                    end
                end
                for i=1:length(CurveX)
                    fprintf(fid2, '%s%s%d%s%d', newCurveID, sps, CurveX(i), sps, CurveY(i));
                    fprintf(fid2, '\r\n');
                end
                fprintf(fid2, '%s', a{u});
                fprintf(fid2, '\r\n');
                nn=1;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
if obj.Bin==1
    Errcode=reloadNetwork(obj);
end
end
function [BinCurveNameID, BinCurveXvalue, BinCurveYvalue, BinCurveAllLines, BinCurveTypes, BinCurveCount, BinCTypes]=CurveInfo(obj)
BinCurveTypes=[];Bintypecode=0;BinCNameID={};BinCurveNameID={};BinCurveCount=0;
BinCurveXvalue=[];BinCurveYvalue=[];BinCurveAllLines={};sect=0;i=1;u=1;BinCTypes=[];
cc=1;uu=1;gg=1;
% Open epanet input file
[~, info] = obj.readInpFile;
for h=1:length(info)
    tline = info{h};
    if ~ischar(tline),   break,   end
    % Get first token in the line
    tok = strtok(tline);
    % Skip blank Clines and comments
    if isempty(tok), continue, end
    ee=regexp(tline, '\w*EFFICIENCY*\w', 'match');
    nn=regexp(tline, '\w*VOLUME*\w', 'match');
    kk=regexp(tline, '\w*HEADLOSS*\w', 'match');
    if strcmp(ee, 'EFFICIENCY'), %typecode=1;   % EFFICIENCY
    elseif strcmp(nn, 'VOLUME'), %typecode=2;   % VOLUME
    elseif strcmp(kk, 'HEADLOSS'), %typecode=3; % HEADLOSS
    else
        if (tok(1) == ';'), continue, end  %typecode=0;
    end
    if (tok(1) == '[')
        % [CURVES] section
        if strcmpi(tok(1:5), '[CURV')
            sect = 1;
            continue;
            % [END]
        elseif strcmpi(tok(1:4), '[END')
            break;
        else
            sect = 0;
            continue;
        end
    end
    if sect == 0
        continue;
        % Curves
    elseif sect == 1
        ee=regexp(tline, '\w*EFFICIENCY*\w', 'match');
        nn=regexp(tline, '\w*VOLUME*\w', 'match');
        kk=regexp(tline, '\w*HEADLOSS*\w', 'match');
        if strcmp(ee, 'EFFICIENCY'), Bintypecode=1;   % EFFICIENCY
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif strcmp(nn, 'VOLUME'), Bintypecode=2;   % VOLUME
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif strcmp(kk, 'HEADLOSS'), Bintypecode=3; % HEADLOSS
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif (~length(strcmp(nn, 'VOLUME')) || ~length(strcmp(ee, 'EFFICIENCY')) || ~length(strcmp(kk, 'HEADLOSS'))) &&  (tok(1)==';'), Bintypecode=0; % HEADLOSS
            BinCurveAllLines{u}=tline;u=u+1;continue;
        else
            a = textscan(tline, '%s %f %f');
            %aa=regexp(tline, '\s', 'split');
            BinCNameID{i}=a{1};
            if i==1
                BinCurveTypes(gg)=Bintypecode;
            elseif ~strcmp(BinCNameID{i-1}, BinCNameID{i})
                if (u-i+1)==length(BinCurveTypes)
                    Bintypecode=0;
                end
                gg=gg+1;
                BinCurveTypes(gg)=Bintypecode;
            end
            BinCTypes(i)=Bintypecode;
        end
        if i==1
            BinCurveXvalue{cc}(uu)=a{2};
            BinCurveYvalue{cc}(uu)=a{3};
        elseif strcmp(BinCNameID{i-1}, BinCNameID{i})
            BinCurveXvalue{cc}(uu)=a{2};
            BinCurveYvalue{cc}(uu)=a{3};
        elseif ~strcmp(BinCNameID{i-1}, BinCNameID{i})
            cc=cc+1;uu=1;
            BinCurveXvalue{cc}(uu)=a{2};
            BinCurveYvalue{cc}(uu)=a{3};
        end
        uu=uu+1;
        BinCurveAllLines{u}=tline;
        i=i+1;u=u+1;
    end
end
if ~isempty(BinCNameID)
    for i=1:length(BinCNameID)
        nn(i)=BinCNameID{i};
    end
    BinCurveNameID=unique(nn);
    BinCurveCount=length(BinCurveNameID);
end
end
function node_index = addBinNode(obj, typeCode, nodeID, coords, varargin)
    if ~iscell(nodeID)
        nodeID = {nodeID};
    end
    nodesInfo = obj.getBinNodesInfo;
    for i = 1:length(nodeID)
        if ismember(nodeID{i}, nodesInfo.BinNodeNameID)
            warning(['Node ', nodeID{i}, ' already exists.'])
            node_index=-1;
            return;
        end
    end
    if typeCode == 1 || typeCode == 2
        if typeCode == 1
            patternID = varargin{3};
        else
            patternID = varargin{2};
        end
        for i = 1:length(patternID)
            if ~isempty(patternID{i})
                if ~ismember(num2str(patternID{i}), obj.getBinPatternsInfo.BinPatternNameID)
                    warning(['Pattern ', patternID{i}, ' does not exist.'])
                    node_index=-1;
                    return;
                end
            end
        end
    end
    fid = fopen(obj.BinTempfile); % Opens the file for read access
    % Creates the string that will be set under the [NODE] section
    if typeCode == 1
        str_junction = str_make(nodeID, varargin{1}, varargin{2}, varargin{3});
        str_demands = str_make(nodeID, varargin{2}, varargin{3}, varargin{4});
        quality = varargin{5};
    elseif typeCode == 2
        str_reserv = str_make(nodeID, varargin{1}, varargin{2});
        quality = varargin{3};
    elseif typeCode == 3
        str_tank = str_make(nodeID, varargin{1}, varargin{3}, varargin{4}, varargin{5}, varargin{2}, varargin{6}, varargin{7});
        quality = varargin{8};
    end
    % Creates the string that will be set under the [QUALITY] section
    str_qual = str_make(nodeID, quality);
    % Creates the string that will be set under the [COORDINATES] section
    str_coords = str_make(nodeID, coords(:, 1), coords(:, 2));
    % Creates the entire text that will replace the .inp file
    texta = char;
    while ~feof(fid)
        aline = fgetl(fid);
        section_checker = regexp(aline,'\s','split','once');
        if length(section_checker)>1
            section_checker = section_checker{1};
        end
        texta = [texta, aline, char(10)];
        if typeCode == 1
            if strcmp(section_checker, '[JUNCTIONS]')
                for i = 1:obj.getBinNodesInfo.BinNodeJunctionCount
                    aline = fgetl(fid);
                    texta = [texta, aline, char(10)];
                end
                texta = [texta, str_junction];
            end
            if strcmp(section_checker, '[DEMANDS]')
                texta = [texta, str_demands];
            end
        elseif typeCode == 2
            if strcmp(section_checker, '[RESERVOIRS]')
                for i = 1:obj.getBinNodesInfo.BinNodeReservoirCount
                    aline = fgetl(fid);
                    texta = [texta, aline, char(10)];
                end
                texta = [texta, str_reserv];
            end
        elseif typeCode == 3
            if strcmp(section_checker, '[TANKS]')
                for i = 1:obj.getBinNodesInfo.BinNodeTankCount
                    aline = fgetl(fid);
                    texta = [texta, aline, char(10)];
                end
                texta = [texta, str_tank];
            end
        end
        if strcmp(section_checker, '[QUALITY]')
            texta = [texta, str_qual];
        end
        if strcmp(section_checker, '[COORDINATES]')
            texta = [texta, str_coords];
        end
    end
    fclose('all');
    fid = fopen(obj.BinTempfile, 'w');   % Opens file for writing and discard existing contents
    fprintf(fid, texta);   % Writes the new text in the .inp file
    fclose('all');
    if obj.Bin, obj.Errcode = reloadNetwork(obj); end
    node_index = zeros(1, length(nodeID));
    for i = 1:length(nodeID)
        node_index(i) = obj.getBinNodeIndex(nodeID{i});
    end
end
function link_index = addBinLink(obj, typeCode, linkID, from, to, varargin)
        if ~iscell(linkID)
            linkID = {linkID};
        end
        if ~iscell(from)
            from = {from};
        end
        if ~iscell(to)
            to = {to};
        end
        LinksInfo = obj.getBinLinksInfo;
        for i = 1:length(linkID)
            if ismember(linkID{i}, LinksInfo.BinLinkNameID)
                warning(['Link ', linkID{i}, ' already exists.'])
                link_index=-1;
                return;
            end
        end

        BinNodeNameID = obj.getBinNodeNameID.BinNodeNameID;
        for i = 1:length(linkID)
            if ~ismember(from{i}, BinNodeNameID)
                warning(['Node ', from{i}, ' does not exist.'])
                link_index=-1;
                return;
            end
            if ~ismember(to{i}, BinNodeNameID)
                warning(['Node ', to{i}, ' does not exist.'])
                link_index=-1;
                return;
            end
        end
        fid = fopen(obj.BinTempfile); % Opens the file for read access
        % Creates the string that will be set under the [NODE] section
        if strcmpi(typeCode, 'PIPE')
            if ~iscell(varargin{5})
             	varargin{5} = {varargin{5}};
            end
            str_pipe = str_make(linkID, from, to, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        elseif strcmpi(typeCode, 'PUMP')
            if ~iscell(varargin{1})
             	varargin{1} = {varargin{1}};
            end
            str_pump = str_make(linkID, from, to, varargin{1});
        elseif strcmpi(typeCode, 'VALVE')
            if ~iscell(varargin{1})
             	varargin{1} = {varargin{1}};
            end
            str_valve = str_make(linkID, from, to, varargin{2}, varargin{1}, varargin{3}, varargin{4});
        end
        % Creates the entire text that will replace the .inp file
        texta = char;
        while ~feof(fid)
            aline = fgetl(fid);
            section_checker = regexp(aline,'\s','split','once');
            if length(section_checker)>1
                section_checker = section_checker{1};
            end
            texta = [texta, aline, char(10)];
            if strcmpi(typeCode, 'PIPE')
                if strcmpi(section_checker, '[PIPES]')
                    for i = 1:LinksInfo.BinLinkPipeCount
                        aline = fgetl(fid);
                        texta = [texta, aline, char(10)];
                    end
                    texta = [texta, str_pipe];
                end
            elseif strcmpi(typeCode, 'PUMP')
                if strcmpi(section_checker, '[PUMPS]')
                    for i = 1:LinksInfo.BinLinkPumpCount
                        aline = fgetl(fid);
                        texta = [texta, aline, char(10)];
                    end
                    texta = [texta, str_pump];
                end
            elseif strcmpi(typeCode, 'VALVE')
                if strcmp(section_checker, '[VALVES]')
                    for i = 1:LinksInfo.BinLinkValveCount
                        aline = fgetl(fid);
                        texta = [texta, aline, char(10)];
                    end
                    texta = [texta, str_valve];
                end
            end
        end
        fclose('all');
        fid = fopen(obj.BinTempfile, 'w');   % Opens file for writing and discard existing contents
        fprintf(fid, texta);   % Writes the new text in the .inp file
        fclose('all');
        if obj.Bin, obj.Errcode = reloadNetwork(obj); end
        link_index = zeros(1, length(linkID));
        for i = 1:length(linkID)
            link_index(i) = obj.getBinLinkIndex(linkID{i});
        end
end
function str = str_make(ID, varargin)
    str = ID{1};
    for i = 1:length(ID)
        if i>1
            str = [str, ID{i}];
        end
        for j = 1:(nargin-1)
            if isnumeric(varargin{j})
                value =  num2str(varargin{j}(i));
            else
                value = varargin{j}{i};
            end
            str = [str, blanks(10), value];
        end
        str = [str, char(10)];
    end
end
function Errcode=addNode(obj, typecode, varargin)
    % addNode - Add node in the network. Node type codes consist of the
    % following constants: EN_JUNCTION 0 Junction node EN_RESERVOIR 1
    % Reservoir node EN_TANK 2 Tank node
    newID=varargin{1};Errcode=0;
    X=varargin{2};
    Y=varargin{3};
    links = obj.getBinLinksInfo;
    nodes = obj.getBinNodesInfo;
    l = unique([links.BinLinkFromNode links.BinLinkToNode]);
    if nodes.BinNodeCount~=length(l)
        cg=ismember(nodes.BinNodeNameID, l);
        if ~(sum(cg)==nodes.BinNodeCount)
            ind=find(cg==0);
            warning('Node %s disconnected.', nodes.BinNodeNameID{ind(1)});
            Errcode=-1; return;
        end
    end
    if sum(typecode==[0, 1]) % junction & reservoir
        if typecode==0
            v=obj.getBinPatternsInfo;
            newidpattern=varargin{6};
            patterns=v.BinPatternNameID;
            if ~sum(strcmp(newidpattern, patterns))
                warning('Invalid argument found.');
                Errcode=-1; return;
            end
            newBaseDemand=varargin{5};
        end
        newElevation=varargin{4};
        initqual=0;
    else
        % Initial TANK
        MaxLevel=varargin{4};
        Diameter=varargin{5};
        Initlevel=varargin{6};
        newElevation=varargin{7};
        initqual=varargin{8};
        MinLevel=varargin{9};
        MinVol=varargin{10};
    end
    % Check if id new already exists
    if isempty(nodes.BinNodeNameID)
        warning('There is no such object in the network.');Errcode=-1; return;
    end
    if ismember(newID, nodes.BinNodeNameID)
        warning('Node "%s" already exists.', newID);
        Errcode=-1;return;
    end
    % check section in inpname, [JUNCTIONS], [RESERVOIRS], [TANKS]
    stank_check=1;
    sreservoir_check=1;
    sjunction_check=1;
    % Open and read inpname
    % Read all file and save in variable info
    [~, info, ~] = obj.readInpFile;
    fid2 = writenewTemp(obj.BinTempfile);
    % Initiality
    qualch=0;qq=0;
    Coordch=0;onetime=1;gg=0;
    sps1=blanks(3);
    for t = 1:length(info)
        c = info{t};
        if ~isempty(c)
            a = regexp(c, '\s*', 'split');
        else
            a='';
        end
        if isempty(a)
            % skip
        elseif isempty(c)
            % skip
        else
            u=1;
            while u < length(a)+1
                % Find [brackets] cnt=2;
                cnt=bracketsCheck(a{u});
                %%%%%%%% Quality Section %%%%%%%%
                if strcmp(a{u}, '[QUALITY]')
                    fprintf(fid2, '[QUALITY]');
                    qualch=1;
                    break;
                end
                if (cnt==2 && qualch==1)
                    fprintf(fid2, '%s%s%d', newID, sps1, initqual);
                    fprintf(fid2, '\r\n');qq=1;
                end
                %%%%%%%% Coordinates Section %%%%%%%%
                if strcmp(a{u}, '[COORDINATES]');
                    fprintf(fid2, '[COORDINATES]');
                    Coordch=1; break;
                end
                if length(strfind(c, ';Node'))==1 && Coordch==1 && cnt~=2
                    break;
                    elseif u==1 && Coordch==1
                    if ((gg==0)) && (typecode==0)
                        fprintf(fid2, '%s%s%d%s%d\n', newID, sps1, X, sps1, Y);
                    end
                    gg=gg+1;
                end
                if isempty(obj.NodeCoordinates) && obj.Bin==1% no bin
                    if strcmp(a{u}, '[END]')
                        fprintf(fid2, '%s', '[COORDINATES]');
                        fprintf(fid2, '\r\n');
                        for qq=1:length(X)
                            fprintf(fid2, '%s%s%d%s%d', char(newID(qq)), sps1, X(qq), sps1, Y(qq));
                            fprintf(fid2, '\r\n');
                        end
                        fprintf(fid2, '%s%s%d%s%d\n', ...
                        newID, sps1, X, sps1, Y);
                        fprintf(fid2, '%s', a{u}); fprintf(fid2, '\r\n');
                    end
                end
                %%%%%%%% Nodes Section %%%%%%%%
                if (cnt==2 && (strcmp(a{u}, '[TANKS]') || strcmp(a{u}, '[JUNCTIONS]') || strcmp(a{u}, '[RESERVOIRS]') || strcmp(a{u}, '[DEMANDS]')))
                    if sjunction_check==0 && typecode==0 && strcmp(a{u}, '[RESERVOIRS]')
                        fprintf(fid2, '[JUNCTIONS]');
                        fprintf(fid2, '\n%s%s%d%s%s\n', newID, sps1, newElevation, sps1, sps1);
                    end
                    if sreservoir_check==0 && typecode==1 && strcmp(a{u}, '[TANKS]')
                        fprintf(fid2, '[RESERVOIRS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s\n', newID, sps1, newElevation, sps1, '', sps1);
                    end
                    if stank_check==0 && typecode==2 && strcmp(a{u}, '[PIPES]')
                        fprintf(fid2, '[TANKS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d\n', newID, sps1, newElevation, sps1, Initlevel, sps1, MinLevel, sps1, ...
                        MaxLevel, sps1, Diameter, sps1, MinVol);
                    end
                    fprintf(fid2, '%s', a{u});
                    %%%%%%%% Jynctions Section %%%%%%%%
                    if typecode==0 && strcmp(a{u}, '[JUNCTIONS]')
                        fprintf(fid2, '\n%s%s%d', newID, sps1, newElevation);
                    end
                    if typecode==0 && strcmp(a{u}, '[DEMANDS]')
                        fprintf(fid2, '\n%s%s%d', newID, sps1, newBaseDemand);
                    end
                    %%%%%%%% Reservoirs Section %%%%%%%%
                    if typecode==1 && strcmp(a{u}, '[RESERVOIRS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s', newID, sps1, newElevation);
                    end
                    %%%%%%%% Tanks Section %%%%%%%%
                    if typecode==2 && strcmp(a{u}, '[TANKS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d', newID, sps1, newElevation, sps1, Initlevel, sps1, MinLevel, sps1, ...
                        MaxLevel, sps1, Diameter, sps1, MinVol);
                    end
                elseif isempty(a{u})
                else
                    fprintf(fid2, '%s%s', a{u}, sps1);
                end
                u=u+1;
            end
            %%%%%%%% Coordinates Section %%%%%%%%
            if gg~=0 && onetime==1
                % Correction Index
                if isempty(char(nodes.BinNodeJunctionNameID))
                    nodes.BinJunctionsID=[];
                end
                if isempty(nodes.BinNodeReservoirNameID)
                    nodes.BinReservoirsID=[];
                end
                if isempty(nodes.BinNodeTankNameID)
                    nodes.BinNodeNameID=[];
                end
                if (gg==length(nodes.BinNodeJunctionNameID)+length(nodes.BinNodeReservoirNameID)-1) && (typecode==1) || (gg==length(nodes.BinNodeNameID)-1) && (typecode==2)
                    fprintf(fid2, '\r\n');
                    fprintf(fid2, '%s%s%d%s%d', newID, sps1, X, sps1, Y);
                    gg=0; onetime=0;
                end
            end
            if qualch==1 && qq==1
                qualch=0;
            end
            fprintf(fid2, '\n');
        end
    end
    fclose(fid2);
    if obj.Bin, Errcode=reloadNetwork(obj); end
end
function Errcode=addLink(obj, typecode, newLink, fromNode, toNode, varargin)
% Link type codes consist of the following constants:
% CVPIPE 0 pipe
% Check Valve 1 pipe
% PUMP 2
% PRV Pressure Reducing Valve 3
% PSV Pressure Sustaining Valve 4
% PBV Pressure Breaker Valve 5
% FCV Flow Control Valve 6
% TCV Throttle Control Valve 7
% GPV General Purpose Valve 8
% Initial PIPE plength, value for length of new pipe pdiameter,
% value for diameter of new pipe proughness,  value for roughness of new pipe
if typecode==1, status='Open';end
if ~typecode
    status='CV';
    typecode=1;
end
if typecode==1 && nargin>5
    plength=varargin{1};
    pdiameter=varargin{2};
    proughness=varargin{3};
elseif typecode==2
    if ~isnumeric(varargin{1})
        curveID=varargin{1};
    else
        power=varargin{1};
        curveID='';
    end
elseif typecode>2
    type_valv = obj.TYPELINK{typecode+1};
    if typecode>2, typecode=3; end
    vdiameter=varargin{1};
    vsetting=varargin{2};
end
[Errcode]=addLinkWarnings(obj, typecode, newLink, toNode);
crvs = obj.getBinCurvesInfo;
% Open and read inpname
% Read all file and save in variable info
[~, info] = obj.readInpFile;
fid2 = writenewTemp(obj.BinTempfile);
% Add pipe
nn=0;sps=blanks(10);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            cnt=bracketsCheck(a{u});
            if (cnt==2 && strcmp(a{u}, '[PIPES]') && nn==0 && typecode==1)
                fprintf(fid2, '%s', a{u});
                fprintf(fid2, '\n%s%s%s%s%s%s%d%s%d%s%d%s%d%s%s', newLink, sps, fromNode, sps, ...
                    toNode, sps, plength, sps, pdiameter, sps, proughness, sps, 0, sps, status);

            elseif (cnt==2 && strcmp(a{u}, '[PUMPS]') && nn==0 && typecode==2)
                if ~isempty(curveID)
                    if isempty(char(crvs.BinCurveNameID))
                        warning('No head curve supplied for pump %s.', newLink);
                        return;
                    end
                    fprintf(fid2, '%s', a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%s%s%s', newLink, sps, fromNode, sps, ...
                        toNode, sps, 'HEAD', sps, curveID);
                else
                    fprintf(fid2, '%s', a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%s%s%.2f', newLink, sps, fromNode, sps, ...
                        toNode, sps, 'POWER', sps, power);
                end
            elseif typecode==3 && strcmp(a{u}, '[VALVES]')
                fprintf(fid2, '%s', a{u});
                fprintf(fid2, '\n%s%s%s%s%s%s%d%s%s%s%s', newLink, sps, fromNode, sps, ...
                    toNode, sps, vdiameter, sps, type_valv, sps, num2str(vsetting));
                nn=1;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
if obj.Bin, Errcode=reloadNetwork(obj); end
end
function [Errcode] = rmNode(obj, NodeID)
% Remove node from the network.
% Check if id new already exists
nodes = obj.getBinNodesInfo;Errcode=0;
if isempty(nodes.BinNodeNameID), return; end
if ~ismember(NodeID, nodes.BinNodeNameID)
    warning('There is no such object in the network.');
    Errcode=-1; return;
end
% if ismember(NodeID, nodes.BinNodeReservoirNameID) || ismember(NodeID, nodes.BinNodeTankNameID)
%     if (nodes.BinNodeReservoirCount+nodes.BinNodeTankCount-1)==0;
%         warning('This tank/reservoir has not removed.');
%         Errcode=-1; return;
%     end
% end
% Get links which delete with function Remove Link
links = obj.getBinLinksInfo;
a=strcmp(links.BinLinkFromNode, NodeID);
linkindex1=find(a);
b=strcmp(links.BinLinkToNode, NodeID);
linkindex2=find(b);
linkindex12=[linkindex1 linkindex2];
checklinks_index=unique(linkindex12);
checklinks=links.BinLinkNameID(checklinks_index);
obj.removeBinControlNodeID(NodeID);% Remove control, code 0(NODE)
obj.removeBinRulesControlNodeID(NodeID); %Remove Rule
[~, info] = obj.readInpFile;
fid2 = writenewTemp(obj.BinTempfile);
out=0; sps=blanks(10);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;x=0;xx=0;q=0;
        while u < length(a)+1
            if isempty(a{u}) && (x==0)
                u=u+1; x=1;xx=1;
                if u==length(a)+1, break; end
            end
            if strcmp(a{u}, '[PIPES]'), out=1; end
            if strcmp(a{u}, '[DEMANDS]'), out=0; end %out=0; delete line
            if strcmp(a{u}, '[PATTERNS]'), out=1; end
            if strcmp(a{u}, '[QUALITY]'), out=0; end
            if strcmp(a{u}, '[SOURCES]'), out=1; end
            if strcmp(a{u}, '[MIXING]'), out=0; end
            if strcmp(a{u}, '[COORDINATES]'), out=0; end
            if strcmp(a{u}, NodeID) && q~=1 && out==0
                if xx==1 || strcmp(a{u}, NodeID)
                    u=length(a)+1;
                end
            else
                q=1;
                fprintf(fid2, '%s%s', a{u}, sps);
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
% Remove links
for i=1:length(checklinks)
    obj.removeBinLinkID(checklinks{i});
end
% Find who other id must be delete
remove_link={''};
remove_link_index = zeros(1, length(links.BinLinkFromNode));
for i=1:length(checklinks)
    remove_link(i)=checklinks(i);
    remove_link_index(i)=i;
    warning('Removed link:%s', char(remove_link(i)));
end

if obj.Bin, Errcode=reloadNetwork(obj); end
end
function Errcode=rmRulesControl(obj, type, id)
% Remove control from the network.
exists=0;Errcode=0;exists1=0;
rulescontrols = obj.getBinRulesControlsInfo;
if type
    if isempty(rulescontrols.BinRulesControlLinksID)
        warning('There is no rule object in the network.');
        Errcode=-1;return;
    end
    for i=length(rulescontrols.BinRulesControlLinksID):-1:1
        exists(i, :) = strcmp(rulescontrols.BinRulesControlLinksID{i}{length(rulescontrols.BinRulesControlLinksID{1})}, char(id));
    end
else
    if isempty(rulescontrols.BinRulesControlNodesID)
        warning('There is no such rule in the network.');
        Errcode=-1;return;
    end
    for i=length(rulescontrols.BinRulesControlNodesID):-1:1
        exists1(i) = strcmp( rulescontrols.BinRulesControlNodesID{i}{length(rulescontrols.BinRulesControlNodesID{1})}, char(id));
    end
end
if ~sum(exists) && ~sum(exists1)
    warning('There is no such rule in the network.');
    Errcode=-1; return;
end
[addSectionCoordinates, addSectionRules] = obj.getBinCoordRuleSections(obj.BinTempfile);
cntRules = cellfun('length', rulescontrols.BinRulesControlsInfo);
endInpIndex=find(~cellfun(@isempty, regexp(addSectionCoordinates, 'END', 'match')));
[~, info] = obj.readInpFile;
info(find(~cellfun(@isempty, regexp(info, 'END', 'match'))))='';
f1 = writenewTemp(obj.BinTempfile);
rulesSectionIndex=find(~cellfun(@isempty, regexp(info, 'RULES', 'match')));
if ~isempty(rulesSectionIndex)
    fprintf(f1, '%s\n', info{1:rulesSectionIndex-1});
    fprintf(f1, '[RULES]\n');
end
if type
    for i=length(rulescontrols.BinRulesControlLinksID):-1:1
        if ~exists(i)
            for j=1:length(rulescontrols.BinRulesControlsInfo{i})
                fprintf(f1, '%s\n', addSectionRules{sum(cntRules(1:i-1))+1+j});
            end
        end
    end
else
    for i=length(rulescontrols.BinRulesControlNodesID):-1:1
        if ~exists1(i)
            for j=1:length(rulescontrols.BinRulesControlsInfo{i})
                fprintf(f1, '%s\n', addSectionRules{sum(cntRules(1:i-1))+1+j});
            end
        end
    end
end
if ~isempty(addSectionCoordinates) % && isempty(coordSectionIndex)
    fprintf(f1, '%s\n', addSectionCoordinates{:});
end
if isempty(endInpIndex), fprintf(f1, '[END]\n'); end
fclose(f1);
if obj.Bin, Errcode=reloadNetwork(obj); end
end
function Errcode=rmControl(obj, type, id)
% Remove control from the network.
Errcode=0;
controls = obj.getBinControlsInfo;
if type
    if isempty(controls.BinControlLinksID)
        warning('There is no such control in the network.');
        Errcode=-1; return;
    end
else
    if isempty(controls.BinControlNodesID)
        warning('There is no such control in the network.');
        Errcode=-1; return;
    end
end
[~, info] = obj.readInpFile;
fid2 = writenewTemp(obj.BinTempfile);
e=0;n=0;kk=1;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            rr = regexp(a, '\w*[\w*]\w*', 'split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets, '[');
            ch2 = strcmp(check_brackets, ']');

            if strcmp(a{u}, '[CONTROLS]')
                fprintf(fid2, '%s', a{u});
                n=1;
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1), break; end
                e=1;
            end
            if strcmp(a{u}, '[END]'),  e=1; fprintf(fid2, '%s', a{u});break;   end

            if n==1 && e==0 && kk==1
                if strcmp(a{u}, '[CONTROLS]'), break; end
                if isempty(a{u})
                elseif strfind(a{u}, ';')
                    break;
                else
                    if type==1
                        tt = strcmp(a{u+1}, id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2, '%s%s', a{u}, sps);
                        end
                    elseif type==0
                        tt = strcmp(a{u+5}, id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2, '%s%s', a{u}, sps);
                        end
                    end
                end
            else
                if isempty(a{u})
                else
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');kk=1;
end
fclose(fid2);
if obj.Bin==1
    Errcode=reloadNetwork(obj);
end
end
function [Errcode] = rmLink(obj, LinkID)
% Remove link from the network.
% Check if id new already exists
links = obj.getBinLinksInfo;Errcode=0;
if isempty(links.BinLinkNameID)
    warning('There is no such object in the network.');
    Errcode=-1; return;
end
if ~ismember(LinkID, links.BinLinkNameID)
    warning('There is no such object in the network.');
    Errcode=-1; return;
else
    index_rmvlink = find(strcmp(LinkID, links.BinLinkNameID));
end
nodes = obj.getBinNodesInfo;
from_node = links.BinLinkFromNode(index_rmvlink);
r = strcmp(nodes.BinNodeNameID, from_node);
if sum(r)==0, from_node=''; end
to_node = links.BinLinkToNode(index_rmvlink);
r = strcmp(nodes.BinNodeNameID, to_node);
if sum(r)==0, to_node=''; end
% Remove control, code 1(LINK)
obj.removeBinControlLinkID(LinkID);
obj.removeBinRulesControlLinkID(LinkID); %Remove Rule

[~, info] = obj.readInpFile;
fid2 = writenewTemp(obj.BinTempfile);

% section [JUNCTIONS]
out=0;YY=0;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;x=0;xx=0;q=0;
        while u < length(a)+1
            if strcmp(a{u}, '[PIPES]'), YY=1;end
            if YY==1
                if isempty(a{u}) && (x==0)
                    u=u+1; x=1;xx=1;
                    if u==length(a)+1
                        break
                    end
                end
                if strcmp(a{u}, '[TAGS]'), out=1; end
                if strcmp(a{u}, '[STATUS]'), out=1; end
                if strcmp(a{u}, '[DEMANDS]'), out=1; end
                if strcmp(a{u}, '[PATTERNS]'), out=1; end

                if strcmp(a{u}, LinkID) && q~=1 && out==0
                    if xx==1 || strcmp(a{1}, LinkID)
                        u=length(a)+1;
                    end
                else
                    q=1;
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            else
                if isempty(a{u})
                    u=u+1;
                    if u==length(a)+1
                        break
                    end
                end
                fprintf(fid2, '%s%s', a{u}, sps);
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
% Get nodes which delete with function Remove Node
links = obj.getBinLinksInfo;
if ~ismember(from_node, [links.BinLinkToNode links.BinLinkFromNode]) && ~isempty(from_node)
    warning('Node %s disconnected.', char(from_node));
end
if ~ismember(to_node, [links.BinLinkToNode links.BinLinkFromNode]) && ~isempty(from_node)
    warning('Node %s disconnected.', char(from_node));
end
if obj.Bin,  Errcode=reloadNetwork(obj); end
end
function [Errcode]=addNewControl(obj, x, status, y_t_c, param, z, varargin)
% syntax
Errcode=0;
if (nargin==6)
    syntax = ['LINK ', x, ' ', status, ' IF NODE ', y_t_c, ' ', param, ' ', num2str(z)];
elseif (nargin==5)
    syntax = ['LINK ', x, ' ', status, ' AT CLOCKTIME ', y_t_c, ' ', param];
elseif (nargin==4)
    syntax = ['LINK ', x, ' ', status, ' AT TIME ', y_t_c];
end
if (nargin==6)
    % Check if id new already exists
    if ~ismember(x, obj.getBinNodesInfo.BinNodeNameID)
        warning('There is no such object in the network.');
        Errcode=-1; return;
    end
end
if (nargin==2)
   controls = x;
else
    % Check if id new already exists
    if ~ismember(x, obj.getBinLinksInfo.BinLinkNameID)
        warning('There is no such object in the network.');
        Errcode=-1; return;
    end
end
type_n='[CONTROLS]';
[~, info] = obj.readInpFile;
m = strfind(info, type_n);
Index = find(not(cellfun('isempty', m)));
fid2 = writenewTemp(obj.BinTempfile);
noo=0;s=0;sps=blanks(15);goOut=0;
for i=1:Index-1
    fprintf(fid2, '%s', info{i});
    fprintf(fid2, '\n');
end
for t = Index:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u}, type_n)
                fprintf(fid2, '[CONTROLS]');
                s=1; break;
            end
            if (s==1) && (noo==0)
                if (nargin==2)
                    for i=1:size(controls, 1)
                        fprintf(fid2, controls(i, :));
                        fprintf(fid2, '\r\n');
                    end
                    for i=t:length(info)
                        fprintf(fid2, '%s', info{i});
                        fprintf(fid2, '\n');
                    end
                    goOut=1;
                end
                if ~goOut
                    fprintf(fid2, '%s', syntax);
                    fprintf(fid2, '\r\n');
                    fprintf(fid2, c);
                    noo=1;
                end
                break;
            elseif isempty(a{u}) && noo==0
            else
                if isempty(a{u}) && noo==1
                else
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    if goOut, break; end
    fprintf(fid2, '\n');
end
fclose(fid2);
if obj.Bin, Errcode=reloadNetwork(obj); end
end
function [Errcode]=rmCurveID(obj, CurveID, varargin)
% Check if id new already exists
Errcode=0;
if ~ismember(CurveID, obj.getBinCurvesInfo.BinCurveNameID)
    warning('There is no such object in the network.');
    Errcode=-1; return;
end
value=obj.getBinLinksInfo;
indCurve = find(strcmp(CurveID, value.BinLinkPumpCurveNameID), 1);
if ~isempty(indCurve)
    warning('Pump %s refers to undefined curve.', value.BinLinkPumpNameID{indCurve});
end
% Open and read inpname
% Read all file and save in variable info
[~, info] = obj.readInpFile;
fid2 = writenewTemp(obj.BinTempfile);
e=0;n=0;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*', 'split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            rr = regexp(a, '\w*[\w*]\w*', 'split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets, '[');
            ch2 = strcmp(check_brackets, ']');
            if strcmp(a{u}, '[CURVES]')
                fprintf(fid2, '%s', a{u});
                n=1;
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1), break; end
                e=1;
            end
            if strcmp(a{u}, '[END]'), e=1; fprintf(fid2, '%s', a{u});break;end
            if n==1 && e==0
                if strcmp(a{u}, '[CURVES]'), break; end
                if isempty(a{u})
                    u=u+1;continue;
                elseif strfind(a{u}, ';')
                    ee=regexp(c, '\w*EFFICIENCY*\w', 'match');
                    nn=regexp(c, '\w*VOLUME*\w', 'match');
                    kk=regexp(c, '\w*HEADLOSS*\w', 'match');
                    if length(strcmp(ee, 'EFFICIENCY')) || length(strcmp(nn, 'VOLUME')) || length(strcmp(kk, 'HEADLOSS')) || length(strcmp(a{1}, ';PUMP:'))
                        fprintf(fid2, '%s%s', a{u}, sps);
                    else
                        break;
                    end
                else
                    tt = strcmp(a{u}, CurveID);
                    if tt==1
                        u = length(a)+1;
                    else
                        fprintf(fid2, '%s%s', a{u}, sps);
                    end
                end
            else
                if isempty(a{u})
                else
                    if strcmp(a{u}, '[CURVES]'), break; end
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
if obj.Bin, Errcode=reloadNetwork(obj); end
end
function [Errcode]=Options(obj, newFlowUnits, headloss, varargin)
% Notes: Flow units codes are as follows:
% CFS cubic feet per second
% GPM gallons per minute
% MGD million gallons per day
% IMGD Imperial mgd
% AFD acre-feet per day
% LPS liters per second
% LPM liters per minute
% MLD million liters per day
% CMH cubic meters per hour
% CMD cubic meters per day
value=obj.getBinOptionsInfo;Errcode=0;
previousFlowUnits=value.BinLinkFlowUnits;
US_Customary=0;
SI_Metric=0;
switch newFlowUnits
    case 'CFS',  US_Customary=1;
    case 'GPM',  US_Customary=1;
    case 'MGD',  US_Customary=1;
    case 'IMGD', US_Customary=1;
    case 'AFD',  US_Customary=1;
    case 'LPS',  SI_Metric=1;
    case 'LPM',  SI_Metric=1;
    case 'MLD',  SI_Metric=1;
    case 'CMH',  SI_Metric=1;
    case 'CMD',  SI_Metric=1;
end
if US_Customary==value.BinUnits_US_Customary
    changes=0; US_Customary=0;
    SI_Metric=0; % feet to feet
elseif SI_Metric==value.BinUnits_SI_Metric
    changes=1; US_Customary=0;
    SI_Metric=0; % meter to meter
elseif value.BinUnits_US_Customary==1 && US_Customary==0
    changes=1; % feet to meter or cubic feet to cubic meter
elseif value.BinUnits_US_Customary==0 && US_Customary==1
    changes=2; % meter to feet or cubic meter to cubic feet
end
Units=US_Customary+SI_Metric;
variables=who;nheadl=0;
if ~sum(strcmp('headloss', variables))
    headloss=value.BinOptionsHeadloss;
    nheadl=1;
end
nodes = obj.getBinNodeNameID;
links = obj.getBinLinksInfo;
controls = obj.getBinControlsInfo;
curves = obj.getBinCurvesInfo;
rules=obj.getBinRulesControlsInfo;

[info] = readAllFile(obj.BinTempfile);
fid2 = writenewTemp(obj.BinTempfile);
sect=0;
nn=0;pp=1;sps=blanks(15);
for t = 1:length(info)
    a = info{t};
    c = cell2mat(a);
    if isempty(a)
        % skip
    elseif isempty(c)
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u}, '[JUNCTIONS]') && Units
                fprintf(fid2, '[JUNCTIONS]');
                sect=1;
                break;
            elseif strcmp(a{u}, '[RESERVOIRS]') && Units
                fprintf(fid2, '[RESERVOIRS]');
                nn=0;pp=1;
                sect=2;
                break;
            elseif strcmp(a{u}, '[TANKS]') && Units
                fprintf(fid2, '[TANKS]');
                nn=0;pp=1;
                sect=3;
                break;
            elseif strcmp(a{u}, '[PIPES]') && Units
                fprintf(fid2, '[PIPES]');
                nn=0;pp=1;
                sect=4;
                break;
            elseif strcmp(a{u}, '[PUMPS]') && Units
                fprintf(fid2, '[PUMPS]');
                nn=0;pp=1;
                sect=5;
                break;
            elseif strcmp(a{u}, '[VALVES]') && Units
                fprintf(fid2, '[VALVES]');
                nn=0;pp=1;
                sect=6;
                break;
            elseif strcmp(a{u}, '[DEMANDS]') && ((Units || ~changes) && nheadl)
                fprintf(fid2, '[DEMANDS]');
                nn=0;pp=1;
                sect=7;
                break;
            elseif strcmp(a{u}, '[EMITTERS]') && ((Units || ~changes) && nheadl)
                fprintf(fid2, '[EMITTERS]');
                nn=0;pp=1;
                sect=8;
                break;
            elseif strcmp(a{u}, '[STATUS]')
                fprintf(fid2, '[STATUS]');
                nn=0;pp=1;
                sect=9;
                break;
            elseif strcmp(a{u}, '[PATTERNS]')
                fprintf(fid2, '[PATTERNS]');
                nn=1;
                sect=10;
                break;
            elseif strcmp(a{u}, '[CURVES]') && ((Units || ~changes) && nheadl)
                fprintf(fid2, '[CURVES]');
                nn=0;pp=1;ww=1;
                sect=11;
                break;
            elseif strcmp(a{u}, '[CONTROLS]') && Units
                fprintf(fid2, '[CONTROLS]');
                nn=0;pp=1;
                sect=12;
                break;
            elseif strcmp(a{u}, '[RULES]') && Units
                fprintf(fid2, '[RULES]');
                nn=0;pp=1;
                sect=13;
                break;
            elseif strcmp(a{u}, '[OPTIONS]')
                fprintf(fid2, '[OPTIONS]');
                sect=14;nn=0;
                break;
            end
                % section [JUNCTIONS]
            if (sect==1) && (nn==0)
                mm=1;
                if pp<length(char(nodes.BinNodeJunctionNameID))+1
                    if strcmp(a{mm}, nodes.BinNodeJunctionNameID{pp})
                        pp=pp+1;
                        fprintf(fid2, '%s%s', char(a{mm}), sps);
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                        end
                        if length(a)>2
                            mm=2;
                            setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                        end
                    end
                else
                    nn=1;
                end
                break;
                % section [RESERVOIRS]
            elseif (sect==2) && (nn==0)
                mm=1;
                if pp<length(char(nodes.BinNodeReservoirNameID))+1
                    if strcmp(a{mm}, nodes.BinNodeReservoirNameID{pp})
                        pp=pp+1;
                        fprintf(fid2, '%s%s', char(a{mm}), sps);
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                        end
                    end
                else
                    nn=1;
                end
                break;
                % section [TANKS]
            elseif (sect==3) && (nn==0)
                mm=1;
                if pp<length(char(nodes.BinNodeTankNameID))+1
                    if strcmp(a{mm}, nodes.BinNodeTankNameID{pp})
                        pp=pp+1;
                        fprintf(fid2, '%s%s', char(a{mm}), sps);
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+2})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+3})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+4})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+5})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+6})*0.02831685, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+2})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+3})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+4})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+5})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+6})*35.3147, sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2, '%s%s', char(a{1}), sps);
                end
                break;
                % section [PIPES]
            elseif (sect==4) && (nn==0)
                mm=1;
                if pp<length(char(links.BinLinkPipeNameID))+1
                    if strcmp(a{mm}, links.BinLinkPipeNameID{pp})
                        pp=pp+1;
                        for mm=mm:mm+2
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+2})*25.4, sps);
                            if nheadl
                                if strcmp('D-W', value.BinOptionsHeadloss)
                                    fprintf(fid2, '%.6f%s', str2double(a{mm+3})*0.3048, sps);
                                    mm=7;
                                else
                                    mm=6;
                                end
                            end
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.03937007874, sps);
                            if nheadl
                                if strcmp('D-W', value.BinOptionsHeadloss)
                                    mm=7;
                                else
                                    mm=6;
                                end
                            end
                        end
                        for mm=mm:length(a)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2, '%s%s', char(a{1}), sps);
                end
                break;
            % section [PUMPS]
            elseif (sect==5) && (nn==0)
                mm=1;
                if pp<length(char(links.BinLinkPumpNameID))+1
                    if strcmp(a{mm}, links.BinLinkPumpNameID{pp})
                        pp=pp+1;
                        for mm=mm:length(a)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        power=regexp(c, 'POWER', 'match');
                        if strcmpi(power, 'POWER')
                            mm=mm-1;
                            if changes==1
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.745699882507324, sps);
                            elseif changes==2
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1})/0.745699882507324, sps);
                            end
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2, '%s%s', char(a{1}), sps);
                end
                break;
            % section [VALVES]
            elseif (sect==6) && (nn==0)
                mm=1;
                if pp<length(char(links.BinLinkValveNameID))+1
                    if strcmp(a{mm}, links.BinLinkValveNameID{pp})
                        pp=pp+1;
                        for mm=mm:(mm+2)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        prv=regexp(c, 'PRV', 'match');if isempty(prv), prv=0; end
                        psv=regexp(c, 'PSV', 'match');if isempty(psv), psv=0; end
                        pbv=regexp(c, 'PBV', 'match');if isempty(pbv), pbv=0; end
                        fcv=regexp(c, 'FCV', 'match');if isempty(fcv), fcv=0; end
%                         tcv=regexp(c, 'TCV', 'match');if isempty(tcv), tcv=0; end
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*25.4, sps);
                            fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                            if strcmpi(prv, 'PRV') || strcmpi(psv, 'PSV') || strcmpi(pbv, 'PBV') %|| strcmpi(tcv, 'TCV')
                                fprintf(fid2, '%s%s', num2str(str2double(a{mm+3})*0.3048), sps);
                            elseif strcmpi(fcv, 'FCV')
                                setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                            end
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.03937007874, sps);
                            fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                            if strcmpi(prv, 'PRV') || strcmpi(psv, 'PSV') || strcmpi(pbv, 'PBV') %|| strcmpi(tcv, 'TCV')
                                fprintf(fid2, '%s%s', num2str(str2double(a{mm+3})/0.3048), sps);
                            elseif strcmpi(fcv, 'FCV')
                                setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                            end
                        end
                        for mm=(mm+4):length(a)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2, '%s%s', char(a{1}), sps);
                end
                break;
                % section [DEMANDS]
            elseif (sect==7) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if pp<length(char(nodes.BinNodeJunctionNameID))+1
                        if strcmp(a{mm}, nodes.BinNodeJunctionNameID{pp})
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                            if length(a)>2
                                fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2, '%s%s', char(a{1}), sps);
                    end
                end
                break;
                % section [EMITTERS]
            elseif (sect==8) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if sum(strcmp(a{mm}, nodes.BinNodeJunctionNameID))
                        fprintf(fid2, '%s%s', char(a{mm}), sps);
                        setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                    end
                end
                break;
                % section [STATUS]
            elseif (sect==9) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if pp<length(char(links.BinLinkInitialStatus))+1
                        if strcmp(a{mm}, links.BinLinkNameID{pp})
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            if length(a)==2
                                fprintf(fid2, '%s%s', char(a{mm+1}), sps);
                            end
                            if length(a)>2
                                fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2, '%s%s', char(a{1}), sps);
                    end
                end
                break;

                % section [CURVES]
            elseif (sect==11) && (nn==0)
                mm=1;
                if strfind(c, ';ID')
                    break;
                end
                if pp<length(curves.BinCurveAllLines)+1 && ~isempty(char(a)) % PUMP % EFFICIENCY % VOLUME
                    if ww<length(curves.BinCTypes)+1
                        if curves.BinCTypes(ww)==0
                            if strfind(c, ';PUMP:')
                                fprintf(fid2, c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm);
                            if changes==1
                                fprintf(fid2, '%.6f%s', str2double(a{mm+2})*0.3048, sps);
                            elseif changes==2
                                fprintf(fid2, '%.6f%s', str2double(a{mm+2})*3.281, sps);
                            else
                                fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                            end
                            break;
                        elseif curves.BintypeCurve(ww)==1
                            ee=regexp(c, '\w*EFFICIENCY*\w', 'match');
                            if length(strcmp(ee, 'EFFICIENCY'))
                                fprintf(fid2, c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                            fprintf(fid2, '%s%s', char(a{mm+2}), sps);
                        elseif curves.BintypeCurve(ww)==2
                            gg=regexp(c, '\w*VOLUME*\w', 'match');
                            if length(strcmp(gg, 'VOLUME'))
                                fprintf(fid2, c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            if changes==1
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1})*2.831685e-02, sps);
                                fprintf(fid2, '%.6f%s', str2double(a{mm+2})*0.3048, sps);
                            else
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
                                fprintf(fid2, '%.6f%s', str2double(a{mm+2}), sps);
                            end
                        elseif curves.BintypeCurve(ww)==3 % HEADLOSS
                            kk=regexp(c, '\w*HEADLOSS*\w', 'match');
                            if length(strcmp(kk, 'HEADLOSS'))
                                fprintf(fid2, c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                            if changes==1
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                            elseif changes==2
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                            else
                                fprintf(fid2, '%s%s', char(a{mm+1}), sps);
                            end
                            mm=mm+1;
                            setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                        end
                        ww=ww+1;
                    end
                    pp=pp+1;
                else
                    if ~(ww<length(curves.BinCTypes)+1), nn=1; end
                end
                if ~isempty(regexp(a{mm}, '[\w]*', 'match'))
                    nn=1;
                    fprintf(fid2, '%s%s', char(c), sps);
                end
                break;
                % section [CONTROLS]
            elseif (sect==12) && (nn==0)
                e=regexp(a, ';', 'match');
                if length(e)>0
                    if ~isempty(e{1})
                        break;
                    end
                end
                mm=1;
                if pp<length(controls.BinControlsInfo)+1
                    pp=pp+1;
                    if length(a)>7
                        if strcmpi(a{mm+6}, 'BELOW') || strcmpi(a{mm+6}, 'ABOVE')
                            for mm=mm:(mm+6)
                                fprintf(fid2, '%s%s', char(a{mm}), sps);
                            end
                            v=find(strcmp(a{2}, obj.BinLinkNameID));
                            index=obj.getBinLinkIndex(obj.BinLinkNameID(v));
                            if ~strcmp(obj.BinLinkType(index), 'TCV') && ~strcmp(obj.BinLinkType(index), 'GPV')
                                if changes==1
                                    if strcmp(obj.BinLinkType(index), 'FCV')
                                        setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                                    else
                                        fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                                    end
                                elseif changes==2
                                    if strcmp(obj.BinLinkType(index), 'FCV')
                                        setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                                    else
                                        fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                                    end
                                end
                            else
                                fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
                            end
                        else
                            for mm=mm:length(a)
                                fprintf(fid2, '%s%s', char(a{mm}), sps);
                            end
                        end
                    else
                        for mm=mm:length(a)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2, '%s%s', char(a{1}), sps);
                end
                break;
                % section [RULES]
            elseif (sect==13) && (nn==0)
                mm=1;
                if pp<rules.BinRulesCount+1
                    pp=pp+1;
                    if strcmpi(regexp(cell2mat(a), '\s*LEVEL*', 'match'), 'LEVEL') %||
                        for mm=mm:(mm+4)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                        end
                    elseif strcmpi(regexp(cell2mat(a), '\s*HEAD*', 'match'), 'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3048, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.281, sps);
                        end
                    elseif strcmpi(regexp(cell2mat(a), '\s*DEMAND*', 'match'), 'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
                    elseif strcmpi(regexp(cell2mat(a), '\s*PRESSURE*', 'match'), 'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                        if changes==1
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})/1.422, sps);
                        elseif changes==2
                            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.422, sps);
                        end
                    else
                        for mm=mm:length(a)
                            fprintf(fid2, '%s%s', char(a{mm}), sps);
                        end
                    end
                else
                    nn=1;
                    for mm=mm:length(a)
                        fprintf(fid2, '%s%s', char(a{mm}), sps);
                    end
                end
                break;
                % section [OPTIONS]
            elseif (sect==14) && (nn==0)
                mm=1;
                if strcmpi(a{mm}, 'UNITS')
                    fprintf(fid2, '%s%s', char(a{mm}), sps);
                    if nheadl
                        fprintf(fid2, '%s%s', char(newFlowUnits), sps);nn=1;
                    else
                        fprintf(fid2, '%s%s', char(previousFlowUnits), sps);
                    end
                elseif strcmpi(a{mm}, 'HEADLOSS')
                    fprintf(fid2, '%s%s', char(a{mm}), sps);
                    fprintf(fid2, '%s%s', char(headloss), sps);
                    nn=1;
                else
                    fprintf(fid2, c);
                end
                break;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2, '%s%s', a{u}, sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2, '\n');
end
fclose(fid2);
if obj.Bin, Errcode=reloadNetwork(obj); end
end
function Errcode=reloadNetwork(obj)
%     obj.closeNetwork;
    Errcode=ENopen([obj.BinTempfile], [obj.BinTempfile(1:end-4), '.txt'], [obj.BinTempfile(1:end-4), '.bin'], obj.LibEPANET);
end
function setflow(previousFlowUnits, newFlowUnits, fid2, a, sps, mm)
if isnan(str2double(a{mm+1}))
    return;
end
if strcmp(previousFlowUnits, 'GPM')
    switch newFlowUnits %(GPM)
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.00222816399286988, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.00144, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.00119905, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.004419191, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0630902, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.785412, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.005450993, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.2271247, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*5.450993, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'CFS')
    switch newFlowUnits
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*448.8312, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.6463169, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.5381711, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.983471, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*28.31685, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1899.011, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*2.446576, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*101.9406, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*2446.576, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'MGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.547229, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*694.4445, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.8326738, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.068883, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*43.81264, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*2628.758, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.785412, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*157.7255, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3785.412, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'IMGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.858145, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*833.9936, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.200951, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.685577, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*52.61681, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3157.008, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*4.546092, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*189.4205, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*4546.092, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'AFD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.5041667, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*226.2857, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.3258514, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.271328, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*14.27641, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*856.5846, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.233482, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*51.39508, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1233.482, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'LPS')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.03531466, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*15.85032, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.02282446, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.01900533, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.07004562, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*60, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0864, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*3.6, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*86.4, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'LPM')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0005885777, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.264172, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0003804078, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0003167556, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0011674272, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.01666667, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.00144, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.06, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1.44, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'MLD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.4087345, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*183.4528, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.264172, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.2199692, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.8107132, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*11.57407, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*694.4445, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*41.66667, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*1000, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'CMH')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.009809635, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*4.402868, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.006340129, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.00527926, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.01945712, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.2777778, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*16.66667, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.024, sps);
        case 'CMD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*24, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
elseif strcmp(previousFlowUnits, 'CMD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0004087345, sps);
        case 'GPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.1834528, sps);
        case 'MGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.000264172, sps);
        case 'IMGD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0002199692, sps);
        case 'AFD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.0008107132, sps);
        case 'LPS'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.01157407, sps);
        case 'LPM'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.6944444, sps);
        case 'MLD'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.001, sps);
        case 'CMH'
            fprintf(fid2, '%.6f%s', str2double(a{mm+1})*0.04166667, sps);
        otherwise
            fprintf(fid2, '%.6f%s', str2double(a{mm+1}), sps);
    end
end
end
function [fid, binfile, rptfile] = runEPANETexe(obj)
    arch = computer('arch');
    [inpfile, rptfile, binfile]= createTempfiles(obj.BinTempfile);
    if strcmp(arch, 'win64') || strcmp(arch, 'win32')
        r = sprintf('"%s%s.exe" "%s" %s %s & exit', obj.LibEPANETpath, obj.LibEPANET, inpfile, rptfile, binfile);
    end
    if isunix
        r = sprintf('%s%s %s %s %s', obj.LibEPANETpath, obj.LibEPANET, obj.BinTempfile, rptfile, binfile);
    elseif ismac
        r = sprintf('%s%s %s %s %s', obj.LibEPANETpath, obj.LibEPANET, obj.BinTempfile, rptfile, binfile);
    end
    if obj.getCMDCODE, [~, ~]=system(r); else system(r); end
    fid = fopen(binfile, 'r');
end
function value = getBinComputedTimeSeries(obj, indParam, varargin)
    [fid, binfile, rptfile] = runEPANETexe(obj);
    value=[];
    if fid~=-1
        data = fread(fid, 'int32');
        BinNodeCount=data(3);
        BinNodeResTankCount=data(4);
        BinLinkCount=data(5);
        BinLinkPumpCount=data(6);
        NumberReportingPeriods = data(end-2);
        if indParam==27
            value = NumberReportingPeriods;
        end
        if indParam==28
            value = data(15); % simulation duration
        end
        clear data;
        % Beginning of file
        fseek(fid, 0, 'bof');
        fread(fid, 15, 'uint32');
        fread(fid, 808, '*char');
        fread(fid, 4, 'uint32');
        fread(fid, 32*BinNodeCount+32*BinLinkCount, '*char');
        fread(fid, BinLinkCount*3, 'uint32');
        fread(fid, BinNodeResTankCount, 'uint32');
        fread(fid, BinNodeResTankCount, 'float');
        switch indParam
            case 1
                value = fread(fid, BinNodeCount, 'float')'; % ElevationEachNode
            case 2
                fread(fid, BinNodeCount, 'float');
                value = fread(fid, BinLinkCount, 'float')'; % LengthEachLink
            case 3
                fread(fid, BinNodeCount+BinLinkCount, 'float');
                value = fread(fid, BinLinkCount, 'float')'; % DiameterEachLink
            case 4
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    value(p) = fread(fid, 1, 'float')'; % PumpIndexListLinks
                    fread(fid, 6, 'float');
                end
            case 5
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 1, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinPumpUtilization
                    fread(fid, 5, 'float');
                end
            case 6
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 2, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinAverageEfficiency
                    fread(fid, 4, 'float');
                end
            case 7
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 3, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinAverageKwattsOrMillionGallons
                    fread(fid, 3, 'float');
                end
            case 8
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 4, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinAverageKwatts
                    fread(fid, 2, 'float');
                end
            case 9
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 5, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinPeakKwatts
                    fread(fid, 1, 'float');
                end
            case 10
                fread(fid, BinNodeCount+BinLinkCount*2, 'float');
                for p=1:BinLinkPumpCount
                    fread(fid, 6, 'float');
                    value(p) = fread(fid, 1, 'float')';  % BinAverageCostPerDay
                end
        end
        if indParam>10
            fread(fid, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*7+1, 'float');
            for i=1:NumberReportingPeriods
                switch indParam
                    case 11
                        value(i, :) = fread(fid, BinNodeCount, 'float')'; % nodeDemand
                        fread(fid, BinNodeCount*3, 'float');
                        fread(fid, BinLinkCount*8, 'float');
                    case 12
                        fread(fid, BinNodeCount, 'float');
                        value(i, :) = fread(fid, BinNodeCount, 'float')'; % nodeHead
                        fread(fid, BinNodeCount*2, 'float');
                        fread(fid, BinLinkCount*8, 'float');
                    case 13
                        fread(fid, BinNodeCount*2, 'float');
                        value(i, :) = fread(fid, BinNodeCount, 'float')'; % nodePressure
                        fread(fid, BinNodeCount, 'float');
                        fread(fid, BinLinkCount*8, 'float');
                    case 14
                        fread(fid, BinNodeCount*3, 'float');
                        value(i, :) = fread(fid, BinNodeCount, 'float')'; % nodeQuality
                        fread(fid, BinLinkCount*8, 'float');
                    case 15
                        fread(fid, BinNodeCount*4, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkFlow
                        fread(fid, BinLinkCount*7, 'float');
                    case 16
                        fread(fid, BinNodeCount*4+BinLinkCount, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkVelocity
                        fread(fid, BinLinkCount*6, 'float');
                    case 17
                        fread(fid, BinNodeCount*4+BinLinkCount*2, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkHeadloss
                        fread(fid, BinLinkCount*5, 'float');
                    case 18
                        fread(fid, BinNodeCount*4+BinLinkCount*3, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkQuality
                        fread(fid, BinLinkCount*4, 'float');
                    case 19
                        % Status Code for Each Link
                        % 0 = closed (max. head exceeded)
                        % 1 = temporarily closed
                        % 2 = closed
                        % 3 = open
                        % 4 = active (partially open)
                        % 5 = open (max. flow exceeded)
                        % 6 = open (flow setting not met)
                        % 7 = open (pressure setting not met)
                        fread(fid, BinNodeCount*4+BinLinkCount*4, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkStatus
                        fread(fid, BinLinkCount*3, 'float');
                    case 20
                        fread(fid, BinNodeCount*4+BinLinkCount*5, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkSetting
                        fread(fid, BinLinkCount*2, 'float');
                    case 21
                        fread(fid, BinNodeCount*4+BinLinkCount*6, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkReactionRate
                        fread(fid, BinLinkCount, 'float');
                    case 22
                        fread(fid, BinNodeCount*4+BinLinkCount*7, 'float');
                        value(i, :) = fread(fid, BinLinkCount, 'float')'; % linkFrictionFactor
                end
                if indParam>22
                    fread(fid, BinNodeCount*4+BinLinkCount*8, 'float');
                end
            end
        end
        switch indParam
            case 23
                value =fread(fid, 1, 'float')'; % AverageBulkReactionRate
            case 24
                fread(fid, 1, 'float');
                value =fread(fid, 1, 'float')'; % AverageWallReactionRate
            case 25
                fread(fid, 2, 'float');
                value =fread(fid, 1, 'float')'; % AverageTankReactionRate
            case 26
                fread(fid, 3, 'float');
                value =fread(fid, 1, 'float')'; % AverageSourceInflowRate
        end
    end
    warning('off'); try fclose(fid); catch, end; try delete(binfile); catch, end
    try delete(rptfile); catch, end; warning('on');
end
function Errcode=addLinkWarnings(obj, typecode, newLink, toNode)
% Check if id new already exists
Nodes = obj.getBinNodesInfo;
Errcode=0;
if isempty(Nodes.BinNodeNameID), Errcode=-1; return; end
if ~ismember(toNode, Nodes.BinNodeNameID)
    warning('There is no node "%s" in the network.', toNode);
    Errcode=-1; return;
end
if ~ismember(0:8, typecode)
    warning('Wrong constant type.');
    Errcode=-1; return;
else
%     type_valv = obj.TYPELINK{typecode+1};
    if typecode>2, typecode=3; end
end
% Valve illegally connected to a tank or reservoir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if typecode==3
    if ismember(toNode, Nodes.BinNodeReservoirNameID) || ismember(toNode, Nodes.BinNodeTankNameID)
        Errcode=-1; warning('Valve "%s" illegally connected to a tank.', newLink);
        return;
    end
end
% Check if newLink already exists
Links = obj.getBinLinksInfo;
if ismember(newLink, Links.BinLinkNameID)
    Errcode=-1; warning('Link %s already exists.', newLink); return;
end

if typecode==2
    crvs = obj.getBinCurvesInfo;
    if isempty(char(crvs.BinCurveNameID))
        Errcode=-1; warning('No head curve supplied for pump %s.', newLink); return;
    end
end
end
function cnt=bracketsCheck(v)
    t =  regexp(v, '[(\w*)]', 'split');
    cnt=length(t(~cellfun('isempty', t)));
end
function setMSXOptions(obj, varargin)
solver=obj.getMSXSolver;
areaunits=obj.getMSXAreaUnits;
rateunits=obj.getMSXRateUnits;
rtol=obj.getMSXRtol;
atol=obj.getMSXAtol;
timestep=obj.getMSXTimeStep;
coupling=obj.getMSXCoupling;
compiler=obj.getMSXCompiler;

for i=1:(nargin/2)
    argument =lower(varargin{2*(i-1)+1});
    switch argument
        case 'areaunits'
            areaunits=varargin{2*i};
        case 'rateunits'
            rateunits=varargin{2*i};
        case 'solver'
            solver=varargin{2*i};
        case 'timestep'
            timestep=varargin{2*i};
        case 'atol'
            atol=varargin{2*i};
        case 'rtol'
            rtol=varargin{2*i};
        case 'coupling'
            coupling=varargin{2*i};
        case 'compiler'
            compiler=varargin{2*i};
        otherwise
            warning('Invalid property found.');
            return
    end
end

[info, tline] = readAllFile(obj.MSXFile);
fid2 = writenewTemp(obj.MSXTempFile);
sect=0;
for t = 1:length(info)
    a = info{t};
    c = cell2mat(a);
    if isempty(a)
        % skip
    elseif isempty(c)
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u}, '[OPTIONS]')
                fprintf(fid2, '[OPTIONS]');
                sect=1;
                break;
            elseif strcmp(a{u}, '[SPECIES]')
                fprintf(fid2, '[SPECIES]');
                sect=0;
                break;
            end
            % section [OPTIONS]
            if (sect==1)
                fprintf(fid2, ['AREA_UNITS', blanks(5), '%s\n'], areaunits);
                fprintf(fid2, ['RATE_UNITS', blanks(5), '%s\n'], rateunits);
                fprintf(fid2, ['SOLVER', blanks(5), '%s\n'], solver);
                fprintf(fid2, ['COUPLING', blanks(5), '%s\n'], coupling);
                fprintf(fid2, ['COMPILER', blanks(5), '%s\n'], compiler);
                fprintf(fid2, ['TIMESTEP', blanks(5), '%d\n'], timestep);
                fprintf(fid2, ['RTOL', blanks(5), '%s\n'], rtol);
                fprintf(fid2, ['ATOL', blanks(5), '%s\n\n'], atol);
                sect=2;
                break;
            elseif sect==0
                fprintf(fid2, tline{t});
                break;
            end
            u=u+1;
        end
    end
    if sect~=2
        fprintf(fid2, '\n');
    end
end
obj.unloadMSX;
obj.loadMSXFile(obj.MSXTempFile, 1);
end
function [indices, value] = getLinkIndices(obj, varargin)
    indices = getIndices(obj.getLinkCount, varargin{1});
    value = zeros(1, length(indices));
end
function indices = getControlIndices(obj, varargin)
    indices = getIndices(obj.getControlRulesCount, varargin{1});
end
function [indices, value] = getNodeIndices(obj, varargin)
    indices = getIndices(obj.getNodeCount, varargin{1});
    value = zeros(1, length(indices));
end
function [indices, value] = getNodeJunctionIndices(obj, varargin)
    % EPANET Version 2.2
    indices = getIndices(obj.getNodeJunctionCount, varargin{1});
    value = zeros(1, length(indices));
end
function [indices, value] = getCurveIndices(obj, varargin)
    indices = getIndices(obj.getCurveCount, varargin{1});
    value = zeros(1, length(indices));
end
function [indices, value] = getPatternIndices(obj, varargin)
    indices = getIndices(obj.getPatternCount, varargin{1});
    value = zeros(1, length(indices));
end
function indices = getIndices(cnt, varargin)
    if isempty(varargin{1})
        indices=1:cnt;
    elseif isempty(varargin{1}{1})
        indices=1:cnt;
    else
        indices=varargin{1}{1};
    end
end
function value = getTimes(obj, r, atline, value)
    tmpt=[0 0];
    if sum(r==':')
        r=regexp(r, ':', 'split');
        tmpt(1)=str2double(r{1});
        tmpt(2)=str2double(r{2});
        secnd=tmpt(1)*3600+tmpt(2)*60;
        if length(r)>2
            secnd=secnd+str2double(r{3});
        end
    else
        tmp=3600;
        switch upper(atline{end})
            case 'MIN'
                tmp=60;
            case 'MINUTES'
                tmp=60;
            case 'SECONDS'
                tmp=1/60;
            case 'SEC'
                tmp=1/60;
            case 'DAYS'
                tmp=24;
        end
        secnd=single(str2double(r)*tmp);
    end
    switch upper(atline{1})
        case 'DURATION'
            value.BinTimeSimulationDuration=secnd;
        case 'HYDRAULIC'
            value.BinTimeHydraulicStep=secnd;
        case 'QUALITY'
            value.BinTimeQualityStep=secnd;
        case 'PATTERN'
            if strcmpi(atline{2}, 'TIMESTEP')
                value.BinTimePatternStep=secnd;
            elseif strcmpi(atline{2}, 'START')
                value.BinTimePatternStart=secnd;
            end
        case 'REPORT'
            if strcmpi(atline{2}, 'TIMESTEP')
                value.BinTimeReportingStep=secnd;
            elseif strcmpi(atline{2}, 'START')
                value.BinTimeReportingStart=secnd;
            end
        case 'STATISTIC'
            value.BinTimeStatisticsIndex=find((strcmpi(obj.TYPESTATS, atline{2})-1)>-1)-1;
            value.BinTimeStatistics=obj.TYPESTATS{value.BinTimeStatisticsIndex+1};
    end
end
function value = getOptionsValues(obj, value, atline)
    switch upper(atline{1})
        case 'UNITS'
            value.BinLinkFlowUnits=atline{2};
        case 'HEADLOSS'
            value.BinOptionsHeadloss=atline{2};
        case 'PRESSURE'
            value.BinNodePressureUnits=atline{2};
        case 'SPECIFIC'
            value.BinOptionsSpecificGravity=str2double(atline{3});
        case 'VISCOSITY'
            value.BinOptionsViscosity=str2double(atline{2});
        case 'TRIALS'
            value.BinOptionsMaxTrials=str2double(atline{2});
        case 'ACCURACY'
            value.BinOptionsAccuracyValue=single(str2double(atline{2}));
        case 'UNBALANCED'
            value.BinOptionsUnbalanced= atline(2:end);
        case 'PATTERN'
            value.BinOptionsPattern=str2double(atline{2});
        case 'DEMAND'
            value.BinOptionsPatternDemandMultiplier=str2double(atline{3});
        case 'EMITTER'
            value.BinOptionsEmitterExponent=str2double(atline{3});
        case 'QUALITY'
            value.BinQualityType=atline{2};% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
            value.BinQualityCode=find((strcmpi(obj.TYPEQUALITY, atline{2})-1)>-1)-1;
            if isempty(value.BinQualityCode)
                value.BinQualityCode=1;
            end
            if length(atline)>2
                if value.BinQualityCode==3
                    value.BinQualityTraceNodeIndex=obj.getBinNodeIndex(atline{3});
                    value.BinQualityTraceNodeID=atline{3};
                else
                    value.BinQualityUnits=atline{3};
                end
            end
        case 'DIFFUSIVITY'
            value.BinOptionsDiffusivity=str2double(atline{2});
        case 'TOLERANCE'
            value.BinOptionsQualityTolerance=single(str2double(atline{2}));
    end
end
function [value, cont, sect, i, t, q, d] = getLV(tok, value, sect, tline, i, t, q, d)
    cont=1;
    if (tok(1) == '[')
           % [PIPES] section
        if strcmpi(tok(1:5), '[PIPE')
            sect=1;
            value.BinLinkPipeNameID={};
            value.BinLinkPipeIndex=[];
            value.BinLinkFromNode={};
            value.BinLinkToNode={};
            value.BinLinkPipeLengths=[];
            value.BinLinkPipeDiameters=[];
            value.BinLinkPipeRoughness=[];
            value.BinLinkPipeMinorLoss=[];
            cont=1;return;
            % [PUMPS] section
        elseif strcmpi(tok(1:5), '[PUMP')
            sect=2;
            value.BinLinkPumpNameID={};
            value.BinLinkPumpIndex=[];
            value.BinLinkPumpPatterns={};
            value.BinLinkPumpCurveNameID={};
            value.BinLinkPumpPower=[];
            value.BinLinkPumpNameIDPower={};
            cont=1;return;
            % [VALVES] section
        elseif strcmpi(tok(1:5), '[VALV')
            sect=3;
            value.BinLinkValveNameID={};
            value.BinLinkValveIndex=[];
            value.BinLinkValveDiameters=[];
            value.BinLinkValveType={};
            value.BinLinkValveSetting=[];
            value.BinLinkValveMinorLoss=[];
            cont=1;return;
            % [STATUS] section
        elseif strcmpi(tok(1:5), '[STAT')
            sect=4;
            value.BinLinkInitialStatus={};
            value.BinLinkInitialStatusNameID={};
            value.BinCountStatuslines=0;
            cont=1;return;
        elseif strcmpi(tok(1:5), '[REAC')
            sect=5;
            value.BinLinkGlobalBulkReactionCoeff=[];
            value.BinLinkGlobalWallReactionCoeff=[];
            value.BinLinkBulkReactionCoeff=[];
            value.BinLinkWallReactionCoeff=[];
            value.BinCountReactionlines=0;
            value.BinLinkPipeCount = length(value.BinLinkPipeNameID);
            value.BinLinkPumpCount = length(value.BinLinkPumpNameID);
            value.BinLinkValveCount = length(value.BinLinkValveNameID);
            value.BinLinkCount = value.BinLinkPipeCount+value.BinLinkPumpCount+value.BinLinkValveCount;
            value.BinLinkNameID = [value.BinLinkPipeNameID value.BinLinkPumpNameID value.BinLinkValveNameID];
            cont=1;return;
            % [END]
        elseif strcmpi(tok(1:4), '[END')
            cont=0;return;
        else
            sect = 0;
            cont=1;return;
        end
    end
    clear atline;
    atline = checktlines(tline);
    if sect==0
        return;
        % Links
    elseif sect==1
        value.BinLinkPipeNameID{t}=atline{1};
        value.BinLinkPipeIndex(t)=t;
        value.BinLinkFromNode{t}=atline{2};
        value.BinLinkToNode{t}=atline{3};
        value.BinLinkPipeLengths(t)=str2double(atline{4});
        value.BinLinkPipeDiameters(t)=str2double(atline{5});
        value.BinLinkPipeRoughness(t)=str2double(atline{6});
        if length(atline)>6
            value.BinLinkPipeMinorLoss(t)=str2double(atline{7});
        end
        value.BinLinkType={};
        if length(atline)>7
            value.BinLinkPipeStatus{t}=atline{8};
        else
            value.BinLinkPipeStatus{t}='Open';
        end
        t=t+1;
    elseif sect==2
        value.BinLinkPumpNameID{q}=atline{1};
        value.BinLinkPumpIndex(q)=t;
        value.BinLinkFromNode{t}=atline{2};
        value.BinLinkToNode{t}=atline{3};
        if strcmp(regexp(tline, '\w*HEAD*\w', 'match'), 'HEAD')
            value.BinLinkPumpCurveNameID{q}=atline{5};
        elseif strcmp(regexp(tline, '\w*POWER*\w', 'match'), 'POWER')
            value.BinLinkPumpPower(q)=str2double(atline{5});
            value.BinLinkPumpNameIDPower{q}=atline{1};
        end
        if length(atline)>6
            value.BinLinkPumpPatterns{q}=atline{7};
        end
        t=t+1;
        q=q+1;
    elseif sect==3
        value.BinLinkValveNameID{i}=atline{1};
        value.BinLinkValveIndex(i)=t;
        value.BinLinkFromNode{t}=atline{2};
        value.BinLinkToNode{t}=atline{3};
        value.BinLinkValveDiameters(i)=str2double(atline{4});
        value.BinLinkValveType{i}=atline{5};
        value.BinLinkValveSetting(i)=str2double(atline{6});
        if length(atline)>6
            value.BinLinkValveMinorLoss(i)=str2double(atline{7});
        end
        t=t+1;
        i=i+1;
        % Status
    elseif sect==4
        value.BinLinkInitialStatus{d}=atline{2};
        value.BinLinkInitialStatusNameID{d}=atline{1};
        value.BinCountStatuslines=d;
        d=d+1;
        % Reactions
    elseif sect==5
        if strcmpi(atline{1}, 'GLOBAL') && strcmpi(atline{2}, 'BULK')
            value.BinLinkGlobalBulkReactionCoeff=str2double(atline{3});
        elseif strcmpi(atline{1}, 'GLOBAL') && strcmpi(atline{2}, 'WALL')
            value.BinLinkGlobalWallReactionCoeff=str2double(atline{3});
            value.BinLinkWallReactionCoeff=value.BinLinkGlobalWallReactionCoeff*ones(1, value.BinLinkCount);
            value.BinLinkBulkReactionCoeff=value.BinLinkGlobalBulkReactionCoeff*ones(1, value.BinLinkCount);
            value.BinLinkWallReactionCoeff(value.BinLinkPumpIndex)=0;
            value.BinLinkWallReactionCoeff(value.BinLinkValveIndex)=0;
            value.BinLinkBulkReactionCoeff(value.BinLinkPumpIndex)=0;
            value.BinLinkBulkReactionCoeff(value.BinLinkValveIndex)=0;
        end
        if strcmpi(atline{1}, 'BULK')
            LinkIndex = find(strcmpi(value.BinLinkNameID, atline{2}));
            value.BinLinkBulkReactionCoeff(LinkIndex)=str2double(atline{3});
        elseif strcmpi(atline{1}, 'WALL')
            LinkIndex = find(strcmpi(value.BinLinkNameID, atline{2}));
            value.BinLinkWallReactionCoeff(LinkIndex)=str2double(atline{3});
        end
        value.countReactionlines=d;
        d=d+1;
    end
end
function [status, result] = runMSX(obj, rptfile, varargin)
    inpfile=obj.BinTempfile;
    arch=computer('arch');
    if strcmp(arch,'win64') || strcmp(arch,'win32')
        [~,lpwd]=system(['cmd /c for %A in ("',obj.MSXLibEPANETPath,'") do @echo %~sA']);
        libtmp=regexp(lpwd,'\s','split');
        libPwd=libtmp{1};
        if nargin<3
            r = sprintf('%s\\epanetmsx.exe "%s" "%s" "%s"',libPwd,inpfile,obj.MSXTempFile,rptfile);
        else
            binfile=varargin{1};
            r = sprintf('%s\\epanetmsx.exe "%s" "%s" "%s" "%s"',libPwd,inpfile,obj.MSXTempFile,rptfile,binfile);
        end
    elseif isunix
        libPwd=obj.MSXLibEPANETPath;
        if nargin<3
            r = sprintf('%s\\epanetmsx "%s" "%s" "%s"',libPwd,inpfile,obj.MSXTempFile,rptfile);
        else
            binfile=varargin{1};
            r = sprintf('%s\\epanetmsx "%s" "%s" "%s" "%s"',libPwd,inpfile,obj.MSXTempFile,rptfile,binfile);
        end
    end
    [status, result] = system(r);
end
% function value = readMSXBinaryFile(binfile)
%     fid = fopen(binfile, 'r');
%     if fid~=-1
%         data = fread(fid, 'int32');
%         fclose(fid);
%         value.BinMSXNumberReportingPeriods = data(end-2);
%         clear data;
%         fid1 = fopen(binfile, 'r');
%
%         % Seek to the 10th byte ('J'), read 5
%         fseek(fid1, 0, 'bof');
%         value.BinMSXmagicnumber=fread(fid1, 1, 'uint32');
%         value.BinMSXLibMSX=fread(fid1, 1, 'uint32');
%         value.BinMSXNumberNodes=fread(fid1, 1, 'uint32');
%         value.BinMSXNumberLinks=fread(fid1, 1, 'uint32');
%         value.BinMSXSpeciesCount=fread(fid1, 1, 'uint32');
%         value.BinMSXReportingTimeStepSec=fread(fid1, 1, 'uint32');
%
%         for i=1:value.BinMSXSpeciesCount
%             value.BinMSXSpeciesNumberChar(i)=fread(fid1, 1, 'uint32');
%             value.BinMSXSpeciesNameID{i}=fread(fid1, value.BinMSXSpeciesNumberChar(i), '*char')';
%         end
%         for i=1:value.BinMSXSpeciesCount
%             value.BinMSXSpeciesUnits{i}=fread(fid1, 15, '*char')';
%         end
%         value.BinMSXSpeciesUnits = regexprep(value.BinMSXSpeciesUnits, '[^\w'']', '');
%
%         fread(fid1, 32, 'float');
%         for i=1:value.BinMSXReportingTimeStepSec/3600:value.BinMSXNumberReportingPeriods
%             for s=1:value.BinMSXSpeciesCount
%                 for u=1:value.BinMSXNumberNodes
%                     value.BinMSXNodesQuality{s}(i, u) = fread(fid1, 1, 'float')'; %%%% edit here
%                 end
%             end
%         end
%         for i=1:value.BinMSXReportingTimeStepSec/3600:value.BinMSXNumberReportingPeriods
%             for s=1:value.BinMSXSpeciesCount
%                 for u=1:value.BinMSXNumberLinks
%                     value.BinMSXLinksQuality{s}(i, :) = fread(fid1, 1, 'float')';
%                 end
%             end
%         end
%     end
% end
function value = readEpanetBin(fid, binfile, varargin)
    value=[];
    if fid~=-1
        data = fread(fid, 'int32');
        value.BinNumberReportingPeriods = data(end-2);
        clear data;
        % Beginning of file
        fseek(fid, 0, 'bof');
        value.Binmagicnumber=fread(fid, 1, 'uint32');
        value.BinLibEPANET=fread(fid, 1, 'uint32');
        value.BinNumberNodes=fread(fid, 1, 'uint32');
        value.BinNumberReservoirsTanks=fread(fid, 1, 'uint32');
        value.BinNumberLinks=fread(fid, 1, 'uint32');
        value.BinNumberPumps=fread(fid, 1, 'uint32');
        value.BinNumberValves=fread(fid, 1, 'uint32');
        value.BinWaterQualityOption=fread(fid, 1, 'uint32');
        value.BinIndexNodeSourceTracing=fread(fid, 1, 'uint32');
        value.BinFlowUnitsOption=fread(fid, 1, 'uint32');
        value.BinPressureUnitsOption=fread(fid, 1, 'uint32');
        value.BinTimeStatisticsFlag=fread(fid, 1, 'uint32');
        value.BinReportingStartTimeSec=fread(fid, 1, 'uint32');
        value.BinReportingTimeStepSec=fread(fid, 1, 'uint32');
        value.BinSimulationDurationSec=fread(fid, 1, 'uint32');
        value.BinProblemTitle1=fread(fid, 80, '*char')';
        value.BinProblemTitle2=fread(fid, 80, '*char')';
        value.BinProblemTitle3=fread(fid, 80, '*char')';
        value.BinNameInputFile=fread(fid, 260, '*char')';
        value.BinNameReportFile=fread(fid, 260, '*char')';
        value.BinNameChemical=regexprep(fread(fid, 16, '*char')', '[^\w'']', '');
        value.BinChemicalConcentrationUnits=regexprep(fread(fid, 32, '*char')', '[^\w'']', '');
        fread(fid, 4, 'uint32');
        for i=1:value.BinNumberNodes
            value.BinIDLabelEachNode{i}=fread(fid, 32, '*char')'; % error NODES*32
            value.BinIDLabelEachNode{i}=value.BinIDLabelEachNode{i}(find(value.BinIDLabelEachNode{i}, 32));
        end
        for i=1:value.BinNumberLinks
            value.BinIDLabelEachLink{i}=fread(fid, 32, '*char')';  % error LINKS*32
            value.BinIDLabelEachLink{i}=value.BinIDLabelEachLink{i}(find(value.BinIDLabelEachLink{i}, 32));
        end
        value.BinIndexStartNodeEachLink=fread(fid, value.BinNumberLinks, 'uint32')';
        value.BinIndexEndNodeEachLink=fread(fid, value.BinNumberLinks, 'uint32')';
        value.BinTypeCodeEachLink=fread(fid, value.BinNumberLinks, 'uint32')';
        value.BinNodeIndexEachReservoirsTank=fread(fid, value.BinNumberReservoirsTanks, 'uint32')'; % error
        value.BinCrossSectionalAreaEachTank=fread(fid, value.BinNumberReservoirsTanks, 'float')';
        value.BinElevationEachNode=fread(fid, value.BinNumberNodes, 'float')';
        value.BinLengthEachLink=fread(fid, value.BinNumberLinks, 'float')';
        value.BinDiameterEachLink=fread(fid, value.BinNumberLinks, 'float')';

        for p=1:value.BinNumberPumps
            value.BinPumpIndexListLinks(p)=fread(fid, 1, 'float')';
            value.BinPumpUtilization(p)=fread(fid, 1, 'float')';
            value.BinAverageEfficiency(p)=fread(fid, 1, 'float')';
            value.BinAverageKwattsOrMillionGallons(p)=fread(fid, 1, 'float')';
            value.BinAverageKwatts(p)=fread(fid, 1, 'float')';
            value.BinPeakKwatts(p)=fread(fid, 1, 'float')';
            value.BinAverageCostPerDay(p)=fread(fid, 1, 'float')';
        end

        fread(fid, 1, 'float');

        for i=1:value.BinNumberReportingPeriods
            value.BinNodeDemand(i,:)         = fread(fid, value.BinNumberNodes, 'float')';
            value.BinNodeHead(i,:)           = fread(fid, value.BinNumberNodes, 'float')';
            value.BinNodePressure(i,:)       = fread(fid, value.BinNumberNodes, 'float')';
            value.BinNodeQuality(i,:)        = fread(fid, value.BinNumberNodes, 'float')';
            value.BinLinkFlow(i,:)           = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkVelocity(i,:)       = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkHeadloss(i,:)       = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkQuality(i,:)        = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkStatus(i,:)         = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkSetting(i,:)        = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkReactionRate(i,:)   = fread(fid, value.BinNumberLinks, 'float')';
            value.BinLinkFrictionFactor(i,:) = fread(fid, value.BinNumberLinks, 'float')';
        end

        value.BinAverageBulkReactionRate=fread(fid, 1, 'float')';
        value.BinAverageWallReactionRate=fread(fid, 1, 'float')';
        value.BinAverageTankReactionRate=fread(fid, 1, 'float')';
        value.BinAverageSourceInflowRate=fread(fid, 1, 'float')';
        value.BinNumberReportingPeriods2=fread(fid, 1, 'uint32')';
        value.BinWarningFlag=fread(fid, 1, 'uint32')';
        value.BinMagicNumber=fread(fid, 1, 'uint32')';
    end
    if ~isempty(varargin)
        v = struct();
        v.Time = (0:value.BinReportingTimeStepSec:value.BinSimulationDurationSec)';

        fields_param = {'BinNodePressure', 'BinNodeDemand', 'BinNodeHead', 'BinNodeQuality',...
                        'BinLinkFlow', 'BinLinkVelocity', 'BinLinkHeadloss', 'BinLinkStatus', 'BinLinkSetting', ...
                        'BinLinkReactionRate', 'BinLinkFrictionFactor', 'BinLinkQuality'};
        fields_new = {'Pressure', 'Demand', 'Head', 'NodeQuality',...
                    'Flow', 'Velocity', 'HeadLoss', 'Status', 'Setting', ...
                    'ReactionRate', 'FrictionFactor', 'LinkQuality'};
        for i=1:length(fields_param)
            v.(fields_new{i}) = eval(['value.', fields_param{i}]);
        end
        clear value;
        value=v;
    end
    warning('off'); try fclose(fid); catch, end; try delete(binfile); catch, end
    warning('on');
end
function [inpfile, rptfile, binfile]= createTempfiles(BinTempfile)
    inpfile=BinTempfile;
    uuID = char(java.util.UUID.randomUUID);
    rptfile=['@#', uuID, '.txt'];
    binfile=['@#', uuID, '.bin'];
end
function typecode = getTypeLink(type)
    ch_ = find(type=='_');
    if ~isempty(ch_)
        type=type(ch_+1:end);
    end
    switch upper(type)
        case 'CVPIPE'
            typecode=0;
        case 'PIPE'
            typecode=1;
        case 'PUMP'
            typecode=2;
        case 'PRV'
            typecode=3;
        case 'PSV'
            typecode=4;
        case 'PBV'
            typecode=5;
        case 'FCV'
            typecode=6;
        case 'TCV'
            typecode=7;
        case 'GPV'
            typecode=8;
        otherwise
            typecode=-1;
    end
end
function fid = writenewTemp(Tempfile)
    fid=fopen(Tempfile, 'w');
    while fid==-1, fid=fopen(Tempfile, 'w'); end
end
function value = BinNodeCoords(obj, vertices)
    BinNodeName = obj.getBinNodeNameID;
    BinLinkName = obj.getBinLinkNameID;
    linkcount =  length(BinLinkName.BinLinkNameID);
    nodecount = length(BinNodeName.BinNodeNameID);

    vx = NaN(nodecount, 1);
    vy = NaN(nodecount, 1);
    vertx = cell(linkcount, 1);
    verty = cell(linkcount, 1);
    nvert = zeros(linkcount, 1);
    % Open epanet input file
    [~, info] = obj.readInpFile;

    sect = 0;
    IndexC = strfind(info,'[COORDINATES]');
    Index = find(not(cellfun('isempty',IndexC)));
    for h=Index:length(info)
        tline = info{h};
        if ~ischar(tline)
            break;
        end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank Clines and comments
        if isempty(tok)
            continue;
        end
        if (tok(1) == ';')
            continue;
        end
        if (tok(1) == '[')

            % [COORDINATES] section
            if strcmpi(tok(1:5), '[COOR')
                sect = 17;
                continue;
            % [VERTICES] section
            elseif strcmpi(tok(1:5), '[VERT')
                sect = 18;
                continue;
            % [END]
            elseif strcmpi(tok(1:4), '[END')
                break;
            else
                sect = 0;
                continue;
            end
        end
        if sect==0
            continue;

        % Coordinates
        elseif sect==17
            if ~vertices
                A = textscan(tline, '%s %f %f');
                mm = strcmp(A{1}, BinNodeName.BinNodeNameID);
                index=strfind(mm, 1);
                if isempty(index)
                    continue;
                end
                vx(index) = A{2};
                vy(index) = A{3};
            end

        % Vertices
        elseif sect==18
            A = textscan(tline, '%s %f %f');
            index =  find(strcmp(BinLinkName.BinLinkNameID, A{1}));
            if isempty(index)
                continue;
            end
            nvert(index) = nvert(index) + 1;
            vertx{index}(nvert(index)) = A{2};
            verty{index}(nvert(index)) = A{3};
        end
    end
    value{1} = vx;
    value{2} = vy;
    value{3} = vertx;
    value{4} = verty;
end
function atline = checktlines(tline)
    atline='';
    a = regexp(tline, '\s*', 'split');uu=1;
    for tt=1:length(a)
        if isempty(a{tt})
            %skip
        elseif sum(a{tt}==';')
            %skip
            if tt>1,  break; end
        else
            atline{uu}=a{tt}; uu=uu+1;
        end
    end
end
function setControlFunction(obj, index, value)
    controlRuleIndex = index;
    [controlTypeIndex, linkIndex,controlSettingValue,...
    nodeIndex, controlLevel] = controlSettings(obj, value);
    [obj.Errcode] = ENsetcontrol(controlRuleIndex, ...
        controlTypeIndex, linkIndex, controlSettingValue, nodeIndex, controlLevel, obj.LibEPANET);
    error(obj.getError(obj.Errcode));
end
function controlRuleIndex = addControlFunction(obj, value)
    if isstruct(value)
        for c=1:length(value)
            [controlTypeIndex, linkIndex, controlSettingValue,...
            nodeIndex, controlLevel] = controlSettings(obj, value(c).Control);
            [obj.Errcode, controlRuleIndex(c)] = ENaddcontrol(controlTypeIndex, linkIndex,...
            controlSettingValue, nodeIndex, controlLevel, obj.LibEPANET);
            error(obj.getError(obj.Errcode));
        end
    else
        [controlTypeIndex, linkIndex, controlSettingValue,...
            nodeIndex, controlLevel] = controlSettings(obj, value);
        [obj.Errcode, controlRuleIndex] = ENaddcontrol(controlTypeIndex, linkIndex,...
            controlSettingValue, nodeIndex, controlLevel, obj.LibEPANET);
        error(obj.getError(obj.Errcode));
    end
end
function [controlTypeIndex, linkIndex,controlSettingValue,...
    nodeIndex, controlLevel] = controlSettings(obj, value)
    splitControl = strsplit(value);
    controlSettingValue = find(strcmpi(obj.TYPESTATUS, splitControl(3)))-1;
    if isempty(controlSettingValue)
        if strcmpi(splitControl(3), 'CLOSE')
            controlSettingValue = 0;
        else
            % control setting Value (type should be int) for pump or valve
            controlSettingValue = str2double(splitControl{3});
        end
    end
    linkIndex = obj.getLinkIndex(splitControl(2));
    if ~linkIndex
        warning('Wrong link ID. Please change your control.')
    end
    switch upper(splitControl{4})
        case 'IF'
            %LINK linkID status IF NODE nodeID ABOVE/BELOW value
            nodeIndex = obj.getNodeIndex(splitControl(6));
            controlTypeIndex = 0; % LOWLEVEL
            if strcmpi(splitControl(7), 'ABOVE')
                controlTypeIndex = 1; % HIGHLEVEL
            end
            controlLevel = str2double(splitControl{8});
        case 'AT'
            if strcmpi(splitControl{5}, 'CLOCKTIME')
                %LINK linkID status AT CLOCKTIME clocktime AM/PM
                nodeIndex = 0;
                controlTypeIndex = 3;
            else
                %LINK linkID status AT TIME time
                nodeIndex = 0;
                controlTypeIndex = 2;
            end
            if isempty(strfind(splitControl{6}, ':'))
                controlLevel = str2double(splitControl{6});
            else
                [~, ~, days, H, MN, S] = datevec(splitControl{6});
                controlLevel = (24*(days-1)+H)*3600+MN*60+S;
            end
        otherwise
    end
end
