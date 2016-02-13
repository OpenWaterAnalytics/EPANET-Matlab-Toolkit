classdef epanet <handle
    %epanet EPANET-Matlab Class: A Matlab Class for EPANET and EPANET-MSX
    %
    %   How to run: d=epanet('Net1_Rossman2000.inp');
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
    %   The latest EPANET files can downloaded by the *unofficial* channel
    %   https://github.com/OpenWaterAnalytics/epanet which is the most
    %   updated libepanet of the software
    %
    %   Some part of the EPANET-Matlab Class are based on the epanet-matlab
    %   wrappers prepared by Jim Uber
    %   (https://github.com/OpenWaterAnalytics/epanet-matlab)
    %
    %   EPANET-Matlab Class Licence:
    %
    %   Copyright 2013 KIOS Research Center for Intelligent Systems and
    %   Networks, University of Cyprus (www.kios.org.cy)
    %
    %   Licensed under the EUPL, libepanet 1.1 or - as soon they will be
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
        Msxlibepanet;
        MsxlibepanetPath;
        MsxConstantsNameID;
        MsxConstantsValue;
        MsxConstantsCount;
        MsxConstantsIndex;
        MsxParametersCount;
        MsxPatternsCount;
        MsxSpeciesCount;
        MsxLinkInitqualValue;
        MsxNodeInitqualValue;
        MsxFile;
        MsxTempFile;
        MsxParametersNameID;
        MsxParametersIndex;
        MsxParametersPipesValue;
        MsxParametersTanksValue;
        MsxPatternsNameID;
        MsxPatternsIndex;
        MsxPatternsLengths;
        MsxPattern;
        MsxEquationsPipes;
        MsxSources;
        MsxSourceLevel;
        MsxSourceNodeNameID;
        MsxSourcePatternNameID;
        MsxSourcePatternIndex;
        MsxSourceType;
        MsxSourceTypeCode;
        MsxSpeciesATOL;
        MsxSpeciesIndex;
        MsxSpeciesNameID;
        MsxSpeciesRTOL;
        MsxSpeciesType;
        MsxSpeciesUnits;
        MsxEquationsTanks;
        MsxEquationsTerms;
        NodeCoordinates; % Coordinates for each node (long/lat & intermediate pipe coordinates)
        NodeJunctionCount; %Number of junctions
        NodeCount; % Number of nodes
        NodeReservoirCount; %Number of reservoirs
        NodeTankCount;% Number of tanks
        NodeTankReservoirCount;   %Number of tanks and reservoirs
        NodeTankDiameterUnits; % Units for tank diameters
        NodeElevationUnits; % Units for elevation
        NodeEmitterCoefficientUnits; %Units for emitter coefficient
        NodeBaseDemands; % Base demands of nodes
        NodeDemandPatternIndex; % Index of demand patterns
        NodeElevations;% Elevation of nodes
        NodeEmitterCoeff; %Emmitter Coefficient of nodes
        NodeIndex; %Index of nodes
        NodeInitialQuality; %Initial quality of nodes
        NodeJunctionIndex; %Index of node junctions
        NodeJunctionNameID; %Name ID of node junctions
        NodeNameID; %Name ID of all nodes
        NodeReservoirIndex; %Index of reservoirs
        NodeReservoirNameID;%Name ID of reservoirs
        NodesConnectingLinksID; %Name IDs of nodes which connect links
        NodesConnectingLinksIndex; %Indices of nodes which connect links
        NodeSourcePatternIndex; %Index of pattern for node sources
        NodeSourceQuality; %Quality of node sources
        NodeSourceTypeIndex; %Index of source type
        NodeType; %ID of node type
        NodeTypeIndex; %Index of nodetype
        NodeHeadUnits; % Nodal head units
        NodeTankBulkReactionCoeff; %Bulk reaction coefficients in tanks
        NodeTankDiameter; %Tank Diameters
        NodeTankIndex; %Indices of Tanks
        NodeTankInitialLevel; %Initial water level in tanks
        NodeTankInitialWaterVolume; %Initial water volume in tanks
        NodeTankMaximumWaterLevel; % Maximum water level in tanks
        NodeTankMinimumFraction; % Fraction of the total tank volume devoted to the inlet/outlet compartment
        NodeTankMinimumWaterLevel; % Minimum water level
        NodeTankMinimumWaterVolume; % Minimum water volume
        NodeTankMixingModelCode; % Code of mixing model (MIXED:0, 2COMP:1, FIFO:2, LIFO:3)
        NodeTankMixingModelType; % Type of mixing model (MIXED, 2COMP, FIFO, or LIFO)
        NodeTankMixZoneVolume;  % Mixing zone volume
        NodeTankNameID; % Name ID of Tanks
        NodeTankVolumeCurveIndex; %Index of curve for tank volumes
        NodeTankVolumeUnits; %Units for volume
        NodePressureUnits; % PUnits for Pressure
        LinkCount; % Number of links
        LinkPipeCount;% Number of pipes
        LinkPumpCount;%Number of pumps
        LinkValveCount;% Number of valves
        LinkPipeDiameterUnits; % Units for pipe diameters
        LinkFrictionFactorUnits; %Units for friction factor
        LinkLengthUnits; % Units of length
        LinkBulkReactionCoeff; % Bulk reaction coefficient of each link
        LinkDiameter; % Diameter of each link
        LinkIndex; % Index of links
        LinkInitialSetting; %Initial settings of links
        LinkInitialStatus; %Initial status of links
        LinkLength; % Length of links
        LinkLengthsUnits; % Units of length
        LinkMinorLossCoeff; %Minor loss coefficient of links
        LinkNameID; % Name ID of links
        LinkPipeIndex; % Index of pipe links
        LinkPipeNameID; % Name ID of pipe links
        LinkPumpIndex; % Index of pumps
        LinkPumpNameID; % Name ID of pumps
        LinkRoughnessCoeff; % Roughness coefficient of links
        LinkType; %ID of link type
        LinkTypeIndex; %Index of link type
        LinkWallReactionCoeff; %Wall reaction coefficient of links
        LinkMinorLossCoeffUnits; %Minor loss coefficient units
        LinkPumpPowerUnits; %Units of power
        LinkPipeRoughnessCoeffUnits; %Pipe roughness coefficient units
        LinkFlowUnits; %Units of flow
        LinkValveIndex; % Index of valves
        LinkValveNameID; % ID name of valves
        LinkStatus;
        LinkVelocityUnits; % Units for velocity
        ControlLevelValues; %The control level values
        ControlLinkIndex; % Set of control types in links
        ControlNodeIndex; % Set of control types in nodes
        ControlSettings; % Settings for the controls
        ControlTypes; % Set of control types
        ControlTypesIndex; %Index of the control types
        ControlRules; %Retrieves the parameters of all control statements
        ControlRulesCount;  % Number of controls
        Controls;
        CurveCount; % Number of curves
        TimeRuleControlStep; % Time step for evaluating rule-based controls
        PatternCount;  % Number of patterns
        PatternDemandsUnits; %Units for demands
        PatternNameID; %ID of the patterns
        PatternIndex; %Indices of the patterns
        PatternLengths; %Length of the patterns
        Pattern; % get all patterns
        TimePatternStart; %Pattern start time
        TimePatternStep; %Pattern Step
        EnergyEfficiencyUnits; % Units for efficiency
        EnergyUnits; %Units for energy
        OptionsHeadloss; %Headloss formula (Hazen-Williams, Darcy-Weisbach or Chezy-Manning)
        OptionsHydraulics; %Save or Use hydraulic soltion. *** Not implemented ***
        OptionsViscosity; %*** Not implemented ***
        OptionsSpecificGravity; %*** Not implemented ***
        OptionsMaxTrials; % Maximum number of trials (40 is default)
        OptionsAccuracyValue; %Convergence value (0.001 is default)
        OptionsPattern; % *** Not implemented ***
        OptionsPatternDemandMultiplier; %Multiply demand values (1 is default)
        OptionsEmitterExponent; %Exponent of pressure at an emmiter node (0.5 is default)
        OptionsQualityTolerance; %Tolerance for water quality (0.01 is default)
        OptionsUnbalanced; %*** Not implemented ***
        OptionsUnbalancedContinueN; % *** Not implemented ***
        QualityCode; % Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        QualityTraceNodeIndex; %Index of trace node (0 if QualityCode<3)
        QualityType; % Water quality analysis type (None/Chemical/Age/Trace)
        QualityReactionCoeffBulkUnits; % Bulk reaction coefficient units
        QualityReactionCoeffWallUnits; % Wall reaction coefficient units
        QualitySourceMassInjectionUnits; % Units for source mass injection
        QualityWaterAgeUnits; %Units for water age
        QualityUnits; %Units for quality concentration.
        TimeHydraulicStep; %Hydraulic time step
        TimeQualityStep; %Quality Step
        TimeReportingPeriods; % Reporting periods
        TimeReportingStart; %Start time for reporting
        TimeReportingStep; %Reporting time step
        TimeSimulationDuration; %Simulation duration
        TimeStatisticsIndex; %Index of time series post-processing type ('NONE':0,'AVERAGE':1,'MINIMUM':2,'MAXIMUM':3, 'RANGE':4)
        TimeStatisticsType; %Type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
        
        %%%%% New version dev2.1 %%%%%
        TimeStartTime;
        TimeHTime;
        TimeHaltFlag;
        TimeNextEvent; %find the lesser of the hydraulic time step length, or the time to next fill/empty
        NodeTankMaxVolume;
        CurvesInfo;
        HeadCurveIndex;
        LinkPumpPatternNameID;
        LinkPumpPatternIndex;
        NodeNumDemandCategories;
        NodeDemandPatternNameID;
        NodeDemandPatternsIndex;
        LinkPumpTypeCode;
        LinkPumpType;
        RelativeError;
        Iterations;
        PatternAveragePatternValue;
        QualityChemName;
        QualityChemUnits;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        errcode; %Code for the EPANET error message
        pathfile;   % The path of the input file
        inputfile;  % Name of the input file
        libepanet; %EPANET library dll
        libepanetpath; %EPANET library dll path
        Version; % EPANET version
        
        SImetric;
        UScustomary;
        
        % Bin 
        Bin;
        Bintempfile;
        BinNodeJunctionNameID;
        BinNodeReservoirNameID;
        BinNodeTankNameID;
        BinNodeNameID;
        BinNodeType;
        BinNodeReservoirIndex;
        BinNodeTankIndex;
        BinNodeJunctionIndex;
        BinNodeElevations;
        BinNodeTankElevation;
        BinNodeTankInitLevel;
        BinNodeTankMinLevel;
        BinNodeTankMaxLevel;
        BinNodeTankDiameter;
        BinNodeTankMinVol;
        BinNodeTankMixID;
        BinNodeTankMixModel;
        BinNodeTankMinimumFraction;
        BinNodeReservoirElevation;
        BinNodeJunctionElevation;
        BinNodeResDemandPatternNameID;
        BinNodeCount;  
        BinNodeJunctionCount; 
        BinNodeReservoirCount; 
        BinNodeTankCount;
        BinNodeTankReservoirCount;   
        BinNodeSourceTypeCode;
        BinNodeSourceType;
        BinNodeSourcePatternIndex;
        BinNodeSourcePatternNameID;
        BinNodeSourceQuality;
        BinNodeSourceTypeIndex;
        BinNodeJunctionsBaseDemands;
        BinNodeJunctionsBaseDemandsID;
        BinNodeDemandPatternNameID;
        BinNodeInitialQuality;
        BinNodeCoordinates;
        BinNodeBaseDemands;
            
        BinPatternNameID;
        BinPatternLengths;
        BinPatternValue; 
        BinPatternMatrix;
        BinPatternCount;

        BinLinkNameID;
        BinLinkType;
        BinLinkPipeNameID;
        BinLinkPumpNameID;
        BinLinkValveNameID;
        BinLinkInitialStatusNameID;
        BinLinkInitialStatus;
        BinLinkPumpPatterns;
        BinLinkPumpCurveNameID;
        BinLinkPipeIndex;
        BinLinkPumpIndex;
        BinLinkValveIndex;
        BinLinkFromNode;
        BinLinkToNode;
        BinLinkCount;
        BinLinkPipeCount;
        BinLinkPumpCount;
        BinLinkValveCount;
        BinLinkPipeLengths;
        BinLinkPipeDiameters;
        BinLinkPipeRoughness;
        BinLinkPipeMinorLoss;
        BinLinkValveDiameters;
        BinLinkValveType;
        BinLinkValveSetting;
        BinLinkValveMinorLoss;
        BinLinkValveStatus;
        BinLinkValveStatusNameID;
        BinLinkPipeStatus;
        BinLinkPumpStatusNameID;   
        BinLinkPumpStatus;  
        BinLinkPumpPower;
        BinLinkPumpNameIDPower;
        BinLinkSettings;
        BinLinkDiameters;
        BinLinkLengths;
        BinLinkRoughnessCoeff;
        BinLinkGlobalWallReactionCoeff;
        BinLinkGlobalBulkReactionCoeff;
        BinLinkWallReactionCoeff;
        BinLinkBulkReactionCoeff;
        BinCurveNameID; %ID name of curves
        BinCurveXvalue; %X-value of curves
        BinCurveYvalue; %Y-value of curves
        BinCurveAllLines;
        BinCurveTypes; % Type of curves 
        BinCurveCount; % Number of curves
        
        BinControlsInfo;
        BinControlLinksID;
        BinControlNodesID;
        BinControlRulesCount;

        BinRulesControlsInfo;
        BinRulesControlLinksID;
        BinRulesControlNodesID;
        BinRulesCount;

        BinTimeSimulationDuration;
        BinTimeHydraulicStep;
        BinTimeQualityStep;
        BinTimePatternStep;
        BinTimePatternStart;
        BinTimeReportingStep;
        BinTimeReportingStart;
        BinTimeStatisticsIndex;
        BinTimeStatistics;

        BinOptionsSpecificGravity;
        BinOptionsViscosity;
        BinOptionsMaxTrials;
        BinOptionsAccuracyValue;
        BinOptionsUnbalanced;
        BinOptionsPattern;
        BinOptionsPatternDemandMultiplier;
        BinOptionsEmitterExponent;
        BinQualityType;% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        BinQualityCode;
        BinQualityTraceNodeIndex;
        BinQualityTraceNodeID;
        BinQualityUnits;
        BinOptionsDiffusivity;
        
        BincountStatuslines;
        BincountInitialQualitylines;
        BincountReactionlines;
        BincountPatternlines;
        BinSImetric;
        BinUScustomary;
        BinUnits;
        BinLinkFlowUnits;
        BinOptionsHeadloss;
        BinOptionsQualityTolerance;
        BinNodePressureUnits;
        
        % EMC Version
        emcversion;
    end
    properties (Constant = true)
        TYPECONTROL={'LOWLEVEL','HIGHLEVEL', 'TIMER', 'TIMEOFDAY'}; % Constants for control: 'LOWLEVEL','HILEVEL', 'TIMER', 'TIMEOFDAY'
        TYPECURVE={'PUMP','EFFICIENCY','VOLUME','HEADLOSS'}; % Constants for pump curves: 'PUMP','EFFICIENCY','VOLUME','HEADLOSS'
        TYPELINK={'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'}; % Constants for links: 'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV', 'VALVE'
        TYPEMIXMODEL={'MIX1','MIX2', 'FIFO','LIFO'}; % Constants for mixing models: 'MIX1','MIX2', 'FIFO','LIFO'
        TYPENODE={'JUNCTION','RESERVOIR', 'TANK'}; % Contants for nodes: 'JUNCTION','RESERVOIR', 'TANK'
        TYPEPUMP={'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'}; % Constants for pumps: 'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'
        TYPEQUALITY={'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'}; % Constants for quality: 'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'
        TYPEREPORT={'YES','NO','FULL'}; % Constants for report: 'YES','NO','FULL'
        TYPESOURCE={'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}; % Constants for sources: 'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'
        TYPESOURCEMSX={'NOSOURCE','CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}; % Constants for sources: 'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'
        TYPESTATS={'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'}; % Constants for statistics: 'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
        TYPEUNITS={'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'}; % Constants for units: 'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'
        
        MsxTYPEAREAUNITS={'FT2','M2','CM2'}; % sets the units used to express pipe wall surface area
        MsxTYPERATEUNITS={'SEC','MIN','HR','DAY'}; % is the units in which all reaction rate terms are expressed
        MsxTYPESOLVER={'EUL','RK5','ROS2'}; % is the choice of numerical integration method used to solve the reaction system
        MsxTYPECOUPLING={'FULL','NONE'}; % is the choice of numerical integration method used to solve the reaction system
        MsxTYPECOMPILER={'NONE','VC','GC'}; % is the choice of numerical integration method used to solve the reaction system
    end
    methods
        function obj = epanet(varargin)
            %Constructor of the EPANET Class
            try unloadlibrary('epanet2');catch e; end
            try unloadlibrary('epanetmsx');catch e; end
            % DLLs
            if strcmp(computer('arch'),'win64')% if no DLL is given, select one automatically
                if exist('64bit')==7 % name is a folder.
                    obj.libepanetpath = [pwd,'\64bit\'];
                else
                    warning('Folder "64bit" does not exit.');return;
                end
            elseif strcmp(computer('arch'),'win32')
                if exist('32bit')==7
                    obj.libepanetpath = [pwd,'\32bit\'];
                else
                    warning('Folder "32bit" does not exit.');return;
                end
            end
            warning on;   
            obj.inputfile=varargin{1}; % Get name of INP file
            % Bin functions
            if nargin==2
                if strcmp(upper(varargin{2}),'BIN')
                    obj.Bintempfile=[obj.inputfile(1:end-4),'_temp.inp'];
                    copyfile(obj.inputfile,obj.Bintempfile);
                    value=obj.getBinCurvesInfo;
                    if ~isempty(value.BinCurveNameID), obj.remAddBinCurvesID(obj.Bintempfile);end
                    obj.inputfile=obj.Bintempfile;
                    obj.Bin=0;
                    obj = BinUpdateClass(obj);
                    obj.saveBinInpFile;
                    return;
                end
            end
            obj.Bin=1;
            if nargin==2
                if ~isempty(find(obj.inputfile==' '))
                    warning(['File "', obj.inputfile, '" is not a valid']);return;
                end
                obj.libepanet=varargin{2}; % Get DLL libepanet (e.g. epanet20012x86 for 32-bit)
                warning off;
                try  loadlibrary([obj.libepanetpath,obj.libepanet],[obj.libepanetpath,obj.libepanet,'.h']); 
                catch e
                   warning on; 
                   obj.errcode=-1;
                   warning(['File "', obj.libepanet, '" is not a valid win application.']);return;
                end
                warning on;
                
            elseif nargin==1
                obj.libepanet = 'epanet2';
                [~,inp]=fileparts(obj.inputfile);
                if ~isempty(find(inp==' '))
                    warning(['File "', obj.inputfile, '" is not a valid.']);return;
                end
            end
            if ~exist(obj.inputfile,'file')
                warning(['File "', obj.inputfile, '" does not exist in folder.']);return;
            end
            if strcmp(computer('arch'),'win64') || strcmp(computer('arch'),'win32')
                %Load EPANET Library
                warning off;
                ENLoadLibrary(obj.libepanetpath,obj.libepanet);
                warning on;
                %Open the file
                obj.errcode=ENopen(obj.inputfile,[obj.inputfile(1:end-4),'.txt'],[obj.inputfile(1:end-4),'.bin'],obj.libepanet);
                if obj.errcode~=0
                    warning('Could not open the file, please check INP file.');return;
                end
                %Save the temporary input file
                obj.Bintempfile=[obj.inputfile(1:end-4),'_temp.inp'];
                obj.saveInputFile(obj.Bintempfile,1); %create a new INP file (Working Copy) using the SAVE command of EPANET
                obj.closeNetwork;  %ENclose; %Close input file
                %Load temporary file
                ENopen(obj.Bintempfile,[obj.inputfile(1:end-4),'_temp.txt'], [obj.inputfile(1:end-4),'_temp.bin'],obj.libepanet);
                obj.pathfile='';
            end
            
            obj.emcversion='dev-2.1';
            % Get type of the parameters
            obj.LinkType=obj.getLinkType;
            obj.NodeType=obj.getNodeType;
            %Get all the countable network parameters
            obj.NodeCount = obj.getNodeCount;
            obj.NodeTankReservoirCount = obj.getNodeTankReservoirCount;
            obj.LinkCount = obj.getLinkCount;
            obj.PatternCount = obj.getPatternCount;
            obj.CurveCount = obj.getCurveCount;
            obj.ControlRulesCount = obj.getControlRulesCount;
            obj.NodeReservoirCount = obj.getNodeReservoirCount;
            obj.NodeTankCount = obj.getNodeTankCount;
            obj.NodeJunctionCount = obj.getNodeJunctionCount;
            obj.LinkPipeCount = obj.getLinkPipeCount;
            obj.LinkPumpCount = obj.getLinkPumpCount;
            obj.LinkValveCount = obj.getLinkValveCount;
            %Get all the controls
            obj.Controls = obj.getControls;
            %Get the flow units
            obj.LinkFlowUnits = obj.getFlowUnits;
            %Get all the link data
            obj.LinkNameID = obj.getLinkNameID;
            obj.LinkIndex = obj.getLinkIndex;
            obj.LinkTypeIndex = obj.getLinkTypeIndex;
            obj.LinkPipeIndex = find(strcmp(obj.LinkType,'PIPE'));
            obj.LinkPumpIndex = find(strcmp(obj.LinkType,'PUMP'));
            obj.LinkValveIndex = find(obj.LinkTypeIndex>2);
            obj.LinkDiameter = obj.getLinkDiameter;
            obj.LinkLength = obj.getLinkLength;
            obj.LinkRoughnessCoeff = obj.getLinkRoughnessCoeff;
            obj.LinkMinorLossCoeff = obj.getLinkMinorLossCoeff;
            obj.LinkInitialStatus  = obj.getLinkInitialStatus;
            obj.LinkInitialSetting = obj.getLinkInitialSetting;
            obj.LinkBulkReactionCoeff = obj.getLinkBulkReactionCoeff;
            obj.LinkWallReactionCoeff = obj.getLinkWallReactionCoeff;
            obj.LinkPipeNameID = obj.LinkNameID(obj.LinkPipeIndex);
            obj.LinkPumpNameID = obj.LinkNameID(obj.LinkPumpIndex);
            obj.LinkValveNameID = obj.LinkNameID(obj.LinkValveIndex);
            obj.LinkStatus = obj.getLinkStatus;
            %Get all the node data
            obj.NodesConnectingLinksIndex = obj.getLinkNodesIndex;
            obj.NodesConnectingLinksID = obj.getNodesConnectingLinksID;
            obj.NodeNameID = obj.getNodeNameID;
            obj.NodeIndex = obj.getNodeIndex;
            obj.NodeTypeIndex = obj.getNodeTypeIndex;
            obj.NodeElevations = obj.getNodeElevations;
            obj.NodeDemandPatternIndex = obj.getNodeDemandPatternIndex;
            obj.NodeEmitterCoeff = obj.getNodeEmitterCoeff;
            obj.NodeInitialQuality = obj.getNodeInitialQuality;
            obj.NodeSourceQuality = obj.getNodeSourceQuality;
            obj.NodeSourcePatternIndex = obj.getNodeSourcePatternIndex;
            obj.NodeSourceTypeIndex = obj.NodeSourceTypeIndex;
            obj.NodeReservoirIndex = find(strcmp(obj.NodeType,'RESERVOIR'));
            obj.NodeTankIndex = find(strcmp(obj.NodeType,'TANK'));
            obj.NodeJunctionIndex = find(strcmp(obj.NodeType,'JUNCTION'));
            obj.NodeReservoirNameID=obj.NodeNameID(obj.NodeReservoirIndex);
            obj.NodeTankNameID=obj.NodeNameID(obj.NodeTankIndex);
            obj.NodeJunctionNameID=obj.NodeNameID(obj.NodeJunctionIndex);
            %Get all tank data
            obj.NodeTankInitialLevel = obj.getNodeTankInitialLevel;
            obj.NodeTankMinimumWaterVolume = obj.getNodeTankInitialWaterVolume;
            obj.NodeTankMixingModelCode = obj.getNodeTankMixingModelCode;
            obj.NodeTankMixingModelType = obj.getNodeTankMixingModelType;
            obj.NodeTankMixZoneVolume = obj.getNodeTankMixZoneVolume;
            obj.NodeTankDiameter = obj.getNodeTankDiameter;
            obj.NodeTankMinimumWaterVolume = obj.getNodeTankMinimumWaterVolume;
            obj.NodeTankVolumeCurveIndex = obj.getNodeTankVolumeCurveIndex;
            obj.NodeTankMinimumWaterLevel = obj.getNodeTankMinimumWaterLevel;
            obj.NodeTankMaximumWaterLevel = obj.getNodeTankMaximumWaterLevel;
            obj.NodeTankMinimumFraction = obj.getNodeTankMinimumFraction;
            obj.NodeTankBulkReactionCoeff = obj.getNodeTankBulkReactionCoeff;
            %Get all options
            obj.OptionsMaxTrials = obj.getOptionsMaxTrials;
            obj.OptionsAccuracyValue = obj.getOptionsAccuracyValue;
            obj.OptionsQualityTolerance = obj.getOptionsQualityTolerance;
            obj.OptionsEmitterExponent = obj.getOptionsEmitterExponent;
            obj.OptionsPatternDemandMultiplier = obj.getOptionsPatternDemandMultiplier;
            %Get pattern data
            obj.PatternNameID = obj.getPatternNameID;
            obj.PatternIndex = obj.getPatternIndex;
            obj.PatternLengths = obj.getPatternLengths;
            obj.Pattern = obj.getPattern;
            %Get quality types
            obj.QualityCode = obj.getQualityCode;
            obj.QualityTraceNodeIndex = obj.getQualityTraceNodeIndex;
            %Bug
            obj.QualityType = obj.getQualityType;
%             n = obj.getQualityInfo;
%             obj.QualityChemUnits = n.QualityChemUnits;
%             obj.QualityChemName= n.QualityChemUnits;
            
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
            try%New version dev2.1.dll libepanet
                obj.TimeStartTime = obj.getTimeStartTime;
                obj.TimeHTime = obj.getTimeHTime;
                obj.TimeHaltFlag = obj.getTimeHaltFlag;
                obj.TimeNextEvent = obj.getTimeNextEvent;
                obj.NodeTankMaxVolume = obj.getNodeTankMaxVolume;
                if sum(strcmp(libfunctions(obj.libepanet),'ENgetbasedemand'))
                    obj.NodeBaseDemands = obj.getNodeBaseDemands;
                    obj.NodeNumDemandCategories = obj.getNodeNumDemandCategories;
                else
                    obj.NodeBaseDemands={};
                    for i=1:obj.NodeCount
                        [obj.errcode, obj.NodeBaseDemands{i}] = ENgetnodevalue(i,1,obj.libepanet);
                    end
                end
                obj.PatternAveragePatternValue = obj.getPatternAveragePatternValue;
                n = obj.getStatistic;
                obj.RelativeError = n.RelativeError;
                obj.Iterations = n.Iterations;
                obj.NodeDemandPatternNameID = obj.getNodeDemandPatternNameID;
                obj.NodeDemandPatternsIndex = obj.getNodeDemandPatternsIndex;
                obj.HeadCurveIndex = obj.getHeadCurveIndex;
                obj.LinkPumpPatternNameID = obj.getLinkPumpPatternNameID;
                obj.LinkPumpPatternIndex = obj.getLinkPumpPatternIndex;
                obj.LinkPumpTypeCode = obj.getLinkPumpTypeCode;
                obj.LinkPumpType = obj.getLinkPumpType;
%                 obj.CurvesInfo = obj.getCurvesInfo; % New version dev2.1
                %Get data from raw file (for information which cannot be
                %accessed by the epanet library)
                value=obj.getNodeCoordinates;
                %Get coordinates
                obj.NodeCoordinates{1} = value{1};
                obj.NodeCoordinates{2} = value{2};
                obj.NodeCoordinates{3} = value{3};
                obj.NodeCoordinates{4} = value{4};
            catch e
            end
            
            %     US Customary - SI metric
            switch char(obj.LinkFlowUnits)
                case 'CFS'
                    obj.UScustomary=1;
                case 'GPM'
                    obj.UScustomary=1;
                case 'MGD'
                    obj.UScustomary=1;
                case 'IMGD'
                    obj.UScustomary=1;
                case 'AFD'
                    obj.UScustomary=1;
                case 'LPS'
                    obj.SImetric=1;
                case 'LPM'
                    obj.SImetric=1;
                case 'MLD'
                    obj.SImetric=1;
                case 'CMH'
                    obj.SImetric=1;
                case 'CMD'
                    obj.SImetric=1;
            end

            if obj.UScustomary==1;
                obj.NodePressureUnits='pounds per square inch';
                obj.PatternDemandsUnits=obj.LinkFlowUnits;
                obj.LinkPipeDiameterUnits='inches';
                obj.NodeTankDiameterUnits='feet';
                obj.EnergyEfficiencyUnits='percent';
                obj.NodeElevationUnits='feet';
                obj.NodeEmitterCoefficientUnits='flow units @ 1 psi drop';
                obj.EnergyUnits='kwatt-hours';
                obj.LinkFrictionFactorUnits='unitless';
                obj.NodeHeadUnits='feet';
                obj.LinkLengthsUnits='feet';
                obj.LinkMinorLossCoeffUnits='unitless';
                obj.LinkPumpPowerUnits='horsepower';
                obj.QualityReactionCoeffBulkUnits='1/day (1st-order)';
                obj.QualityReactionCoeffWallUnits='mass/sq-ft/day (0-order), ft/day (1st-order)';
                obj.LinkPipeRoughnessCoeffUnits='millifeet(Darcy-Weisbach), unitless otherwise';
                obj.QualitySourceMassInjectionUnits='mass/minute';
                obj.LinkVelocityUnits='ft/sec';
                obj.NodeTankVolumeUnits='cubic feet';
                obj.QualityWaterAgeUnits='hours';
            else % SI Metric
                obj.NodePressureUnits='meters';
                obj.PatternDemandsUnits=obj.LinkFlowUnits;
                obj.LinkPipeDiameterUnits='millimeters';
                obj.NodeTankDiameterUnits='meters';
                obj.EnergyEfficiencyUnits='percent';
                obj.NodeElevationUnits='meters';
                obj.NodeEmitterCoefficientUnits='flow units @ 1 meter drop';
                obj.EnergyUnits='kwatt-hours';
                obj.LinkFrictionFactorUnits='unitless';
                obj.NodeHeadUnits='meters';
                obj.LinkLengthsUnits='meters';
                obj.LinkMinorLossCoeffUnits='unitless';
                obj.LinkPumpPowerUnits='kwatts';
                obj.QualityReactionCoeffBulkUnits='1/day (1st-order)';
                obj.QualityReactionCoeffWallUnits='mass/sq-m/day(0-order), meters/day (1st-order)';
                obj.LinkPipeRoughnessCoeffUnits='mm(Darcy-Weisbach), unitless otherwise';
                obj.QualitySourceMassInjectionUnits='mass/minute';
                obj.LinkVelocityUnits='meters/sec';
                obj.NodeTankVolumeUnits='cubic meters';
                obj.QualityWaterAgeUnits='hours';
            end
            
        end % End of epanet class constructor
        function errcode = LoadFile(obj,varargin)
           [errcode] = ENopen(varargin{1},[varargin{1}(1:end-4),'.rpt'],[varargin{1}(1:end-4),'.bin'],obj.libepanet); 
        end
        function [value] = plot(obj,varargin)
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
            % d.plot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
            % d.plot('line','no');
            % d.plot('point','no','linksindex','yes');
            % d.plot('linksindex','yes','fontsize',8);
            % d.plot('nodesindex','yes','fontsize',14);
            [value] = ENplot(obj,'bin',0,varargin{:});
        end
        function value = getControls(obj)
            %Retrieves the parameters of all control statements
            cnt=obj.getControlRulesCount;
            if cnt
                obj.ControlTypes{cnt}=[];
                obj.ControlTypesIndex(cnt)=NaN;
                obj.ControlLinkIndex(cnt)=NaN;
                obj.ControlSettings(cnt)=NaN;
                obj.ControlNodeIndex(cnt)=NaN;
                obj.ControlLevelValues(cnt)=NaN;
                for i=1:cnt
                    [obj.errcode, obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)] = ENgetcontrol(i,obj.libepanet);
                    obj.ControlTypes(i)={obj.TYPECONTROL(obj.ControlTypesIndex(i)+1)};
                    value{i}={obj.ControlTypes{i},obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)};
                end
            else
                value=-1;
            end
        end
        function value = getNodeCount(obj)
            % Retrieves the number of nodes
            [obj.errcode, value] = ENgetcount(0,obj.libepanet);
        end
        function value = getNodeTankReservoirCount(obj)
            % Retrieves the number of tanks
            [obj.errcode, value] = ENgetcount(1,obj.libepanet);
        end
        function value = getLinkCount(obj)
            % Retrieves the number of links
            [obj.errcode, value] = ENgetcount(2,obj.libepanet);
        end
        function value = getPatternCount(obj)
            % Retrieves the number of patterns
            [obj.errcode, value] = ENgetcount(3,obj.libepanet);
        end
        function value = getCurveCount(obj)
            % Retrieves the number of curves
            [obj.errcode, value] = ENgetcount(4,obj.libepanet);
        end
        function value = getControlRulesCount(obj)
            % Retrieves the number of controls
            [obj.errcode, value] = ENgetcount(5,obj.libepanet);
        end
        function value = getNodeTankCount(obj)
            % Retrieves the number of Tanks
            value = sum(strcmp(obj.getNodeType,'TANK'));
        end
        function value = getNodeReservoirCount(obj)
            % Retrieves the number of Reservoirs
            value = sum(strcmp(obj.getNodeType,'RESERVOIR'));
        end
        function value = getNodeJunctionCount(obj)
            % Retrieves the number of junction nodes
            value = sum(strcmp(obj.getNodeType,'JUNCTION'));
        end
        function value = getLinkPipeCount(obj)
            % Retrieves the number of pipes
            value = sum(strcmp(obj.getLinkType,'PIPE'))+sum(strcmp(obj.getLinkType,'CVPIPE'));
        end
        function value = getLinkPumpCount(obj)
            % Retrieves the number of pumps
            value = sum(strcmp(obj.getLinkType,'PUMP'));
        end
        function value = getLinkValveCount(obj)
            % Retrieves the number of valves
            value = obj.getLinkCount - (obj.getLinkPipeCount + obj.getLinkPumpCount);
        end
        function errmssg = getError(obj,errcode)
            %Retrieves the text of the message associated with a particular error or warning code.
            [errmssg , errcode] = ENgeterror(errcode,obj.libepanet);
        end
        function value = getFlowUnits(obj)
            %Retrieves flow units used to express all flow rates.
            [obj.errcode, flowunitsindex] = ENgetflowunits(obj.libepanet);
            obj.LinkFlowUnits=obj.TYPEUNITS(flowunitsindex+1);
            value=obj.LinkFlowUnits;
        end
        function value = getLinkNameID(obj,varargin)
            % Retrieves the ID label(s) of all links, or the IDs of an index set of links
            if isempty(varargin)
                value{obj.getLinkCount}=[];
                for i=1:obj.getLinkCount
                    [obj.errcode, value{i}]=ENgetlinkid(i,obj.libepanet);
                end
            else
                k=1;
                if ~isempty(varargin{1})
                    value{length(varargin{1})}=[];
                    for i=varargin{1}
                        [obj.errcode, value{k}]=ENgetlinkid(i,obj.libepanet);
                        k=k+1;
                    end
                else
                    value=[];
                end
            end
        end
        function value = getLinkPipeNameID(obj)
            %Retrieves the pipe id
            value=obj.getLinkNameID(obj.getLinkPipeIndex);
        end
        function value = getLinkPumpNameID(obj)
            %Retrieves the pump id
            value=obj.getLinkNameID(obj.getLinkPumpIndex);
        end
        function value = getLinkValveNameID(obj)
            %Retrieves the valve id
            value=obj.getLinkNameID(obj.getLinkValveIndex);
        end
        function value = getLinkIndex(obj,varargin)
            %Retrieves the indices of all links, or the indices of an ID set of links
            if isempty(varargin)
                value=1:obj.getLinkCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetlinkindex(varargin{1}{j},obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetlinkindex(varargin{1},obj.libepanet);
            end
        end
        function value = getLinkPipeIndex(obj)
            %Retrieves the pipe index
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PIPE'));
%             if isempty(value), value=-1; end
        end
        function value = getLinkPumpIndex(obj)
            %Retrieves the pipe index
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PUMP'));
%             if isempty(value), value=-1; end
        end
        function value = getLinkValveIndex(obj)
            %Retrieves the pipe index
            value = obj.getLinkPipeCount+obj.getLinkPumpCount+1:obj.getLinkCount;
%             if isempty(value), value=-1; end
        end
        function value = getLinkNodesIndex(obj)
            %Retrieves the indexes of the from/to nodes of all links.
            value(obj.getLinkCount,1:2)=[nan nan];
            for i=1:obj.getLinkCount
                [obj.errcode,linkFromNode,linkToNode] = ENgetlinknodes(i,obj.libepanet);
                value(i,:)= [linkFromNode,linkToNode];
            end
        end
        function value = getNodesConnectingLinksID(obj)
            %Retrieves the id of the from/to nodes of all links.
            obj.NodesConnectingLinksIndex=obj.getLinkNodesIndex;
            if obj.getLinkCount
                for i=1:obj.getLinkCount
                    value(i,1)=obj.getNodeNameID(obj.NodesConnectingLinksIndex(i,1));
                    value(i,2) = obj.getNodeNameID(obj.NodesConnectingLinksIndex(i,2));
                end
            else
                value=-1;
            end
        end
        function value = getLinkType(obj)
            %Retrieves the link-type code for all links.
            index=obj.getLinkTypeIndex;
            for i=1:obj.getLinkCount
                value(i)=obj.TYPELINK(index(i)+1);
            end
        end
        function value = getLinkTypeIndex(obj)
            %Retrieves the link-type code for all links.
            for i=1:obj.getLinkCount
                [obj.errcode,value(i)] = ENgetlinktype(i,obj.libepanet);
            end
        end
        function value = getLinkDiameter(obj)
            %Retrieves the value of all link diameters
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,0,obj.libepanet);
            end
        end
        function value = getLinkLength(obj)
            %Retrieves the value of all link lengths
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,1,obj.libepanet);
            end
        end
        function value = getLinkRoughnessCoeff(obj)
            %Retrieves the value of all link roughness
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,2,obj.libepanet);
            end
        end
        function value = getLinkMinorLossCoeff(obj)
            %Retrieves the value of all link minor loss coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,3,obj.libepanet);
            end
        end
        function value = getLinkInitialStatus(obj)
            %Retrieves the value of all link initial status
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,4,obj.libepanet);
            end
        end
        function value = getLinkInitialSetting(obj)
            %Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,5,obj.libepanet);
            end
        end
        function value = getLinkBulkReactionCoeff(obj)
            %Retrieves the value of all link bulk reaction coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,6,obj.libepanet);
            end
        end
        function value = getLinkWallReactionCoeff(obj)
            %Retrieves the value of all link wall reaction coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,7,obj.libepanet);
            end
        end
        function value = getLinkFlows(obj)
            %Retrieves the value of all computed link flow rates
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,8,obj.libepanet);
            end
        end
        function value = getLinkVelocity(obj)
            %Retrieves the value of all computed link velocities
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,9,obj.libepanet);
            end
        end
        function value = getLinkHeadloss(obj)
            %Retrieves the value of all computed link headloss
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,10,obj.libepanet);
            end
        end
        function value = getLinkStatus(obj)
            %Retrieves the value of all computed link status (0 = closed, 1 = open)
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,11,obj.libepanet);
            end
        end
        function value = getLinkSettings(obj)
            %Retrieves the value of all computed link roughness for pipes or actual speed for pumps or actual setting for valves
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,12,obj.libepanet);
            end
        end
        function value = getLinkPumpEnergy(obj)
            %Retrieves the value of all computed energy in kwatts
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,13,obj.libepanet);
            end
        end
        function value = getLinkQuality(obj)
            %New version dev2.1
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,14,obj.libepanet);
            end
        end
        function value = getLinkPumpPatternIndex(obj)
            %New version dev2.1
            value=zeros(1,obj.getLinkPumpCount);v=1;
            for i=obj.getLinkPumpIndex
                [obj.errcode, value(v)] = ENgetlinkvalue(i,15,obj.libepanet);
                v=v+1;
            end
        end
        function value = getLinkPumpPatternNameID(obj)
            %New version dev2.1
            v = obj.getLinkPumpPatternIndex;
            value = obj.getPatternNameID(v);
        end
        function value = getNodeNameID(obj,varargin)
            %Retrieves the ID label of all nodes or some nodes with a specified index.
            if isempty(varargin)
                value{obj.getNodeCount}=[];
                for i=1:obj.getNodeCount
                    [obj.errcode, value{i}]=ENgetnodeid(i,obj.libepanet);
                end
            else
                k=1;
                value{length(varargin{1})}=[];
                for i=varargin{1}
                    [obj.errcode, value{k}]=ENgetnodeid(i,obj.libepanet);
                    k=k+1;
                end
            end
        end
        function value = getNodeReservoirNameID(obj)
            %Retrieves the reservoir id
            value=obj.getNodeNameID(obj.getNodeReservoirIndex);
        end
        function value = getNodeJunctionNameID(obj)
            %Retrieves the junction id
            value=obj.getNodeNameID(obj.getNodeJunctionIndex);
        end
        function value = getNodeIndex(obj,varargin)
            %Retrieves the indices of all nodes or some nodes with a specified ID
            if isempty(varargin)
                value=1:obj.getNodeCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetnodeindex(varargin{1}{j},obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetnodeindex(varargin{1},obj.libepanet);
            end
        end
        function value = getNodeReservoirIndex(obj)
            %Retrieves the reservoir index
            if obj.getNodeReservoirCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'RESERVOIR'));
            else
                value=0;
            end
        end
        function value = getNodeJunctionIndex(obj)
            %Retrieves the junction index
            if obj.getNodeJunctionCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'JUNCTION'));
            else
                value=0;
            end
        end
        function value = getNodeType(obj)
            %Retrieves the node-type code for all nodes
            for i=1:obj.getNodeCount
                [obj.errcode,obj.NodeTypeIndex(i)] = ENgetnodetype(i,obj.libepanet);
                value(i)=obj.TYPENODE(obj.NodeTypeIndex(i)+1);
            end
        end
        function value = getNodeTypeIndex(obj)
            %Retrieves the node-type code for all nodes
            for i=1:obj.getNodeCount
                [obj.errcode,value(i)] = ENgetnodetype(i,obj.libepanet);
            end
        end
        function value = getNodeElevations(obj)
            %Retrieves the value of all node elevations
            ndcount = obj.getNodeCount;
            value=zeros(1,ndcount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,0,obj.libepanet);
            end
        end
        function value = getNodeBaseDemands(obj)
            %New version dev2.1
            chckfunctions=libfunctions(obj.libepanet);
            if sum(strcmp(chckfunctions,'ENgetbasedemand'))
                numdemands = obj.getNodeNumDemandCategories;
                val=zeros(max(numdemands),obj.getNodeCount);
                for i=obj.getNodeIndex
                    v=1;
                    for u=1:numdemands(i)
                        [obj.errcode, val(v,i)] = ENgetbasedemand(i,u,obj.libepanet);v=v+1;
                    end
                end
                for i=1:size(val,1)
                    value{i} = val(i,:);
                end
            else%version epanet20012
                %Retrieves the value of all node base demands
                value=zeros(1,obj.getNodeCount);
                for i=1:obj.getNodeCount
                    [obj.errcode, value(i)] = ENgetnodevalue(i,1,obj.libepanet);
                end
            end
        end
        function value = getNodeNumDemandCategories(obj)
            %New version dev2.1
            value=zeros(1,obj.getNodeCount);
            for i=obj.getNodeIndex
                [obj.errcode, value(i)] = ENgetnumdemands(i,obj.libepanet);
            end
        end
        function value = getNodeDemandPatternsIndex(obj)
            %New version dev2.1
            numdemands = obj.getNodeNumDemandCategories;
            val=zeros(max(numdemands),obj.getNodeCount);
            for i=obj.getNodeIndex
                v=1;
                for u=1:numdemands(i)
                    [obj.errcode, val(v,i)] = ENgetdemandpattern(i,u,obj.libepanet);v=v+1;
                end             
            end
            for i=1:size(val,1)
                value{i} = val(i,:);
            end
        end
        function value = getNodeDemandPatternNameID(obj)
            %New version dev2.1
            value={};
            v = obj.getNodeDemandPatternsIndex;
            m = obj.getPatternNameID;
            numdemands = obj.getNodeNumDemandCategories;
            for i=1:obj.getNodeCount
                for u=1:numdemands(i)
                    if v{u}(i)~=0 
                        value{u,i}= char(m(v{u}(i))); 
                    else
                        value{u,i}= '';
                    end
                end
                if numdemands(i)==0
                    value{1,i}= [];
                end
            end
        end
        function value = getStatistic(obj)
            %New version dev2.1
            %Input: none
            %Output: *iter = # of iterations to reach solution
            %*relerr = convergence error in solution
            %returns error code
            [obj.errcode, value.Iterations] = ENgetstatistic(0,obj.libepanet);
            [obj.errcode, value.RelativeError] = ENgetstatistic(1,obj.libepanet);
        end    
        function value = getNodeDemandPatternIndex(obj)
            %Retrieves the value of all node demand pattern indices
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,2,obj.libepanet);
            end
        end
        function value = getNodeEmitterCoeff(obj)
            %Retrieves the value of all node emmitter coefficients
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,3,obj.libepanet);
            end
        end
        function value = getNodeInitialQuality(obj)
            %Retrieves the value of all node initial quality
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,4,obj.libepanet);
            end
        end
        function value = getNodeSourceQuality(obj)
            %Retrieves the value of all nodes source quality
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,5,obj.libepanet);
            end
        end
        function value = getNodeSourcePatternIndex(obj)
            %Retrieves the value of all node source pattern index
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,6,obj.libepanet);
            end
        end
        function value = getNodeSourceType(obj)
            %Retrieves the value of all node source type
            value=cell(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, temp] = ENgetnodevalue(i,7,obj.libepanet);
                if ~isnan(temp)
                    value(i)=obj.TYPESOURCE(temp+1);
                end
            end
        end
        function value = getNodeTankInitialLevel(obj)
            %Retrieves the value of all tank initial water levels
            value=nan(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,8,obj.libepanet);
            end
        end
        function value = getNodeActualDemand(obj)
            %Retrieves the computed value of all actual demands
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,9,obj.libepanet);
            end
        end
        function value = getNodeActualDemandSensingNodes(obj,varargin)
            %Retrieves the computed demand values at some sensing nodes
            value=zeros(1,length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.errcode, value(v)] = ENgetnodevalue(i,9,obj.libepanet);v=v+1;
            end
        end
        function value = getNodeHydaulicHead(obj)
            %Retrieves the computed values of all hydraulic heads
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,10,obj.libepanet);
            end
        end
        function value = getNodePressure(obj)
            %Retrieves the computed values of all node pressures
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,11,obj.libepanet);
            end
        end
        function value = getNodeActualQuality(obj)
            %Retrieves the computed values of the actual quality for all nodes
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,12,obj.libepanet);
            end
        end
        function value = getNodeMassFlowRate(obj)
            %Retrieves the computed mass flow rates per minute of chemical sources
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,13,obj.libepanet);
            end
        end
        function value = getNodeActualQualitySensingNodes(obj,varargin)
            %Retrieves the computed quality values at some sensing nodes
            value=zeros(1,length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.errcode, value(v)] = ENgetnodevalue(i,12,obj.libepanet);v=v+1;
            end
        end
        function value = getNodeTankInitialWaterVolume(obj)
            %Retrieves the tank initial volume
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i,14,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMixiningModel(obj)
            %Retrieves the tank mixing mode (mix1, mix2, fifo, lifo)
            obj.NodeTankInitialWaterVolume=zeros(1,obj.getNodeCount);
            obj.NodeTankMixingModelCode=nan(1,obj.getNodeCount);
            obj.NodeTankMixingModelType={};
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, obj.NodeTankMixingModelCode(i)] = ENgetnodevalue(i, 15,obj.libepanet);
                    obj.NodeTankMixingModelType(i)=obj.TYPEMIXMODEL(obj.NodeTankMixingModelCode(i)+1);
                end
            end
            value={obj.NodeTankMixingModelCode obj.NodeTankMixingModelType};
        end
        function value = getNodeTankMixingModelCode(obj)
            %Retrieves the tank mixing model code (mix1, mix2, fifo, lifo)
            value = obj.getNodeTankMixiningModel{1};
        end
        function value = getNodeTankMixingModelType(obj)
            %Retrieves the tank mixing model type (mix1, mix2, fifo, lifo)
            value = obj.getNodeTankMixiningModel{2};
        end
        function value = getNodeTankMixZoneVolume(obj)
            %Retrieves the tank mixing zone volume
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i,16,obj.libepanet);
                end
            end
        end
        function value = getNodeTankDiameter(obj)
            %Retrieves the tank diameters
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 17,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMinimumWaterVolume(obj)
            %Retrieves the tank minimum volume
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 18,obj.libepanet);
                end
            end
        end
        function value = getNodeTankVolumeCurveIndex(obj)
            %Retrieves the tank volume curve index
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 19,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMinimumWaterLevel(obj)
            %Retrieves the tank minimum water level
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 20,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMaximumWaterLevel(obj)
            %Retrieves the tank maximum water level
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 21,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMinimumFraction(obj)
            %Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 22,obj.libepanet);
                end
            end
        end
        function value = getNodeTankBulkReactionCoeff(obj)
            %Retrieves the tank bulk rate coefficient
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 23,obj.libepanet);
                end
            end
        end
        function value = getNodeTankVolume(obj)
            %New version dev2.1
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 24,obj.libepanet);
                end
            end
        end
        function value = getNodeTankMaxVolume(obj)
            %New version dev2.1
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 25,obj.libepanet);
                end
            end
        end
        function value = getNodeTankIndex(obj)
            %Retrieves the tank index
            if obj.getNodeTankCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'TANK'));
            else
                value=0;
            end
        end
        function value = getNodeTankNameID(obj)
            %Retrieves the tank id
            value=obj.getNodeNameID(obj.getNodeTankIndex);
        end
        function value = getOptionsMaxTrials(obj)
            % Retrieve maximum number of analysis trials
            [obj.errcode, value] = ENgetoption(0,obj.libepanet);
        end
        function value = getOptionsAccuracyValue(obj)
            % Retrieve the analysis convergence criterion (0.001)
            [obj.errcode, value] = ENgetoption(1,obj.libepanet);
        end
        function value = getOptionsQualityTolerance(obj)
            % Retrieve the water quality analysis tolerance
            [obj.errcode, value] = ENgetoption(2,obj.libepanet);
        end
        function value = getOptionsEmitterExponent(obj)
            % Retrieve power exponent for the emmitters (0.5)
            [obj.errcode, value] = ENgetoption(3,obj.libepanet);
        end
        function value = getOptionsPatternDemandMultiplier(obj)
            % Retrieve the demand multiplier (x1)
            [obj.errcode, value] = ENgetoption(4,obj.libepanet);
        end
        function value = getPatternNameID(obj,varargin)
            %Retrieves the ID label of all or some time patterns indices
            if isempty(varargin)
                value{obj.getPatternCount}=[];
                for i=1:obj.getPatternCount
                    [obj.errcode, value{i}]=ENgetpatternid(i,obj.libepanet);
                end
            else
                if obj.getLinkPumpCount
                    k=1;
                    value{length(varargin{1})}=[];
                    for i=varargin{1}
                        [obj.errcode, value{k}]=ENgetpatternid(i,obj.libepanet);
                        k=k+1;
                    end
                else
                    value={};
                end
            end
        end
        function value = getCurveNameID(obj,varargin)
            %Retrieves ID of a curve with specific index
            %New version dev2.1
            if isempty(varargin)
                value{obj.getCurveCount}=[];
                for i=1:obj.getCurveCount
                    [obj.errcode, value{i}]=ENgetcurveid(i,obj.libepanet);
                end
            else
                if obj.getLinkPumpCount
                    k=1;
                    value{length(varargin{1})}=[];
                    for i=varargin{1}
                        [obj.errcode, value{k}]=ENgetcurveid(i,obj.libepanet);
                        k=k+1;
                    end
                else
                    value={};
                end
            end
        end
        function value = getCurveLengths(obj,varargin)
            %Retrieves number of points in a curve 
            %New version dev2.1
            if isempty(varargin)
                tmpCurves=1:obj.getCurveCount;
                for i=tmpCurves
                    [obj.errcode, value(i)]=ENgetcurvelen(i,obj.libepanet);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetcurvelen(obj.getCurveIndex(varargin{1}{j}),obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetcurvelen(obj.getCurveIndex(varargin{1}),obj.libepanet);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errcode, value(k)]=ENgetcurvelen(i,obj.libepanet);
                    k=k+1;
                end
            end
        end
        function value = getCurveIndex(obj,varargin)
            %Retrieves index of curve with specific ID
            %New version dev2.1
            if isempty(varargin)
                value=1:obj.getCurveCount;
            elseif isa(varargin{1},'cell')
                k=1;
                value{length(varargin{1})}=[];
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetcurveindex(varargin{1}{j},obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetcurveindex(varargin{1},obj.libepanet);
            end
        end
        function setCurve(obj,index,curveVector)
            %Sets x,y values for a specific curve
            %New version dev2.1
            nfactors=size(curveVector,1);%x = number of points in curve
            [obj.errcode] = ENsetcurve(index, curveVector(:,1), curveVector(:,2), nfactors,obj.libepanet);
        end
        function value = getCurveValue(obj,curveIndex, curvePnt)
            %Retrieves x,y point for a specific point number and curve
            %New version dev2.1
            [obj.errcode, x, y] = ENgetcurvevalue(curveIndex, curvePnt,obj.libepanet);
            value = [x y];
        end
        function setCurveValue(obj,index, curvePnt, value)
            %Retrieves x,y point for a specific point number and curve
            %New version dev2.1
            x=value(1); y=value(2);
            [obj.errcode] = ENsetcurvevalue(index, curvePnt, x, y, obj.libepanet);
        end
        function value = getPatternIndex(obj,varargin)
            %Retrieves the index of all or some time patterns IDs
            if isempty(varargin)
                value=1:obj.getPatternCount;
            elseif isa(varargin{1},'cell')
                k=1;
                value{length(varargin{1})}=[];
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetpatternindex(varargin{1}{j},obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetpatternindex(varargin{1},obj.libepanet);
            end
        end
        function value = getPatternLengths(obj,varargin)
            %Retrieves the number of time periods in all or some patterns
            if isempty(varargin)
                tmpPatterns=1:obj.getPatternCount;
                for i=tmpPatterns
                    [obj.errcode, value(i)]=ENgetpatternlen(i,obj.libepanet);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetpatternlen(obj.getPatternIndex(varargin{1}{j}),obj.libepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetpatternlen(obj.getPatternIndex(varargin{1}),obj.libepanet);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errcode, value(k)]=ENgetpatternlen(i,obj.libepanet);
                    k=k+1;
                end
            end
        end
        function value = getPattern(obj)
            %Retrieves the multiplier factor for all patterns and all times
            tmpmaxlen=max(obj.getPatternLengths);
            value=nan(obj.getPatternCount,tmpmaxlen);
            for i=1:obj.getPatternCount
                tmplength=obj.getPatternLengths(i);
                for j=1:tmplength
                    [obj.errcode, value(i,j)] = ENgetpatternvalue(i, j,obj.libepanet);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                
            end
        end
        function value = getPatternValue(obj,patternIndex, patternStep)
            %Retrieves the multiplier factor for a certain pattern and time
            [obj.errcode, value] = ENgetpatternvalue(patternIndex, patternStep,obj.libepanet);
        end
        function value = getQualityType(obj,varargin)
            %Retrieves the type of water quality analysis type
% %             [obj.errcode, obj.QualityCode,obj.QualityTraceNodeIndex] = ENgetqualinfo(obj.libepanet); bug
% %             value=obj.TYPEQUALITY(obj.QualityCode+1);
            if nargin>1
                obj.saveInputFile(obj.Bintempfile,1);
            else
                obj.saveInputFile([obj.pathfile,obj.Bintempfile]);
            end
            value = {obj.getBinOptionsInfo.BinQualityType};
        end 
        function value = getQualityInfo(obj)
            [obj.errcode, ~,value.QualityChemName,value.QualityChemUnits,~] = ENgetqualinfo(obj.libepanet);
        end
        function value = getQualityCode(obj)
            %Retrieves the code of water quality analysis type
            [obj.errcode, value,obj.QualityTraceNodeIndex] = ENgetqualtype(obj.libepanet);
        end
        function value = getQualityTraceNodeIndex(obj)
            %Retrieves the trace node index of water quality analysis type
            [obj.errcode, obj.QualityCode,value] = ENgetqualtype(obj.libepanet);
        end
        function value = getTimeSimulationDuration(obj)
            %Retrieves the value of simulation duration
            [obj.errcode, value] = ENgettimeparam(0,obj.libepanet);
        end
        function value = getTimeHydraulicStep(obj)
            %Retrieves the value of the hydraulic time step
            [obj.errcode, value] = ENgettimeparam(1,obj.libepanet);
        end
        function value = getTimeQualityStep(obj)
            %Retrieves the value of the water quality time step
            [obj.errcode, value] = ENgettimeparam(2,obj.libepanet);
        end
        function value = getTimePatternStep(obj)
            %Retrieves the value of the pattern time step
            [obj.errcode, value] = ENgettimeparam(3,obj.libepanet);
        end
        function value = getTimePatternStart(obj)
            %Retrieves the value of pattern start time
            [obj.errcode, value] = ENgettimeparam(4,obj.libepanet);
        end
        function value = getTimeReportingStep(obj)
            %Retrieves the value of the reporting time step
            [obj.errcode, value] = ENgettimeparam(5,obj.libepanet);
        end
        function value = getTimeReportingStart(obj)
            %Retrieves the value of the reporting start time
            [obj.errcode, value] = ENgettimeparam(6,obj.libepanet);
        end
        function value = getTimeRuleControlStep(obj)
            %Retrieves the time step for evaluating rule-based controls
            [obj.errcode, value] = ENgettimeparam(7,obj.libepanet);
        end
        function value = getTimeStatisticsType(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.errcode, obj.TimeStatisticsIndex] = ENgettimeparam(8,obj.libepanet);
            value=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
        end
        function value = getTimeStatisticsIndex(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.errcode, value] = ENgettimeparam(8,obj.libepanet);
        end
        function value = getTimeReportingPeriods(obj)
            %Retrieves the number of reporting periods saved to the binary
            [obj.errcode, value] = ENgettimeparam(9,obj.libepanet);
        end
        %%%%% New version dev2.1 %%%%%
        function value = getTimeStartTime(obj)
            %Retrieves the number of start time
            [obj.errcode, value] = ENgettimeparam(10,obj.libepanet);
        end
        function value = getTimeHTime(obj)
            %Retrieves the number of htime
            [obj.errcode, value] = ENgettimeparam(11,obj.libepanet);
        end
        function value = getTimeHaltFlag(obj)
            %Retrieves the number of  halt flag
            [obj.errcode, value] = ENgettimeparam(12,obj.libepanet);
        end
        function value = getTimeNextEvent(obj)
            %Retrieves the number of next event 
            [obj.errcode, value] = ENgettimeparam(13,obj.libepanet);
        end
%         function value = getCurvesInfo(obj)
%             %New version dev2.1
%             %Input:   curveIndex = curve index
%             %Output:  *nValues = number of points on curve
%             %         *xValues = values for x
%             %         *yValues = values for y
%             %Returns: error code
%             %Purpose: retrieves end nodes of a specific link
%             for i=1:obj.getCurveCount
%                 [obj.errcode, value.CurveNameID{i}, value.CurveNvalue{i}, value.CurveXvalue{i}, value.CurveYvalue{i}] = ENgetcurve(i,obj.libepanet);
%             end
%         end
        function value = getConnectivityMatrix(obj)
            conn = obj.getNodesConnectingLinksID;
            nodesID = obj.getNodeNameID;
            value = zeros(obj.getNodeCount,obj.getNodeCount);
            for i=1:obj.getNodeCount
                mm = strcmp(nodesID(i),conn);
                mS = mm(:,1)+mm(:,2);       
                chIndex = find(mS);
                connFinal = conn(chIndex,:);
                mmFinal = mm(chIndex,:);
                Ok = connFinal(find(mmFinal==0));
                nodesIndOk = obj.getNodeIndex(Ok);
                value(i,nodesIndOk) = 1;
            end
        end
        function valueIndex = addCurve(obj,varargin)
            %New version dev2.1
            %Adds a new curve appended to the end of the existing curves
            valueIndex=-1;
            if nargin==2
                [obj.errcode] = ENaddcurve(varargin{1},obj.libepanet);
                valueIndex = getCurveIndex(obj,varargin{1});
            elseif nargin==3
                [obj.errcode] = ENaddcurve(varargin{1},obj.libepanet);
                valueIndex = getCurveIndex(obj,varargin{1});
                setCurve(obj,valueIndex,varargin{2});
            end
        end
        function value = getCurveXY(obj,index)
            %New version dev2.1
            %Retrieves x,y point for a specific point number and curve
            tmplen=obj.getCurveLengths;
            value=zeros(tmplen(index),2);
            for i=1:tmplen(index)
                [obj.errcode, value(i,1), value(i,2)] = ENgetcurvevalue(index, i,obj.libepanet);
            end
        end
        function value = getHeadCurveIndex(obj)
            %New version dev2.1
            %Retrieves index of a head curve for specific link index
            v=1;
            for i=obj.getLinkPumpIndex
                [obj.errcode, value(v)] = ENgetheadcurveindex(i,obj.libepanet);
                v=v+1;
            end
        end
        function value = getLinkPumpTypeCode(obj)
            %New version dev2.1
            %Input:   index = index of pump in list of links
            %Output:  type = PumpType
            %Returns: error code                              
            %Purpose: retrieves type of a pump for specific link index
            v=1;
            if obj.getLinkPumpCount
                for i=obj.getLinkPumpIndex
                    [obj.errcode, value(v)] = ENgetpumptype(i,obj.libepanet);
                    v=v+1;
                end
            else
                value=[];
            end
        end
        function value = getLinkPumpType(obj)
            %New version dev2.1
            v = obj.getLinkPumpTypeCode;
            value = obj.TYPEPUMP(v+1);
        end
        function value = getPatternAveragePatternValue(obj)
            %New version dev2.1
            for i=obj.getPatternIndex
                [obj.errcode, value(i)] = ENgetaveragepatternvalue(i,obj.libepanet);
            end
        end
        function value = getENfunctionsImpemented(obj)
            [tline]=regexp(fileread([mfilename,'.m']), '\n', 'split');
            value={};
            u=1;
            for i=1:length(tline)
                if ~isempty(regexp(tline{i},'EN\w','match'))
                    enfunc = regexp(tline{i},'\w*EN\w*','match');
                    checkL = isstrprop(enfunc,'upper');
                    if length(checkL{1})>2
                        if strcmp(enfunc{1}(1:2),'EN') && ~checkL{1}(3) && ~strcmp(enfunc{1}(1:3),'EN_')
                            value(u) = enfunc;
                            u=u+1;
                        end
                    end
                end
            end
            value = unique(value)';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = getVersion(obj)
            % Retrieve the current EPANET libepanet
            [obj.errcode, value] = ENgetversion(obj.libepanet);
        end
        function value = getComputedHydraulicTimeSeries(obj,varargin)
            % Compute hydraulic simulation and retrieve all time-series
            obj.openHydraulicAnalysis;
            obj.initializeHydraulicAnalysis
            tstep=1;
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep+1;
            initnodematrix=zeros(totalsteps, obj.getNodeCount);
            initlinkmatrix=zeros(totalsteps, obj.getLinkCount);
            if size(varargin,2)==0
                varargin={'time','pressure','demand','head','tankvolume','flow','velocity','headloss','status','setting','energy'};
            else
                for i=1:length(varargin)
                    if isnumeric(varargin{i})
                        sensingnodes=i;
                    end
                end
            end
            if find(strcmpi(varargin,'time'))
                value.Time=zeros(totalsteps,1);
            end
            if find(strcmpi(varargin,'pressure'))
                value.Pressure=initnodematrix;
            end
            if find(strcmpi(varargin,'demand'))
                value.Demand=initnodematrix;
            end
            if find(strcmpi(varargin,'demandSensingNodes'))
                value.DemandSensingNodes=zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices=varargin{sensingnodes};
            end
            if find(strcmpi(varargin,'head'))
                value.Head=initnodematrix;
            end
            if find(strcmpi(varargin,'tankvolume'))
                value.TankVolume=initnodematrix;
            end
            if find(strcmpi(varargin,'flow'))
                value.Flow=initlinkmatrix;
            end
            if find(strcmpi(varargin,'velocity'))
                value.Velocity=initlinkmatrix;
            end
            if find(strcmpi(varargin,'headloss'))
                value.HeadLoss=initlinkmatrix;
            end
            if find(strcmpi(varargin,'status'))
                value.Status=initlinkmatrix;
            end
            if find(strcmpi(varargin,'setting'))
                value.Setting=initlinkmatrix;
            end
            if find(strcmpi(varargin,'energy'))
                value.Energy=initlinkmatrix;
            end
            k=1;
            while (tstep>0)
                t=obj.runHydraulicAnalysis;
                if find(strcmpi(varargin,'time'))
                    value.Time(k,:)=t;
                end
                if find(strcmpi(varargin,'pressure'))
                    value.Pressure(k,:)=obj.getNodePressure;
                end
                if find(strcmpi(varargin,'demand'))
                    value.Demand(k,:)=obj.getNodeActualDemand;
                end
                if find(strcmpi(varargin,'demandSensingNodes'))
                    value.DemandSensingNodes(k,:)=obj.getNodeActualDemandSensingNodes(varargin{sensingnodes});
                end
                if find(strcmpi(varargin,'head'))
                    value.Head(k,:)=obj.getNodeHydaulicHead;
                end
                if find(strcmpi(varargin,'tankvolume'))
                    value.TankVolume(k,:)=obj.getNodeTankVolume;
                end
                if find(strcmpi(varargin,'flow'))
                    value.Flow(k,:)=obj.getLinkFlows;
                end
                if find(strcmpi(varargin,'velocity'))
                    value.Velocity(k,:)=obj.getLinkVelocity;
                end
                if find(strcmpi(varargin,'headloss'))
                    value.HeadLoss(k,:)=obj.getLinkHeadloss;
                end
                if find(strcmpi(varargin,'status'))
                    value.Status(k,:)=obj.getLinkStatus;
                end
                if find(strcmpi(varargin,'setting'))
                    value.Setting(k,:)=obj.getLinkSettings;
                end
                if find(strcmpi(varargin,'energy'))
                    value.Energy(k,:)=obj.getLinkPumpEnergy;
                end
                tstep = obj.nextHydraulicAnalysisStep;
                k=k+1;
            end
            obj.closeHydraulicAnalysis;
        end
        function value = getComputedQualityTimeSeries(obj,varargin)
            % Compute Quality simulation and retrieve all or some time-series
            obj.openQualityAnalysis
            obj.initializeQualityAnalysis
            %tleft=obj.nextQualityAnalysisStep;
            tleft=1;
            totalsteps=0;%obj.getTimeSimulationDuration%;/obj.getTimeQualityStep;
            initnodematrix=zeros(totalsteps, obj.getNodeCount);
            if size(varargin,2)==0
                varargin={'time', 'quality', 'mass'};
            else
                for i=1:length(varargin)
                    if isnumeric(varargin{i})
                        sensingnodes=i;
                    end
                end
            end
            if find(strcmpi(varargin,'time'))
                value.Time=zeros(totalsteps,1);
            end
            if find(strcmpi(varargin,'quality'))
                value.Quality=initnodematrix;
            end
            if find(strcmpi(varargin,'qualitySensingNodes'))
                value.QualitySensingNodes=zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices=varargin{sensingnodes};
            end
            if find(strcmpi(varargin,'demandSensingNodes'))
                value.DemandSensingNodes=zeros(totalsteps, length(varargin{sensingnodes}));
                value.SensingNodesIndices=varargin{sensingnodes};
            end
            if find(strcmpi(varargin,'mass'))
                value.MassFlowRate=initnodematrix;
            end
            if find(strcmpi(varargin,'demand'))
                value.Demand=initnodematrix;
            end
            k=1;t=1;
            while (tleft>0)||(t<obj.getTimeSimulationDuration)
                t=obj.runQualityAnalysis;
                if find(strcmpi(varargin,'time'))
                    value.Time(k,:)=t;
                end
                if find(strcmpi(varargin,'quality'))
                    value.Quality(k,:)=obj.getNodeActualQuality;
                end
                if find(strcmpi(varargin,'mass'))
                    value.MassFlowRate(k,:)=obj.getNodeMassFlowRate;
                end
                if find(strcmpi(varargin,'demand'))
                    value.Demand(k,:)=obj.getNodeActualDemand;
                end
                if find(strcmpi(varargin,'qualitySensingNodes'))
                    value.QualitySensingNodes(k,:)=obj.getNodeActualQualitySensingNodes(varargin{2});
                end
                if find(strcmpi(varargin,'demandSensingNodes'))
                    value.DemandSensingNodes(k,:)=obj.getNodeActualDemandSensingNodes(varargin{sensingnodes});
                end
                tleft = obj.nextQualityAnalysisStep;
                k=k+1;
                if t==obj.getTimeSimulationDuration
                    t=obj.getTimeSimulationDuration+1;
                end
            end
            obj.closeQualityAnalysis;
        end
        function solveCompleteHydraulics(obj)
            [obj.errcode] = ENsolveH(obj.libepanet);
        end
        function solveCompleteQuality(obj)
            [obj.errcode] = ENsolveQ(obj.libepanet);
        end
        function valueIndex = addPattern(obj,varargin)
            valueIndex=-1;
            if nargin==2
                [obj.errcode] = ENaddpattern(varargin{1},obj.libepanet);
                valueIndex = getPatternIndex(obj,varargin{1});
            elseif nargin==3
                [obj.errcode] = ENaddpattern(varargin{1},obj.libepanet);
                valueIndex = getPatternIndex(obj,varargin{1});
                setPattern(obj,valueIndex,varargin{2});
            end
        end
        function setControl(obj,controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel)
            % Example: d.setControl(1,1,13,1,11,150)
            if controlRuleIndex<=obj.getControlRulesCount
                [obj.errcode] = ENsetcontrol(controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel,obj.libepanet);
            else
                disp('New rules cannot be added in this libepanet')
            end
        end
        function setLinkDiameter(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 0, value(i),obj.libepanet);
            end
        end
        function setLinkLength(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 1, value(i),obj.libepanet);
            end
        end
        function setLinkRoughnessCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 2, value(i),obj.libepanet);
            end
        end
        function setLinkMinorLossCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 3, value(i),obj.libepanet);
            end
        end
        function setLinkInitialStatus(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 4, value(i),obj.libepanet);
                e=obj.getError(obj.errcode); % Cannot set status for a check valve
                if ~isempty(e)
                   disp(e);
                end
            end
        end
        function setLinkInitialSetting(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 5, value(i),obj.libepanet);
            end
        end
        function setLinkBulkReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 6, value(i),obj.libepanet);
            end
        end
        function setLinkWallReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 7, value(i),obj.libepanet);
            end
        end
        function setLinkStatus(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 11, value(i),obj.libepanet);
            end
        end
        function setLinkSettings(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 12, value(i),obj.libepanet);
            end
        end
        function setNodeElevations(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 0, value(i),obj.libepanet);
            end
        end
        function setNodeBaseDemands(obj, value)
            %New version dev2.1
            chckfunctions=libfunctions(obj.libepanet);
            if sum(strcmp(chckfunctions,'ENsetbasedemand'))
                NodeNumDemandC=obj.getNodeNumDemandCategories;
                for u=1:obj.getNodeJunctionCount
                    [obj.errcode] = ENsetbasedemand(u, NodeNumDemandC(u), value{NodeNumDemandC(u)}(u),obj.libepanet);
                end
            else %version epanet20012
                for i=1:length(value)
                    [obj.errcode] = ENsetnodevalue(i, 1, value(i),obj.libepanet);
                end
            end
        end
        function setNodeCoordinates(obj, value)
            for i=1:length(value)
                x=value{1}(i);
                y=value{2}(i);
                [obj.errcode] = ENsetcoord(i,x,y,obj.libepanet);
            end
        end
        function setNodeDemandPatternIndex(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 2, value(i),obj.libepanet);
            end
        end
        function setNodeEmitterCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 3, value(i),obj.libepanet);
            end
        end
        function setNodeInitialQuality(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 4, value(i),obj.libepanet);
            end
        end
        function setNodeTankInitialLevel(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 8, value(i),obj.libepanet);
            end
        end
        function setNodeTankMixingModelType(obj, value)
            for i=obj.getNodeTankIndex
                code=strfind(strcmpi(value(i),obj.TYPEMIXMODEL),1)-1;
                [obj.errcode] = ENsetnodevalue(i, 15, code,obj.libepanet);
            end
        end
        function setNodeTankDiameter(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 17, value(i),obj.libepanet);
            end
        end
        function setNodeTankMinimumWaterLevel(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 20, value(i),obj.libepanet);
            end
        end
        function setNodeTankMinimumWaterVolume(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 18, value(i),obj.libepanet);
            end
        end
        function setNodeTankMaximumWaterLevel(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 21, value(i),obj.libepanet);
            end
        end
        function setNodeTankMinimumFraction(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 22, value(i),obj.libepanet);
            end
        end
        function setNodeTankBulkReactionCoeff(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 23, value(i),obj.libepanet);
            end
        end
        function setNodeSourceQuality(obj, value)
            value(find(isnan(value)))='';
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 5, value(i),obj.libepanet);
            end
        end
        function setNodeSourcePatternIndex(obj, value)
            value(find(isnan(value)))='';
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 6, value(i),obj.libepanet);
            end
        end
        function setNodeSourceType(obj, index, value)
            value=find(strcmpi(obj.TYPESOURCE,value)==1)-1;
            [obj.errcode] = ENsetnodevalue(index, 7, value,obj.libepanet);
        end
        function setOptionsMaxTrials(obj,value)
            [obj.errcode] = ENsetoption(0,value,obj.libepanet);
        end
        function setOptionsAccuracyValue(obj,value)
            [obj.errcode] = ENsetoption(1,value,obj.libepanet);
        end
        function setOptionsQualityTolerance(obj,value)
            [obj.errcode] = ENsetoption(2,value,obj.libepanet);
        end
        function setOptionsEmitterExponent(obj,value)
            [obj.errcode] = ENsetoption(3,value,obj.libepanet);
        end
        function setOptionsPatternDemandMultiplier(obj,value)
            [obj.errcode] = ENsetoption(4,value,obj.libepanet);
        end
        function setTimeSimulationDuration(obj,value)
            [obj.errcode] = ENsettimeparam(0,value,obj.libepanet);
        end
        function setTimeHydraulicStep(obj,value)
            % Hstep = min(Pstep,Hstep)
            % Hstep = min(Rstep,Hstep)
            % Hstep = min(Qstep,Hstep)
            [obj.errcode] = ENsettimeparam(1,value,obj.libepanet);
        end
        function setTimeQualityStep(obj,value)
            % Qstep = min(Qstep,Hstep)
            [obj.errcode] = ENsettimeparam(2,value,obj.libepanet);
        end
        function setTimePatternStep(obj,value)
            [obj.errcode] = ENsettimeparam(3,value,obj.libepanet);
        end
        function setTimePatternStart(obj,value)
            [obj.errcode] = ENsettimeparam(4,value,obj.libepanet);
        end
        function setTimeReportingStep(obj,value)
            [obj.errcode] = ENsettimeparam(5,value,obj.libepanet);
        end
        function setTimeReportingStart(obj,value)
            [obj.errcode] = ENsettimeparam(6,value,obj.libepanet);
        end
        function setTimeRuleControlStep(obj,value)
            [obj.errcode] = ENsettimeparam(7,value,obj.libepanet);
        end
        function setTimeStatisticsType(obj,value)
            %'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
            tmpindex=find(strcmpi(obj.TYPESTATS,value)==1)-1;
            [obj.errcode] = ENsettimeparam(8,tmpindex,obj.libepanet);
        end
        function setTimeHTime(obj,value)
            % epanet20100
            [obj.errcode] = ENsettimeparam(11,value,obj.libepanet);
        end
        function setTimeHaltFlag(obj,value)
            % epanet20100
            [obj.errcode] = ENsettimeparam(12,value,obj.libepanet);
        end
%         function value = setTimeReportingPeriods(obj)
%             [obj.errcode, value] = ENgettimeparam(9,obj.libepanet);
%         function value = setTimeStartTime(obj)
%             [obj.errcode, value] = ENgettimeparam(10,obj.libepanet);
%         function value = setTimeNextEvent(obj)
%             [obj.errcode, value] = ENgettimeparam(13,obj.libepanet);
%         end
        function setPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errcode] = ENsetpattern(index, patternVector, nfactors,obj.libepanet);
        end
        function setPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errcode] = ENsetpattern(i, patternMatrix(i,:), nfactors,obj.libepanet);
            end
        end
        function setPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.errcode] = ENsetpatternvalue(index, patternTimeStep, patternFactor,obj.libepanet);
        end
        function setQualityType(obj,varargin)
            qualcode=0;chemname='';chemunits='';tracenode='';
            if find(strcmpi(varargin,'none')==1)
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.libepanet);
            elseif find(strcmpi(varargin,'age')==1)
                qualcode=2;
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.libepanet);
            elseif find(strcmpi(varargin,'chem')==1)
                qualcode=1;
                chemname=varargin{1};
                if nargin<3
                    chemunits='mg/L';
                else
                    chemunits=varargin{2};
                end
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.libepanet);
            elseif find(strcmpi(varargin,'trace')==1)
                qualcode=3;
                tracenode=varargin{2};
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.libepanet);
            else
                qualcode=1;
                chemname=varargin{1};
                if nargin<3
                    chemunits='mg/L';
                else
                    chemunits=varargin{2};
                end
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.libepanet);                
            end
        end
        function setReportFormatReset(obj)
            [obj.errcode]=ENresetreport(obj.libepanet);
        end
        function setReportStatus(obj,value)
            %'yes','no','full'
            statuslevel=find(strcmpi(obj.TYPEREPORT,value)==1)-1;
            [obj.errcode] = ENsetstatusreport(statuslevel,obj.libepanet);
        end
        function setReport(obj,value)
            [obj.errcode] = ENsetreport(value,obj.libepanet);
        end
        function closeNetwork(obj)
            [obj.errcode] = ENclose(obj.libepanet);
        end
        function closeHydraulicAnalysis(obj)
            [obj.errcode] = ENcloseH(obj.libepanet);
        end
        function closeQualityAnalysis(obj)
            [obj.errcode] = ENcloseQ(obj.libepanet);
        end
        function saveHydraulicFile(obj,hydname)
            [obj.errcode]=ENsavehydfile(hydname,obj.libepanet);
        end
        function useHydraulicFile(obj,hydname)
            [obj.errcode]=ENusehydfile(hydname,obj.libepanet);
        end
        function initializeHydraulicAnalysis(obj)
            [obj.errcode] = ENinitH(1,obj.libepanet);
        end
        function initializeQualityAnalysis(obj)
            [obj.errcode] = ENinitQ(1,obj.libepanet);
        end
        function tstep = nextHydraulicAnalysisStep(obj)
            [obj.errcode, tstep] = ENnextH(obj.libepanet);
        end
        function tstep = nextQualityAnalysisStep(obj)
            [obj.errcode, tstep] = ENnextQ(obj.libepanet);
        end
        function openHydraulicAnalysis(obj)
            [obj.errcode] = ENopenH(obj.libepanet);
        end
        function openQualityAnalysis(obj)
            [obj.errcode] = ENopenQ(obj.libepanet);
        end
        function tstep = runHydraulicAnalysis(obj)
            [obj.errcode, tstep] = ENrunH(obj.libepanet);
        end
        function tstep = runQualityAnalysis(obj)
            [obj.errcode, tstep] = ENrunQ(obj.libepanet);
        end
        function saveHydraulicsOutputReportingFile(obj)
            [obj.errcode] = ENsaveH(obj.libepanet);
        end
        function tleft=stepQualityAnalysisTimeLeft(obj)
            [obj.errcode, tleft] = ENstepQ(obj.libepanet);
        end
        function errcode = saveInputFile(obj,inpname,varargin)
            if strcmp(inpname,obj.Bintempfile) && nargin<3 %&& ~isempty(varargin)
                addSectionCoordinates=obj.getBinCoordinatesSection;
                addSectionRules = obj.getBinRulesSection;
                [obj.errcode] = ENsaveinpfile(inpname,obj.libepanet);
                % Open epanet input file
                [~,info] = obj.readInpFile;
                endSectionIndex=find(~cellfun(@isempty,regexp(info,'END','match')));
                info(endSectionIndex)='';
                f1=fopen([obj.pathfile,obj.Bintempfile],'w');
                fprintf(f1, '%s\n', info{:});
                if ~isempty(addSectionRules)
                    fprintf(f1, '%s\n', addSectionRules{:});
                end
                if ~isempty(addSectionCoordinates)
                    fprintf(f1, '%s\n', addSectionCoordinates{:});
                end
                fclose(f1);return;
            end
            [errcode] = ENsaveinpfile(inpname,obj.libepanet);
            % The code below is because of a bug in EPANET 2.00.12
            % When saving using ENsaveinpfile, it does not save the type of the curves.
            obj.remAddBinCurvesID(inpname);
        end
        function writeLineInReportFile(obj, line)
            [obj.errcode] = ENwriteline (line,obj.libepanet);
        end
        function writeReport(obj)
            %Writes a formatted text report on simulation results to the Report file
            [obj.errcode]=ENreport(obj.libepanet);
        end
        function unload(obj)
            ENclose(obj.libepanet);
            ENMatlabCleanup(obj.libepanet);
            if exist([obj.Bintempfile(1:end-4),'.bin'])==2
                delete([obj.Bintempfile(1:end-4),'.bin']);
            end
            delete(obj.Bintempfile);
            if exist([obj.Bintempfile(1:end-4),'.txt'])==2
                delete([obj.Bintempfile(1:end-4),'.txt']);
            end
            [p,f]=fileparts(obj.inputfile);
            if exist([p,'/',f,'.txt'])==2
                delete([p,'/',f,'.txt']);
            end
            if exist(obj.MsxTempFile)==2
                delete(obj.MsxTempFile);
            end
            disp('EPANET Class is unloaded')
        end
        function msx(obj,msxname,varargin)
            if isempty(varargin)
                MSXMatlabSetup(obj,msxname);
            else
                MSXMatlabSetup(obj,msxname,varargin);
            end
        end
        function value = getMsxEquationsTerms(obj)
            [value,~,~] = getEquations(obj.MsxFile);
        end
        function value = getMsxEquationsPipes(obj)
            [~,value,~] = getEquations(obj.MsxFile);
        end
        function value = getMsxEquationsTanks(obj)
            [~,~,value] = getEquations(obj.MsxFile);
        end
        function value = getMsxTimeStep(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.timestep;
        end
        function value = getMsxSolver(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.solver;
        end
        function value = getMsxAreaUnits(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.areaunits;
        end        
        function value = getMsxRateUnits(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.rateunits;
        end 
        function value = getMsxRtol(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.rtol;
        end 
        function value = getMsxAtol(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.atol;
        end
        function value = getMsxCoupling(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.coupling;
        end
        function value = getMsxCompiler(obj)
            [value] = getMsxOptions(obj.MsxFile);
            value = value.compiler;
        end
        function value = getMsxSpeciesCount(obj)
            % Species, Constants, Parameters, Patterns
            [obj.errcode, value] = MSXgetcount(3,obj.Msxlibepanet);
        end
        function value = getMsxConstantsCount(obj)
            [obj.errcode, value] = MSXgetcount(6,obj.Msxlibepanet);
        end
        function value = getMsxParametersCount(obj)
            [obj.errcode, value] = MSXgetcount(5,obj.Msxlibepanet);
        end
        function value = getMsxPatternsCount(obj)
            [obj.errcode, value] = MSXgetcount(7,obj.Msxlibepanet);
        end
        function value = getMsxSpeciesNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode, len] = MSXgetIDlen(3,i,obj.Msxlibepanet);
                    [obj.errcode, value{i}]=MSXgetID(3,i,len,obj.Msxlibepanet);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(3,i,obj.Msxlibepanet);
                    [obj.errcode, value{k}]=MSXgetID(3,i,len,obj.Msxlibepanet);
                    k=k+1;
                end
            end
        end
        function value = getMsxSpeciesType(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,value{i},~,~,~] = MSXgetspecies(i,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesUnits(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,value{i},~,~] = MSXgetspecies(i,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesATOL(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,value,~] = MSXgetspecies(i,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesRTOL(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,~,value] = MSXgetspecies(i,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMsxSpeciesCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = MSXgetindex(3,varargin{1}{j},obj.Msxlibepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(obj.Msxlibepanet,3,varargin{1});
            end
        end
        function value = getMsxConstantsNameID(obj)
            if obj.getMsxConstantsCount
                for i=1:obj.getMsxConstantsCount
                    [obj.errcode, len] = MSXgetIDlen(6,i,obj.Msxlibepanet);
                    [obj.errcode, value{i}] = MSXgetID(6,i,len,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxConstantsValue(obj)
            if obj.getMsxConstantsCount
                for i=1:obj.getMsxConstantsCount
                    [obj.errcode, value(i)] = MSXgetconstant(i,obj.Msxlibepanet);
                end
            else
                value=-1;
            end
        end
        function value = getMsxConstantsIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMsxConstantsCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = MSXgetindex(6,varargin{1}{j},obj.Msxlibepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(obj.Msxlibepanet,6,varargin{1});
            end
        end
        function value = getMsxParametersNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxParametersCount
                    [obj.errcode, len] = MSXgetIDlen(5,i,obj.Msxlibepanet);
                    [obj.errcode, value{i}]=MSXgetID(5,i,len,obj.Msxlibepanet);
                end
                if ~obj.getMsxParametersCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(5,i,obj.Msxlibepanet);
                    [obj.errcode, value{k}]=MSXgetID(5,i,len,obj.Msxlibepanet);
                    k=k+1;
                end
            end
        end
        function value = getMsxParametersIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMsxParametersCount;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = MSXgetindex(5,varargin{1}{j},obj.Msxlibepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(5,varargin{1},obj.Msxlibepanet);
            end
        end
        function value = getMsxParametersTanksValue(obj)
            value={};
            if ~obj.getMsxParametersCount, value=0;return;end
            if ~length(obj.NodeTankIndex), value=0;return;end
            for i=1:length(obj.NodeTankIndex)
                for j=1:obj.MsxParametersCount%
                    [obj.errcode, value{obj.NodeTankIndex(i)}(j)] = MSXgetparameter(0,obj.NodeTankIndex(i),j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxParametersPipesValue(obj)
            if ~obj.getMsxParametersCount
                value=0;return;
            end
            for i=1:obj.getLinkPipeCount
                for j=1:obj.getMsxParametersCount
                    [obj.errcode, value{i}(j)] = MSXgetparameter(1,i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxPatternsNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxPatternsCount
                    [obj.errcode, len] = MSXgetIDlen(7,i,obj.Msxlibepanet);
                    [obj.errcode, value{i}]=MSXgetID(7,i,len,obj.Msxlibepanet);
                end
                if ~obj.getMsxPatternsCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(7,i,obj.Msxlibepanet);
                    [obj.errcode, value{k}]=MSXgetID(7,i,len,obj.Msxlibepanet);
                    k=k+1;
                end
            end
        end
        function value = getMsxPatternsIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMsxPatternsCount;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value] = MSXgetindex(obj.Msxlibepanet,7,varargin{1});
                    if obj.errcode
                        value{k}=0;
                    end
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(obj.Msxlibepanet,7,varargin{1});
                if obj.errcode
                    value=0;
                end
            end
        end
        function value = getMsxPatternsLengths(obj,varargin)
            if isempty(varargin)
                if obj.getMsxPatternsCount
                    for i=obj.getMsxPatternsIndex
                        [obj.errcode, value(i)]=MSXgetpatternlen(i,obj.Msxlibepanet);
                    end
                else
                    value=-1;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = MSXgetpatternlen(obj.getMsxPatternsIndex(varargin{1}{j}),obj.Msxlibepanet);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetpatternlen(obj.getMsxPatternsIndex(varargin{1}),obj.Msxlibepanet);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errcode, value(k)]=MSXgetpatternlen(i,obj.Msxlibepanet);
                    k=k+1;
                end
            end
        end
        function value = getMsxNodeInitqualValue(obj)
            if obj.getMsxSpeciesCount==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, value{i}(j)] = MSXgetinitqual(0,i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxLinkInitqualValue(obj)
            if obj.getMsxSpeciesCount==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.getLinkCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, value{i}(j)] = MSXgetinitqual(1,i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxSources(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, obj.MsxSourceType{i}{j},obj.MsxSourceLevel{i}(j),obj.MsxSourcePatternIndex{i}(j)] = MSXgetsource(i,j,obj.Msxlibepanet);
                    obj.MsxSourceTypeCode{i}(j)=find(strcmp(obj.TYPESOURCEMSX,obj.MsxSourceType{i}{j}))-2;
                end
            end
            SnodeID=obj.getMsxSourceNodeNameID;
            value={obj.MsxSourceType,obj.MsxSourceTypeCode,obj.MsxSourceLevel,obj.MsxSourcePatternIndex,SnodeID};
            warning off; value.MsxSourceType=obj.MsxSourceType; warning on;
            value.MsxSourceTypeCode=obj.MsxSourceTypeCode;
            value.MsxSourceLevel=obj.MsxSourceLevel;
            value.MsxSourcePatternIndex=obj.MsxSourcePatternIndex;
            value.MsxSourceNodeNameID=SnodeID;
        end
        function value = getMsxSourceNodeNameID(obj)
            value = obj.getNodeNameID;
        end
        function value = getMsxSourceType(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, value{i}{j},~,~] = MSXgetsource(i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxSourceLevel(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,value{i}(j),~] = MSXgetsource(i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxSourcePatternIndex(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,value{i}(j)] = MSXgetsource(i,j,obj.Msxlibepanet);
                end
            end
        end
        function value = getMsxPattern(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getMsxPatternsLengths);
            value=nan(obj.getMsxPatternsCount,tmpmaxlen);
            for i=1:obj.getMsxPatternsCount
                tmplength=obj.getMsxPatternsLengths(i);
                for j=1:tmplength
                    [obj.errcode, value(i,j)] = MSXgetpatternvalue(i,j,obj.Msxlibepanet);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                
            end
        end
        function value = getMsxPatternValue(obj,patternIndex, patternStep) %Mass flow rate per minute of a chemical source
            [obj.errcode, value] = MSXgetpatternvalue(patternIndex, patternStep, obj.Msxlibepanet);
        end
        function value = getMsxSpeciesConcentration(obj, type, index, species)
            [obj.errcode, value] = MSXgetqual(type, index, species, obj.Msxlibepanet);
        end
        function value = getMsxComputedQualityNode(obj,varargin)
            if obj.getMsxSpeciesCount==0
                value=0;
                return;
            end
            if ~isempty(varargin)
                if length(varargin)==1
                    ss=varargin{1};%index node
                    uu=1:obj.getMsxSpeciesCount;
                elseif length(varargin)==2
                    ss=varargin{1};%index node
                    uu=varargin{2};%index species
                end
            else
                ss=1:obj.getNodeCount;%index node
                uu=1:obj.getMsxSpeciesCount;
            end
            % Obtain a hydraulic solution
            obj.MsxSolveCompleteHydraulics(obj.Msxlibepanet);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.MsxInitializeQualityAnalysis(0);
            % Retrieve species concentration at node
            k=2; tleft=1;t=0;i=1;nl=1;
            value.Time(k,:)=0;
            timeSmle=obj.getTimeSimulationDuration;%bug at time
            while(tleft>0 && obj.errcode==0 && timeSmle~=t)
                [t, tleft]=obj.MsxStepQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<3600 || t==3600
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMsxNodeInitqualValue{ss}(j);
                        end
                    else
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMsxSpeciesConcentration(0, ss, j);%node code0
                        end
                    end
                else
                    if t<3600 || t==3600
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMsxNodeInitqualValue{(nl)}(j);
                            end
                        end
                    else
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMsxSpeciesConcentration(0, (nl), j);%node code0
                            end
                        end
                    end
                    nl=nl+1;
                end
                value.Time(k,:)=t;
                k=k+1;
            end
        end
        function value = getMsxComputedQualityLink(obj,varargin)
            if obj.getMsxSpeciesCount==0
                value=0;
                return;
            end
            if ~isempty(varargin)
                if length(varargin)==1
                    ss=varargin{1};%index link
                    uu=1:obj.getMsxSpeciesCount;
                elseif length(varargin)==2
                    ss=varargin{1};%index link
                    uu=varargin{2};%index species
                end
            else
                ss=1:obj.getLinkCount;%index link
                uu=1:obj.getMsxSpeciesCount;
            end
            % Obtain a hydraulic solution
            obj.MsxSolveCompleteHydraulics(obj.Msxlibepanet);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.MsxInitializeQualityAnalysis(0);
            
            % Retrieve species concentration at node
            k=2;tleft=1;i=1;
            value.Time(k,:)=0;
            while(tleft>0 && obj.errcode==0)
                [t, tleft]=obj.MsxStepQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<3600 || t==3600
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMsxLinkInitqualValue{ss}(j);
                        end
                    else
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMsxSpeciesConcentration(1, ss, j);
                        end
                    end
                else
                    if t<3600 || t==3600
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMsxLinkInitqualValue{(nl)}(j);
                            end
                        end
                    else
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMsxSpeciesConcentration(1, (nl), j);%link code1
                            end
                        end
                    end
                    nl=nl+1;
                end
                value.Time(k,:)=t;
                k=k+1;
            end
        end
        function MsxPlotConcentrationSpeciesOfNodes(obj,varargin)
            s=obj.getMsxComputedQualityNode(varargin{1},varargin{2});
            nodesID=obj.getNodeNameID;
            SpeciesNameID=obj.getMsxSpeciesNameID;
            for l=varargin{1}
                nodeID=nodesID(l);
                figure('Name',['NODE ',char(nodeID)]);
                for i=varargin{2}
                    specie(:,i)=s.Quality{i,1};
                    time(:,i)=s.Time;
                end
                plot(time,specie);
                title(['NODE ',char(nodeID)]);
                ylabel('Quantity');
                xlabel('Time(s)');
                legend(SpeciesNameID(varargin{2}));
            end
        end
        function MsxPlotConcentrationSpeciesOfLinks(obj,varargin)
            s=obj.getMsxComputedQualityLink(varargin{1},varargin{2});
            linksID=obj.getLinkNameID;
            SpeciesNameID=obj.getMsxSpeciesNameID;
            for l=varargin{1}
                linkID=linksID(l);
                figure('Name',['LINK ',char(linkID)]);
                for i=varargin{2}
                    specie(:,i)=s.Quality{i,1};
                    time(:,i)=s.Time;
                end
                plot(time,specie);
                title(['LINK ',char(linkID)]);
                ylabel('Quantity');
                xlabel('Time(s)');
                legend(SpeciesNameID(varargin{2}));
            end
        end
        function value = getMsxError(obj,errcode)
            [obj.errcode, value] = MSXgeterror(errcode,obj.Msxlibepanet);
        end
        function MsxSolveCompleteHydraulics(obj,varargin)
            [obj.errcode] = MSXsolveH(obj.Msxlibepanet);
        end
        function MsxSolveCompleteQuality(obj,varargin)
            [obj.errcode] = MSXsolveQ(obj.Msxlibepanet);
        end
        function MsxWriteReport(obj,varargin)
            [obj.errcode]=MSXreport(obj.Msxlibepanet);
        end
        function index = MsxAddPattern(obj,varargin)
            index=-1;
            if nargin==2
                [obj.errcode] = MSXaddpattern(varargin{1},obj.Msxlibepanet);
                [obj.errcode, index] = MSXgetindex(obj.Msxlibepanet,7,varargin{1});
            elseif nargin==3
                [obj.errcode] = MSXaddpattern(varargin{1},obj.Msxlibepanet);
                [obj.errcode, index] = MSXgetindex(obj.Msxlibepanet,7,varargin{1});
                setMsxPattern(obj,index,varargin{2});
            end
        end
        function setMsxSources(obj, node, species, type, level, pat)
            MSXsetsource(node, species, type, level, pat, obj.Msxlibepanet);
        end
        function setMsxConstantsValue(obj, value)
            for i=1:length(value)
                [obj.errcode] = MSXsetconstant(i, value(i), obj.Msxlibepanet);
            end
        end
        function setMsxParametersTanksValue(obj, NodeTankIndex, paramindex, value)
            if ~sum(NodeTankIndex==obj.NodeTankIndex)
                fprintf('>> Invalid Tank Index <<\n');obj.NodeTankIndex
                return;
            end
            [obj.errcode] = MSXsetparameter(0, NodeTankIndex, paramindex, value, obj.Msxlibepanet);
        end
        function setMsxParametersPipesValue(obj, pipeIndex, value)
            for i=1:length(value)
                [obj.errcode] = MSXsetparameter(1, pipeIndex, i, value(i), obj.Msxlibepanet);
            end
        end
        function setMsxNodeInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.errcode] = MSXsetinitqual(0, i, j, value{i}(j), obj.Msxlibepanet);
                end
            end
        end
        function setMsxLinkInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.errcode] = MSXsetinitqual(1, i, j, value{i}(j), obj.Msxlibepanet);
                end
            end
        end
        function setMsxPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errcode] = MSXsetpattern(index, patternVector, nfactors, obj.Msxlibepanet);
        end
        function setMsxPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errcode] = MSXsetpattern(i, patternMatrix(i,:), nfactors, obj.Msxlibepanet);
            end
        end
        function setMsxPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.errcode] = MSXsetpatternvalue(index, patternTimeStep, patternFactor, obj.Msxlibepanet);
        end
        function setMsxTimeStep(obj,timestep)
            setMsxOptions(obj,'timestep',timestep);
        end
        function setMsxAreaUnitsFT2(obj)
            setMsxOptions(obj,'areaunits','FT2');
        end
        function setMsxAreaUnitsM2(obj)
            setMsxOptions(obj,'areaunits','M2');
        end
        function setMsxAreaUnitsCM2(obj)
            setMsxOptions(obj,'areaunits','CM2');
        end
        function setMsxRateUnitsSEC(obj)
            setMsxOptions(obj,'rateunits','SEC');
        end
        function setMsxRateUnitsMIN(obj)
            setMsxOptions(obj,'rateunits','MIN');
        end
        function setMsxRateUnitsHR(obj)
            setMsxOptions(obj,'rateunits','HR');
        end        
        function setMsxRateUnitsDAY(obj)
            setMsxOptions(obj,'rateunits','DAY');
        end
        function setMsxSolverEUL(obj)
            setMsxOptions(obj,'solver','EUL');
        end
        function setMsxSolverRK5(obj)
            setMsxOptions(obj,'solver','RK5');
        end        
        function setMsxSolverROS2(obj)
            setMsxOptions(obj,'solver','ROS2');
        end
        function setMsxCouplingFULL(obj)
            setMsxOptions(obj,'coupling','FULL');
        end
        function setMsxCouplingNONE(obj)
            setMsxOptions(obj,'coupling','NONE');
        end
        function setMsxCompilerNONE(obj)
            setMsxOptions(obj,'compiler','NONE');
        end
        function setMsxCompilerVC(obj)
            setMsxOptions(obj,'compiler','VC');
        end
        function setMsxCompilerGC(obj)
            setMsxOptions(obj,'compiler','GC');
        end
        function setMsxAtol(obj,atol)
            setMsxOptions(obj,'atol',atol);
        end
        function setMsxRtol(obj,rtol)
            setMsxOptions(obj,'rtol',rtol);
        end        
        function MsxSaveQualityFile(obj,outfname)
            [obj.errcode]=MSXsaveoutfile(outfname,obj.Msxlibepanet);
        end
        function MsxUseHydraulicFile(obj,hydname)
            [obj.errcode]=MSXusehydfile(hydname,obj.Msxlibepanet);
        end
        function MsxInitializeQualityAnalysis(obj,flag)
            [obj.errcode] = MSXinit(flag,obj.Msxlibepanet);
        end
        function [t, tleft]=MsxStepQualityAnalysisTimeLeft(obj)
            [obj.errcode, t, tleft] = MSXstep(obj.Msxlibepanet);
        end
        function MsxSaveFile(obj,msxname)
            [obj.errcode] = MSXsavemsxfile(msxname,obj.Msxlibepanet);
        end
        function MsxUnload(obj)
            MSXclose(obj);
            MSXMatlabCleanup(obj);
        end
        function obj = BinUpdateClass(obj)
            sect=0;i=1;t=1;q=1;
            typecode=0;x=1;b=1;d=1;
            if obj.Bin
                obj.saveInputFile([obj.pathfile,obj.Bintempfile]);
            end
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[JUNC')
                        sect=1;k=1;
                        obj.BinNodeJunctionNameID={};
                        obj.BinNodeJunctionIndex=[];
                        obj.BinNodeJunctionElevation=[];
                        obj.BinNodeType={};
                        obj.BinNodeJunctionsBaseDemands=[];
                        obj.BinNodeJunctionsBaseDemandsID={};
                        obj.BinNodeDemandPatternNameID={};
                        continue;
                        % [RESERVOIRS] section
                    elseif strcmpi(tok(1:5),'[RESE')
                        sect=2;r=1;
                        obj.BinNodeReservoirNameID={};
                        obj.BinNodeReservoirIndex=[];
                        obj.BinNodeReservoirElevation=[];
                        obj.BinNodeResDemandPatternNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5),'[TANK')
                        sect=3;p=1;
                        obj.BinNodeTankNameID={};
                        obj.BinNodeTankIndex=[];
                        obj.BinNodeTankElevation=[];
                        obj.BinNodeTankInitLevel=[]; 
                        obj.BinNodeTankMinLevel=[];  
                        obj.BinNodeTankMaxLevel=[]; 
                        obj.BinNodeTankDiameter=[];  
                        obj.BinNodeTankMinVol=[];
                        continue;
                        % [PIPES] section
                    elseif strcmpi(tok(1:5),'[PIPE')
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
                        obj.BinNodeInitialQuality=zeros(1,obj.BinNodeCount);  
                        continue;
                        % [PUMPS] section
                    elseif strcmpi(tok(1:5),'[PUMP')
                        sect=5;
                        obj.BinLinkPumpNameID={};
                        obj.BinLinkPumpIndex=[];
                        obj.BinLinkPumpPatterns={};
                        obj.BinLinkPumpPower=[];
                        obj.BinLinkPumpNameIDPower={};
                        continue;
                        % [VALVES] section
                    elseif strcmpi(tok(1:5),'[VALV')
                        sect=6;
                        obj.BinLinkValveNameID={};
                        obj.BinLinkValveIndex=[];
                        obj.BinLinkValveDiameters=[];
                        obj.BinLinkValveType={};
                        obj.BinLinkValveSetting=[];
                        obj.BinLinkValveMinorLoss=[];                        
                        continue;
                        % [CURVES] section
                    elseif strcmpi(tok(1:5),'[CURV')
                        sect=7;
                        obj.BinCurveAllLines={};
                        obj.BinCurveXvalue=[];
                        obj.BinCurveYvalue=[];
                        obj.BinCurveAllLines={};
                        continue;
                        % [DEMANDS] section
                    elseif strcmpi(tok(1:5),'[DEMA') %&& max(obj.BinNodeJunctionsBaseDemands)==0
                        sect=8;d=1;pd=1;pk=1;                       
                        continue;
                        % [STATUS] section
                    elseif strcmpi(tok(1:5),'[STAT')
                        sect=9;d=1;
                        obj.BinLinkInitialStatus={};
                        obj.BinLinkInitialStatusNameID={};
                        obj.BincountStatuslines=0;
                        continue;
                        % [PATTERNS]
                    elseif strcmpi(tok(1:5),'[PATT')
                        sect=10;
                        obj.BinPatternLengths=[];
                        obj.BinPatternNameID={};
                        obj.BinPatternValue={};
                        obj.BincountPatternlines=0;d=1;h=1;
                        continue;                    
                        % [CONTROLS] section
                    elseif strcmpi(tok(1:5),'[CONT')
                        sect=11;d=1;
                        obj.BinControlsInfo={};
                        obj.BinControlLinksID={};
                        obj.BinControlNodesID={};
                        obj.BinPatternCount=length(obj.BinPatternNameID);
                        tmpmaxlen=max(obj.BinPatternLengths);
                        obj.BinPatternMatrix=nan(obj.BinPatternCount,tmpmaxlen);
                        for i=1:obj.BinPatternCount
                            tmplength=obj.BinPatternLengths(i);
                            for j=1:tmplength
                                obj.BinPatternMatrix(i,j)=single(obj.BinPatternValue{i}(j));
                            end
                            if tmplength<tmpmaxlen
                                for j=(tmplength+1):tmpmaxlen
                                    obj.BinPatternMatrix(i,j)=obj.BinPatternMatrix(i,j-tmplength);
                                end
                            end
                        end
                        continue;
                        % [RULES] section
                    elseif strcmpi(tok(1:5),'[RULE')
                        sect=20;d=1;dd=0; 
                        obj.BinRulesControlsInfo={};
                        obj.BinRulesControlLinksID={};
                        obj.BinRulesControlNodesID={};
                        obj.BinRulesCount=1;
                        continue;
                        % [QUALITY] section
                    elseif strcmpi(tok(1:5),'[QUAL')
                        sect=12;d=1;
                        obj.BincountInitialQualitylines=0;
                        continue;
                        % [SOURCES] section
                    elseif strcmpi(tok(1:5),'[SOUR')
                        sect=13;
                        obj.BinNodeSourcePatternIndex = nan(1,obj.BinNodeCount);
                        obj.BinNodeSourceQuality = nan(1,obj.BinNodeCount);
                        obj.BinNodeSourceTypeCode = nan(1,obj.BinNodeCount);
                        obj.BinNodeSourceType = cell(1,obj.BinNodeCount);
                        obj.BinNodeSourcePatternNameID = cell(1,obj.BinNodeCount);
                        continue;
                        % [MIXING] section
                    elseif strcmpi(tok(1:5),'[MIXI')
                        sect=14;d=1;
                        obj.BinNodeTankMixModel={};
                        obj.BinNodeTankMixID={};
                        obj.BinNodeTankMinimumFraction=[];
                        continue;
                        % [REACTIONS] section
                    elseif strcmpi(tok(1:5),'[REAC')
                        sect=15;d=1;
                        obj.BinLinkGlobalBulkReactionCoeff=[];
                        obj.BinLinkGlobalWallReactionCoeff=[];
                        obj.BinLinkBulkReactionCoeff=[];
                        obj.BinLinkWallReactionCoeff=[];
                        obj.BincountReactionlines=0;                        
                        continue;
                        % [TIMES] section
                    elseif strcmpi(tok(1:5),'[TIME')
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
                    elseif strcmpi(tok(1:5),'[OPTI')
                        sect=17;
                        vx = NaN(obj.BinNodeCount,1);
                        vy = NaN(obj.BinNodeCount,1);
                        vertx = cell(obj.BinLinkCount,1);
                        verty = cell(obj.BinLinkCount,1);
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
                    elseif strcmpi(tok(1:5),'[COOR')
                        sect=18;
                        continue;
                        % [VERTICES] section
                    elseif strcmpi(tok(1:5),'[VERT')
                        sect=19;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline, '\s*','split');uu=1;
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
                if sect==0
                    continue;
                    % Nodes
                elseif sect==1                 
                    obj.BinNodeJunctionNameID{k}=atline{1};
                    obj.BinNodeJunctionIndex(k)=k;
                    obj.BinNodeJunctionElevation(k)=str2num(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            obj.BinNodeJunctionsBaseDemands(k)=str2num(atline{3});
                            if length(atline)>3
                                if ~sum(atline{4}==';')
                                    obj.BinNodeDemandPatternNameID{k}=atline{4};
                                else
                                    obj.BinNodeDemandPatternNameID{k}='';
                                end
                            else
                                obj.BinNodeDemandPatternNameID{k}='';
                            end
                        end
                    end
                    k=k+1;
                elseif sect==2
                    obj.BinNodeReservoirNameID{r}=atline{1};
                    obj.BinNodeReservoirIndex(r)=k;
                    obj.BinNodeReservoirElevation(r)=str2num(atline{2});
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
                    obj.BinNodeTankElevation(p)=str2num(atline{2});
                    obj.BinNodeTankInitLevel(p)=single(str2num(atline{3}));
                    obj.BinNodeTankMinLevel(p)=str2num(atline{4});
                    obj.BinNodeTankMaxLevel(p)=single(str2num(atline{5}));
                    obj.BinNodeTankDiameter(p)=str2num(atline{6});
                    obj.BinNodeTankMinVol(p)=single(str2num(atline{7}));
                    k=k+1;
                    p=p+1;
                    % Links
                elseif sect==4
                    obj.BinLinkPipeNameID{t}=atline{1};
                    obj.BinLinkPipeIndex(t)=t;
                    obj.BinLinkFromNode{t}=atline{2};
                    obj.BinLinkToNode{t}=atline{3};
                    obj.BinLinkPipeLengths(t)=str2num(atline{4});
                    obj.BinLinkPipeDiameters(t)=str2num(atline{5});
                    obj.BinLinkPipeRoughness(t)=str2num(atline{6});
                    obj.BinLinkPipeMinorLoss(t)=str2num(atline{7});
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
                    if strcmp(regexp(tline,'\w*HEAD*\w','match'),'HEAD')
                        obj.BinLinkPumpCurveNameID{q}=atline{5};
                    elseif strcmp(regexp(tline,'\w*POWER*\w','match'),'POWER')
                        obj.BinLinkPumpPower(q)=str2num(atline{5});
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
                    obj.BinLinkValveDiameters(i)=str2num(atline{4});
                    obj.BinLinkValveType{i}=atline{5};
                    obj.BinLinkValveSetting(i)=str2num(atline{6});
                    if length(atline)>6
                        obj.BinLinkValveMinorLoss(i)=str2num(atline{7});
                    end
                    t=t+1;
                    i=i+1;
                    % Curves
                elseif sect==7
                    ee=regexp(tline,'\w*EFFICIENCY*\w','match');
                    nn=regexp(tline,'\w*VOLUME*\w','match');
                    kk=regexp(tline,'\w*HEADLOSS*\w','match');

                    if strcmp(ee,'EFFICIENCY'), typecode=1;   % EFFICIENCY
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif strcmp(nn,'VOLUME'), typecode=2;   % VOLUME
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif strcmp(kk,'HEADLOSS'), typecode=3; % HEADLOSS
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    elseif (~length(strcmp(nn,'VOLUME')) || ~length(strcmp(ee,'EFFICIENCY')) || ~length(strcmp(kk,'HEADLOSS'))) &&  (tok(1)==';'), typecode=0; % HEADLOSS
                        obj.BinCurveAllLines{b}=tline;b=b+1;continue;
                    else
                        obj.BinCurveTypes(x)=typecode;
                    end
                    a = textscan(tline,'%s %f %f');l=1;
                    aa = regexp(tline, '\s*','split');
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
                    indd=find(strcmpi(obj.BinNodeNameID,atline{1}));
                    if ~isempty(obj.BinNodeJunctionsBaseDemandsID)
                        if strcmp(obj.BinNodeJunctionsBaseDemandsID{end},atline{1})
                            pd=pd+1;obj.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2})); 
                        else
                            obj.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2})); pd=1;
                        end
                    else
                        obj.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2})); pd=1;
                    end                    
                    if length(atline)>2
                        if ~isempty(obj.BinNodeJunctionsBaseDemandsID)
                            if strcmp(obj.BinNodeJunctionsBaseDemandsID{end},atline{1})
                                obj.BinNodeDemandPatternNameID{indd}=atline{3};
                            else
                                obj.BinNodeDemandPatternNameID{indd}=atline{3}; 
                            end
                        else
                            obj.BinNodeDemandPatternNameID{indd}=atline{3};
                        end
                    else
                        obj.BinNodeDemandPatternNameID{indd}='';
                    end
                    obj.BinNodeJunctionsBaseDemandsID{d}=atline{1};                       
                    d=d+1;
                    % Status
                elseif sect==9
                    obj.BinLinkInitialStatus{d}=atline{2};
                    obj.BinLinkInitialStatusNameID{d}=atline{1};
                    obj.BincountStatuslines=d;
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
                    obj.BincountPatternlines=h;
                    d=d+1;h=h+1;                
                    % Controls
                elseif sect==11
                    obj.BinControlsInfo{d}=atline;
                    obj.BinControlLinksID{d}=atline{2};
                    t = regexp(tline, '\w*TIME\w*','match');
                    if length(t)==0
                        obj.BinControlNodesID{d}=atline{6};
                    end
                    d=d+1;
                    % Rules
                elseif sect==20
                    if strcmp(upper(atline{1}),{'RULE'})
                        dd=dd+1;d=1;
                    end
                    obj.BinRulesControlsInfo{dd}{d}=atline;
                    
                    if sum(strcmp(upper(atline{2}),{'LINK','PIPE','PUMP','VALVE'}))
                        obj.BinRulesControlLinksID{dd}{d}=atline{3}; 
                    elseif sum(strcmp(upper(atline{2}),{'NODE','JUNCTION','RESERVOIR','TANK'}))
                        obj.BinRulesControlNodesID{dd}{d}=atline{3}; 
                    end
                    d=d+1;     
                    obj.BinRulesCount=obj.BinRulesCount+1;
                    % Quality
                elseif sect==12
                    h=find(strcmpi(obj.BinNodeNameID,atline{1}));
                    obj.BinNodeInitialQuality(h)=str2num(atline{2});
                    obj.BincountInitialQualitylines=d;
                    d=d+1;
                    % Sources
                elseif sect==13
                    indexPat=find(strcmpi(obj.BinPatternNameID,atline{1}));
                    indexNode=find(strcmpi(obj.BinNodeNameID,atline{1}));
                    if length(atline)>3
                        obj.BinNodeSourcePatternIndex(indexPat)= find(strcmpi(obj.BinPatternNameID,atline{4}));
                        obj.BinNodeSourcePatternNameID{indexNode}=atline{4};
                    end
                    obj.BinNodeSourceQuality(indexNode)=str2num(atline{3});
                    obj.BinNodeSourceTypeCode(indexNode)=find((strcmpi(obj.TYPESOURCE,atline{2})-1)>-1)-1;
                    obj.BinNodeSourceType{indexNode}=obj.TYPESOURCE{obj.BinNodeSourceTypeCode(indexNode)+1};
                    % Mixing
                elseif sect==14
                    obj.BinNodeTankMixID{d}=atline{1};
                    obj.BinNodeTankMixModel{d}=atline{2};
                    obj.BinNodeTankMinimumFraction(d)=str2num(atline{3});
                    d=d+1;
                    % Reactions
                elseif sect==15
                    if strcmpi(upper(atline{1}),'GLOBAL') && strcmpi(upper(atline{2}),'BULK')
                        obj.BinLinkGlobalBulkReactionCoeff=str2num(atline{3});
                    elseif strcmpi(upper(atline{1}),'GLOBAL') && strcmpi(upper(atline{2}),'WALL')
                        obj.BinLinkGlobalWallReactionCoeff=str2num(atline{3});
                        obj.BinLinkWallReactionCoeff=obj.BinLinkGlobalWallReactionCoeff*ones(1,obj.BinLinkCount);
                        obj.BinLinkBulkReactionCoeff=obj.BinLinkGlobalBulkReactionCoeff*ones(1,obj.BinLinkCount);
                        obj.BinLinkWallReactionCoeff(obj.BinLinkPumpIndex)=0;
                        obj.BinLinkWallReactionCoeff(obj.BinLinkValveIndex)=0;
                        obj.BinLinkBulkReactionCoeff(obj.BinLinkPumpIndex)=0;
                        obj.BinLinkBulkReactionCoeff(obj.BinLinkValveIndex)=0;
                    end
                    if strcmpi(upper(atline{1}),'BULK')
                        LinkIndex = find(strcmpi(obj.BinLinkNameID,atline{2}));
                        obj.BinLinkBulkReactionCoeff(LinkIndex)=str2num(atline{3});
                    elseif strcmpi(upper(atline{1}),'WALL')
                        LinkIndex = find(strcmpi(obj.BinLinkNameID,atline{2}));
                        obj.BinLinkWallReactionCoeff(LinkIndex)=str2num(atline{3});
                    end
                    obj.BincountReactionlines=d;
                    d=d+1;
                    % Times
                elseif sect==16
                    if strcmp(upper(atline{1}),'DURATION')
                        r=atline{2};
                    elseif length(atline)>2
                        r=atline{3};
                    end
                    if sum(r==':')
                        hrs=str2num(r(1:find(r==':')-1));
                        min=str2num(r((find(r==':')+1):end-find(r==':')-1));
                        if isempty(min), min=0; end
                        if isempty(hrs), hrs=0; end
                        secnd=hrs*3600+min*60;
                    elseif sum(r=='.')
                        min=str2num(r(1:find(r=='.')-1));
                        secnd1=str2num(r((find(r=='.')+1):end));
                        if isempty(min), min=0; end
                        if isempty(secnd1), secnd1=0; end
                        secnd=min*60+secnd1;
                    elseif ~sum(r==':') && ~sum(r=='.')
                        secnd=str2num(r)*3600;
                    end
                    if strcmp(upper(atline{1}),'DURATION')
                        obj.BinTimeSimulationDuration=secnd;
                    elseif strcmp(upper(atline{1}),'HYDRAULIC')
                        obj.BinTimeHydraulicStep=secnd;
                    elseif strcmp(upper(atline{1}),'QUALITY')
                        obj.BinTimeQualityStep=secnd;
                    elseif strcmp(upper(atline{1}),'PATTERN') && strcmp(upper(atline{2}),'TIMESTEP')
                        obj.BinTimePatternStep=secnd;
                    elseif strcmp(upper(atline{1}),'PATTERN') && strcmp(upper(atline{2}),'START')
                        obj.BinTimePatternStart=secnd;
                    elseif strcmp(upper(atline{1}),'REPORT') && strcmp(upper(atline{2}),'TIMESTEP')
                        obj.BinTimeReportingStep=secnd;
                    elseif strcmp(upper(atline{1}),'REPORT') && strcmp(upper(atline{2}),'START')
                        obj.BinTimeReportingStart=secnd;
                    elseif strcmp(upper(atline{1}),'STATISTIC') 
                        obj.BinTimeStatisticsIndex=find((strcmpi(obj.TYPESTATS,atline{2})-1)>-1)-1;
                        obj.BinTimeStatistics=obj.TYPESTATS{obj.BinTimeStatisticsIndex+1};
                    end                  
                    % Options
                elseif sect==17
                    if strcmp(upper(atline{1}),'UNITS')
                        obj.BinLinkFlowUnits=atline{2};
                    elseif strcmp(upper(atline{1}),'HEADLOSS')
                        obj.BinOptionsHeadloss=atline{2};
                    elseif strcmp(upper(atline{1}),'PRESSURE')
                        obj.BinNodePressureUnits=atline{2};
                    elseif strcmp(upper(atline{1}),'SPECIFIC')
                        obj.BinOptionsSpecificGravity=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'VISCOSITY')
                        obj.BinOptionsViscosity=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'TRIALS')
                        obj.BinOptionsMaxTrials=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'ACCURACY')
                        obj.BinOptionsAccuracyValue=single(str2num(atline{2}));
                    elseif strcmp(upper(atline{1}),'UNBALANCED')
                        obj.BinOptionsUnbalanced= atline(2:end);
                    elseif strcmp(upper(atline{1}),'PATTERN')
                        obj.BinOptionsPattern=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'DEMAND')
                        obj.BinOptionsPatternDemandMultiplier=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'EMITTER')
                        obj.BinOptionsEmitterExponent=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'QUALITY')
                        obj.BinQualityType=atline{2};% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                        obj.BinQualityCode=find((strcmpi(obj.TYPEQUALITY,atline{2})-1)>-1)-1;
                        if isempty(obj.BinQualityCode)
                            obj.BinQualityCode=1;
                        end
                        if obj.BinQualityCode==3
                            obj.BinQualityTraceNodeIndex=find(strcmpi(obj.BinNodeNameID,atline{3}));
                            obj.BinQualityTraceNodeID=atline{3};
                        else
                            if length(atline)>2
                                obj.BinQualityUnits=atline{3};
                            end
                        end
                    elseif strcmp(upper(atline{1}),'DIFFUSIVITY')
                        obj.BinOptionsDiffusivity=str2num(atline{2});% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                    elseif strcmp(upper(atline{1}),'TOLERANCE')
                        obj.BinOptionsQualityTolerance=single(str2num(atline{2}));% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                    end
                    % Coordinates
                elseif sect==18
                    A = textscan(tline,'%s %f %f');
                    % get the node index
                    a=strcmp(A{1},obj.BinNodeNameID);
                    index=strfind(a,1);
                    if length(index)==0
                        return;
                    end
                    vx(index) = A{2};
                    vy(index) = A{3};
                    % Vertices
                elseif sect==19
                    A = textscan(tline,'%s %f %f');
                    vertx(indexV) = A(2);
                    verty(indexV) = A(3);
                    indexV = indexV + 1;
                end
            end
            if ~sum(obj.BinLinkBulkReactionCoeff)
                obj.BinLinkBulkReactionCoeff=obj.BinLinkGlobalBulkReactionCoeff*ones(1,obj.BinLinkCount);
                obj.BinLinkBulkReactionCoeff(obj.BinLinkPumpIndex)=0;
                obj.BinLinkBulkReactionCoeff(obj.BinLinkValveIndex)=0;
            end
            if ~sum(obj.BinLinkWallReactionCoeff)
                obj.BinLinkWallReactionCoeff=obj.BinLinkGlobalWallReactionCoeff*ones(1,obj.BinLinkCount);
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
            
            obj.BinLinkSettings = [obj.BinLinkPipeRoughness zeros(1,obj.BinLinkPumpCount) obj.BinLinkValveSetting]';
            obj.BinNodeElevations = single([obj.BinNodeJunctionElevation obj.BinNodeReservoirElevation obj.BinNodeTankElevation]);
            obj.BinLinkDiameters = single([obj.BinLinkPipeDiameters zeros(1,obj.BinLinkPumpCount) obj.BinLinkValveDiameters]);
            obj.BinLinkLengths = single([obj.BinLinkPipeLengths zeros(1,obj.BinLinkPumpCount) zeros(1,obj.BinLinkValveCount)]);
            obj.BinLinkRoughnessCoeff = [obj.BinLinkPipeRoughness zeros(1,obj.BinLinkPumpCount) zeros(1,obj.BinLinkValveCount)];
%             obj.BinNodeJunctionsBaseDemands(length(obj.BinNodeJunctionsBaseDemands):obj.BinNodeJunctionCount)=0;
            obj.BinNodeBaseDemands = single([obj.BinNodeJunctionsBaseDemands zeros(1,obj.BinNodeReservoirCount) zeros(1,obj.BinNodeTankCount)]);
            obj.BinNodeDemandPatternNameID=[obj.BinNodeDemandPatternNameID obj.BinNodeResDemandPatternNameID];
            for i=obj.BinNodeTankIndex
               obj.BinNodeDemandPatternNameID{i}=''; 
            end
            
            b={};
            for i=1:obj.BinLinkCount
                ind=find((strcmp(obj.BinLinkInitialStatusNameID,obj.BinLinkNameID{i}))==1);
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
                obj.BinLinkPumpStatusNameID=obj.BinLinkInitialStatusNameID(1:obj.BinLinkPumpCount);
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
        end
        function [errcode]=setBinNodeInitialQuality(obj,varargin)
            parameter=varargin{1};
            zz=obj.BinNodeCount-obj.BincountInitialQualitylines+1;
            sections={'[QUALITY]','[SOURCES]'};
            [errcode]=setBinParam2(obj,parameter,sections,zz);   
        end
        function [errcode]=setBinLinkReactionCoeff(obj,varargin)
            wall=[];errcode=0;
            bulk=[];
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'wall' 
                        wall=varargin{2*i};
                    case 'bulk' 
                        bulk=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end   
           links=obj.getBinLinkNameID;
           BinLinkCount=length(links.BinLinkNameID);
           [tlines]=regexp( fileread(obj.Bintempfile), '\n', 'split');
           fid = fopen([obj.pathfile,obj.Bintempfile],'w');start=0;
           for i=1:length(tlines)
               tt=regexp(tlines{i}, '\s*', 'split');
               tok = strtok(tlines{i});m=1;
               % Skip blank Clines and comments
               if isempty(tok), continue;
               elseif isempty(tt{m})
                   m=m+1;
               end
               if strcmp(upper(tt{m}),'GLOBAL') && strcmp(upper(tt{m+1}),'WALL')
                   start=i;
               end
               if strcmp(tt{m},'[MIXING]')
                   stop1=i; 
               end
               if strcmp(tt{m},'[ENERGY]')
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
               if isempty(tok)
                   tok='1';
               end
               if strcmp(tok(1),';')
               else
                   if ~isempty(wall) 
                       for e=1:BinLinkCount
                           clear atlines;
                           atlines{1} = 'WALL';  
                           atlines{2} = links.BinLinkNameID{e};  
                           atlines{3} = num2str(wall(e));  
                           newlines=[];
                           for pp=1:length(atlines)
                               newlines = [newlines, atlines{pp},blanks(12)];
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
                               newlines = [newlines, atlines{pp},blanks(12)];
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
               errcode=closeOpenNetwork(obj);
           end
        end
        function [errcode]=setBinQualityChem(obj,varargin)
            sections={'[OPTIONS]','[REPORT]'};
            indexParameter=1;
            parameter='Quality            	chem   mg/L';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [errcode]=setBinQualityNone(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	None';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [errcode]=setBinQualityAge(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	Age';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [errcode]=setBinQualityTrace(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1; 
            if ~sum(strcmp(varargin{1},obj.BinNodeNameID))
                warning('Invalid property found.');
                return
            end
            parameter=['Quality            	Trace ',varargin{1}];
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [errcode]=setBinTimeSimulationDuration(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=1;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeHydraulicStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=2;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeQualityStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=3;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimePatternStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=4;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimePatternStart(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=5;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeReportingStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=6;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeReportingStart(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=7;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeStatisticsNone(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='None';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeStatisticsAverage(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='AVERAGE';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeStatisticsMinimum(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='MINIMUM';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeStatisticsMaximum(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='MAXIMUM';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinTimeStatisticsRange(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='RANGE';
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=2;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankInitLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=3;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankMinLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=4;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankMaxLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=5;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankDiameter(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=6;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeTankMinVol(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=7;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinLinkGlobalWallReactionCoeff(obj,varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]','[MIXING]'};
            [errcode]=setBinParam(obj,3,parameter,sections);        
        end
        function [errcode]=setBinLinkGlobalBulkReactionCoeff(obj,varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]','[MIXING]'};
            [errcode]=setBinParam(obj,1,parameter,sections);        
        end
        function BinClose(obj)
            if exist([obj.Bintempfile(1:end-4),'.bin'])==2
                delete([obj.Bintempfile(1:end-4),'.bin'])
            end
            delete(obj.Bintempfile)
            delete([obj.Bintempfile(1:end-4),'.txt'])
            if exist([obj.Bintempfile(1:end-4),'.msx'])==2
                delete([obj.Bintempfile(1:end-4),'.msx'])
            end
        end
        function [errcode]=setBinLinkValvesParameters(obj,varargin)
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
                            aa(u)=sum(strcmp(Type{u},obj.TYPELINK));
                        end
                        if ~sum(aa)
                            warning('Invalid property found.');errcode=-1;
                            return;
                        end
                    case 'setting' 
                        Setting=varargin{2*i};
                    case 'minorloss' 
                        MinorLoss=varargin{2*i};
                    case 'status'  
                        if sum(strcmp(lower(varargin{2*i}),'closed')+strcmp(lower(varargin{2*i}),'open')+strcmp(lower(varargin{2*i}),'none')+strcmp(lower(varargin{2*i}),'nonestatus'))==obj.LinkValveCount
                            Status=varargin{2*i};
                        else
                            warning('Invalid argument found.');errcode=-1;
                            return;
                        end   
                        zz=abs(obj.BinLinkPumpCount+obj.BinLinkValveCount-obj.BincountStatuslines);
                        sections={'[STATUS]','[PATTERNS]','valve'};
                        setBinParam2(obj,Status,sections,zz);   
                    otherwise
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp( fileread(obj.Bintempfile), '\n', 'split');
            fid = fopen([obj.pathfile,obj.Bintempfile],'w');
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1},'[VALVES]')
                    start=i;
                end
                if strcmp(tt{1},'[DEMANDS]')
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
                   if isempty(tok), tok(1)='1'; end
                   % Skip blank Clines and comments
                   if strcmp(tok(1),';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       a = regexp(tlines{i}, '\s*','split');uu=1;
                       clear atlines;
                       for tt=1:length(a)
                            if isempty(a{tt})
                                %skip
                            elseif sum(a{tt}==';')
                                %skip
                                if tt>1,  break; end
                            else
                                atlines{uu}=a{tt}; uu=uu+1;
                            end
                       end
                       if ll<obj.BinLinkValveCount+1
                           atlines(1)=obj.BinLinkNameID(obj.BinLinkValveIndex(ll));
                           atlines(2)=obj.BinLinkFromNode(obj.BinLinkValveIndex(ll));
                           atlines(3)=obj.BinLinkToNode(obj.BinLinkValveIndex(ll));
                           if ~isempty(Diameter)%Diameters
                               atlines{4} = num2str(Diameter(ll));  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(Type)%Type
                               atlines{5} = num2str(Type{ll});  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(Setting)%Setting
                               atlines{6} = num2str(Setting(ll));  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                           if ~isempty(MinorLoss)%MinorLoss
                               atlines{7} = num2str(MinorLoss(ll));  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
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
                errcode=closeOpenNetwork(obj);
            end      
        end
        function [errcode]=setBinNodeResDemandPatternNameID(obj,varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]','[TANKS]'};
            indexParameter=3;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeReservoirElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]','[TANKS]'};
            indexParameter=2;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [errcode]=setBinNodeJunctionElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=2;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinNodeJunctionsBaseDemands(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=3;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinNodeDemandPatternNameID(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=4;
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinPattern(obj,varargin)
            idpattern=varargin{1};
            %ind=obj.getPatternIndex(idpattern);
            values=varargin{2}; %obj.getPatternValue{ind};
            sections={'[PATTERNS]','[CURVES]'};
            v=obj.getBinPatternsInfo;
            patterns=v.BinPatternNameID;
            if sum(strcmp(idpattern,patterns))
                errcode=setBinParam(obj,idpattern,values,sections);
            else
                warning('Invalid argument found.');errcode=-1;
                return
            end
        end
        function [errcode]=addBinPattern(obj,varargin)
            newidpattern=varargin{1};
            values=varargin{2};
            sections={'[PATTERNS]','[CURVES]'};
            v=obj.getBinPatternsInfo;
            patterns=v.BinPatternNameID;zz=0;
            if ~sum(strcmp(newidpattern,patterns))
                for i=1:length(patterns)
                    if mod(length(v.BinPatternValue{i}),6)==0
                        zz=zz+length(v.BinPatternValue{i})/6;
                    else
                        zz=zz+1;
                        if mod(length(values),6)
                            zz=zz+1;
                        end
                    end
                end
                errcode=setBinParam2(obj,values,sections,zz,newidpattern);
            else
                warning('Invalid argument found.');errcode=-1;
                return;
            end
        end
        function [errcode]=setBinNodeSourceQuality(obj,varargin)
            sections={'[SOURCES]','[MIXING]'};
            values=varargin{1};
            [errcode]=setBinParam(obj,11,values,sections);
        end
        function saveBinInpFile(obj)
            [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
            f = fopen([obj.pathfile,obj.Bintempfile],'w');
            % /*Write [TITLE] section */
               for i=1:length(tlines)
                   tok = strtok(tlines{i});
                   if sum(tlines{i}=='[') && ~strcmp(tok,'[TITLE]')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [JUNCTIONS] section */
               fprintf(f,'\n[JUNCTIONS]');
               sps=blanks(10);
               for u=obj.BinNodeJunctionIndex
                   fprintf(f,'\n%s%s%f',obj.BinNodeNameID{u},sps,obj.BinNodeElevations(u));
               end
            % % /*Write [RESERVOIRS] section */
               fprintf(f,'\n[RESERVOIRS]');
               for u=obj.BinNodeReservoirIndex
                   fprintf(f,'\n%s%s%f',obj.BinNodeNameID{u},sps,obj.BinNodeElevations(u));
               end
            % % /*Write [TANKS] section */
               fprintf(f,'\n[TANKS]');b=1;
               for u=obj.BinNodeTankIndex
%                    InitLevel   	MinLevel    	MaxLevel    	Diameter    	MinVol      	VolCurve
                   fprintf(f,'\n%s%s%f%s%f%s%f%s%f%s%f%s%f',obj.BinNodeNameID{u},sps,obj.BinNodeElevations(u),sps,obj.BinNodeTankInitLevel(b),sps,obj.BinNodeTankMinLevel(b),...
                       sps,obj.BinNodeTankMaxLevel(b),sps,obj.BinNodeTankDiameter(b),sps,obj.BinNodeTankMinVol(b));
                   b=b+1;
               end
            % % /*Write [PIPES] section */
               fprintf(f,'\n[PIPES]');
               for u=obj.BinLinkPipeIndex
%                   ;ID              	Node1           	Node2           	Length      	Diameter    	Roughness   	MinorLoss   	Status
                   fprintf(f,'\n%s%s%s%s%s%s%f%s%f%s%f%s%f%s%s',obj.BinLinkNameID{u},sps,obj.BinLinkFromNode{u},sps,obj.BinLinkToNode{u},sps,obj.BinLinkLengths(u),sps,obj.BinLinkDiameters(u),...
                       sps,obj.BinLinkPipeRoughness(u),sps,obj.BinLinkPipeMinorLoss(u),sps,obj.BinLinkInitialStatus{u});
               end
            % % /*Write [PUMPS] section */
               fprintf(f,'\n[PUMPS]');
               par={';ID',';Node',';Junction',';Demand',';Type',';Tank',';Link'};
               for pp=1:length(par)
                   if find(strcmp(strtok(tlines),par{pp}))
                       for i=find(strcmp(strtok(tlines),par{pp}))
                          tlines{i}='';
                       end
                   end
               end
               for i=find(strcmp(strtok(tlines),'[PUMPS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [VALVES] section */
               fprintf(f,'\n[VALVES]');
               for i=find(strcmp(strtok(tlines),'[VALVES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [DEMANDS] section */
               fprintf(f,'\n[DEMANDS]');
               for u=1:obj.BinNodeJunctionCount
                   if isempty(obj.BinNodeDemandPatternNameID{u})
                       fprintf(f,'\n%s%s%f%s%s',obj.BinNodeNameID{u},sps,obj.BinNodeBaseDemands(u),sps,obj.BinPatternNameID{1});
                   else
                       fprintf(f,'\n%s%s%f%s%s',obj.BinNodeNameID{u},sps,obj.BinNodeBaseDemands(u),sps,obj.BinNodeDemandPatternNameID{u});
                   end
               end
               % % /*Write [EMITTERS] section */
               fprintf(f,'\n[EMITTERS]');
               for i=find(strcmp(strtok(tlines),'[EMITTERS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [STATUS] section */
                fprintf(f,'\n[STATUS]');
                for i=find(strcmp(strtok(tlines),'[STATUS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [PATTERNS] section */
                fprintf(f,'\n[PATTERNS]');
                for i=find(strcmp(strtok(tlines),'[PATTERNS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [CURVES] section */
                fprintf(f,'\n[CURVES]');
                for i=find(strcmp(strtok(tlines),'[CURVES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [CONTROLS] section */
                fprintf(f,'\n[CONTROLS]');
                for i=find(strcmp(strtok(tlines),'[CONTROLS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                      break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [QUALITY] section */
                fprintf(f,'\n[QUALITY]');
                for i=find(strcmp(strtok(tlines),'[QUALITY]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
                end
            % % /*Write [SOURCES] section */
               fprintf(f,'\n[SOURCES]');
               for i=find(strcmp(strtok(tlines),'[SOURCES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [MIXING] section */
               fprintf(f,'\n[MIXING]');
               for i=find(strcmp(strtok(tlines),'[MIXING]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   for u=1:obj.BinNodeTankCount
                       if isempty(obj.BinNodeTankMixModel)
                           obj.BinNodeTankMixModel{u}='MIXED';
                           obj.BinNodeTankMinimumFraction(u)=0;
                           obj.BinNodeTankMixID{u}=obj.BinNodeTankNameID{u};
                           fprintf(f,'\n%s%s%s%s%f',obj.BinNodeTankNameID{u},sps,obj.BinNodeTankMixModel{u},sps,obj.BinNodeTankMinimumFraction);
                       end
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [REACTIONS] section */
               fprintf(f,'\n[REACTIONS]');
               ff=find(strcmp(strtok(tlines),'[REACTIONS]'));
               for u=1:length(ff)
                   for i=ff(u)+1:length(tlines)
                       if sum(tlines{i}=='[')
                           break;
                       end
                       fprintf(f,'\n%s',tlines{i});
                   end
               end
            % % /*Write [ENERGY] section */
               fprintf(f,'\n[ENERGY]');
               for i=find(strcmp(strtok(tlines),'[ENERGY]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [TIMES] section */
               fprintf(f,'\n[TIMES]');
               for i=find(strcmp(strtok(tlines),'[TIMES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [OPTIONS] section */
            fprintf(f,'\n[OPTIONS]');
               for i=find(strcmp(strtok(tlines),'[OPTIONS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [REPORT] section */
            fprintf(f,'\n[REPORT]');
               for i=find(strcmp(strtok(tlines),'[REPORT]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [TAGS] section */
            fprintf(f,'\n[TAGS]');
               for i=find(strcmp(strtok(tlines),'[TAGS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [RULES] section */
            fprintf(f,'\n[RULES]');
               for i=find(strcmp(strtok(tlines),'[RULES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end         
            % % /*Write [COORDINATES] section */
            fprintf(f,'\n[COORDINATES]');
               for i=find(strcmp(strtok(tlines),'[COORDINATES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end     
            % % /*Write [VERTICES] section */
            fprintf(f,'\n[VERTICES]');
               for i=find(strcmp(strtok(tlines),'[VERTICES]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end
            % % /*Write [LABELS] section */
            fprintf(f,'\n[LABELS]');
               for i=find(strcmp(strtok(tlines),'[LABELS]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end    
            % % /*Write [BACKDROP] section */
            fprintf(f,'\n[BACKDROP]');
               for i=find(strcmp(strtok(tlines),'[BACKDROP]'))+1:length(tlines)
                   if sum(tlines{i}=='[')
                       break;
                   end
                   fprintf(f,'\n%s',tlines{i});
               end     
               fprintf(f, '\n[END]');
        end
        function [errcode]=addBinJunction(obj,varargin)
            newID=varargin{1};
            X=varargin{2};
            Y=varargin{3};
            newElevation=varargin{4}; 
            newBaseDemand=varargin{5};
            newDemandPattern=varargin{6};
            toNode=varargin{8}; 
            fromNode=varargin{1};
            switch upper(varargin{end})
                case 'PIPE' % typecode PIPE=1;PUMP=2 PRV=3, PSV=4, PBV=5, FCV=6, TCV=7, GPV=8
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
            end                    
            errcode=addLinkWarnings(obj,typecode,varargin{7},toNode); 
            if errcode==-1, return; end
            [errcode]=addNode(obj,0,newID,X,Y,newElevation,newBaseDemand,newDemandPattern);
            if errcode~=200 && obj.Bin==1
                return;
            end
            if typecode==1
                newPipeID=varargin{7}; 
                newLength=varargin{9};
                newDiameter=varargin{10};
                newRoughness=varargin{11};
                [errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{7}; 
                newCurveIDofPump=varargin{9}; 
                newCurveXvalue=varargin{10}; 
                newCurveYvalue=varargin{11}; 
                newCurveType=varargin{12};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{7}; 
                newValveDiameter=varargin{9}; 
                newValveSetting=varargin{10}; 
                [errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [errcode]=addBinReservoir(obj,varargin)
            newID=varargin{1};
            X=varargin{2};
            Y=varargin{3};
            newElevation=varargin{4}; 
            fromNode=varargin{1};
            toNode=varargin{6};
            switch upper(varargin{end})
                case 'PIPE' % typecode PIPE=1;PUMP=2 PRV=3, PSV=4, PBV=5, FCV=6, TCV=7, GPV=8
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
            end 
            errcode=addLinkWarnings(obj,typecode,varargin{5},toNode); 
            if errcode==-1, return; end            
            [errcode]=addNode(obj,1,newID,X,Y,newElevation);
            if errcode
                return;
            end
            if strcmp(upper(varargin{end}),'PIPE')
                newPipeID=varargin{5}; 
                newLength=varargin{7};
                newDiameter=varargin{8};
                newRoughness=varargin{9};
                [errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{5}; 
                newCurveIDofPump=varargin{7}; 
                newCurveXvalue=varargin{8}; 
                newCurveYvalue=varargin{9}; 
                newCurveType=varargin{10};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{5}; 
                newValveDiameter=varargin{7}; 
                newValveSetting=varargin{8}; 
                [errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [errcode]=addBinTank(obj,varargin)
            newID=varargin{1};
            X=varargin{2};
            Y=varargin{3};
            MaxLevel=varargin{4};
            Diameter=varargin{5};
            Initlevel=varargin{6};
            newElevation=varargin{7};
            initqual=varargin{8};
            MinLevel=varargin{9};
            MinVol=varargin{10};
            fromNode=varargin{1};
            toNode=varargin{12};
            switch upper(varargin{end})
                case 'PIPE' % typecode PIPE=1;PUMP=2 PRV=3, PSV=4, PBV=5, FCV=6, TCV=7, GPV=8
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
            end 
            errcode=addLinkWarnings(obj,typecode,varargin{11},toNode); 
            if errcode==-1, return; end            
            [errcode]=addNode(obj,2,newID,X,Y,MaxLevel,Diameter,Initlevel,newElevation,initqual,MinLevel,MinVol);
            if errcode
                return;
            end
            if strcmp(upper(varargin{end}),'PIPE')
                newPipeID=varargin{11}; 
                newLength=varargin{13};
                newDiameter=varargin{14};
                newRoughness=varargin{15};
                [errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{11}; 
                newCurveIDofPump=varargin{13}; 
                newCurveXvalue=varargin{14}; 
                newCurveYvalue=varargin{15}; 
                newCurveType=varargin{16};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{11}; 
                newValveDiameter=varargin{13}; 
                newValveSetting=varargin{14}; 
                [errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [errcode]=addBinPipe(obj,newLink,fromNode,toNode,newLength,newDiameter,newRoughness)
            [errcode]=addLink(obj,1,newLink,fromNode,toNode,newLength,newDiameter,newRoughness);
        end
        function [errcode]=addBinPump(obj,newPumpID,fromNode,toNode,varargin)
            if length(varargin)==4
                newCurveIDofPump=varargin{1};
                newCurveXvalue=varargin{2};
                newCurveYvalue=varargin{3};
                newCurveType=varargin{4};
                if strcmp(upper(newCurveType),'PUMP')
                    errcode=addBinCurvePump(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Head
                elseif strcmp(upper(newCurveType),'EFFICIENCY')
                    errcode=addBinCurveEfficiency(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Efficiency
                elseif strcmp(upper(newCurveType),'VOLUME')
                    errcode=addBinCurveVolume(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Heigh-Volume
                elseif strcmp(upper(newCurveType),'HEADLOSS')
                    errcode=addBinCurveHeadloss(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Headloss
                end
                [errcode]=addLink(obj,2,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            elseif length(varargin)==1
                power=varargin{1};
                [errcode]=addLink(obj,2,newPumpID,fromNode,toNode,power);
            end
        end
        function [errcode]=addBinValvePRV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,3,newLink,fromNode,toNode,diameter,setting); % Pressure Reducing Valve
        end
        function [errcode]=addBinValvePSV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,4,newLink,fromNode,toNode,diameter,setting); % Pressure Sustaining Valve
        end
        function [errcode]=addBinValvePBV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,5,newLink,fromNode,toNode,diameter,setting); % Pressure Breaker Valve
        end
        function [errcode]=addBinValveFCV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,6,newLink,fromNode,toNode,diameter,setting); % Flow Control Valve
        end
        function [errcode]=addBinValveTCV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,7,newLink,fromNode,toNode,diameter,setting); % Throttle Control Valve
        end
        function [errcode]=addBinValveGPV(obj,newLink,fromNode,toNode,diameter,setting)
            [errcode]=addLink(obj,8,newLink,fromNode,toNode,diameter,setting); % General Purpose Valve
        end
        function [errcode]=addBinCurvePump(obj,newCurveID,varargin)
            CurveX=varargin{1};
            CurveY=varargin{2};
            [errcode]=addCurve(obj,newCurveID,CurveX,CurveY,0);  %ID Flow-OptionsHeadloss
        end
        function [errcode]=addBinCurveEfficiency(obj,newCurveID,CurveX,CurveY)
            [errcode]=addCurve(obj,newCurveID,CurveX,CurveY,1);  %ID Flow-Efficiency
        end
        function [errcode]=addBinCurveVolume(obj,newCurveID,CurveX,CurveY)
            [errcode]=addCurve(obj,newCurveID,CurveX,CurveY,2);  %ID Heigh-Volume
        end
        function [errcode]=addBinCurveHeadloss(obj,newCurveID,CurveX,CurveY)
            [errcode]=addCurve(obj,newCurveID,CurveX,CurveY,3);  %ID Flow-OptionsHeadloss
        end
        function [errcode]=addBinControl(obj,x,status,y_t_c,param,z,varargin)
            if nargin==6
                [errcode]=addNewControl(obj,x,status,y_t_c,param,z);
            elseif nargin==5
                [errcode]=addNewControl(obj,x,status,y_t_c,param);
            else
                [errcode]=addNewControl(obj,x,status,y_t_c);
            end
        end
        function [errcode]=removeBinNodeID(obj,NodeID)
            [errcode]=rmNode(obj,NodeID);
        end
        function [errcode]=removeBinCurveID(obj,CurveID)
            [errcode]=rmCurveID(obj,CurveID);
        end
        function [errcode]=removeBinLinkID(obj,LinkID)
            [errcode]=rmLink(obj,LinkID);
        end
        function [errcode]=removeBinControlLinkID(obj,ID)
            [errcode]=rmControl(obj,1,ID);
        end
        function [errcode]=removeBinRulesControlLinkID(obj,ID)
            [errcode]=rmRulesControl(obj,1,ID);
        end
        function [errcode]=removeBinRulesControlNodeID(obj,ID)
            [errcode]=rmRulesControl(obj,1,ID);
        end
        function [errcode]=removeBinControlNodeID(obj,ID)
            [errcode]=rmControl(obj,0,ID);
        end
        function [errcode]=setBinFlowUnitsGPM(obj)
            [errcode]=Options(obj,'GPM'); %gallons per minute
        end
        function [errcode]=setBinFlowUnitsLPS(obj)
            [errcode]=Options(obj,'LPS'); %liters per second
        end
        function [errcode]=setBinFlowUnitsMGD(obj)
            [errcode]=Options(obj,'MGD'); %million gallons per day
        end
        function [errcode]=setBinFlowUnitsIMGD(obj)
            [errcode]=Options(obj,'IMGD'); %Imperial mgd
        end
        function [errcode]=setBinFlowUnitsCFS(obj)
            [errcode]=Options(obj,'CFS'); %cubic feet per second
        end
        function [errcode]=setBinFlowUnitsAFD(obj)
            [errcode]=Options(obj,'AFD'); %acre-feet per day
        end
        function [errcode]=setBinFlowUnitsLPM(obj)
            [errcode]=Options(obj,'LPM'); %liters per minute
        end
        function [errcode]=setBinFlowUnitsMLD(obj)
            [errcode]=Options(obj,'MLD'); %million liters per day
        end
        function [errcode]=setBinFlowUnitsCMH(obj)
            [errcode]=Options(obj,'CMH'); %cubic meters per hour
        end
        function [errcode]=setBinFlowUnitsCMD(obj)
            [errcode]=Options(obj,'CMD'); %cubic meters per day
        end
        function [errcode]=setBinHeadlossHW(obj)
            [errcode]=Options(obj,'','H-W');  %Hazen-Wiliams
        end
        function [errcode]=setBinHeadlossDW(obj)
            [errcode]=Options(obj,'','D-W');  %Darcy-Weisbach
        end
        function [errcode]=setBinHeadlossCM(obj)
            [errcode]=Options(obj,'','C-M');  %Chezy-Manning
        end        
        function [errcode]=setBinLinkPipeLengths(obj,varargin)
            parameter=varargin{1};
            indexParameter=4;
            sections={'[PIPES]', '[PUMPS]'};
            [errcode]=setBinParam(obj,indexParameter,parameter,sections);
        end
        function [errcode]=setBinLinkPipeDiameters(obj,varargin)
            parameter=varargin{1};
            indexParameter=5;
            sections={'[PIPES]', '[PUMPS]'};
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinLinkPipeRoughness(obj,varargin)
            parameter=varargin{1};
            indexParameter=6;
            sections={'[PIPES]', '[PUMPS]'};
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinLinkPipeMinorLoss(obj,varargin)
            parameter=varargin{1};
            indexParameter=7;
            sections={'[PIPES]', '[PUMPS]'};
            [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinLinkPipeStatus(obj,varargin)
           indexParameter=8;
           if sum(strcmp(lower(varargin{1}),'closed')+strcmp(lower(varargin{1}),'open')+strcmp(lower(varargin{1}),'cv'))==obj.BinLinkPipeCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');errcode=-1;
                return;
           end
           sections={'[PIPES]', '[PUMPS]'};
           [errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [errcode]=setBinLinkPumpStatus(obj,varargin)
            if sum(strcmp(lower(varargin{1}),'closed')+strcmp(lower(varargin{1}),'open'))==obj.BinLinkPumpCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');errcode=-1;
                return;
            end   
            zz=abs(obj.BinLinkPumpCount+1+obj.BinLinkValveCount-obj.BincountStatuslines);
            sections={'[STATUS]','[PATTERNS]','pump'};
            [errcode]=setBinParam2(obj,parameter,sections,zz);   
        end
        function [errcode]=setBinLinkPipesParameters(obj,varargin)
            % Initiality
            Lengths=[];errcode=0;
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
                        if sum(strcmp(lower(varargin{2*i}),'closed')+strcmp(lower(varargin{2*i}),'open'))==obj.BinLinkPipeCount
                            Status=varargin{2*i};
                        else
                            warning('Invalid argument found.');errcode=-1;
                            return;
                        end
                    otherwise
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
            fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1},'[PIPES]')
                    pipes=i;
                end
                if strcmp(tt{1},'[PUMPS]')
                    pumps=i;
                end
            end
            ll=1;clear atlines;
           for i=pipes:pumps
               % Get first token in the line
               tok = strtok(tlines{i});
%                if isempty(tok), break; end
               % Skip blank Clines and comments
               if isempty(tok)
               elseif strcmp(tok(1),';')
               elseif sum(tlines{i}=='[')
               % skip
               else
                   a = regexp(tlines{i}, '\s*','split');uu=1;
                   clear atlines;
                   for tt=1:length(a)
                        if isempty(a{tt})
                            %skip
                        elseif (a{tt}==';')
                            %skip
                        else
                            atlines{uu}=a{tt}; uu=uu+1;
                        end
                   end
                   if ~isempty(Diameters)%Diameters
                       atlines{5} = num2str(Diameters(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Lengths)%Lengths
                       atlines{4} = num2str(Lengths(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Roughness)%Roughness
                       atlines{6} = num2str(Roughness(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Minorloss)%Minorloss
                       atlines{7} = num2str(Minorloss(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Status)%Status
                       atlines{8} = Status{ll};  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid);
            if obj.Bin==1
                errcode=closeOpenNetwork(obj);
            end
        end
        function [errcode]=setBinNodeJunctionsParameters(obj,varargin)
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
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp(fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
            fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1},'[JUNCTIONS]')
                    junctions=i;
                end
                if strcmp(tt{1},'[RESERVOIRS]')
                    reservoirs=i;
                end
                if strcmp(tt{1},'[DEMANDS]')
                    demands=i;
                end
                if strcmp(tt{1},'[STATUS]')
                    status=i;
                end
            end
            ll=1;
            if junctions
               for i=junctions:reservoirs
                   % Get first token in the line
                   tok = strtok(tlines{i});
                   % Skip blank Clines and comments
                   if isempty(tok), break; end
                   if strcmp(tok(1),';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       a = regexp(tlines{i}, '\s*','split');uu=1;
                       clear atlines;
                       for tt=1:length(a)
                            if isempty(a{tt})
                                %skip
                            elseif (a{tt}==';')
                                %skip
                            else
                                atlines{uu}=a{tt}; uu=uu+1;
                            end
                       end
                       if ~isempty(Elevations)%Elevations
                           atlines{2} = num2str(Elevations(ll));  
                           newlines=[];
                           for pp=1:length(atlines)
                               newlines = [newlines, atlines{pp},blanks(10)];
                           end
                           tlines{i}=newlines;
                       end
                       if ~isempty(BaseDemands) && length(atlines)>2
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{3} = num2str(BaseDemands(ll));  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       if ~isempty(patterns) && length(atlines)>3
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{4} = num2str(patterns{ll});  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
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
                   if isempty(tok), break; end
                   if strcmp(tok(1),';')
                   elseif sum(tlines{i}=='[')
                   elseif isempty(tok)
                   % skip
                   else
                       a = regexp(tlines{i}, '\s*','split');uu=1;
                       clear atlines;
                       for tt=1:length(a)
                            if isempty(a{tt})
                                %skip
                            elseif (a{tt}==';')
                                %skip
                            else
                                atlines{uu}=a{tt}; uu=uu+1;
                            end
                       end
                       if ~isempty(BaseDemands) && length(atlines)>2%BaseDemands
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{2} = num2str(BaseDemands(ll));  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
                               end
                               tlines{i}=newlines;
                           end
                       end
                       if ~isempty(patterns) && length(atlines)>2
                           if ~isempty(atlines{3}) && ~sum(atlines{3}==';')
                               atlines{3} = num2str(patterns{ll});  
                               newlines=[];
                               for pp=1:length(atlines)
                                   newlines = [newlines, atlines{pp},blanks(10)];
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
                errcode=closeOpenNetwork(obj);
            end
        end
        function [errcode]=setBinNodeTankParameters(obj,varargin)
            % Initiality
            elevations=[];errcode=0;
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
%                         v=obj.getBinNodesInfo;
                        mixm=upper(MixModel(find(cellfun('isempty',(varargin{2*i}))==0)));
                        for u=1:length(mixm)
                            if ~sum(strcmp(mixm(u),{'MIX1', 'FIFO','LIFO','MIXED','2COMP'}))
                              warning('Invalid argument found.');errcode=-1;
                              return;
                            end
                        end
                    case 'mixfraction'
                        MixFraction=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
            fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1},'[TANKS]')
                    tanks=i;
                end
                if strcmp(tt{1},'[PIPES]')
                    pipes=i;
                end
                if strcmp(tt{1},'[MIXING]')
                    start1=i;
                end
                if strcmp(tt{1},'[REACTIONS]')
                    stop1=i;
                end
            end
            ll=1;clear atlines;
           for i=tanks:pipes
               % Get first token in the line
               tok = strtok(tlines{i});
               if isempty(tok), break; end
               % Skip blank Clines and comments
               if strcmp(tok(1),';')
               elseif sum(tlines{i}=='[')
               elseif isempty(tok)
               % skip
               else
                   a = regexp(tlines{i}, '\s*','split');uu=1;
                   clear atlines;
                   for tt=1:length(a)
                        if isempty(a{tt})
                            %skip
                        elseif sum(a{tt}==';')
                            %skip
                            if tt>1,  break; end
                        else
                            atlines{uu}=a{tt}; uu=uu+1;
                        end
                   end
                   if ~isempty(elevations)%elevations
                       atlines{2} = num2str(elevations(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(InitLevel) 
                       atlines{3} = num2str(InitLevel(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MinLevel) 
                       atlines{4} = num2str(MinLevel(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MaxLevel) 
                       atlines{5} = num2str(MaxLevel(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(Diameter) 
                       atlines{6} = num2str(Diameter(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MinVol) 
                       atlines{7} = num2str(MinVol(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
           if ~isempty(MixModel)
               ll=find(cellfun('isempty',(MixModel))==0);
               ll=ll(1);
           else
               ll=1;
           end
          for i=start1+1:stop1-1
               % Get first token in the line
               tok = strtok(tlines{i});
               if isempty(tok)
                   tok='1';
               end
               if isempty(tlines{i})
               elseif strcmp(tok(1),';')
               else
                   a = regexp(tlines{i}, '\s*','split');uu=1;
                   clear atlines;
                   for tt=1:length(a)
                        if isempty(a{tt})
                            %skip
                        elseif sum(a{tt}==';')
                            %skip
                            if tt>1,  break; end
                        else
                            atlines{uu}=a{tt}; uu=uu+1;
                        end
                   end
                   if ~isempty(MixModel) && ll<obj.BinNodeTankIndex
                       atlines{2} = num2str(MixModel{ll});  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},'              	'];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(MixFraction) && ll<obj.BinNodeTankIndex
                       atlines{3} = num2str(MixFraction(ll)); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},'              	'];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
          end
          fprintf(fid, '%s\n', tlines{:});
          fclose(fid);  
           if obj.Bin==1
               errcode=closeOpenNetwork(obj);
           end
        end
        function [errcode]=setBinNodeReservoirParameters(obj,varargin)
            % Initiality
            elevations=[];errcode=0;
            patterns={};
            for i=1:(nargin/2)
                argument =lower(varargin{2*(i-1)+1});
                switch argument
                    case 'elevation' 
                        elevations=varargin{2*i};
                    case 'pattern' 
                        patterns=varargin{2*i};
                    otherwise
                        warning('Invalid property found.');errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
            fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
            for i=1:length(tlines)
                tt=regexp(tlines{i}, '\s*', 'split');
                if strcmp(tt{1},'[RESERVOIRS]')
                    reservoirs=i;
                end
                if strcmp(tt{1},'[TANKS]')
                    tanks=i;
                end
            end
           ll=1;clear atlines;
           for i=reservoirs:tanks
               % Get first token in the line
               tok = strtok(tlines{i});
               % Skip blank Clines and comments
               if isempty(tok)
               elseif strcmp(tok(1),';')
               elseif sum(tlines{i}=='[')
               % skip
               else
                   a = regexp(tlines{i}, '\s*','split');uu=1;
                   clear atlines;
                   for tt=1:length(a)
                        if isempty(a{tt})
                            %skip
                        elseif sum(a{tt}==';')
                            %skip
                            if tt>1,  break; end
                        else
                            atlines{uu}=a{tt}; uu=uu+1;
                        end
                   end
                   if ~isempty(elevations)%elevations
                       atlines{2} = num2str(elevations(ll));  
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   if ~isempty(patterns) 
                       atlines{3} = num2str(patterns{ll}); 
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
                   ll=ll+1;
               end
           end
            fprintf(fid, '%s\n', tlines{:});
            fclose(fid); 
            if obj.Bin==1
                errcode=closeOpenNetwork(obj);
            end
        end
        function value = getBinNodeSourceInfo(obj,varargin)
            sections = {'[SOURCES]' '[MIXING]'};
            value = getBinParam(obj,sections);
        end
        function value = getBinPatternsInfo(obj,varargin)
            sections = {'[PATTERNS]' '[CURVES]'};
            value = getBinParam(obj,sections);
        end
        function value = getBinNodeIndex(obj,varargin)
            v=obj.getBinNodeNameID;
            if isempty(varargin)
                value=1:length(v.BinNodeNameID);
            else
                value=find(strcmpi(v.BinNodeNameID,varargin{1}));
            end            
        end
        function value = getBinLinkIndex(obj,varargin)
            v=obj.getBinLinkNameID;
            if isempty(varargin)
                value=1:length(v.BinLinkNameID);
            else
                value=find(strcmpi(v.BinLinkNameID,varargin{1}));
            end
        end
        function value = getBinPatternIndex(obj, varargin)
            v=obj.getBinPatternsInfo;
            value = find(strcmpi(v.BinPatternNameID,varargin{1}));
        end
        function value = getBinCurvesInfo(obj)
            [value.BinCurveNameID,value.BinCurveXvalue,value.BinCurveYvalue,value.BinCurveAllLines,value.BinCurveTypes,value.BinCurveCount,value.BinCTypes]=CurveInfo(obj);
        end
        function remAddBinCurvesID(obj,newinp)
            remAddCurve(obj,newinp);
        end
        function [value] = Binplot(obj,varargin)
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
            % d.Binplot('nodes','yes','links','yes','highlightnode',{'10','11'},'highlightlink',{'10'},'fontsize',8);
            % d.Binplot('line','no');
            % d.Binplot('point','no','linksindex','yes');
            % d.Binplot('linksindex','yes','fontsize',8);
            % d.Binplot('nodesindex','yes','fontsize',14);
            [value] = ENplot(obj,'bin',1,varargin{:});
        end
        function value = getBinNumberReportingPeriods(obj,varargin)
            value = getBinComputedTimeSeries(obj,27);
        end
        function value = getBinSimulationDuration(obj,varargin)
            value = getBinComputedTimeSeries(obj,28);
        end
        function value = getBinElevationEachNode(obj,varargin)
            value = getBinComputedTimeSeries(obj,1);
        end
        function value = getBinLengthEachLink(obj,varargin)
            value = getBinComputedTimeSeries(obj,2);
        end
        function value = getBinDiameterEachLink(obj,varargin)
            value = getBinComputedTimeSeries(obj,3);
        end
        function value = getBinComputedPumpIndexListLinks(obj,varargin)
            value = getBinComputedTimeSeries(obj,4);
        end
        function value = getBinComputedPumpUtilization(obj,varargin)
            value = getBinComputedTimeSeries(obj,5);
        end
        function value = getBinComputedAverageEfficiency(obj,varargin)
            value = getBinComputedTimeSeries(obj,6);
        end
        function value = getBinComputedAverageKwattsOrMillionGallons(obj,varargin)
            value = getBinComputedTimeSeries(obj,7);
        end
        function value = getBinComputedAverageKwatts(obj,varargin)
            value = getBinComputedTimeSeries(obj,8);
        end
        function value = getBinComputedPeakKwatts(obj,varargin)
            value = getBinComputedTimeSeries(obj,9);
        end
        function value = getBinComputedAverageCostPerDay(obj,varargin)
            value = getBinComputedTimeSeries(obj,10);
        end
        function value = getBinComputedNodeDemand(obj,varargin)
            value = getBinComputedTimeSeries(obj,11);
        end
        function value = getBinComputedNodeHead(obj,varargin)
            value = getBinComputedTimeSeries(obj,12);
        end
        function value = getBinComputedNodePressure(obj,varargin)
            value = getBinComputedTimeSeries(obj,13);
        end
        function value = getBinComputedNodeQuality(obj,varargin)
            value = getBinComputedTimeSeries(obj,14);
        end
        function value = getBinComputedLinkFlow(obj,varargin)
            value = getBinComputedTimeSeries(obj,15);
        end
        function value = getBinComputedLinkVelocity(obj,varargin)
            value = getBinComputedTimeSeries(obj,16);
        end
        function value = getBinComputedLinkHeadloss(obj,varargin)
            value = getBinComputedTimeSeries(obj,17);
        end
        function value = getBinComputedLinkQuality(obj,varargin)
            value = getBinComputedTimeSeries(obj,18);
        end
        function value = getBinComputedLinkStatus(obj,varargin)
            value = getBinComputedTimeSeries(obj,19);
        end
        function value = getBinComputedLinkSetting(obj,varargin)
            value = getBinComputedTimeSeries(obj,20);
        end
        function value = getBinComputedLinkReactionRate(obj,varargin)
            value = getBinComputedTimeSeries(obj,21);
        end
        function value = getBinComputedLinkFrictionFactor(obj,varargin)
            value = getBinComputedTimeSeries(obj,22);
        end
        function value = getBinComputedAverageBulkReactionRate(obj,varargin)
            value = getBinComputedTimeSeries(obj,23);
        end
        function value = getBinComputedAverageWallReactionRate(obj,varargin)
            value = getBinComputedTimeSeries(obj,24);
        end
        function value = getBinComputedAverageTankReactionRate(obj,varargin)
            value = getBinComputedTimeSeries(obj,25);
        end
        function value = getBinComputedAverageSourceInflow(obj,varargin)
            value = getBinComputedTimeSeries(obj,26);
        end
        function value = getBinComputedAllParameters(obj)
            value=[];
            [fid,binfile,msg] = makebatfile(obj);
            if fid~=-1
                data = fread(fid,'int32');
                fclose(fid);
                value.BinNumberReportingPeriods = data(end-2);

                fid1 = fopen(binfile, 'r');

                % Seek to the 10th byte ('J'), read 5
                fseek(fid1, 0, 'bof');
                value.Binmagicnumber=fread(fid1, 1, 'uint32');
                value.Binlibepanet=fread(fid1, 1, 'uint32');
                value.BinNumberNodes=fread(fid1, 1, 'uint32');
                value.BinNumberReservoirsTanks=fread(fid1, 1, 'uint32');
                value.BinNumberLinks=fread(fid1, 1, 'uint32');
                value.BinNumberPumps=fread(fid1, 1, 'uint32');
                value.BinNumberValves=fread(fid1, 1, 'uint32');
                value.BinWaterQualityOption=fread(fid1, 1, 'uint32');
                value.BinIndexNodeSourceTracing=fread(fid1, 1, 'uint32');
                value.BinFlowUnitsOption=fread(fid1, 1, 'uint32');
                value.BinPressureUnitsOption=fread(fid1, 1, 'uint32');
                value.BinTimeStatisticsFlag=fread(fid1, 1, 'uint32');
                value.BinReportingStartTimeSec=fread(fid1, 1, 'uint32');
                value.BinReportingTimeStepSec=fread(fid1, 1, 'uint32');
                value.BinSimulationDurationSec=fread(fid1, 1, 'uint32');
                value.BinProblemTitle1=fread(fid1, 80, '*char')';
                value.BinProblemTitle2=fread(fid1, 80, '*char')';
                value.BinProblemTitle3=fread(fid1, 80, '*char')';
                value.BinNameInputFile=fread(fid1, 260, '*char')';
                value.BinNameReportFile=fread(fid1, 260, '*char')';
                value.BinNameChemical=fread(fid1, 16, '*char')';
                value.BinChemicalConcentrationUnits=fread(fid1, 32, '*char')'; 
                fread(fid1, 4, 'uint32');
                for i=1:value.BinNumberNodes
                    value.BinIDLabelEachNode{i}=fread(fid1, 32, '*char')'; % error NODES*32
                    value.BinIDLabelEachNode{i}=value.BinIDLabelEachNode{i}(find(value.BinIDLabelEachNode{i}));
                end
                for i=1:value.BinNumberLinks
                    value.BinIDLabelEachLink{i}=fread(fid1, 32, '*char')';  % error LINKS*32
                    value.BinIDLabelEachLink{i}=value.BinIDLabelEachLink{i}(find(value.BinIDLabelEachLink{i}));
                end
                value.BinIndexStartNodeEachLink=fread(fid1, value.BinNumberLinks, 'uint32')';
                value.BinIndexEndNodeEachLink=fread(fid1, value.BinNumberLinks, 'uint32')';
                value.BinTypeCodeEachLink=fread(fid1, value.BinNumberLinks, 'uint32')';
                value.BinNodeIndexEachReservoirsTank=fread(fid1, value.BinNumberReservoirsTanks, 'uint32')'; % error
                value.BinCrossSectionalAreaEachTank=fread(fid1, value.BinNumberReservoirsTanks, 'float')';
                value.BinElevationEachNode=fread(fid1, value.BinNumberNodes, 'float')';
                value.BinLengthEachLink=fread(fid1, value.BinNumberLinks, 'float')';
                value.BinDiameterEachLink=fread(fid1, value.BinNumberLinks, 'float')';

                value.BinPumpIndexListLinks=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinPumpUtilization=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinAverageEfficiency=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinAverageKwattsOrMillionGallons=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinAverageKwatts=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinPeakKwatts=fread(fid1, value.BinNumberPumps, 'float')';
                value.BinAverageCostPerDay=fread(fid1, value.BinNumberPumps, 'float')';

                fread(fid1, 1, 'float');
                
                for i=1:value.BinNumberReportingPeriods
                    value.BinnodeDemand(:,i)         = fread(fid1, value.BinNumberNodes, 'float')';
                    value.BinnodeHead(:,i)           = fread(fid1, value.BinNumberNodes, 'float')';
                    value.BinnodePressure(:,i)       = fread(fid1, value.BinNumberNodes, 'float')';
                    value.BinnodeQuality(:,i)        = fread(fid1, value.BinNumberNodes, 'float')';
                    value.BinlinkFlow(:,i)           = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkVelocity(:,i)       = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkHeadloss(:,i)       = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkQuality(:,i)        = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkStatus(:,i)         = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkSetting(:,i)        = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkReactionRate(:,i)   = fread(fid1, value.BinNumberLinks, 'float')';
                    value.BinlinkFrictionFactor(:,i) = fread(fid1, value.BinNumberLinks, 'float')';
                end
                value.BinAverageBulkReactionRate=fread(fid1, 1, 'float')';
                value.BinAverageWallReactionRate=fread(fid1, 1, 'float')';
                value.BinAverageTankReactionRate=fread(fid1, 1, 'float')';
                value.BinAverageSourceInflowRate=fread(fid1, 1, 'float')';
                value.BinNumberReportingPeriods2=fread(fid1, 1, 'uint32')';
                value.BinWarningFlag=fread(fid1, 1, 'uint32')';
                value.BinMagicNumber=fread(fid1, 1, 'uint32')';
                fclose(fid1);
            end
            if sum(strcmp(regexp(msg,'\s','split'),'errors.'))
                fprintf('"Run was unsuccessful."\n');
            else
                fprintf('"Run was successful."\n');
            end
        end
        function [info,tline,allines] = readInpFile(obj,varargin)
            if ~sum(strcmp(who,'varargin'))
                [info,tline,allines] = readAllFile(obj.Bintempfile);
            else
                [info,tline,allines] = readAllFile(obj.Bintempfile);
            end
        end
        function value = getBinNodesInfo(obj)
            valueL = obj.getBinLinksInfo;
            % Open epanet input file
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[JUNC')
                        sect=1;
                        value.BinNodeJunctionNameID={};
                        value.BinNodeType={};
                        value.BinNodeJunctionIndex=[];
                        value.BinNodeJunctionElevation=[];
                        value.BinNodeJunctionsBaseDemands=[];
                        value.BinNodeDemandPatternNameID={};
                        value.BinNodeJunctionsBaseDemandsID={};k=1;
                        continue;
                    elseif strcmpi(tok(1:5),'[RESE')
                        sect=2;r=1;
                        value.BinNodeReservoirNameID={};
                        value.BinNodeReservoirIndex=[];
                        value.BinNodeReservoirElevation=[];
                        value.BinNodeResDemandPatternNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5),'[TANK')
                        sect=3;p=1;
                        value.BinNodeTankNameID={};
                        value.BinNodeTankIndex=[];
                        value.BinNodeTankElevation=[];
                        value.BinNodeTankInitLevel=[]; 
                        value.BinNodeTankMinLevel=[];  
                        value.BinNodeTankMaxLevel=[]; 
                        value.BinNodeTankDiameter=[];  
                        value.BinNodeTankMinVol=[];
                        continue;
                    elseif strcmpi(tok(1:5),'[PIPE')
                        value.BinNodeJunctionCount = length(value.BinNodeJunctionNameID); 
                        value.BinNodeReservoirCount = length(value.BinNodeReservoirNameID); 
                        value.BinNodeTankCount = length(value.BinNodeTankNameID);
                        value.BinNodeTankReservoirCount = value.BinNodeReservoirCount + value.BinNodeTankCount;
                        value.BinNodeCount = value.BinNodeJunctionCount+value.BinNodeTankReservoirCount;
                        value.BinNodeNameID=[value.BinNodeJunctionNameID value.BinNodeReservoirNameID value.BinNodeTankNameID];
                        vx = NaN(value.BinNodeCount,1);
                        vy = NaN(value.BinNodeCount,1);
                        vertx = cell(valueL.BinLinkCount,1);
                        verty = cell(valueL.BinLinkCount,1);
                        nvert = zeros(valueL.BinLinkCount,1);                        
                        value.BinNodeInitialQuality=zeros(1,value.BinNodeCount);
                        sect=4;
                        continue;
                    elseif strcmpi(tok(1:5),'[DEMA') %&& max(value.BinNodeJunctionsBaseDemands)==0
                        sect=8;d=1;
%                         value.BinNodeDemandPatternNameID={};
                        continue;
                        % [QUALITY] section
                    elseif strcmpi(tok(1:5),'[QUAL')
                        sect=12;d=1;
                        value.BincountInitialQualitylines=0;
                        continue;                        
                        % [MIXING] section
                    elseif strcmpi(tok(1:5),'[MIXI')
                        sect=14;d=1;
                        value.BinNodeTankMixModel={};
                        value.BinNodeTankMixID={};
                        value.BinNodeTankMinimumFraction=[];
                        continue;
                        % [COORDINATES] section
                    elseif strcmpi(tok(1:5),'[COOR')
                        sect=17;
                        continue;
                        % [VERTICES] section
                    elseif strcmpi(tok(1:5),'[VERT')
                        sect=18;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline,'\s*','split');uu=1;
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
                if sect==0
                    continue;
                    % Nodes
                elseif sect==1
                    value.BinNodeJunctionNameID{k}=atline{1};
                    value.BinNodeJunctionIndex(k)=k;
                    value.BinNodeJunctionElevation(k)=str2num(atline{2});
                    if length(atline)>2
                        if ~isempty(atline{3}) && ~sum(atline{3}==';')
                            value.BinNodeJunctionsBaseDemands(k)=single(str2num(atline{3}));
                            if length(atline)>3
                                if ~sum(atline{4}==';')
                                    value.BinNodeDemandPatternNameID{k}=atline{4};
                                else
                                    value.BinNodeDemandPatternNameID{k}='';
                                end
                            else
                                value.BinNodeDemandPatternNameID{k}='';
                            end
                        end
                    end
                    k=k+1;
                elseif sect==2
                    value.BinNodeReservoirNameID{r}=atline{1};
                    value.BinNodeReservoirIndex(r)=k;
                    value.BinNodeReservoirElevation(r)=str2num(atline{2});
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
                    value.BinNodeTankElevation(p)=str2num(atline{2});
                    value.BinNodeTankInitLevel(p)=single(str2num(atline{3}));
                    value.BinNodeTankMinLevel(p)=str2num(atline{4});
                    value.BinNodeTankMaxLevel(p)=single(str2num(atline{5}));
                    value.BinNodeTankDiameter(p)=str2num(atline{6});
                    value.BinNodeTankMinVol(p)=single(str2num(atline{7}));
                    k=k+1;
                    p=p+1;
                    % Demands
                elseif sect==8
                    indd=find(strcmpi(value.BinNodeNameID,atline{1}));
                    if ~isempty(value.BinNodeJunctionsBaseDemandsID)
                        if strcmp(value.BinNodeJunctionsBaseDemandsID{end},atline{1})
                            value.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2})); 
                        else
                            value.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2}));
                        end
                    else
                        value.BinNodeJunctionsBaseDemands(indd)=single(str2num(atline{2}));
                    end
                    value.BinNodeJunctionsBaseDemandsID{d}=atline{1};
                    if length(atline)>2
                        value.BinNodeDemandPatternNameID{indd}=atline{3};
                    else
                        value.BinNodeDemandPatternNameID{indd}='';
                    end
                    d=d+1;   
                    % Quality
                elseif sect==12
                    Hh=find(strcmpi(value.BinNodeNameID,atline{1}));
                    value.BinNodeInitialQuality(Hh)=str2num(atline{2});
                    value.BincountInitialQualitylines=d;
                    d=d+1;
                elseif sect==14
                    value.BinNodeTankMixID{d}=atline{1};
                    value.BinNodeTankMixModel{d}=atline{2};
                    value.BinNodeTankMinimumFraction(d)=str2num(atline{3});
                    d=d+1;
                elseif sect==17
                    A = textscan(tline,'%s %f %f');
                    % get the node index
                    a=strcmp(A{1},value.BinNodeNameID);
                    index=strfind(a,1);
                    if length(index)==0
                        return;
                    end
                    vx(index) = A{2};
                    vy(index) = A{3};
                    % Vertices
                elseif sect==18
                    A = textscan(tline,'%s %f %f');
                    index =  find(strcmp(valueL.BinLinkNameID,A{1}));
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
%             value.BinNodeJunctionsBaseDemands(length(value.BinNodeJunctionsBaseDemands):value.BinNodeJunctionCount)=0;
            value.BinNodeDemandPatternNameID=[value.BinNodeDemandPatternNameID value.BinNodeResDemandPatternNameID];
            for i=value.BinNodeTankIndex
               value.BinNodeDemandPatternNameID{i}=''; 
            end
            value.BinNodeBaseDemands = single([value.BinNodeJunctionsBaseDemands zeros(1,value.BinNodeReservoirCount) zeros(1,value.BinNodeTankCount)]);
            value.BinNodeIndex=obj.getBinNodeIndex;
        end
        function value = getBinCoordinatesSection(obj)
            % Open epanet input file
            [info]=regexp( fileread(obj.Bintempfile), '\n', 'split');%[~,info] = obj.readInpFile;
            sect=0;d=1;value={};
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
                    if strcmpi(tok(1:5),'[COOR')
                        sect=17;
                        value{d}=tline;
                        continue;
                    end
                end
                if sect==0
                    continue;
                elseif sect==17
                    d=d+1;
                    value{d}=tline;
                end
            end
        end
        function value = getBinRulesSection(obj)
            % Open epanet input file
            [~,info] = obj.readInpFile;
            sect=0;d=1;value={};
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
                    if strcmpi(tok(1:5),'[RULE')
                        sect=17;
                        value{d}=tline;
                        continue;
                    end
                end
                if sect==0
                    continue;
                elseif sect==17
                    d=d+1;
                    value{d}=tline;
                end
            end
        end
        function value = getBinNodeCoordinates(obj)
            BinNodeName = obj.getBinNodeNameID;
            BinLinkName = obj.getBinLinkNameID;
            linkcount =  length(BinLinkName.BinLinkNameID);
            nodecount = length(BinNodeName.BinNodeNameID);
            vx = NaN(nodecount,1);
            vy = NaN(nodecount,1);
            vertx = cell(linkcount,1);
            verty = cell(linkcount,1);
            nvert = zeros(linkcount,1);
            % Open epanet input file
            [~,info] = obj.readInpFile;
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
                        % [COORDINATES] section
                    if strcmpi(tok(1:5),'[COOR')
                        sect=17;
                        continue;
                        % [VERTICES] section
                    elseif strcmpi(tok(1:5),'[VERT')
                        sect=18;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline,'\s*','split');uu=1;
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
                if sect==0
                    continue;
                elseif sect==17
                    A = textscan(tline,'%s %f %f');
                    % get the node index
                    a=strcmp(A{1},BinNodeName.BinNodeNameID);
                    index=strfind(a,1);
                    if length(index)==0
                        return;
                    end
                    vx(index) = A{2};
                    vy(index) = A{3};
                    % Vertices
                elseif sect==18
                    A = textscan(tline,'%s %f %f');
                    index =  find(strcmp(BinLinkName.BinLinkNameID,A{1}));
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
        function value = getNodeCoordinates(obj)
            vx = NaN(obj.getNodeCount,1);
            vy = NaN(obj.getNodeCount,1);
            vertx = cell(obj.getLinkCount,1);
            verty = cell(obj.getLinkCount,1);
            nvert = zeros(obj.getLinkCount,1);
            % Open epanet input file
            [~,info] = obj.readInpFile;
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
                        % [VERTICES] section
                    if strcmpi(tok(1:5),'[VERT')
                        sect=18;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline,'\s*','split');uu=1;
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
                if sect==0
                    continue;
                    % Vertices
                elseif sect==18
                    A = textscan(tline,'%s %f %f');
                    [~,index] = ENgetlinkindex(char(A{1}),obj.libepanet);
                    nvert(index) = nvert(index) + 1;
                    vertx{index}(nvert(index)) = A{2};
                    verty{index}(nvert(index)) = A{3};
                end
            end
            for i=1:obj.getNodeCount
                [~,vx(i),vy(i)]=ENgetcoord(i,obj.libepanet);
            end
            value{1} = vx;
            value{2} = vy;
            value{3} = vertx;
            value{4} = verty;
        end
        function value = getBinNodeNameID(obj)
            % Open epanet input file
            % Open epanet input file
            [~,info] = obj.readInpFile;sect=0;
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
                    if strcmpi(tok(1:5),'[JUNC')
                        sect=1;k=1;
                        value.BinNodeJunctionNameID={};
                        continue;
                    elseif strcmpi(tok(1:5),'[RESE')
                        sect=2;r=1;
                        value.BinNodeReservoirNameID={};
                        continue;
                        % [TANKS] section
                    elseif strcmpi(tok(1:5),'[TANK')
                        sect=3;p=1;
                        value.BinNodeTankNameID={};
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        value.BinNodeNameID=[value.BinNodeJunctionNameID value.BinNodeReservoirNameID value.BinNodeTankNameID];
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline, '\s*','split');uu=1;
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
        function value = getBinLinkNameID(obj)
            sect=0;i=1;t=1;q=1;
            % Open epanet input file
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[PIPE')
                        sect=1;
                        value.BinLinkPipeNameID={};
                        continue;
                        % [PUMPS] section
                    elseif strcmpi(tok(1:5),'[PUMP')
                        sect=2;
                        value.BinLinkPumpNameID={};
                        continue;
                        % [VALVES] section
                    elseif strcmpi(tok(1:5),'[VALV')
                        sect=3;
                        value.BinLinkValveNameID={};                      
                        continue;
                    elseif strcmpi(tok(1:5),'[REAC')
                        sect=4;
                        continue;                        
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline, '\s*','split');uu=1;
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
            sect=0;i=1;t=1;q=1;
            % Open epanet input file
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[PIPE')
                        sect=1;
                        value.BinLinkPipeNameID={};
                        value.BinLinkPipeIndex=[];
                        value.BinLinkFromNode={};
                        value.BinLinkToNode={};
                        value.BinLinkPipeLengths=[];
                        value.BinLinkPipeDiameters=[];
                        value.BinLinkPipeRoughness=[];
                        value.BinLinkPipeMinorLoss=[];
                        continue;
                        % [PUMPS] section
                    elseif strcmpi(tok(1:5),'[PUMP')
                        sect=2;
                        value.BinLinkPumpNameID={};
                        value.BinLinkPumpIndex=[];
                        value.BinLinkPumpPatterns={};
                        value.BinLinkPumpCurveNameID={};
                        value.BinLinkPumpPower=[];
                        value.BinLinkPumpNameIDPower={};
                        continue;
                        % [VALVES] section
                    elseif strcmpi(tok(1:5),'[VALV')
                        sect=3;
                        value.BinLinkValveNameID={};
                        value.BinLinkValveIndex=[];
                        value.BinLinkValveDiameters=[];
                        value.BinLinkValveType={};
                        value.BinLinkValveSetting=[];
                        value.BinLinkValveMinorLoss=[];                        
                        continue;
                        % [STATUS] section
                    elseif strcmpi(tok(1:5),'[STAT')
                        sect=4;d=1;
                        value.BinLinkInitialStatus={};
                        value.BinLinkInitialStatusNameID={};
                        value.BincountStatuslines=0;
                        continue;                        
                    elseif strcmpi(tok(1:5),'[REAC')
                        sect=5;d=1;
                        value.BinLinkGlobalBulkReactionCoeff=[];
                        value.BinLinkGlobalWallReactionCoeff=[];
                        value.BinLinkBulkReactionCoeff=[];
                        value.BinLinkWallReactionCoeff=[];
                        value.BincountReactionlines=0;         
                        value.BinLinkPipeCount = length(value.BinLinkPipeNameID);
                        value.BinLinkPumpCount = length(value.BinLinkPumpNameID);
                        value.BinLinkValveCount = length(value.BinLinkValveNameID);
                        value.BinLinkCount = value.BinLinkPipeCount+value.BinLinkPumpCount+value.BinLinkValveCount;
                        value.BinLinkNameID = [value.BinLinkPipeNameID value.BinLinkPumpNameID value.BinLinkValveNameID];
                        continue;                        
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
                        break;
                    else
                        sect = 0;
                        continue;
                    end
                end
                clear atline;
                a = regexp(tline, '\s*','split');uu=1;
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
                if sect==0
                    continue;
                    % Links
                elseif sect==1
                    value.BinLinkPipeNameID{t}=atline{1};
                    value.BinLinkPipeIndex(t)=t;
                    value.BinLinkFromNode{t}=atline{2};
                    value.BinLinkToNode{t}=atline{3};
                    value.BinLinkPipeLengths(t)=str2num(atline{4});
                    value.BinLinkPipeDiameters(t)=str2num(atline{5});
                    value.BinLinkPipeRoughness(t)=str2num(atline{6});
                    if length(atline)>6
                        value.BinLinkPipeMinorLoss(t)=str2num(atline{7});
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
                    if strcmp(regexp(tline,'\w*HEAD*\w','match'),'HEAD')
                        value.BinLinkPumpCurveNameID{q}=atline{5};
                    elseif strcmp(regexp(tline,'\w*POWER*\w','match'),'POWER')
                        value.BinLinkPumpPower(q)=str2num(atline{5});
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
                    value.BinLinkValveDiameters(i)=str2num(atline{4});
                    value.BinLinkValveType{i}=atline{5};
                    value.BinLinkValveSetting(i)=str2num(atline{6});
                    if length(atline)>6
                        value.BinLinkValveMinorLoss(i)=str2num(atline{7});
                    end
                    t=t+1;
                    i=i+1;
                    % Status
                elseif sect==4
                    value.BinLinkInitialStatus{d}=atline{2};
                    value.BinLinkInitialStatusNameID{d}=atline{1};
                    value.BincountStatuslines=d;
                    d=d+1;                    
                    % Reactions
                elseif sect==5
                    if strcmpi(upper(atline{1}),'GLOBAL') && strcmpi(upper(atline{2}),'BULK')
                        value.BinLinkGlobalBulkReactionCoeff=str2num(atline{3});
                    elseif strcmpi(upper(atline{1}),'GLOBAL') && strcmpi(upper(atline{2}),'WALL')
                        value.BinLinkGlobalWallReactionCoeff=str2num(atline{3});
                        value.BinLinkWallReactionCoeff=value.BinLinkGlobalWallReactionCoeff*ones(1,value.BinLinkCount);
                        value.BinLinkBulkReactionCoeff=value.BinLinkGlobalBulkReactionCoeff*ones(1,value.BinLinkCount);
                        value.BinLinkWallReactionCoeff(value.BinLinkPumpIndex)=0;
                        value.BinLinkWallReactionCoeff(value.BinLinkValveIndex)=0;
                        value.BinLinkBulkReactionCoeff(value.BinLinkPumpIndex)=0;
                        value.BinLinkBulkReactionCoeff(value.BinLinkValveIndex)=0;
                    end
                    if strcmpi(upper(atline{1}),'BULK')
                        LinkIndex = find(strcmpi(value.BinLinkNameID,atline{2}));
                        value.BinLinkBulkReactionCoeff(LinkIndex)=str2num(atline{3});
                    elseif strcmpi(upper(atline{1}),'WALL')
                        LinkIndex = find(strcmpi(value.BinLinkNameID,atline{2}));
                        value.BinLinkWallReactionCoeff(LinkIndex)=str2num(atline{3});
                    end
                    value.countReactionlines=d;
                    d=d+1;
                end
            end
            if ~sum(value.BinLinkBulkReactionCoeff)
                value.BinLinkBulkReactionCoeff=value.BinLinkGlobalBulkReactionCoeff*ones(1,value.BinLinkCount);
                value.BinLinkBulkReactionCoeff(value.BinLinkPumpIndex)=0;
                value.BinLinkBulkReactionCoeff(value.BinLinkValveIndex)=0;
            end
            if ~sum(value.BinLinkWallReactionCoeff)
                value.BinLinkWallReactionCoeff=value.BinLinkGlobalWallReactionCoeff*ones(1,value.BinLinkCount);
                value.BinLinkWallReactionCoeff(value.BinLinkPumpIndex)=0;
                value.BinLinkWallReactionCoeff(value.BinLinkValveIndex)=0;
            end
            value.BinLinkSettings = [value.BinLinkPipeRoughness zeros(1,value.BinLinkPumpCount) value.BinLinkValveSetting]';
            value.BinLinkDiameters = single([value.BinLinkPipeDiameters zeros(1,value.BinLinkPumpCount) value.BinLinkValveDiameters]);
            value.BinLinkLengths = single([value.BinLinkPipeLengths zeros(1,value.BinLinkPumpCount) zeros(1,value.BinLinkValveCount)]);
            value.BinLinkRoughnessCoeff = [value.BinLinkPipeRoughness zeros(1,value.BinLinkPumpCount) zeros(1,value.BinLinkValveCount)];
            value.BinLinkType(value.BinLinkPipeIndex)=obj.TYPELINK(2);
            value.BinLinkType(value.BinLinkPumpIndex)=obj.TYPELINK(3);
            value.BinLinkType(value.BinLinkValveIndex)=value.BinLinkValveType;
            
            b={};
            for i=1:value.BinLinkCount
                ind=find((strcmp(value.BinLinkInitialStatusNameID,value.BinLinkNameID{i}))==1);
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
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[CONT')
                        sect=1;d=1;
                        value.BinControlsInfo={};
                        value.BinControlLinksID={};
                        value.BinControlNodesID={};
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
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
                    a = regexp(tline,'\s*','split');uu=1;
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
                    value.BinControlsInfo{d}=atline;
                    if length(atline)>1
                        value.BinControlLinksID{d}=atline{2};
                        t = regexp(tline, '\w*TIME\w*','match');
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
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[RULE')
                        sect=1;d=1;dd=0;
                        value.BinRulesControlsInfo={};
                        value.BinRulesControlLinksID={};
                        value.BinRulesControlNodesID={};
                        value.BinRulesCount=1;
                        continue;
                        % [END]
                    elseif strcmpi(tok(1:4),'[END')
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
                    a = regexp(tline,'\s*','split');uu=1;
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
                    if strcmp(upper(atline{1}),{'RULE'})
                        dd=dd+1;d=1;
                    end
                    value.BinRulesControlsInfo{dd}{d}=atline;
                    if sum(strcmp(upper(atline{2}),{'LINK','PIPE','PUMP','VALVE'}))
                        value.BinRulesControlLinksID{dd}{d}=atline{3};
                    elseif sum(strcmp(upper(atline{2}),{'NODE','JUNCTION','RESERVOIR','TANK'}))
                        value.BinRulesControlNodesID{dd}{d}=atline{3};
                    end
                    d=d+1;
                    value.BinRulesCount=value.BinRulesCount+1;
                end
            end
        end
        function value = getBinOptionsInfo(obj)
            % Open epanet input file
            [~,info] = obj.readInpFile;
            value.BinSImetric=0; value.BinUScustomary=0;
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
                    if strcmpi(tok(1:5),'[OPTI')
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
                    elseif strcmpi(tok(1:4),'[END')
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
                    a = regexp(tline,'\s*','split');uu=1;
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
                    if strcmp(upper(atline{1}),'UNITS')
                        value.BinLinkFlowUnits=atline{2};
                    elseif strcmp(upper(atline{1}),'HEADLOSS')
                        value.BinOptionsHeadloss=atline{2};
                    elseif strcmp(upper(atline{1}),'PRESSURE')
                        value.BinNodePressureUnits=atline{2};
                    elseif strcmp(upper(atline{1}),'SPECIFIC')
                        value.BinOptionsSpecificGravity=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'VISCOSITY')
                        value.BinOptionsViscosity=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'TRIALS')
                        value.BinOptionsMaxTrials=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'ACCURACY')
                        value.BinOptionsAccuracyValue=single(str2num(atline{2}));
                    elseif strcmp(upper(atline{1}),'UNBALANCED')
                        value.BinOptionsUnbalanced= atline(2:end);
                    elseif strcmp(upper(atline{1}),'PATTERN')
                        value.BinOptionsPattern=str2num(atline{2});
                    elseif strcmp(upper(atline{1}),'DEMAND')
                        value.BinOptionsPatternDemandMultiplier=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'EMITTER')
                        value.BinOptionsEmitterExponent=str2num(atline{3});
                    elseif strcmp(upper(atline{1}),'QUALITY')
                        value.BinQualityType=atline{2};% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                        value.BinQualityCode=find((strcmpi(obj.TYPEQUALITY,atline{2})-1)>-1)-1;
                        if isempty(value.BinQualityCode)
                            value.BinQualityCode=1;
                        end
                        if value.BinQualityCode==3 && length(atline)>2
                            value.BinQualityTraceNodeIndex=getBinNodeIndex(obj, atline{3});
                            value.BinQualityTraceNodeID=atline{3};
                        else
                            if length(atline)>2
                                value.BinQualityUnits=atline{3};
                            end
                        end
                    elseif strcmp(upper(atline{1}),'DIFFUSIVITY')
                        value.BinOptionsDiffusivity=str2num(atline{2});% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                    elseif strcmp(upper(atline{1}),'TOLERANCE')
                        value.BinOptionsQualityTolerance=single(str2num(atline{2}));% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
                    end
                end
            end
        %     US Customary - SI metric
            switch char(value.BinLinkFlowUnits)
                case 'CFS'
                    value.BinUScustomary=1;
                case 'GPM'
                    value.BinUScustomary=1;
                case 'MGD'
                    value.BinUScustomary=1;
                case 'IMGD'
                    value.BinUScustomary=1;
                case 'AFD'
                    value.BinUScustomary=1;
                case 'LPS'
                    value.BinSImetric=1;
                case 'LPM'
                    value.BinSImetric=1;
                case 'MLD'
                    value.BinSImetric=1;
                case 'CMH'
                    value.BinSImetric=1;
                case 'CMD'
                    value.BinSImetric=1;
            end
        end
        function value = getBinUnits(obj)
            value.BinLinkFlowUnits=obj.getBinOptionsInfo.BinLinkFlowUnits;
            if obj.getBinOptionsInfo.BinUScustomary==1;
                value.BinQualityUnits='mg/L or ug/L';
                value.BinNodePressureUnits='psi';
                value.BinPatternDemandsUnits=value.BinLinkFlowUnits;
                value.BinLinkPipeDiameterUnits='inches';
                value.BinNodeTankDiameterUnits='feet';
                value.BinEnergyEfficiencyUnits='percent';
                value.BinNodeElevationUnits='feet';
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
            [~,info] = obj.readInpFile;
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
                    if strcmpi(tok(1:5),'[TIME')
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
                    elseif strcmpi(tok(1:4),'[END')
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
                    a = regexp(tline,'\s*','split');uu=1;
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
                    if strcmp(upper(atline{1}),'DURATION')
                        r=atline{2};
                    elseif length(atline)>2
                        r=atline{3};
                    end
                    if sum(r==':')
                        hrs=str2num(r(1:find(r==':')-1));
                        min=str2num(r((find(r==':')+1):end-find(r==':')-1));
                        if isempty(min), min=0; end
                        if isempty(hrs), hrs=0; end
                        if ~min && ~hrs
                            g=find(r==':'); 
                            secnd=str2num(r(g(end)+1:end)); 
                        else
                        secnd=hrs*3600+min*60;
                        end
                    elseif sum(r=='.')
%                         min=str2num(r(1:find(r=='.')-1));
%                         secnd1=str2num(r((find(r=='.')+1):end));
                        if length(atline)>2; 
                            if strcmp(upper(atline{3}),'HOURS')
                                secnd=str2num(r)*3600;
                            end
                        else
                            if isempty(min), min=0; end
    %                         if isempty(secnd1), secnd1=0; end
                            min=str2num(r);
                            secnd=single(min*60);%+secnd1;
                        end
                    elseif ~sum(r==':') && ~sum(r=='.')
                        secnd=single(str2num(r)*3600);
                    end
                    if strcmp(upper(atline{1}),'DURATION')
                        value.BinTimeSimulationDuration=secnd;
                    elseif strcmp(upper(atline{1}),'HYDRAULIC')
                        value.BinTimeHydraulicStep=secnd;
                    elseif strcmp(upper(atline{1}),'QUALITY')
                        value.BinTimeQualityStep=secnd;
                    elseif strcmp(upper(atline{1}),'PATTERN') && strcmp(upper(atline{2}),'TIMESTEP')
                        value.BinTimePatternStep=secnd;
                    elseif strcmp(upper(atline{1}),'PATTERN') && strcmp(upper(atline{2}),'START')
                        value.BinTimePatternStart=secnd;
                    elseif strcmp(upper(atline{1}),'REPORT') && strcmp(upper(atline{2}),'TIMESTEP')
                        value.BinTimeReportingStep=secnd;
                    elseif strcmp(upper(atline{1}),'REPORT') && strcmp(upper(atline{2}),'START')
                        value.BinTimeReportingStart=secnd;
                    elseif strcmp(upper(atline{1}),'STATISTIC') 
                        value.BinTimeStatisticsIndex=find((strcmpi(obj.TYPESTATS,atline{2})-1)>-1)-1;
                        value.BinTimeStatistics=obj.TYPESTATS{value.BinTimeStatisticsIndex+1};
                    end       
                end
            end
        end
    end
end
function [errcode] = ENwriteline (line,libepanet)
    [errcode]=calllib(libepanet,'ENwriteline',line);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode] = ENaddpattern(patid,libepanet)
    errcode=calllib(libepanet,'ENaddpattern',patid);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode] = ENclose(libepanet)
    [errcode]=calllib(libepanet,'ENclose');
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode] = ENcloseH(libepanet)
    [errcode]=calllib(libepanet,'ENcloseH');
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode, value] = ENgetnodevalue(index, paramcode,libepanet)
    value=single(0);
    %p=libpointer('singlePtr',value);
    index=int32(index);
    paramcode=int32(paramcode);
    [errcode, value]=calllib(libepanet,'ENgetnodevalue',index, paramcode,value);
    if errcode==240
        value=NaN;
    end
end
function [errcode, value] = ENgetbasedemand(index,numdemands,libepanet)
    %epanet20100
    [errcode,value]=calllib(libepanet,'ENgetbasedemand',index,numdemands,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode, value] = ENgetnumdemands(index,libepanet)
    %epanet20100
    [errcode,value]=calllib(libepanet,'ENgetnumdemands',index,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode, value] = ENgetdemandpattern(index,numdemands,libepanet)
    %epanet20100
    [errcode,value]=calllib(libepanet,'ENgetdemandpattern',index,numdemands,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode, value] = ENgetstatistic(code,libepanet)
    %epanet20100
    [errcode,value]=calllib(libepanet,'ENgetstatistic',code,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode] = ENcloseQ(libepanet)
    [errcode]=calllib(libepanet,'ENcloseQ');
    if errcode
        ENgeterror(errcode,libepanet);
    end
end

function [errcode, ctype,lindex,setting,nindex,level] = ENgetcontrol(cindex,libepanet)
    [errcode, ctype,lindex,setting,nindex,level]=calllib(libepanet,'ENgetcontrol',cindex,0,0,0,0,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errcode, count] = ENgetcount(countcode,libepanet)
    [errcode,count]=calllib(libepanet,'ENgetcount',countcode,0);
    if errcode
        ENgeterror(errcode,libepanet);
    end
end
function [errmsg, e] = ENgeterror(errcode,libepanet)
    if errcode
        errmsg = char(32*ones(1,80));
        [e,errmsg] = calllib(libepanet,'ENgeterror',errcode,errmsg,80);
    else
        e=0;
        errmsg='';
    end
end
function [errcode,flowunitsindex] = ENgetflowunits(libepanet)
[errcode, flowunitsindex]=calllib(libepanet,'ENgetflowunits',0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,id] = ENgetlinkid(index,libepanet)
id=char(32*ones(1,17));
[errcode,id]=calllib(libepanet,'ENgetlinkid',index,id);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,index] = ENgetlinkindex(id,libepanet)
[errcode,~,index]=calllib(libepanet,'ENgetlinkindex',id,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,from,to] = ENgetlinknodes(index,libepanet)
[errcode,from,to]=calllib(libepanet,'ENgetlinknodes',index,0,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, type] = ENgetlinktype(index,libepanet)
[errcode,type]=calllib(libepanet,'ENgetlinktype',index,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetlinkvalue(index, paramcode,libepanet)
[errcode,value]=calllib(libepanet,'ENgetlinkvalue',index, paramcode, 0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,id] = ENgetnodeid(index,libepanet)
id=char(32*ones(1,17));
[errcode,id]=calllib(libepanet,'ENgetnodeid',index,id);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,index] = ENgetnodeindex(id,libepanet)
[errcode, ~, index]=calllib(libepanet,'ENgetnodeindex',id,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, type] = ENgetnodetype(index,libepanet)
[errcode,type]=calllib(libepanet,'ENgetnodetype',index,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetoption(optioncode,libepanet)
[errcode,value]=calllib(libepanet,'ENgetoption',optioncode,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, id] = ENgetpatternid(index,libepanet)
id=char(32*ones(1,31));
[errcode,id]=calllib(libepanet,'ENgetpatternid',index,id);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, id] = ENgetcurveid(index,libepanet)
%New version dev2.1
id=char(32*ones(1,31));
[errcode,id]=calllib(libepanet,'ENgetcurveid',index,id);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, index] = ENgetpatternindex(id,libepanet)
[errcode,~, index]=calllib(libepanet,'ENgetpatternindex',id,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, len] = ENgetpatternlen(index,libepanet)
[errcode,len]=calllib(libepanet,'ENgetpatternlen',index,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetpatternvalue(index, period,libepanet)
[errcode,value]=calllib(libepanet,'ENgetpatternvalue',index, period, 0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode,qualcode,tracenode] = ENgetqualtype(libepanet)
[errcode,qualcode,tracenode]=calllib(libepanet,'ENgetqualtype',0,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, timevalue] = ENgettimeparam(paramcode,libepanet)
[errcode,timevalue]=calllib(libepanet,'ENgettimeparam',paramcode,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, libepanet] = ENgetversion(libepanet)
[errcode,libepanet]=calllib(libepanet,'ENgetversion',0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENinitH(flag,libepanet)
[errcode]=calllib(libepanet,'ENinitH',flag);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENinitQ(saveflag,libepanet)
[errcode]=calllib(libepanet,'ENinitQ',saveflag);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function ENMatlabCleanup(libepanet)
% Load library
if libisloaded(libepanet)
    unloadlibrary(libepanet);
else
    errstring =['Library ', libepanet, '.dll was not loaded.'];
    disp(errstring);
end
end
function ENLoadLibrary(libepanetpath,libepanet)
if ~libisloaded(libepanet)
    loadlibrary([libepanetpath,libepanet],[libepanetpath,libepanet,'.h'])
end
if libisloaded(libepanet)
    libepanetString = 'EPANET loaded sucessfuly.';
    disp(libepanetString);
else
    warning('There was an error loading the EPANET library (DLL).')
end
end
function [errcode, tstep] = ENnextH(libepanet)
[errcode,tstep]=calllib(libepanet,'ENnextH',int32(0));
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, tstep] = ENnextQ(libepanet)
[errcode,tstep]=calllib(libepanet,'ENnextQ',int32(0));
if errcode
    ENgeterror(errcode,libepanet);
end
tstep = double(tstep);
end
function [errcode] = ENopen(inpname,repname,binname,libepanet) %DE
    errcode=calllib(libepanet,'ENopen',inpname,repname,binname);
    if errcode
       [~,errmsg] = calllib(libepanet,'ENgeterror',errcode,char(32*ones(1,80)),80);
       warning(errmsg);
    end
end
function [errcode] = ENopenH(libepanet)
[errcode]=calllib(libepanet,'ENopenH');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENopenQ(libepanet)
[errcode]=calllib(libepanet,'ENopenQ');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENreport(libepanet)
[errcode]=calllib(libepanet,'ENreport');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENresetreport(libepanet)
[errcode]=calllib(libepanet,'ENresetreport');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, t] = ENrunH(libepanet)
[errcode,t]=calllib(libepanet,'ENrunH',int32(0));
t = double(t);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, t] = ENrunQ(libepanet)
t=int32(0);
[errcode,t]=calllib(libepanet,'ENrunQ',t);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsaveH(libepanet)
[errcode]=calllib(libepanet,'ENsaveH');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsavehydfile(fname,libepanet)
[errcode]=calllib(libepanet,'ENsavehydfile',fname);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsaveinpfile(inpname,libepanet)
errcode=calllib(libepanet,'ENsaveinpfile',inpname);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetcontrol(cindex,ctype,lindex,setting,nindex,level,libepanet)
[errcode]=calllib(libepanet,'ENsetcontrol',cindex,ctype,lindex,setting,nindex,level);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetlinkvalue(index, paramcode, value,libepanet)
[errcode]=calllib(libepanet,'ENsetlinkvalue',index, paramcode, value);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetnodevalue(index, paramcode, value,libepanet)
[errcode]=calllib(libepanet,'ENsetnodevalue',index, paramcode, value);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetoption(optioncode,value,libepanet)
[errcode]=calllib(libepanet,'ENsetoption',optioncode,value);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetpattern(index, factors, nfactors,libepanet)
[errcode]=calllib(libepanet,'ENsetpattern',index,factors,nfactors);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetpatternvalue(index, period, value,libepanet)
[errcode]=calllib(libepanet,'ENsetpatternvalue',index, period, value);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,libepanet)
[errcode]=calllib(libepanet,'ENsetqualtype',qualcode,chemname,chemunits,tracenode);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetreport(command,libepanet)
[errcode]=calllib(libepanet,'ENsetreport',command);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetstatusreport(statuslevel,libepanet)
[errcode]=calllib(libepanet,'ENsetstatusreport',statuslevel);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsettimeparam(paramcode, timevalue,libepanet)
paramcode=int32(paramcode);
timevalue=int32(timevalue);
[errcode]=calllib(libepanet,'ENsettimeparam',paramcode,timevalue);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsolveH(libepanet)
[errcode]=calllib(libepanet,'ENsolveH');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsolveQ(libepanet)
[errcode]=calllib(libepanet,'ENsolveQ');
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, tleft] = ENstepQ(libepanet)
tleft=int32(0);
[errcode,tleft]=calllib(libepanet,'ENstepQ',tleft);
if errcode
    ENgeterror(errcode,libepanet);
end
tleft=double(tleft);
end
function [errcode] = ENusehydfile(hydfname,libepanet)
[errcode]=calllib(libepanet,'ENusehydfile',hydfname);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetcurve(index, x, y, nfactors,libepanet)
% New version dev2.1
[errcode]=calllib(libepanet,'ENsetcurve',index,x,y,nfactors);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, x, y] = ENgetcurvevalue(index, period,libepanet)
% New version dev2.1
[errcode,x, y]=calllib(libepanet,'ENgetcurvevalue',index, period, 0, 0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, x, y] = ENsetcurvevalue(index, pnt, x, y, libepanet)
% New version dev2.1
% index  = curve index
% pnt    = curve's point number
% x      = curve x value
% y      = curve y value                            
% sets x,y point for a specific point and curve 
[errcode]=calllib(libepanet,'ENsetcurvevalue',index, pnt, x, y);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, index] = ENgetcurveindex(id,libepanet)
% New version dev2.1
[errcode,~, index]=calllib(libepanet,'ENgetcurveindex',id,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENaddcurve(cid,libepanet)
% New version dev2.1
errcode=calllib(libepanet,'ENaddcurve',cid);
if errcode
    ENgeterror(errcode,libepanet);
end
end
% function [errcode, ids, nvalue, xvalue, yvalue] = ENgetcurve(value,libepanet)
% % New version dev2.1 bug
% [~,~,nvalue,~,~]=calllib(libepanet,'ENgetcurve',value,'',0,0,0);
% [errcode,ids,~, xvalue, yvalue]=calllib(libepanet,'ENgetcurve',value,char(32*ones(1,17)),0,zeros(1,nvalue),zeros(1,nvalue));
% if errcode
%     ENgeterror(errcode,libepanet);
% end
% end
function [errcode, len] = ENgetcurvelen(index,libepanet)
% New version dev2.1
[errcode,len]=calllib(libepanet,'ENgetcurvelen',index,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetheadcurveindex(pumpindex,libepanet)
% New version dev2.1
[errcode,value]=calllib(libepanet,'ENgetheadcurveindex',pumpindex,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetpumptype(pumpindex,libepanet)
% New version dev2.1
[errcode,value]=calllib(libepanet,'ENgetpumptype',pumpindex,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, value] = ENgetaveragepatternvalue(index,libepanet) 
% return  average pattern value
% New version dev2.1
[errcode,value]=calllib(libepanet,'ENgetaveragepatternvalue',index,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode, x, y] = ENgetcoord(index,libepanet)
% New version dev2.1
[errcode,x,y]=calllib(libepanet,'ENgetcoord',index,0,0);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetcoord(index,x,y,libepanet)
% New version dev2.1
[errcode]=calllib(libepanet,'ENsetcoord',index,x,y);
if errcode
    ENgeterror(errcode,libepanet);
end
end
function [errcode] = ENsetbasedemand(index, demandIdx, value, libepanet)
% New version dev2.1
[errcode]=calllib(libepanet,'ENsetbasedemand',index, demandIdx, value);
if errcode
    ENgeterror(errcode,libepanet);
end
end
% function [errcode,qualcode,chemname,chemunits,tracenode] = ENgetqualinfo(libepanet)
% % New version dev2.1 bug
% chm=libpointer('cstring', char(32*ones(1,16)));
% [errcode,qualcode,chemname,chemunits,tracenode]=calllib(libepanet,'ENgetqualinfo',0,chm,chm,0);
% if errcode
%     ENgeterror(errcode,libepanet);
% end
% end
function [obj] = MSXMatlabSetup(obj,msxname,varargin)
if ~isempty(varargin)
    if varargin{1}{1}~=1
        if nargin==3 
            obj.MsxlibepanetPath=char(varargin{1});
            obj.MsxlibepanetPath=[fileparts(obj.MsxlibepanetPath),'\'];
            if isempty(varargin{1})
                obj.MsxlibepanetPath='';
            end
        end
    end
else
    obj.MsxlibepanetPath=obj.libepanetpath;
end
obj.Msxlibepanet='epanetmsx'; % Get DLL libepanet (e.g. epanet20012x86 for 32-bit)
if ~libisloaded(obj.Msxlibepanet)
    loadlibrary([obj.MsxlibepanetPath,obj.Msxlibepanet],[obj.MsxlibepanetPath,[obj.Msxlibepanet,'.h']]);
end

obj.MsxFile = char(msxname);
%Save the temporary msx file
if ~isempty(varargin)
    if varargin{1}{1}==1
        if ~iscell(varargin)
            obj.MsxTempFile=obj.MsxFile;
        end
    end
else
    obj.MsxTempFile=[obj.MsxFile(1:end-4),'_temp.msx'];
    copyfile(obj.MsxFile,obj.MsxTempFile);
end
%Open the file
[obj.errcode] = MSXopen(obj);
obj.MsxEquationsTerms = obj.getMsxEquationsTerms;
obj.MsxEquationsPipes = obj.getMsxEquationsPipes;
obj.MsxEquationsTanks = obj.getMsxEquationsTanks;
obj.MsxSpeciesCount = obj.getMsxSpeciesCount;
obj.MsxConstantsCount = obj.getMsxConstantsCount;
obj.MsxParametersCount = obj.getMsxParametersCount;
obj.MsxPatternsCount = obj.getMsxPatternsCount;
obj.MsxSpeciesIndex = obj.getMsxSpeciesIndex;
obj.MsxSpeciesNameID = obj.getMsxSpeciesNameID;
obj.MsxSpeciesType = obj.getMsxSpeciesType;
obj.MsxSpeciesUnits = obj.getMsxSpeciesUnits;
obj.MsxSpeciesATOL = obj.getMsxSpeciesATOL;
obj.MsxSpeciesRTOL = obj.getMsxSpeciesRTOL;
obj.MsxConstantsNameID = obj.getMsxConstantsNameID;
obj.MsxConstantsValue  = obj.getMsxConstantsValue;
obj.MsxConstantsIndex = obj.getMsxConstantsIndex;
obj.MsxParametersNameID = obj.getMsxParametersNameID;
obj.MsxParametersIndex = obj.getMsxParametersIndex;
obj.MsxParametersTanksValue = obj.getMsxParametersTanksValue;
obj.MsxParametersPipesValue = obj.getMsxParametersPipesValue;
obj.MsxPatternsNameID = obj.getMsxPatternsNameID;
obj.MsxPatternsIndex = obj.getMsxPatternsIndex;
obj.MsxPatternsLengths = obj.getMsxPatternsLengths;
obj.MsxNodeInitqualValue = obj.getMsxNodeInitqualValue;
obj.MsxLinkInitqualValue = obj.getMsxLinkInitqualValue;
obj.MsxSources = obj.getMsxSources;
obj.MsxSourceType = obj.getMsxSourceType;
obj.MsxSourceLevel = obj.getMsxSourceLevel;
obj.MsxSourcePatternIndex = obj.getMsxSourcePatternIndex;
obj.MsxSourceNodeNameID = obj.getMsxSourceNodeNameID;
obj.MsxPattern = obj.getMsxPattern;
end
function [errcode] = MSXopen(obj,varargin)
[errcode] = calllib(obj.Msxlibepanet,'MSXopen',obj.MsxTempFile);
if errcode
    MSXerror(errcode,obj.Msxlibepanet);
end
if (errcode == 520)
    disp('current MSX project will be closed and the new project will be opened');
    [errcode] = MSXclose(obj);
    if errcode
        MSXerror(errcode,obj.Msxlibepanet);
    else
        [errcode] = calllib(obj.Msxlibepanet,'MSXopen',obj.MsxTempFile);
        if errcode
            MSXerror(errcode,obj.Msxlibepanet);
        end
    end
end
end
function [errcode] = MSXclose(obj)
[errcode] = calllib(obj.Msxlibepanet,'MSXclose');
if errcode
    MSXerror(errcode,obj.Msxlibepanet);
end
end
function [e] = MSXerror(errcode,Msxlibepanet)
len=80;
errstring=char(32*ones(1,len+1));
[e,errstring] = calllib(Msxlibepanet,'MSXgeterror',errcode,errstring,len);
disp(errstring);
end
function [errcode, count] = MSXgetcount(code,Msxlibepanet)
count=0;
[errcode,count] = calllib(Msxlibepanet,'MSXgetcount',code,count);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, index] = MSXgetindex(Msxlibepanet,varargin)
index =0;
if ~isnumeric(varargin{1})
    varargin{1}=varargin{2};
    varargin{2}=varargin{3};
end
[errcode,~,index]=calllib(Msxlibepanet,'MSXgetindex',varargin{1},varargin{2},index);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, id] = MSXgetID(type,index,len,Msxlibepanet)
id=char(32*ones(1,len+1));
[errcode,id]=calllib(Msxlibepanet,'MSXgetID',type,index,id,len);
id=id(1:len);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, len] = MSXgetIDlen(type,index,Msxlibepanet)
len=0;
[errcode,len]=calllib(Msxlibepanet,'MSXgetIDlen',type,index,len);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, type, units, atol, rtol] = MSXgetspecies(index,Msxlibepanet)
type=0; rtol=0; atol=0;
units=char(32*ones(1,16));
[errcode,type,units,atol,rtol]=calllib(Msxlibepanet,'MSXgetspecies',index,type,units,atol,rtol);
switch type
    case 0
        type='BULK';   % for a bulk water species
    case 1
        type='WALL';   % for a pipe wall surface species
end
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, value] = MSXgetconstant(index,Msxlibepanet)
value=0;
[errcode,value]=calllib(Msxlibepanet,'MSXgetconstant',index,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, value] = MSXgetparameter(type,index,param,Msxlibepanet)
value=0;
[errcode,value]=calllib(Msxlibepanet,'MSXgetparameter',type,index,param,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, patlen] = MSXgetpatternlen(patindex,Msxlibepanet)
patlen=0;
[errcode,patlen]=calllib(Msxlibepanet,'MSXgetpatternlen',patindex,patlen);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, value] = MSXgetpatternvalue(patindex,period,Msxlibepanet)
value=0;
[errcode,value]=calllib(Msxlibepanet,'MSXgetpatternvalue',patindex,period,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, value] = MSXgetinitqual(obj,index,species,Msxlibepanet)
value=0;
[errcode,value]=calllib(Msxlibepanet,'MSXgetinitqual',obj,index,species,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, type, level, pat] = MSXgetsource(node,species,Msxlibepanet)
type=0;
level=0;
pat=0;
[errcode,type,level,pat]=calllib(Msxlibepanet,'MSXgetsource',node,species,type,level,pat);
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
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function MSXMatlabCleanup(obj)
% Unload library
if libisloaded(obj.Msxlibepanet)
    unloadlibrary(obj.Msxlibepanet);
else
    errstring =['Library ', obj.Msxlibepanet, '.dll was not loaded'];
    disp(errstring);
end
end
function [errcode] = MSXsaveoutfile(outfname,Msxlibepanet)
[errcode] = calllib(Msxlibepanet,'MSXsaveoutfile',outfname);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsavemsxfile(msxname,Msxlibepanet)
[errcode] = calllib(Msxlibepanet,'MSXsavemsxfile',msxname);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetconstant(index, value, Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetconstant',index,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetparameter(type, index, param, value, Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetparameter',type,index,param,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetinitqual(type,index,species,value,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetinitqual',type,index,species,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetpattern(index, factors, nfactors, Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetpattern',index,factors,nfactors);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetpatternvalue(pat, period, value,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetpatternvalue',pat,period,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsolveQ(Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsolveQ');
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsolveH(Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsolveH');
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXaddpattern(patid,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXaddpattern',patid);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXusehydfile(hydfname,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXusehydfile',hydfname);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode, t, tleft] = MSXstep(Msxlibepanet)
t=int32(0);
tleft=int32(0);
[errcode,t,tleft]=calllib(Msxlibepanet,'MSXstep',t,tleft);
t = double(t);
tleft = double(tleft);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXinit(flag,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXinit',flag);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXreport(Msxlibepanet)
[errcode] = calllib(Msxlibepanet,'MSXreport');
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [e, errmsg] = MSXgeterror(errcode,Msxlibepanet)
errmsg = char(32*ones(1,80));
[e,errmsg] = calllib(Msxlibepanet,'MSXgeterror',errcode,errmsg,80);
if e
    MSXerror(e);
end
end
function [errcode, value] = MSXgetqual(type, index, species,Msxlibepanet)
value=0;
[errcode,value]=calllib(Msxlibepanet,'MSXgetqual',type,index,species,value);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [errcode] = MSXsetsource(node,species,type,level,pat,Msxlibepanet)
[errcode]=calllib(Msxlibepanet,'MSXsetsource',node,species,type,level,pat);
if errcode
    MSXerror(errcode,Msxlibepanet);
end
end
function [Terms,Pipes,Tanks] = getEquations(msxname)
% Open epanet input file
[fid,message]=fopen(msxname,'rt');
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
        if strcmpi(tok(1:5),'[TERM')
            sect = 1;
            continue;
            % [PIPES] section
        elseif strcmpi(tok(1:5),'[PIPE')
            sect = 2;
            continue;
            % [TANKS]
        elseif strcmpi(tok(1:5),'[TANK')
            sect = 3;
            continue;
            % [END]
        elseif strcmpi(tok(1:4),'[END')
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
end
function value=getMsxOptions(msxname)

if isempty(msxname)
    warning('Please load Msx File.');
    return;
end
% Open epanet input file
[fid,message] = fopen(msxname,'rt');
if fid < 0
    disp(message)
    return
end
% DEFAULT OPTIONS
value.areaunits='FT2'; 
value.rateunits='HR';
value.solver='EUL';
value.timestep=300;
value.atol=0.01;
value.rtol=0.001;
value.coupling='NONE';
value.compiler='NONE';
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
        if strcmpi(tok(1:5),'[OPTI')
            sect=1;
            continue;
            % [END]
        elseif strcmpi(tok(1:4),'[REP')
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
        atline={};
        a = regexp(tline,'\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        if strcmp(upper(atline{1}),'TIMESTEP')
            value.timestep=str2num(atline{2});return;
        elseif strcmp(upper(atline{1}),'AREA_UNITS')
            value.areaunits=atline{2};return;
        elseif strcmp(upper(atline{1}),'RATE_UNITS')
            value.rateunits=atline{2};return;
        elseif strcmp(upper(atline{1}),'SOLVER')
            value.solver=atline{2};return;
        elseif strcmp(upper(atline{1}),'RTOL')
            value.rtol=str2num(atline{2});return;      
        elseif strcmp(upper(atline{1}),'ATOL')
            value.atol=str2num(atline{2});return;      
        elseif strcmp(upper(atline{1}),'COUPLING')
            value.coupling=atline{2};return;     
        elseif strcmp(upper(atline{1}),'COMPILER')
            value.compiler=atline{2};return;  
        end
    end
end
end
function [axesid] = ENplot(obj,varargin)
% Initiality
highlightnode=0;
highlightlink=0;
highlightnodeindex=[];
highlightlinkindex=[];
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
for i=1:(nargin/2)
    argument =lower(varargin{2*(i-1)+1});
    switch argument
        case 'nodes' % Nodes
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            Node=varargin{2*i};
        case 'links' % Nodes
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            Link=varargin{2*i};
        case 'nodesindex' % Nodes
            if ~strcmp(lower(varargin{2*i}),'yes') 
                warning('Invalid argument.');
                return
            end
            NodeInd=varargin{2*i};
        case 'linksindex' % Links
            if ~strcmp(lower(varargin{2*i}),'yes')
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
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            npoint=varargin{2*i};
        case 'line' % color
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            lline=varargin{2*i};
        case 'axes' % color
            try
                axesid=axes('Parent',varargin{2*i});
            catch e
                axesid=varargin{2*i};
            end
        case 'bin' % color
            bin=varargin{2*i};
        otherwise
            warning('Invalid property founobj.');
            return
    end
end

if axesid==0
   g=figure;
   axesid=axes('Parent',g);
end

if cellfun('isempty',selectColorNode)==1
    init={'r'};
    for i=1:length(highlightnode)
        selectColorNode=[init selectColorNode];
    end
end
if cellfun('isempty',selectColorLink)==1
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
    v.nodesconnlinks=obj.getNodesConnectingLinksID;
    chckfunctions=libfunctions(obj.libepanet);
    if sum(strcmp(chckfunctions,'ENgetcoord'))
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

% Get node names and x, y coordiantes
if isa(highlightnode,'cell')
    for i=1:length(highlightnode)
        n = strcmp(v.nodenameid,highlightnode{i});
        if sum(n)==0
            warning('Undefined node with id "%s" in function call therefore the index is zero.', char(highlightnode{i}));
        else
            highlightnodeindex(i) = strfind(n,1);
        end
    end
end

if isa(highlightlink,'cell')
    for i=1:length(highlightlink)
        n = strcmp(v.linknameid,highlightlink{i});
        if sum(n)==0
            warning('Undefined link with id "%s" in function call therefore the index is zero.', char(highlightlink{i}));
        else
            highlightlinkindex(i) = strfind(n,1);
        end
    end
end

if (strcmp(lower(lline),'yes'))
    for i=1:v.linkcount
        FromNode=strfind(strcmp(v.nodesconnlinks(i,1),v.nodenameid),1);
        ToNode=strfind(strcmp(v.nodesconnlinks(i,2),v.nodenameid),1);

        if FromNode
            x1 = double(v.nodecoords{1}(FromNode));
            y1 = double(v.nodecoords{2}(FromNode));
        end
        if ToNode
            x2 = double(v.nodecoords{1}(ToNode));
            y2 = double(v.nodecoords{2}(ToNode));
        end

        hh=strfind(highlightlinkindex,i);

        if length(hh) && ~isempty(selectColorLink)
            line([x1 v.nodecoords{3}{i} x2],[y1 v.nodecoords{4}{i} y2],'LineWidth',1,'Color',[.5 .5 .5],'Parent',axesid);
        end
        if ~length(hh)
            h(:,4)=line([x1 v.nodecoords{3}{i} x2],[y1 v.nodecoords{4}{i} y2],'LineWidth',1,'Parent',axesid);
        end
            
        legendString{4} = char('Pipes');
        % Plot Pumps
        if sum(strfind(v.pumpindex,i))
            colornode = 'm';
            if length(hh) && isempty(selectColorLink)
                colornode = 'r';
            end
            h(:,5)=plot((x1+x2)/2,(y1+y2)/2,'mv','LineWidth',2,'MarkerEdgeColor','m',...
                'MarkerFaceColor','m',...
                'MarkerSize',5,'Parent',axesid);
            plot((x1+x2)/2,(y1+y2)/2,'mv','LineWidth',2,'MarkerEdgeColor',colornode,...
                'MarkerFaceColor',colornode,...
                'MarkerSize',5,'Parent',axesid);

            legendString{5} = char('Pumps');
        end

        % Plot Valves
        if sum(strfind(v.valveindex,i))
            colornode = 'k';
            if length(hh) && isempty(selectColorLink)
                colornode = 'r';
            end
            h(:,6)=plot((x1+x2)/2,(y1+y2)/2,'k*','LineWidth',2,'MarkerEdgeColor',colornode,...
                'MarkerFaceColor',colornode,'MarkerSize',7,'Parent',axesid);
            legendString{6} = char('Valves');
        end

        % Show Link id
        if (strcmp(lower(Link),'yes') && ~length(hh))
            text((x1+x2)/2,(y1+y2)/2,v.linknameid(i),'Fontsize',fontsize,'Parent',axesid);
        end
        % Show Link Index
        if (strcmp(lower(LinkInd),'yes') && ~length(hh))
            text((x1+x2)/2,(y1+y2)/2,num2str(v.linkindex(i)),'Fontsize',fontsize,'Parent',axesid);
        end

        if length(hh) && isempty(selectColorLink)
            line([x1,x2],[y1,y2],'LineWidth',1,'Color','r','Parent',axesid);
            text((x1+x2)/2,(y1+y2)/2,v.linknameid(i),'Fontsize',fontsize,'Parent',axesid);
        elseif length(hh) && ~isempty(selectColorLink)
            try 
                tt=length(selectColorLink{hh});
            catch err
                tt=2;
            end
           if tt>1
                if length(selectColorLink(hh))==1
                    nm{1}=selectColorLink(hh);
                else
                    nm=selectColorLink(hh);
                end
                if iscell(nm{1}) 
                    line([x1 v.nodecoords{3}{i} x2],[y1 v.nodecoords{4}{i} y2],'LineWidth',1,'Color',nm{1}{1},'Parent',axesid);
                else
                    line([x1 v.nodecoords{3}{i} x2],[y1 v.nodecoords{4}{i} y2],'LineWidth',1,'Color',nm{1},'Parent',axesid);
                end
            else
                line([x1 v.nodecoords{3}{i} x2],[y1 v.nodecoords{4}{i} y2],'LineWidth',1,'Color',char(selectColorLink(hh)),'Parent',axesid);
            end
        end
        hold on
    end
end

if (strcmp(lower(npoint),'yes'))
    if (strcmp(lower(npoint),'yes'))
        % Coordinates for node FROM
        for i=1:v.nodecount
            [x] = double(v.nodecoords{1}(i));
            [y] = double(v.nodecoords{2}(i));

            hh=strfind(highlightnodeindex,i);
            if ~length(hh)
                h(:,1)=plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','b',...
                'MarkerFaceColor','b',...
                'MarkerSize',5,'Parent',axesid);
                legendString{1}= char('Junctions');
            end

            % Plot Reservoirs
            if sum(strfind(v.resindex,i))
                colornode = 'g';
                if length(hh) && isempty(selectColorNode)
                    colornode = 'r';
                end
                h(:,2)=plot(x,y,'s','LineWidth',2,'MarkerEdgeColor','g',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',13,'Parent',axesid);
                plot(x,y,'s','LineWidth',2,'MarkerEdgeColor', colornode,...
                    'MarkerFaceColor', colornode,...
                    'MarkerSize',13,'Parent',axesid);
                legendString{2} = char('Reservoirs');
            end
            % Plot Tanks
            if sum(strfind(v.tankindex,i))
                colornode='c';
                if length(hh) && isempty(selectColorNode)
                    colornode='r';
                elseif length(hh) && ~isempty(selectColorNode)
                    colornode= 'c';
                end
                h(:,3)=plot(x,y,'p','LineWidth',2,'MarkerEdgeColor','c',...
                    'MarkerFaceColor','c',...
                    'MarkerSize',16,'Parent',axesid);

                plot(x,y,'p','LineWidth',2,'MarkerEdgeColor',colornode,...
                    'MarkerFaceColor',colornode,...
                    'MarkerSize',16,'Parent',axesid);

                legendString{3} = char('Tanks');
            end

            % Show Node id
            if (strcmp(lower(Node),'yes') && ~length(hh))
                text(x,y,v.nodenameid(i),'Fontsize',fontsize);%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
            end
            % Show Node index
            if (strcmp(lower(NodeInd),'yes') && ~length(hh))
                text(x,y,num2str(v.nodeindex(i)),'Fontsize',fontsize,'Parent',axesid);%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
            end

            if length(hh) && isempty(selectColorNode)
                plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'MarkerSize',5,'Parent',axesid);
                text(x,y,v.nodenameid(i),'Fontsize',fontsize,'Parent',axesid)%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
            elseif length(hh) && ~isempty(selectColorNode)
                try 
                    tt=length(selectColorNode{hh});
                catch err
                    tt=2;
                end
               if tt>1
                    if length(selectColorNode(hh))==1
                        nm{1}=selectColorNode(hh);
                        nmplot=nm{1}{1};
                    else
                        nm=selectColorNode(hh);
                        nmplot=nm{1};
                    end
                    if iscell(nm{1}) 
                        plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',nmplot,'MarkerFaceColor',nmplot,'MarkerSize',5,'Parent',axesid);
                    else
                        plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',nmplot,'MarkerFaceColor',nmplot,'MarkerSize',5,'Parent',axesid);
                    end
                    if sum(find(i==v.resindex))
                       plot(x,y,'s','LineWidth',2,'MarkerEdgeColor', nmplot,...
                       'MarkerFaceColor', nmplot,...
                       'MarkerSize',13,'Parent',axesid);
                    end
                    if sum(find(i==v.tankindex))
                       plot(x,y,'p','LineWidth',2,'MarkerEdgeColor',nmplot,...
                       'MarkerFaceColor',nmplot,...
                       'MarkerSize',16,'Parent',axesid);
                    end
               else
                    nmplot=char(selectColorNode(hh));
                    plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',nmplot,'MarkerFaceColor',nmplot,...
                        'MarkerSize',5,'Parent',axesid);
                    if sum(find(i==v.resindex))
                       plot(x,y,'s','LineWidth',2,'MarkerEdgeColor', nmplot,...
                       'MarkerFaceColor', nmplot,...
                       'MarkerSize',13,'Parent',axesid);
                    end
                    if sum(find(i==v.tankindex))
                       plot(x,y,'p','LineWidth',2,'MarkerEdgeColor',nmplot,...
                       'MarkerFaceColor',nmplot,...
                       'MarkerSize',16,'Parent',axesid);
                    end
               end
            end
            hold on
        end
    end
end
% Legend Plots
u=1;
for i=1:length(h)
    if h(i)~=0
        String{u} = legendString{i};
        hh(:,u) = h(i);
        u=u+1;
    end
end

legend(hh,String);
% Axis OFF and se Background
[xmax,~]=max(v.nodecoords{1});
[xmin,~]=min(v.nodecoords{1});
[ymax,~]=max(v.nodecoords{2});
[ymin,~]=min(v.nodecoords{2});

if ~isnan(ymax)
    if ymax==ymin
        xlim([xmin-((xmax-xmin)*.1),xmax+((xmax-xmin)*.1)]);
        ylim([ymin-.1,ymax+.1]);
    elseif xmax==xmin
        xlim([xmin-.1,xmax+.1]);
        ylim([ymin-(ymax-ymin)*.1,ymax+(ymax-ymin)*.1]);
    else
        xlim([xmin-((xmax-xmin)*.1),xmax+((xmax-xmin)*.1)]);
        ylim([ymin-(ymax-ymin)*.1,ymax+(ymax-ymin)*.1]);
    end
else
    warning('Undefined coordinates.');
end
axis off
whitebg('w');
set(axesid,'position',[0 0 1 1],'units','normalized');
end
function [info,tline,allines] = readAllFile(inpname)
    fid = fopen(inpname, 'rt');%or msxname
    allines = textscan(fid, '%s', 'delimiter', '\n');
    [tline]=regexp( fileread(inpname), '\n', 'split');
    for i=1:length(tline)
        str=regexp( tline{i}, '\s', 'split');
        info{i} = str(~cellfun('isempty', str));
    end
    fclose(fid);
end
function [errcode]=setBinParam(obj,indexParameter,parameter,sections,varargin)
    ok=0;errcode=0;
    if ~isempty(parameter) && (strcmpi(sections{1},'[SOURCES]')) && indexParameter==11
        indices=find(parameter.BinNodeSourceQuality>-1);
        sources=obj.getBinNodeNameID.BinNodeNameID(indices);ok=1;
    end
    if strcmp(sections{1},'[PATTERNS]')
        if ischar(indexParameter)
            param=indexParameter;
            indexParameter=21;
            pat=obj.getBinPatternsInfo;
        end
    end
    [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
    cntIDpat=0;start=0;stop1=0;stop11=0;itsOkQual=0;stop2=0;stop22=0;
    for i=1:length(tlines)
        tt=regexp(tlines{i}, '\s*', 'split');
        if strcmp(tt{1},sections{1})
            start=i;
        end
        if ~strcmp(sections{1},'[REACTIONS]')
            if strcmp(tt{1},sections{2})
                stop=i;
            end
        else
            if strcmp(tt{1},'[MIXING]')
                stop1=i;
            end
            if strcmp(tt{1},'[ENERGY]')
                stop11=i;
            end
            stop=max([stop1 stop11]);
        end
        cnts=obj.BinLinkPipeCount;
        if length(sections)>2
            if strcmp(tt{1},sections{3})
                start2=i;
            end
            if strcmp(tt{1},sections{5})
                stop2=i;
            elseif strcmp(tt{1},sections{4})
                stop22=i;
            end
            stop_2=max([stop2 stop22]);
            cnts=obj.BinNodeJunctionCount;
        end
    end
    if strcmpi(sections{1},'[RESERVOIRS]')
        cnts=obj.BinNodeReservoirCount;
    elseif strcmpi(sections{1},'[TANKS]')
        cnts=obj.BinNodeTankCount;
    end
    fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
    ll=1;clear atlines;
    for i=start:stop
       % Get first token in the line
       tok = strtok(tlines{i});
       if isempty(tok), tok(1)=';'; end
       % Skip blank Clines and comments
       if strcmp(tok(1),';') && ok==0
       elseif sum(tlines{i}=='[') && ok==0
       elseif isempty(tok) && ok==0
       % skip
       else
           a = regexp(tlines{i}, '\s*','split');uu=1;
           clear atlines;
           for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif sum(a{tt}==';')
                    %skip
                    if tt>1,  break; end
                else
                    atlines{uu}=a{tt}; uu=uu+1;
                end
           end
           if ll<cnts+1
               if (~isempty(parameter) && indexParameter~=3 && length(sections)<3 || indexParameter==2) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[REACTIONS]')) && (~strcmpi(sections{1},'[TIMES]')) && (~strcmpi(sections{1},'[OPTIONS]')) && (~strcmpi(sections{1},'[PATTERNS]')) && ~strcmpi(sections{1},'[TANKS]')
                   if indexParameter ~= 8  
                       atlines{indexParameter} = num2str(parameter(ll));  
                   else
                       atlines{indexParameter} = num2str(parameter{ll});  
                   end
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
                   end
                   tlines{i}=newlines;
               end
               if (~isempty(parameter) && length(atlines)>2 && indexParameter ~= 8) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[PATTERNS]')) && (~strcmpi(sections{1},'[REACTIONS]')) && (~strcmpi(sections{1},'[TIMES]')) && (~strcmpi(sections{1},'[OPTIONS]')) && ~strcmpi(sections{1},'[TANKS]')
                   if length(atlines)<indexParameter, atlines{indexParameter}={''}; end
                   if ~isempty(atlines{indexParameter}) && ~sum(atlines{3}==';')
                       if indexParameter==4 && length(sections)>2
                          atlines{indexParameter} = num2str(parameter{ll});  
                       else
                          atlines{indexParameter} = num2str(parameter(ll));  
                       end
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
               end
               if ~isempty(parameter) && (strcmpi(sections{1},'[RESERVOIRS]') || strcmpi(sections{1},'[TANKS]')) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[TIMES]')) && (~strcmpi(sections{1},'[OPTIONS]')) && (~strcmpi(sections{1},'[PATTERNS]')) 
                   if (indexParameter==2 && strcmpi(sections{1},'[RESERVOIRS]')) || strcmpi(sections{1},'[TANKS]')
                      atlines{indexParameter} = num2str(parameter(ll)); 
                   elseif strcmpi(sections{1},'[RESERVOIRS]')
                      atlines{indexParameter} = num2str(parameter{ll}); 
                   end
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
                   end
                   tlines{i}=newlines;
               end
           end
           if ~isempty(parameter) && (strcmpi(sections{1},'[REACTIONS]')) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[TIMES]')) && (~strcmpi(sections{1},'[OPTIONS]')) && (~strcmpi(sections{1},'[PATTERNS]')) 
               if strcmp(lower(atlines{1}),'global') 
                  if strcmp(lower(atlines{2}),'bulk') && indexParameter==1
                     atlines{3} = num2str(parameter);
                  end
                  if strcmpi(lower(atlines{2}),'wall') && indexParameter==3
                     atlines{3} = num2str(parameter);
                  end
               end
               newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp},blanks(10)];
               end
               tlines{i}=newlines;
           end
           mins=1;
           if ~isempty(parameter) && (strcmpi(sections{1},'[TIMES]')) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[OPTIONS]')) && (~strcmpi(sections{1},'[PATTERNS]')) 
               if strcmp(upper(atlines{1}),'DURATION') && indexParameter==1 
                    [mm,mins]=sec2hrs(parameter);
                    atlines{2} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'HYDRAULIC') && indexParameter==2
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'QUALITY') && indexParameter==3 && ~strcmp(upper(atlines{2}),'TRACE')
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'PATTERN') && strcmp(upper(atlines{2}),'TIMESTEP') && indexParameter==4
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'PATTERN') && strcmp(upper(atlines{2}),'START') && indexParameter==5
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'REPORT') && strcmp(upper(atlines{2}),'TIMESTEP') && indexParameter==6
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'REPORT') && strcmp(upper(atlines{2}),'START') && indexParameter==7
                    [mm,mins]=sec2hrs(parameter);
                    atlines{3} = num2str(mm);
               elseif strcmp(upper(atlines{1}),'STATISTIC') && indexParameter==8
                   atlines{2} = parameter;
               end
               if mins==0 && length(atlines)>3
                   atlines{4}='';
               end
                   newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp},blanks(10)];
               end
               tlines{i}=newlines;
           end    
           if ~isempty(parameter) && (strcmpi(sections{1},'[OPTIONS]')) && (~strcmpi(sections{1},'[SOURCES]')) && (~strcmpi(sections{1},'[PATTERNS]')) 
               if strcmp(upper(atlines{1}),'QUALITY') && indexParameter==1 && itsOkQual==0
                   clear atlines;
                   atlines{1}=parameter;itsOkQual=1;
               end
               newlines=[];
               for pp=1:length(atlines)
                   newlines = [newlines, atlines{pp},blanks(10)];
               end
               tlines{i}=newlines;
           end
           if ~isempty(parameter) && (strcmpi(sections{1},'[PATTERNS]')) && (~strcmpi(sections{1},'[SOURCES]'))
               idpat=param;
               if strcmp(atlines{1},idpat) && cntIDpat==0
                   atlines=[idpat blanks(12)];
                   cntIDpat=1;
                   newlines=atlines;
                   zz=0;lll=0;
                   lengthparam=length(parameter);
                   mlen=pat.BinPatternValue(find(strcmp(idpat,pat.BinPatternNameID)));
                   if lengthparam<length(mlen{1})
                       for j=(lengthparam+1):length(mlen{1})
                           parameter(1,j)=parameter(1,j-lengthparam);
                       end
                   end
                   for pp=1:size(parameter,1)
                       if mod(lengthparam,6)==0
                           zz=zz+lengthparam/6;
                       else
                           zz=zz+1;
                           if mod(lengthparam,6)
                               zz=zz+1;
                           end
                       end
                       m=1;
                       for k=lll+1:zz
                           if mod(lengthparam,6) && k==zz
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
                   if strcmp(atlines{1},idpat)
                       newlines='';
                   else
                   newlines=[atlines{1} blanks(12)];
                   cntIDpat=1;
                   for pp=2:length(atlines)
                       newlines = [newlines, num2str(atlines{pp}),blanks(12)];
                   end
                   end
               else 
                   newlines=tlines{i};
               end
               if isempty(parameter) && (~strcmpi(sections{1},'[PATTERNS]')) && (strcmpi(sections{1},'[SOURCES]'))
                  tlines{i}=newlines;
               end
           end
           if ~isempty(parameter) && (strcmpi(sections{1},'[SOURCES]')) && indexParameter==11
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
                   if sum(strcmp(sources,atlines{1})) || ok==1
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
                           newlines = [newlines, atlines{pp},blanks(10)];
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
           if isempty(tok), break; end
           if strcmp(tok(1),';')
           elseif sum(tlines{i}=='[')
           elseif isempty(tok)
           % skip
           else
               a = regexp(tlines{i}, '\s*','split');uu=1;
               clear atlines;
               for tt=1:length(a)
                    if isempty(a{tt})
                        %skip
                    elseif (a{tt}==';')
                        %skip
                    else
                        atlines{uu}=a{tt}; uu=uu+1;
                    end
               end
               if ~isempty(parameter) && length(atlines)>1 && indexParameter~=2%BaseDemands
                   if ~length(strfind(cell2mat(atlines),';')) %~isempty(atlines{3}) && 
                       if newInd==3
                          atlines{newInd} = num2str(parameter{ll});  
                       else
                          atlines{newInd} = num2str(parameter(ll));  
                       end
                       newlines=[];
                       for pp=1:length(atlines)
                           newlines = [newlines, atlines{pp},blanks(10)];
                       end
                       tlines{i}=newlines;
                   end
               end
               ll=ll+1;
           end
        end
    end
    clear patameter;
    fprintf(fid, '%s\n', tlines{:});
    fclose(fid);
    if obj.Bin==1
        errcode=closeOpenNetwork(obj);
    end
end
function [mm,mins] = sec2hrs(parameter)
   mm='';hrs=0;mins=0;
   if parameter >= 3600
       hrs=floor(parameter/3600);
       mm=[num2str(hrs), ':'];
   end
   if parameter >= 60
       mins=((parameter - 3600*hrs)/60);
   end
   if hrs
       mm=[mm sprintf('%d',(parameter-3600*hrs-60*mins))];
   elseif hrs==0 && mins==0
       mm=parameter;
       if mm<10
           mm=['00:00:0' num2str(mm)];
       else
           mm=['00:00:' num2str(mm)];
       end
   else
       mm=[sprintf('%.20f',mins) '       min'];
   end
end
function [errcode]=setBinParam2(obj,parameter,sections,zz,varargin)
    errcode=0;
    if strcmp(sections{1},'[STATUS]')
        value =obj.getBinLinksInfo;
        if strcmp(sections{3},'pump')
            nameID=value.BinLinkPumpStatusNameID;
            cntlv=obj.BinLinkPumpCount;
        elseif strcmp(sections{3},'valve')
            nameID=value.BinLinkValveStatusNameID;
            cntlv=obj.BinLinkValveCount;
            if strcmp(upper(parameter),'NONE'), errcode=-1;return;end
        end
    elseif strcmp(sections{1},'[PATTERNS]')
        value=obj.getBinPatternsInfo;
        if ~isempty(value.BinPatternValue)
            paramAll=[value.BinPatternValue parameter];
        else
            paramAll{1}=[value.BinPatternValue parameter];
        end
        patternsid=[value.BinPatternNameID varargin];
    end
    [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
    fid = fopen([obj.pathfile,obj.Bintempfile], 'w');
    for i=1:length(tlines)
        tt=regexp(tlines{i}, '\s*', 'split');
        tok = strtok(tlines{i});m=1;
        % Skip blank Clines and comments
        if isempty(tok), continue;
        elseif isempty(tt{m})
            m=m+1;
        end
        if strcmp(tt{m},sections{1})
            start=i;
        end
        if strcmp(tt{m},sections{2})
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
       if isempty(tok)
           tok='1';
       end
       if strcmp(tok(1),';')
       else
           clear atlines;
           if ~isempty(parameter) && strcmp(sections{1},'[STATUS]')
               for ee=1:cntlv
                   atlines{1} = nameID{ee};  
                   atlines{2} = num2str(parameter{ee});  
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
                   end
                   tlines{i+ee}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1},'[QUALITY]')
               for e=1:obj.BinNodeCount
                   atlines{1} = obj.BinNodeNameID{e};  
                   atlines{2} = num2str(parameter(e));  
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
                   end
                   tlines{i+e}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1},'[PATTERNS]')
               zz=0;ll=0;
               for e=1:length(patternsid)
                   if e<length(patternsid)
                       if mod(length(value.BinPatternValue{e}),6)==0
                           zz=zz+length(value.BinPatternValue{e})/6;
                       else 
                           zz=zz+1;
                       end
                   else
                       zz=zz+1;
                       if mod(length(paramAll{e}),6)
                           zz=zz+1;
                       end
                   end
                   m=1;
                   for k=ll+1:zz
                       if mod(length(paramAll{e}),6) && k==zz
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
           if ~isempty(parameter) && strcmp(sections{1},'Global') && varargin{1}==1
               for e=1:obj.BinLinkCount
                   atlines{1} = 'WALL';  
                   atlines{2} = obj.BinLinkNameID{e};  
                   atlines{3} = num2str(parameter(e));  
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
                   end
                   tlines{i+e}=newlines;
               end
           end
           if ~isempty(parameter) && strcmp(sections{1},'Global') && varargin{1}==2
               for e=1:obj.BinLinkCount
                   atlines{1} = 'BULK';  
                   atlines{2} = obj.BinLinkNameID{e};  
                   atlines{3} = num2str(parameter(e));  
                   newlines=[];
                   for pp=1:length(atlines)
                       newlines = [newlines, atlines{pp},blanks(10)];
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
        errcode=closeOpenNetwork(obj);
    end
end
function value = getBinParam(obj,sections,varargin)  
    warning off;
    [tlines]=regexp( fileread([obj.pathfile,obj.Bintempfile]), '\n', 'split');
    if strcmp(sections{1},'[SOURCES]')
        value.BinNodeSourcePatternIndex = nan(1,obj.BinNodeCount);
        value.BinNodeSourceQuality = nan(1,obj.BinNodeCount);
        value.BinNodeSourceTypeCode = nan(1,obj.BinNodeCount);
        value.BinNodeSourceType = cell(1,obj.BinNodeCount);
        value.BinNodeSourcePatternNameID = cell(1,obj.BinNodeCount);
    else 
        value=0;
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
        if strcmp(tt{m},sections{1})
            start=i;
        end
        if strcmp(tt{m},sections{2})
            stop=i;
        end
    end
    d=1;
   for i=start+1:stop-1
       % Get first token in the line
       tok = strtok(tlines{i});
       if isempty(tok)
           tok='1';
       end
       if strcmp(tok(1),';')
       else
           clear atline;
            a = regexp(tlines{i}, '\s*','split');uu=1;
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
            if strcmp(sections{1},'[STATUS]')
                if sum(strcmp(who,'atline'))
                    value.BinLinkInitialStatus{d}=atline{2};
                    value.BinLinkInitialStatusNameID{d}=atline{1};
                    d=d+1;
                end
           end
           if strcmp(sections{1},'[PATTERNS]')
                if sum(strcmp(who,'atline'))
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
                    value.BincountPatternlines=d;
                end
           end
           if strcmp(sections{1},'[SOURCES]')
               if sum(strcmp(who,'atline'))
                   if length(atline)>2
                       indexPat=getBinNodeIndex(obj, atline{1});
                       indexNode=getBinNodeIndex(obj, atline{1});
                       if length(atline)>3
                           value.BinNodeSourcePatternIndex(indexPat)=getBinPatternIndex(obj, atline{4});
                           value.BinNodeSourcePatternNameID{indexNode}=atline{4};
                       end
                       value.BinNodeSourceQuality(indexNode)=str2num(atline{3});
                       value.BinNodeSourceTypeCode(indexNode)=find((strcmpi(obj.TYPESOURCE,atline{2})-1)>-1)-1;
                       value.BinNodeSourceType{indexNode}=obj.TYPESOURCE{value.BinNodeSourceTypeCode(indexNode)+1};
                   end
               end
           end
       end
   end
   if strcmp(sections{1},'[STATUS]')
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
function [errcode]=addCurve(obj,newCurveID,varargin)
v=obj.getBinCurvesInfo;errcode=0;
CurveX=varargin{1};
CurveY=varargin{2};
typecode=varargin{3};
% PUMP 0 EFFICIENCY 1 VOLUME 2 HEADLOSS 3
for i=1:length(CurveX)
    if i+1<length(CurveX)+1
        if CurveX(i)>=CurveX(i+1)
            if strfind([0 1 3],typecode)
                warning('Flow values are not in ascending order.');
                errcode=-1;
                return;
            elseif typecode==2
                warning('Heigh values are not in ascending order.');
                errcode=-1;
                return;
            end
        end
    end
end

% Check if new ID already exists
i=1;exists=0;
while i<length(v.BinCurveNameID)+1
    exists(i) = strcmp(newCurveID,v.BinCurveNameID{i});
    i=i+1;
end
if sum(exists)>0
    s = sprintf('Curve "%s" already exists.',newCurveID);
    warning(s);errcode=-1;
    return
end
sect=0;
% Open and read inpname
% Read all file and save in variable info
[~,info] = obj.readInpFile;
% write
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
sps=blanks(18);
nn=0;yy=0;
for t = 1:length(info)
    a = regexp(info{t}, '\s*','split');
    if isempty(a)
        % skip
    elseif isempty(info{t})
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u},'[CURVES]')
                fprintf(fid2,'[CURVES]');
                sect=1; break;
            end
            if (sum(info{t}=='[') && nn==0)
                if yy==0
                    if sect==0
                        fprintf(fid2,'[CURVES]\n;ID                X-Value            Y-Value\n');
                    end
                    if typecode==0
                        fprintf(fid2,';PUMP: PUMP:%sX-Value%sY-Value\n',sps,sps); yy=1;
                    elseif typecode==1
                        fprintf(fid2,';PUMP: EFFICIENCY:\n'); yy=1;
                    elseif typecode==2
                        fprintf(fid2,';PUMP: VOLUME:\n'); yy=1;
                    elseif typecode==3
                        fprintf(fid2,';PUMP: HEADLOSS:\n'); yy=1;
                    end
                end
                for i=1:length(CurveX)
                    fprintf(fid2, '%s%s%d%s%d', newCurveID,sps,CurveX(i),sps,CurveY(i));
                    fprintf(fid2,'\r\n');
                end
                fprintf(fid2,'%s',a{u});
                fprintf(fid2,'\r\n');
                nn=1;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [BinCurveNameID,BinCurveXvalue,BinCurveYvalue,BinCurveAllLines,BinCurveTypes,BinCurveCount,BinCTypes]=CurveInfo(obj)
BinCurveTypes=[];Bintypecode=0;BinCNameID={};BinCurveNameID={};BinCurveCount=0;
BinCurveXvalue=[];BinCurveYvalue=[];BinCurveAllLines={};sect=0;i=1;u=1;BinCTypes=[];
cc=1;uu=1;gg=1;
% Open epanet input file
[~,info] = obj.readInpFile;
for h=1:length(info)
    tline = info{h};
    if ~ischar(tline),   break,   end
    % Get first token in the line
    tok = strtok(tline);
    % Skip blank Clines and comments
    if isempty(tok), continue, end
    ee=regexp(tline,'\w*EFFICIENCY*\w','match');
    nn=regexp(tline,'\w*VOLUME*\w','match');
    kk=regexp(tline,'\w*HEADLOSS*\w','match');
    if strcmp(ee,'EFFICIENCY'), %typecode=1;   % EFFICIENCY
    elseif strcmp(nn,'VOLUME'), %typecode=2;   % VOLUME
    elseif strcmp(kk,'HEADLOSS'), %typecode=3; % HEADLOSS
    else
        if (tok(1) == ';'), continue, end  %typecode=0;
    end
    if (tok(1) == '[')
        % [CURVES] section
        if strcmpi(tok(1:5),'[CURV')
            sect = 1;
            continue;
            % [END]
        elseif strcmpi(tok(1:4),'[END')
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
        ee=regexp(tline,'\w*EFFICIENCY*\w','match');
        nn=regexp(tline,'\w*VOLUME*\w','match');
        kk=regexp(tline,'\w*HEADLOSS*\w','match');
        if strcmp(ee,'EFFICIENCY'), Bintypecode=1;   % EFFICIENCY
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif strcmp(nn,'VOLUME'), Bintypecode=2;   % VOLUME
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif strcmp(kk,'HEADLOSS'), Bintypecode=3; % HEADLOSS
            BinCurveAllLines{u}=tline;u=u+1;continue;
        elseif (~length(strcmp(nn,'VOLUME')) || ~length(strcmp(ee,'EFFICIENCY')) || ~length(strcmp(kk,'HEADLOSS'))) &&  (tok(1)==';'), Bintypecode=0; % HEADLOSS
            BinCurveAllLines{u}=tline;u=u+1;continue;
        else
            a = textscan(tline,'%s %f %f');
            %aa=regexp(tline,'\s','split');
            BinCNameID{i}=a{1};
            if i==1
                BinCurveTypes(gg)=Bintypecode;
            elseif ~strcmp(BinCNameID{i-1},BinCNameID{i})
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
        elseif strcmp(BinCNameID{i-1},BinCNameID{i})
            BinCurveXvalue{cc}(uu)=a{2};
            BinCurveYvalue{cc}(uu)=a{3};
        elseif ~strcmp(BinCNameID{i-1},BinCNameID{i})
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
function remAddCurve(obj,newinp,varargin)
% Open and read inpname
% Read all file and save in variable info
[info,~] = readAllFile(newinp);
fid2 = fopen(newinp,'w');
sect=0;nn=0;sps=blanks(5);
for t=1:length(info)
    a = info{t};
    if isempty(a)
        % skip
    elseif isempty(info{t})
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u},'[CURVES]')
                fprintf(fid2,'[CURVES]');
                nn=0;sect=1;
                fprintf(fid2,'\n');break;
            elseif strcmp(a{u},'[CONTROLS]')
                for g=1:length(a)
                    fprintf(fid2, '%s%s', a{g},sps);
                end
                fprintf(fid2,'\n');
                nn=1;sect=0;
                break;
            elseif sect==0
                for g=1:length(a)
                    fprintf(fid2, '%s%s', a{g},sps);
                end
                fprintf(fid2,'\n');
                break;
            end
            % section [CURVES]
            if (sect==1) && (nn==0)
                for g=1:length(a)
                    fprintf(fid2, '%s%s', a{g},sps);
                end
                fprintf(fid2,'\n');
                break;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
end
fclose(fid2);
end
function errcode=addNode(obj,typecode,varargin)
    % addNode - Add node in the network. Node type codes consist of the
    % following constants: EN_JUNCTION 0 Junction node EN_RESERVOIR 1
    % Reservoir node EN_TANK 2 Tank node
    newID=varargin{1};errcode=0;
    X=varargin{2};
    Y=varargin{3};
    links = obj.getBinLinksInfo;
    nodes = obj.getBinNodesInfo;
    l = unique([links.BinLinkFromNode links.BinLinkToNode]);
    if nodes.BinNodeCount~=length(l)
        for i=1:nodes.BinNodeCount
            t = strcmpi(nodes.BinNodeNameID(i),l);
            if sum(t)==0
                s = sprintf('Node %s disconnected.',char(nodes.BinNodeNameID(i)));
                warning(s);errcode=-1;
                return;
            end
        end
    end     
    if typecode==1 || typecode==0 % junction & reservoir
        if typecode==0
            v=obj.getBinPatternsInfo;
            newidpattern=varargin{6};
            patterns=v.BinPatternNameID;
            if ~sum(strcmp(newidpattern,patterns))
                warning('Invalid argument found.');
                errcode=-1;
                return;
            end
            newBaseDemand=varargin{5};
        end
        newElevation=varargin{4};
        initqual=0;
    end
    if typecode==2
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
    if length(nodes.BinNodeNameID)==0
        warning('There is no such object in the network.');errcode=-1;
        return;
    end
    for i=1:length(nodes.BinNodeNameID)
        exists(i) = strcmp(newID,char(nodes.BinNodeNameID(i)));
    end
    if (sum(exists)==1)
        s = sprintf('Node "%s" already exists.',newID);
        warning(s);errcode=-1;
        return;
    end
    % Get type of node
    A = [0 1 2];
    code=strfind(A,typecode);
    if length(code)==0
        warning('There is no such typecode(0-2)!');errcode=-1;
        return;
    end
    % check section in inpname, [JUNCTIONS], [RESERVOIRS], [TANKS]
    stank_check=1;
    sreservoir_check=1;
    sjunction_check=1;
    % Open and read inpname
    % Read all file and save in variable info
    [~,info,~] = obj.readInpFile;
    fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
    % Initiality
    qualch=0;qq=0;
    Coordch=0;onetime=1;gg=0;
    sps1=blanks(3);
    for t = 1:length(info)
        c = info{t};
        if ~isempty(c)
            a = regexp(c, '\s*','split');
        else
            a='';
        end
        if isempty(a)
            % skip
        elseif isempty(c)
            % skip
        else
            i=i+1;
            u=1;
            while u < length(a)+1
                % Find [brackets] cnt=2;
                cnt=bracketsCheck(a{u});
                %%%%%%%% Quality Section %%%%%%%%
                if strcmp(a{u},'[QUALITY]')
                    fprintf(fid2,'[QUALITY]');
                    qualch=1;
                    break;
                end
                if (cnt==2 && qualch==1)
                    fprintf(fid2, '%s%s%d', newID,sps1,initqual);
                    fprintf(fid2,'\r\n');qq=1;
                end
                %%%%%%%% Coordinates Section %%%%%%%%
                if strcmp(a{u},'[COORDINATES]');
                    fprintf(fid2,'[COORDINATES]');
                    Coordch=1; break;
                end
                if length(strfind(c,';Node'))==1 && Coordch==1 && cnt~=2
                    break;
                    elseif u==1 && Coordch==1
                    if ((gg==0)) && (typecode==0)
                        fprintf(fid2, '%s%s%d%s%d\n', newID,sps1,X,sps1,Y);
                    end
                    gg=gg+1;
                end
                if isempty(obj.NodeCoordinates) && obj.Bin==1% no bin
                    if strcmp(a{u},'[END]')
                        fprintf(fid2,'%s','[COORDINATES]');
                        fprintf(fid2,'\r\n');
                        for qq=1:length(X)
                            fprintf(fid2,'%s%s%d%s%d',char(newID(qq)),sps1,X(qq),sps1,Y(qq));
                            fprintf(fid2,'\r\n');
                        end
                        fprintf(fid2, '%s%s%d%s%d\n',...
                        newID,sps1,X,sps1,Y);
                        fprintf(fid2,'%s',a{u}); fprintf(fid2,'\r\n');
                    end
                end
                %%%%%%%% Nodes Section %%%%%%%%
                if (cnt==2 && (strcmp(a{u},'[TANKS]') || strcmp(a{u},'[JUNCTIONS]') || strcmp(a{u},'[RESERVOIRS]') || strcmp(a{u},'[DEMANDS]')))
                    if sjunction_check==0 && typecode==0 && strcmp(a{u},'[RESERVOIRS]')
                        fprintf(fid2,'[JUNCTIONS]');
                        fprintf(fid2, '\n%s%s%d%s%s\n',newID,sps1,newElevation,sps1,sps1);
                    end
                    if sreservoir_check==0 && typecode==1 && strcmp(a{u},'[TANKS]')
                        fprintf(fid2,'[RESERVOIRS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s\n', newID,sps1,newElevation,sps1,'',sps1);
                    end
                    if stank_check==0 && typecode==2 && strcmp(a{u},'[PIPES]')
                        fprintf(fid2,'[TANKS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d\n', newID,sps1,newElevation,sps1,Initlevel,sps1,MinLevel,sps1,...
                        MaxLevel,sps1,Diameter,sps1,MinVol);
                    end
                    fprintf(fid2,'%s',a{u});
                    %%%%%%%% Jynctions Section %%%%%%%%
                    if typecode==0 && strcmp(a{u},'[JUNCTIONS]')
                        fprintf(fid2, '\n%s%s%d', newID,sps1,newElevation);
                    end
                    if typecode==0 && strcmp(a{u},'[DEMANDS]')
                        fprintf(fid2, '\n%s%s%d', newID,sps1,newBaseDemand);
                    end
                    %%%%%%%% Reservoirs Section %%%%%%%%
                    if typecode==1 && strcmp(a{u},'[RESERVOIRS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s', newID,sps1,newElevation);
                    end
                    %%%%%%%% Tanks Section %%%%%%%%
                    if typecode==2 && strcmp(a{u},'[TANKS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d', newID,sps1,newElevation,sps1,Initlevel,sps1,MinLevel,sps1,...
                        MaxLevel,sps1,Diameter,sps1,MinVol);
                    end
                elseif isempty(a{u})
                else
                    fprintf(fid2,'%s%s',a{u},sps1);
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
                    fprintf(fid2,'\r\n');
                    fprintf(fid2, '%s%s%d%s%d', newID,sps1,X,sps1,Y);
                    gg=0; onetime=0;
                end
            end
            if qualch==1 && qq==1
                qualch=0;
            end
            fprintf(fid2,'\n');
        end
    end
    fclose(fid2);
    if obj.Bin==1
        errcode=closeOpenNetwork(obj);
    end
end
function errcode=addLink(obj,typecode,newLink,fromNode,toNode,varargin)
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
if typecode==1 && nargin>5
    plength=varargin{1};
    pdiameter=varargin{2};
    proughness=varargin{3};
elseif typecode==2
    if ~isnumeric(varargin{1})
        curveID=varargin{1};
    else
        power=varargin{1};
        curveID={};
    end
elseif typecode>2
    if typecode==0, type_valv = 'CVPIPE';  end
    if typecode==3, type_valv = 'PRV';     end
    if typecode==4, type_valv = 'PSV';     end
    if typecode==5, type_valv = 'PBV';     end
    if typecode==6, type_valv = 'FCV';     end
    if typecode==7, type_valv = 'TCV';     end
    if typecode==8, type_valv = 'GPV';     end
    if typecode~=1 && typecode~=2
        typecode=3;
    end
    vdiameter=varargin{1};
    vsetting=varargin{2};    
end
[errcode]=addLinkWarnings(obj,typecode,newLink,toNode);
crvs = obj.getBinCurvesInfo;
% Open and read inpname
% Read all file and save in variable info
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
% Add pipe
nn=0;sps=blanks(10);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            cnt=bracketsCheck(a{u});
            if (cnt==2 && strcmp(a{u},'[PIPES]') && nn==0 && typecode==1)
                fprintf(fid2,'%s',a{u});
                fprintf(fid2, '\n%s%s%s%s%s%s%d%s%d%s%d%s%d',newLink,sps,fromNode,sps,...
                    toNode,sps,plength,sps,pdiameter,sps,proughness,sps,0);
                
            elseif (cnt==2 && strcmp(a{u},'[PUMPS]') && nn==0 && typecode==2)
                if ~isempty(curveID)
                    if isempty(char(crvs.BinCurveNameID))
                        s = sprintf('No head curve supplied for pump %s.',newLink);
                        warning(s);
                        return;
                    end
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%s%s%s',newLink,sps,fromNode,sps,...
                        toNode,sps,'HEAD',sps,curveID);
                else
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%s%s%.2f',newLink,sps,fromNode,sps,...
                        toNode,sps,'POWER',sps,power);
                end
            elseif typecode==3 && strcmp(a{u},'[VALVES]')
                fprintf(fid2,'%s',a{u});
                fprintf(fid2, '\n%s%s%s%s%s%s%d%s%s%s%d',newLink,sps,fromNode,sps,...
                    toNode,sps,vdiameter,sps,type_valv,sps,vsetting);
                nn=1;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [errcode] = rmNode(obj,NodeID)
% Remove node from the network.
% Check if id new already exists
nodes = obj.getBinNodesInfo;errcode=0;
if length(nodes.BinNodeNameID)==0
    return;
end
i=1;
while i<length(nodes.BinNodeNameID)+1
    exists(i) = strcmp(NodeID,char(nodes.BinNodeNameID(i)));
    i=i+1;
end
if (sum(exists)==0)
    s = sprintf('There is no such object in the network.');
    warning(s);
    errcode=-1;
    return;
end
if sum(strcmp(NodeID,nodes.BinNodeReservoirNameID)) || sum(strcmp(NodeID,nodes.BinNodeTankNameID))
    if (length(char(nodes.BinNodeReservoirNameID))+length(char(nodes.BinNodeTankNameID))-1)==0;
        warning('This tank/reservoir has not removed.');
        errcode=-1;
        return;
    end
end
% Get links which delete with function Remove Link
links = obj.getBinLinksInfo;
a=strcmp(links.BinLinkFromNode,NodeID);
linkindex1=find(a);
b=strcmp(links.BinLinkToNode,NodeID);
linkindex2=find(b);
linkindex12=[linkindex1 linkindex2];
checklinks_index=unique(linkindex12);
checklinks=links.BinLinkNameID(checklinks_index);
obj.removeBinControlNodeID(NodeID);% Remove control, code 0(NODE)
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
out=0; sps=blanks(10);
for t = 1:length(info)
    c = info{t};
    a = regexp(c,'\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;x=0;xx=0;q=0;
        while u < length(a)+1
            if isempty(a{u}) && (x==0)
                u=u+1; x=1;xx=1;
                if u==length(a)+1
                    break;
                end
            end
            if strcmp(a{u},'[PIPES]'), out=1; end
            if strcmp(a{u},'[DEMANDS]'), out=0; end %out=0; delete line
            if strcmp(a{u},'[PATTERNS]'), out=1; end 
            if strcmp(a{u},'[QUALITY]'), out=0; end
            if strcmp(a{u},'[SOURCES]'), out=1; end
            if strcmp(a{u},'[MIXING]'), out=0; end
            if strcmp(a{u},'[COORDINATES]'), out=0; end
            if strcmp(a{u},NodeID) && q~=1 && out==0
                if xx==1 || strcmp(a{u},NodeID)
                    u=length(a)+1;
                end
            else
                q=1;
                fprintf(fid2,'%s%s',a{u},sps);
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
% Remove links
if length(checklinks)
    for i=1:length(checklinks)
        warn1(i)=obj.removeBinLinkID(checklinks{i});
    end
else
    warn1=1;
end
% Find who other id must be delete
remove_link={''};
remove_link_index = zeros(1,length(links.BinLinkFromNode));
if length(checklinks)
    for i=1:length(checklinks)
        remove_link(i)=checklinks(i);
        remove_link_index(i)=i;
        s=sprintf('Removed link:%s',char(remove_link(i)));
        warning(s);
    end
else
    warn1=sum(warn1)+1;
end
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function errcode=rmRulesControl(obj,type,id)
% Remove control from the network.
exists=0;errcode=0;
exists1=0;
rulescontrols = obj.getBinRulesControlsInfo;
if type==1
    if length(rulescontrols.BinRulesControlLinksID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);errcode=-1;
        return;
    end
    for i=1:length(rulescontrols.BinRulesControlLinksID) 
        if type==1
            exists(i) = strcmp(rulescontrols.BinRulesControlLinksID{i}{length(rulescontrols.BinRulesControlLinksID{1})},char(id));
        end
    end
end
if type==0
    if length(rulescontrols.BinRulesControlNodesID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);errcode=-1;
        return;
    end
    for i=1:length(rulescontrols.BinRulesControlNodesID)
        if type==0
            exists1(i) = strcmp( rulescontrols.BinRulesControlNodesID{i}{length(rulescontrols.BinRulesControlNodesID{1})},char(id));
        else
            warning('Type is NODE(0) or LINK(1)');errcode=-1;
            return;
        end
    end
end
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
e=0;n=0;kk=1;sps=blanks(15);tt=0;
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            rr = regexp(a,'\w*[\w*]\w*','split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets,'[');
            ch2 = strcmp(check_brackets,']');
            
            if strcmp(a{u},'[RULES]')
                fprintf(fid2,'%s',a{u});
                n=1;
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1), break; end
                e=1;
            end
            if strcmp(a{u},'[END]'),  e=1; fprintf(fid2,'%s',a{u});break;   end
            
            if n==1 && e==0 && kk==1
                if strcmp(a{u},'[RULES]'), break; end
                if isempty(a{u})
                elseif strfind(a{u},';')
                    break;
                else
                    if type==1
                        if tt==1
                            if strcmp(upper(a{1}),'PRIORITY') && type==1
                                break;
                            end
                        end
                        tt = strcmp(a{u+2},id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2,'%s%s',a{u},sps);
                        end
                    elseif type==0
                        tt = strcmp(a{u+5},id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2,'%s%s',a{u},sps);
                        end
                    end
                end
            else
                if isempty(a{u})
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');kk=1;
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function errcode=rmControl(obj,type,id)
% Remove control from the network.
exists=0;errcode=0;
exists1=0;
controls = obj.getBinControlsInfo;
if type==1
    if length(controls.BinControlLinksID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);errcode=-1;
        return;
    end
    for i=1:length(controls.BinControlLinksID) 
        if type==1
            exists(i) = strcmp(controls.BinControlLinksID{i},char(id));
        end
    end
end
if type==0
    if length(controls.BinControlNodesID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);errcode=-1;
        return;
    end
    for i=1:length(controls.BinControlNodesID)
        if type==0
            exists1(i) = strcmp(controls.BinControlNodesID{i},char(id));
        else
            warning('Type is NODE(0) or LINK(1)');errcode=-1;
            return;
        end
    end
end
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
e=0;n=0;kk=1;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            rr = regexp(a,'\w*[\w*]\w*','split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets,'[');
            ch2 = strcmp(check_brackets,']');
            
            if strcmp(a{u},'[CONTROLS]')
                fprintf(fid2,'%s',a{u});
                n=1;
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1), break; end
                e=1;
            end
            if strcmp(a{u},'[END]'),  e=1; fprintf(fid2,'%s',a{u});break;   end
            
            if n==1 && e==0 && kk==1
                if strcmp(a{u},'[CONTROLS]'), break; end
                if isempty(a{u})
                elseif strfind(a{u},';')
                    break;
                else
                    if type==1
                        tt = strcmp(a{u+1},id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2,'%s%s',a{u},sps);
                        end
                    elseif type==0
                        tt = strcmp(a{u+5},id); kk=0;
                        if tt==1
                            break;
                        else
                            fprintf(fid2,'%s%s',a{u},sps);
                        end
                    end
                end
            else
                if isempty(a{u})
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');kk=1;
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [errcode] = rmLink(obj,LinkID)
% Remove link from the network.
% Check if id new already exists
links = obj.getBinLinksInfo;errcode=0;
if length(links.BinLinkNameID)==0
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
countLinks=length(links.BinLinkNameID);
i=1;
while i<countLinks+1
    exists(i) = strcmp(LinkID,char(links.BinLinkNameID(i)));
    if exists(i)==1
        index_rmvlink=i;
    end
    i=i+1;
end
if (sum(exists)~=1)
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
nodes = obj.getBinNodesInfo;
from_node = links.BinLinkFromNode(index_rmvlink);
r = strcmp(nodes.BinNodeNameID,from_node);
if sum(r)==0, from_node=''; end
to_node = links.BinLinkToNode(index_rmvlink);
r = strcmp(nodes.BinNodeNameID,to_node);
if sum(r)==0, to_node=''; end
% Remove control, code 1(LINK)
obj.removeBinControlLinkID(LinkID);

[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');

% section [JUNCTIONS]
out=0;YY=0;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;x=0;xx=0;q=0;
        while u < length(a)+1
            if strcmp(a{u},'[PIPES]'), YY=1;end
            if YY==1
                if isempty(a{u}) && (x==0)
                    u=u+1; x=1;xx=1;
                    if u==length(a)+1
                        break
                    end
                end
                if strcmp(a{u},'[TAGS]'), out=1; end
                if strcmp(a{u},'[STATUS]'), out=1; end
                if strcmp(a{u},'[DEMANDS]'), out=1; end
                if strcmp(a{u},'[PATTERNS]'), out=1; end
                
                if strcmp(a{u},LinkID) && q~=1 && out==0
                    if xx==1 || strcmp(a{1},LinkID)
                        u=length(a)+1;
                    end
                else
                    q=1;
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            else
                if isempty(a{u})
                    u=u+1;
                    if u==length(a)+1
                        break
                    end
                end
                fprintf(fid2,'%s%s',a{u},sps);
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
% Get nodes which delete with function Remove Node
links = obj.getBinLinksInfo;
for i=1:length(links.BinLinkNameID)
    t(i) = strcmp(from_node,char(links.BinLinkFromNode(i)));
    tt(i) = strcmp(to_node,char(links.BinLinkToNode(i)));
    ttt(i) = strcmp(to_node,char(links.BinLinkFromNode(i)));
    tttt(i) = strcmp(from_node,char(links.BinLinkToNode(i)));
    i=i+1;
end
if sum(t)+sum(tttt)==0 || sum(tt)+sum(ttt)==0
    if isempty(char(from_node)) && isempty(char(to_node))
        warning('Call function Removenode or Addlink.');
    end
end
if sum(t)+sum(tttt)==0
    if ~isempty(char(from_node))
        s = sprintf('Node %s disconnected.',char(from_node));
        warning(s);
    end
end
if sum(tt)+sum(ttt)==0
    if ~isempty(char(to_node))
        s = sprintf('Node %s disconnected.',char(to_node));
        warning(s);
    end
end
if sum(t)+sum(tttt)==0 || sum(tt)+sum(ttt)==0
    if ~isempty(char(from_node)) || ~isempty(char(to_node))
        if ~sum(strcmp(from_node,nodes.BinNodeReservoirNameID)) || ~sum(strcmp(from_node,nodes.BinNodeTankNameID))...
                || ~sum(strcmp(to_node,nodes.BinNodeReservoirNameID)) || ~sum(strcmp(to_node,nodes.BinNodeTankNameID))
            errcode=0;
        else
            errcode=-1;
        end
    end
end
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [errcode]=addNewControl(obj,x,status,y_t_c,param,z,varargin)
% syntax
if (nargin==6)
    syntax = sprintf('LINK     %s     %s     IF     NODE     %s     %s     %d',x,status,y_t_c,param,z);
elseif (nargin==5)
    syntax = sprintf('LINK     %s     %s     AT     CLOCKTIME     %s     %s',x,status,y_t_c,param);
elseif (nargin==4)
    syntax = sprintf('LINK     %s     %s     AT     TIME     %s',x,status,y_t_c);
end
if (nargin==6)
    % Check if id new already exists
    nodes = obj.getBinNodesInfo;
    if length(char(nodes.BinNodeNameID))==0
        return
    end
    i=1;
    while i<length(char(nodes.BinNodeNameID))+1
        exists(i) = strcmp(y_t_c,char(nodes.BinNodeNameID(i)));
        i=i+1;
    end
    
    if (sum(exists)~=1)
        s = sprintf('There is no such object in the network.');
        warning(s);errcode=-1;
        return;
    end
end
% Check if id new already exists
links = obj.getBinLinksInfo; errcode=0;
if length(char(links.BinLinkNameID))==0
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
i=1;
while i<length(char(links.BinLinkNameID))+1
    exists(i) = strcmp(x,char(links.BinLinkNameID(i)));
    i=i+1;
end
if (sum(exists)~=1)
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
type_n='[CONTROLS]';
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
noo=0;s=0;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u},type_n);
                fprintf(fid2,'[CONTROLS]');
                s=1; break;
            end
            rr = regexp(a,'\w*[\w*]\w*','split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets,'[');
            ch2 = strcmp(check_brackets,']');
            if (ch1(1)==1 && ch2(2)==1 && (s==1) && noo==0)
                fprintf(fid2, '%s',syntax);
                fprintf(fid2,'\r\n');
                fprintf(fid2,'\r\n');
                fprintf(fid2,'%s',a{u});
                fprintf(fid2,'\r\n');
                noo=1;
            elseif isempty(a{u}) && noo==0
            else
                if isempty(a{u}) && noo==1
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [errcode]=rmCurveID(obj,CurveID,varargin)
% Check if id new already exists
valueC=obj.getBinCurvesInfo;errcode=0;
if length(valueC.BinCurveNameID)==0
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
i=1;
while i<length(valueC.BinCurveNameID)+1
    exists(i) = strcmp(CurveID,char(valueC.BinCurveNameID(i)));
    i=i+1;
end
if (sum(exists)==0)
    s = sprintf('There is no such object in the network.');
    warning(s);errcode=-1;
    return;
end
value=obj.getBinLinksInfo;
i=1;
while (i<length(value.BinLinkPumpCurveNameID)+1)
    p = strcmp(CurveID,value.BinLinkPumpCurveNameID{i});
    if p==1
        s = sprintf('Pump %s refers to undefined curve.',value.BinLinkPumpNameID{i});
        warning(s);
    end
    i=i+1;
end
% Open and read inpname
% Read all file and save in variable info
[~,info] = obj.readInpFile;
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
e=0;n=0;sps=blanks(15);
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
    elseif isempty(c)
    else
        u=1;
        while u < length(a)+1
            rr = regexp(a,'\w*[\w*]\w*','split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets,'[');
            ch2 = strcmp(check_brackets,']');
            if strcmp(a{u},'[CURVES]')
                fprintf(fid2,'%s',a{u});
                n=1;
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1), break; end
                e=1;
            end
            if strcmp(a{u},'[END]'), e=1; fprintf(fid2,'%s',a{u});break;end
            if n==1 && e==0
                if strcmp(a{u},'[CURVES]'), break; end
                if isempty(a{u})
                    u=u+1;continue;
                elseif strfind(a{u},';') 
                    ee=regexp(c,'\w*EFFICIENCY*\w','match');
                    nn=regexp(c,'\w*VOLUME*\w','match');
                    kk=regexp(c,'\w*HEADLOSS*\w','match');
                    if length(strcmp(ee,'EFFICIENCY')) || length(strcmp(nn,'VOLUME')) || length(strcmp(kk,'HEADLOSS')) || length(strcmp(a{1},';PUMP:'))
                        fprintf(fid2,'%s%s',a{u},sps);
                    else
                        break;
                    end
                else
                    tt = strcmp(a{u},CurveID);
                    if tt==1
                        u = length(a)+1;
                    else
                        fprintf(fid2,'%s%s',a{u},sps);
                    end
                end
            else
                if isempty(a{u})
                else
                    if strcmp(a{u},'[CURVES]'), break; end
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function [errcode]=Options(obj,newFlowUnits,headloss,varargin)
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
value=obj.getBinOptionsInfo;errcode=0;
previousFlowUnits=value.BinLinkFlowUnits;
newUScustomary=0;
newSImetric=0;
switch newFlowUnits
    case 'CFS'
        newUScustomary=1;
    case 'GPM'
        newUScustomary=1;
    case 'MGD'
        newUScustomary=1;
    case 'IMGD'
        newUScustomary=1;
    case 'AFD'
        newUScustomary=1;
    case 'LPS'
        newSImetric=1;
    case 'LPM'
        newSImetric=1;
    case 'MLD'
        newSImetric=1;
    case 'CMH'
        newSImetric=1;
    case 'CMD'
        newSImetric=1;
end
if newUScustomary==value.BinUScustomary
    changes=0; newUScustomary=0;
    newSImetric=0; % feet to feet
elseif newSImetric==value.BinSImetric
    changes=1; newUScustomary=0;
    newSImetric=0; % meter to meter
elseif value.BinUScustomary==1 && newUScustomary==0
    changes=1; % feet to meter or cubic feet to cubic meter
elseif value.BinUScustomary==0 && newUScustomary==1
    changes=2; % meter to feet or cubic meter to cubic feet
end
Units=newUScustomary+newSImetric;
variables=who;nheadl=0;
if ~sum(strcmp('headloss',variables))
    headloss=value.BinOptionsHeadloss;
    nheadl=1;
end
nodes = obj.getBinNodeNameID;
links = obj.getBinLinkNameID;
controls = obj.getBinControlsInfo;
curves = obj.getBinCurvesInfo;
rules=obj.getBinRulesControlsInfo;

[info] = readAllFile(obj.Bintempfile);
fid2 = fopen([obj.pathfile,obj.Bintempfile],'w');
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
            if strcmp(a{u},'[JUNCTIONS]') && Units
                fprintf(fid2,'[JUNCTIONS]');
                sect=1;
                break;
            elseif strcmp(a{u},'[RESERVOIRS]') && Units
                fprintf(fid2,'[RESERVOIRS]');
                nn=0;pp=1;
                sect=2;
                break;
            elseif strcmp(a{u},'[TANKS]') && Units
                fprintf(fid2,'[TANKS]');
                nn=0;pp=1;
                sect=3;
                break;
            elseif strcmp(a{u},'[PIPES]') && Units
                fprintf(fid2,'[PIPES]');
                nn=0;pp=1;
                sect=4;
                break;
            elseif strcmp(a{u},'[PUMPS]') && Units
                fprintf(fid2,'[PUMPS]');
                nn=0;pp=1;
                sect=5;
                break;
            elseif strcmp(a{u},'[VALVES]') && Units
                fprintf(fid2,'[VALVES]');
                nn=0;pp=1;
                sect=6;
                break;
            elseif strcmp(a{u},'[DEMANDS]') && ((Units || ~changes) && nheadl)
                fprintf(fid2,'[DEMANDS]');
                nn=0;pp=1;
                sect=7;
                break;            
            elseif strcmp(a{u},'[EMITTERS]') && ((Units || ~changes) && nheadl)
                fprintf(fid2,'[EMITTERS]');
                nn=0;pp=1;
                sect=8;
                break;
            elseif strcmp(a{u},'[STATUS]') 
                fprintf(fid2,'[STATUS]');
                nn=0;pp=1;
                sect=9;
                break;
            elseif strcmp(a{u},'[PATTERNS]') 
                fprintf(fid2,'[PATTERNS]');
                nn=1;
                sect=10;
                break;
            elseif strcmp(a{u},'[CURVES]') && ((Units || ~changes) && nheadl)
                fprintf(fid2,'[CURVES]');
                nn=0;pp=1;ww=1;
                sect=11;
                break;
            elseif strcmp(a{u},'[CONTROLS]') && Units
                fprintf(fid2,'[CONTROLS]');
                nn=0;pp=1;
                sect=12;
                break;
            elseif strcmp(a{u},'[RULES]') && Units
                fprintf(fid2,'[RULES]');
                nn=0;pp=1;
                sect=13;
                break;                
            elseif strcmp(a{u},'[OPTIONS]')
                fprintf(fid2,'[OPTIONS]');
                sect=14;nn=0;
                break;
            end
                % section [JUNCTIONS]
            if (sect==1) && (nn==0)
                mm=1;
                if pp<length(char(nodes.BinNodeJunctionNameID))+1
                    if strcmp(a{mm},nodes.BinNodeJunctionNameID{pp})
                        pp=pp+1;
                        fprintf(fid2,'%s%s',char(a{mm}),sps);
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                        end
                        if length(a)>2
                            mm=2;
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
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
                    if strcmp(a{mm},nodes.BinNodeReservoirNameID{pp})
                        pp=pp+1;
                        fprintf(fid2,'%s%s',char(a{mm}),sps);
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
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
                    if strcmp(a{mm},nodes.BinNodeTankNameID{pp})
                        pp=pp+1;
                        fprintf(fid2,'%s%s',char(a{mm}),sps);
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+4})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+5})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+6})*0.02831685),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+4})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+5})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+6})*35.3147),sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2,'%s%s',char(a{1}),sps);
                end
                break;
                % section [PIPES]
            elseif (sect==4) && (nn==0)
                mm=1;
                if pp<length(char(links.BinLinkPipeNameID))+1
                    if strcmp(a{mm},links.BinLinkPipeNameID{pp})
                        pp=pp+1;
                        for mm=mm:mm+2
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*25.4),sps);
                            if nheadl
                                if strcmp('D-W',value.BinOptionsHeadloss)
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*0.3048),sps);
                                    mm=7;
                                else
                                    mm=6;
                                end
                            end
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.03937007874),sps);
                            if nheadl
                                if strcmp('D-W',value.BinOptionsHeadloss)
                                    mm=7;
                                else
                                    mm=6;
                                end
                            end
                        end
                        for mm=mm:length(a)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2,'%s%s',char(a{1}),sps);
                end
                break;
            % section [PUMPS]
            elseif (sect==5) && (nn==0)
                mm=1; 
                if pp<length(char(links.BinLinkPumpNameID))+1
                    if strcmp(a{mm},links.BinLinkPumpNameID{pp})
                        pp=pp+1;
                        for mm=mm:length(a)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        power=regexp(c,'POWER','match');
                        if strcmpi(upper(power),'POWER')
                            mm=mm-1;
                            if changes==1 
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.745699882507324),sps);
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})/0.745699882507324),sps);
                            end
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2,'%s%s',char(a{1}),sps);
                end
                break;
            % section [VALVES]
            elseif (sect==6) && (nn==0)
                mm=1;
                if pp<length(char(links.BinLinkValveNameID))+1
                    if strcmp(a{mm},links.BinLinkValveNameID{pp})
                        pp=pp+1;
                        for mm=mm:(mm+2)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        prv=regexp(c,'PRV','match');if isempty(prv), prv=0; end
                        psv=regexp(c,'PSV','match');if isempty(psv), psv=0; end
                        pbv=regexp(c,'PBV','match');if isempty(pbv), pbv=0; end
                        fcv=regexp(c,'FCV','match');if isempty(fcv), fcv=0; end
%                         tcv=regexp(c,'TCV','match');if isempty(tcv), tcv=0; end
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*25.4),sps);
                            fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                            if strcmpi(upper(prv),'PRV') || strcmpi(upper(psv),'PSV') || strcmpi(upper(pbv),'PBV') %|| strcmpi(upper(tcv),'TCV')
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*0.3048),sps);
                            elseif strcmpi(upper(fcv),'FCV')
                                setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            end
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.03937007874),sps);
                            fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                            if strcmpi(upper(prv),'PRV') || strcmpi(upper(psv),'PSV') || strcmpi(upper(pbv),'PBV') %|| strcmpi(upper(tcv),'TCV')
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})/0.3048),sps);
                            elseif strcmpi(upper(fcv),'FCV')
                                setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            end
                        end
                        for mm=(mm+4):length(a)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2,'%s%s',char(a{1}),sps);
                end
                break;
                % section [DEMANDS]
            elseif (sect==7) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if pp<length(char(nodes.BinNodeJunctionNameID))+1
                        if strcmp(a{mm},nodes.BinNodeJunctionNameID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            if length(a)>2
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps);
                    end
                end
                break;
                % section [EMITTERS]
            elseif (sect==8) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if sum(strcmp(a{mm},nodes.BinNodeJunctionNameID))
                        fprintf(fid2,'%s%s',char(a{mm}),sps);
                        setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                    end
                end
                break;
                % section [STATUS]
            elseif (sect==9) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if pp<length(char(nodes.BinNodeJunctionNameID))+1
                        if strcmp(a{mm},nodes.BinNodeJunctionNameID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            if length(a)==2
                                fprintf(fid2,'%s%s',char(a{mm+1}),sps);
                            end
                            if length(a)>2
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps);
                    end
                end
                break;
                
                % section [CURVES]
            elseif (sect==11) && (nn==0)
                mm=1;
                if strfind(c,';ID')
                    break;
                end
                if pp<length(curves.BinCurveAllLines)+1 && ~isempty(char(a)) % PUMP % EFFICIENCY % VOLUME
                    if ww<length(curves.BinCTypes)+1
                        if curves.BinCTypes(ww)==0
                            if strfind(c,';PUMP:')
                                fprintf(fid2,c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm);
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps);
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*3.281),sps);
                            else
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                            end
                            break;
                        elseif curves.BintypeCurve(ww)==1
                            ee=regexp(c,'\w*EFFICIENCY*\w','match');
                            if length(strcmp(ee,'EFFICIENCY'))
                                fprintf(fid2,c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            fprintf(fid2,'%s%s',char(a{mm+2}),sps);
                        elseif curves.BintypeCurve(ww)==2
                            gg=regexp(c,'\w*VOLUME*\w','match');
                            if length(strcmp(gg,'VOLUME'))
                                fprintf(fid2,c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.831685e-02),sps);
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps);
                            else
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})),sps);
                            end
                        elseif curves.BintypeCurve(ww)==3 % HEADLOSS
                            kk=regexp(c,'\w*HEADLOSS*\w','match');
                            if length(strcmp(kk,'HEADLOSS'))
                                fprintf(fid2,c);break;
                            end
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                            else
                                fprintf(fid2,'%s%s',char(a{mm+1}),sps);
                            end
                            mm=mm+1;
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                        end
                        ww=ww+1;
                    end
                    pp=pp+1;
                else
                    if ~(ww<length(curves.BinCTypes)+1), nn=1; end
                end
                if ~isempty(regexp(a{mm},'[\w]*','match'))
                    nn=1;
                    fprintf(fid2,'%s%s',char(c),sps);
                end
                break;
                % section [CONTROLS]
            elseif (sect==12) && (nn==0)
                e=regexp(a,';','match');
                if length(e)>0
                    if ~isempty(e{1})
                        break;
                    end
                end
                mm=1;
                if pp<length(controls.BinControlsInfo)+1
                    pp=pp+1;
                    if length(a)>7
                        if strcmp(upper(a{mm+6}),'BELOW') || strcmp(upper(a{mm+6}),'ABOVE')
                            for mm=mm:(mm+6)
                                fprintf(fid2,'%s%s',char(a{mm}),sps);
                            end
                            v=find(strcmp(a{2},obj.BinLinkNameID));
                            index=obj.getBinLinkIndex(obj.BinLinkNameID(v));
                            if ~strcmp(obj.BinLinkType(index),'TCV') && ~strcmp(obj.BinLinkType(index),'GPV')
                                if changes==1
                                    if strcmp(obj.BinLinkType(index),'FCV')
                                        setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                    else
                                        fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                                    end
                                elseif changes==2
                                    if strcmp(obj.BinLinkType(index),'FCV')
                                        setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                    else
                                        fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                                    end
                                end
                            else
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
                            end
                        else
                            for mm=mm:length(a)
                                fprintf(fid2,'%s%s',char(a{mm}),sps);
                            end
                        end
                    else
                        for mm=mm:length(a)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                    end
                else
                    nn=1;
                    fprintf(fid2,'%s%s',char(a{1}),sps);
                end
                break;
                % section [RULES]
            elseif (sect==13) && (nn==0)
                mm=1;
                if pp<rules.BinRulesCount+1
                    pp=pp+1;
                    if strcmp(upper(regexp(cell2mat(a),'\s*LEVEL*','match')),'LEVEL') %|| 
                        for mm=mm:(mm+4)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                        end
                    elseif strcmp(upper(regexp(cell2mat(a),'\s*HEAD*','match')),'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps);
                        end
                    elseif strcmp(upper(regexp(cell2mat(a),'\s*DEMAND*','match')),'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                    elseif strcmp(upper(regexp(cell2mat(a),'\s*PRESSURE*','match')),'HEAD')
                        for mm=mm:(mm+4)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                        if changes==1
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})/1.422),sps);
                        elseif changes==2
                            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.422),sps);
                        end
                    else
                        for mm=mm:length(a)
                            fprintf(fid2,'%s%s',char(a{mm}),sps);
                        end
                    end
                else
                    nn=1;
                    for mm=mm:length(a)
                        fprintf(fid2,'%s%s',char(a{mm}),sps);
                    end
                end
                break;                
                % section [OPTIONS]
            elseif (sect==14) && (nn==0)
                mm=1;
                if strcmp(upper(a{mm}),'UNITS')
                    fprintf(fid2,'%s%s',char(a{mm}),sps);
                    if nheadl
                        fprintf(fid2,'%s%s',char(newFlowUnits),sps);nn=1;
                    else
                        fprintf(fid2,'%s%s',char(previousFlowUnits),sps);
                    end
                elseif strcmp(upper(a{mm}),'HEADLOSS')
                    fprintf(fid2,'%s%s',char(a{mm}),sps);
                    fprintf(fid2,'%s%s',char(headloss),sps);
                    nn=1;
                else
                    fprintf(fid2,c);
                end
                break;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps);
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose(fid2);
if obj.Bin==1
    errcode=closeOpenNetwork(obj);
end
end
function errcode=closeOpenNetwork(obj)
    obj.closeNetwork;  %Close input file 
    errcode=ENopen([obj.pathfile,obj.Bintempfile],[obj.pathfile,[obj.Bintempfile(1:end-4),'.txt']],[obj.pathfile,[obj.Bintempfile(1:end-4),'.bin']],obj.libepanet);
end
function setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
if strcmp(previousFlowUnits,'GPM')
    switch newFlowUnits %(GPM)
        case 'CFS' 
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00222816399286988),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00144),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00119905),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.004419191),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0630902),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.785412),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.005450993),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2271247),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*5.450993),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'CFS')
    switch newFlowUnits
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*448.8312),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.6463169),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.5381711),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.983471),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*28.31685),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1899.011),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.446576),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*101.9406),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2446.576),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'MGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.547229),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*694.4445),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.8326738),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.068883),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*43.81264),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2628.758),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.785412),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*157.7255),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3785.412),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'IMGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.858145),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*833.9936),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.200951),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.685577),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*52.61681),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3157.008),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.546092),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*189.4205),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4546.092),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'AFD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.5041667),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*226.2857),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3258514),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.271328),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*14.27641),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*856.5846),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.233482),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*51.39508),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1233.482),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'LPS')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.03531466),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*15.85032),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.02282446),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01900533),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.07004562),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*60),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0864),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.6),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*86.4),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'LPM')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0005885777),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.264172),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0003804078),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0003167556),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0011674272),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01666667),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00144),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.06),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.44),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'MLD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.4087345),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*183.4528),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.264172),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2199692),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.8107132),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*11.57407),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*694.4445),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*41.66667),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1000),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'CMH')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.009809635),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.402868),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.006340129),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00527926),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01945712),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2777778),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*16.66667),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.024),sps);
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*24),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
elseif strcmp(previousFlowUnits,'CMD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0004087345),sps);
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.1834528),sps);
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.000264172),sps);
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0002199692),sps);
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0008107132),sps);
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01157407),sps);
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.6944444),sps);
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.001),sps);
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.04166667),sps);
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps);
    end
end
end
function [fid,binfile,msg] = makebatfile(obj)
    binfile=[obj.Bintempfile(1:end-4),'.bin'];
    if exist(binfile)==2, fclose all; delete(binfile); end
    if strcmp(computer('arch'),'win64')
            folder='64bit';
        r = sprintf('%s\\%s\\epanet2d.exe %s %s %s',pwd,folder,obj.Bintempfile,[obj.Bintempfile(1:end-4),'.txt'],binfile);
    elseif strcmp(computer('arch'),'win32')
            folder='32bit';
        r = sprintf('%s\\%s\\epanet2d.exe %s %s %s',pwd,folder,obj.Bintempfile,[obj.Bintempfile(1:end-4),'.txt'],binfile);
    else
        r = sprintf('./runcode2 %s %s %s',obj.Bintempfile,[obj.Bintempfile(1:end-4),'.txt'],binfile);
    end
    f=fopen([obj.pathfile,'Simulate.bat'],'w');
    try fprintf(f,'%s \n',r); fclose(f); catch e; end
    [~,msg]=system([obj.pathfile,'Simulate.bat']);
    delete([obj.pathfile,'Simulate.bat']);
    binfile=[obj.pathfile,obj.Bintempfile(1:end-4),'.bin'];
    fid = fopen(binfile,'r');
end
function value = getBinComputedTimeSeries(obj,indParam,varargin)
    [fid,binfile,msg] = makebatfile(obj);
    value=[];
    if fid~=-1
        data = fread(fid,'int32');
        fclose(fid);
        BinNodeCount=data(3);
        BinNodeResTankCount=data(4);
        BinLinkCount=data(5);
        BinLinkPumpCount=data(6);
        NumberReportingPeriods = data(end-2);
        if indParam==27
            value = NumberReportingPeriods; return;
        end
        if indParam==28
            value = data(15); return; % simulation duration
        end
        fid1 = fopen(binfile, 'r');
        fread(fid1, 15, 'uint32');
        fread(fid1, 808, '*char');
        fread(fid1, 4, 'uint32');
        fread(fid1, 32*BinNodeCount+32*BinLinkCount, '*char'); % error NODES*32
        fread(fid1, BinLinkCount*3, 'uint32');
        fread(fid1, BinNodeResTankCount*2, 'uint32'); % error
        if indParam==1
            value =fread(fid1, BinNodeCount, 'float')';return; % ElevationEachNode
        elseif indParam==2
            fread(fid1, BinNodeCount, 'float');
            value =fread(fid1, BinLinkCount, 'float')';return; % LengthEachLink
        elseif indParam==3
            fread(fid1, BinNodeCount+BinLinkCount, 'float');
            value =fread(fid1, BinLinkCount, 'float')';return; % DiameterEachLink
        elseif indParam==4
            fread(fid1, BinNodeCount+BinLinkCount*2, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % PumpIndexListLinks
        elseif indParam==5
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % PumpUtilization
        elseif indParam==6
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*2, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % AverageEfficiency
        elseif indParam==7
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*3, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % AverageKwattsOrMillionGallons
        elseif indParam==8
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*4, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % AverageKwatts
        elseif indParam==9
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*5, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % PeakKwatts
        elseif indParam==10
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*6, 'float');
            value =fread(fid1, BinLinkPumpCount, 'float')';return; % AverageCostPerDay
        end
        if indParam>10
            fread(fid1, BinNodeCount+BinLinkCount*2+BinLinkPumpCount*7+1, 'float');
        end
        for i=1:NumberReportingPeriods
            if indParam==11
                value(:,i) = fread(fid1, BinNodeCount, 'float')'; % nodeDemand
                fread(fid1, BinNodeCount*3, 'float');
                fread(fid1, BinLinkCount*8, 'float');
            elseif indParam==12
                fread(fid1, BinNodeCount, 'float');
                value(:,i) = fread(fid1, BinNodeCount, 'float')'; % nodeHead
                fread(fid1, BinNodeCount*2, 'float');
                fread(fid1, BinLinkCount*8, 'float');
            elseif indParam==13
                fread(fid1, BinNodeCount*2, 'float');
                value(:,i) = fread(fid1, BinNodeCount, 'float')'; % nodePressure
                fread(fid1, BinNodeCount, 'float');
                fread(fid1, BinLinkCount*8, 'float');
            elseif indParam==14
                fread(fid1, BinNodeCount*3, 'float');
                value(:,i) = fread(fid1, BinNodeCount, 'float')'; % nodeQuality
                fread(fid1, BinLinkCount*8, 'float');
            elseif indParam==15
                fread(fid1, BinNodeCount*4, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkFlow
                fread(fid1, BinLinkCount*7, 'float');
            elseif indParam==16
                fread(fid1, BinNodeCount*4+BinLinkCount, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkVelocity
                fread(fid1, BinLinkCount*6, 'float');
            elseif indParam==17
                fread(fid1, BinNodeCount*4+BinLinkCount*2, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkHeadloss
                fread(fid1, BinLinkCount*5, 'float');
            elseif indParam==18
                fread(fid1, BinNodeCount*4+BinLinkCount*3, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkQuality
                fread(fid1, BinLinkCount*4, 'float');
            elseif indParam==19
                fread(fid1, BinNodeCount*4+BinLinkCount*4, 'float')
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkStatus
                fread(fid1, BinLinkCount*3, 'float');
            elseif indParam==20
                fread(fid1, BinNodeCount*4+BinLinkCount*5, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkSetting
                fread(fid1, BinLinkCount*2, 'float');
            elseif indParam==21
                fread(fid1, BinNodeCount*4+BinLinkCount*6, 'float')
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkReactionRate
                fread(fid1, BinLinkCount, 'float');
            elseif indParam==22
                fread(fid1, BinNodeCount*4+BinLinkCount*7, 'float');
                value(:,i) = fread(fid1, BinLinkCount, 'float')'; % linkFrictionFactor
            elseif indParam>22
                fread(fid1, BinNodeCount*4+BinLinkCount*8, 'float');
            end
        end
        if indParam==23
            value =fread(fid1, 1, 'float')'; % AverageBulkReactionRate
        elseif indParam==24
            fread(fid1, 1, 'float');
            value =fread(fid1, 1, 'float')'; % AverageWallReactionRate
        elseif indParam==25
            fread(fid1, 2, 'float');
            value =fread(fid1, 1, 'float')'; % AverageTankReactionRate
        elseif indParam==26
            fread(fid1, 3, 'float');
            value =fread(fid1, 1, 'float')'; % AverageSourceInflowRate
        end
        fclose(fid1);
    end
    if sum(strcmp(regexp(msg,'\s','split'),'errors.'))
        fprintf('"Run was unsuccessful."\n');
    else
        fprintf('"Run was successful."\n');
    end
end
function errcode=addLinkWarnings(obj,typecode,newLink,toNode)
% Check if id new already exists
Nodes = obj.getBinNodesInfo;
errcode=0;
if length(Nodes.BinNodeNameID)==0
    errcode=-1;
    return;
end
existsTo=0;
for i=1:length(Nodes.BinNodeNameID) 
    existsTo(i) = strcmp(toNode,char(Nodes.BinNodeNameID{i}));
end
if sum(existsTo)~=1
    s = sprintf('There is no node "%s" in the network.',toNode);
    warning(s);errcode=-1;
    return;
end
A = [0 1 2 3 4 5 6 7 8];
code = strfind(A,typecode);
if length(code)==0
    warning('There is no such typecode(0-8)');
    errcode=-1;
    return;
else
    if typecode==0, type_valv = 'CVPIPE';  end
    if typecode==3, type_valv = 'PRV';     end
    if typecode==4, type_valv = 'PSV';     end
    if typecode==5, type_valv = 'PBV';     end
    if typecode==6, type_valv = 'FCV';     end
    if typecode==7, type_valv = 'TCV';     end
    if typecode==8, type_valv = 'GPV';     end
    if typecode~=1 && typecode~=2
        typecode=3;
    end
end
% Valve illegally connected to a tank or reservoir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if typecode==3
    i=1;ifToReservoir=0;
    while i<length(Nodes.BinNodeReservoirNameID)+1
        ifToReservoir(i) = strcmp(toNode,char(Nodes.BinNodeReservoirNameID{i}));
        i=i+1;
    end
    i=1;ifFromTank=0;ifToTank=0;
    while i<length(Nodes.BinNodeTankNameID)+1
        ifToTank(i) = strcmp(toNode,char(Nodes.BinNodeTankNameID{i}));
        i=i+1;
    end
    if sum(ifToReservoir)==1 || sum(ifToTank)==1
        s = sprintf('Valve "%s" illegally connected to a tank.',newLink);
        errcode=-1;
        warning(s);
        return;
    end
end
% Check if newLink already exists
Links = obj.getBinLinksInfo;
for i=1:length(Links.BinLinkNameID)
    exists_link = strcmp(newLink,char(Links.BinLinkNameID(i)));
    if exists_link==1
        s = sprintf('Link %s already exists.',newLink);
        errcode=-1;
        warning(s);
        return;
    end
end

if typecode==2
    crvs = obj.getBinCurvesInfo;
    if isempty(char(crvs.BinCurveNameID))
        s = sprintf('No head curve supplied for pump %s.',newLink);
        errcode=-1;
        warning(s);
        return;
    end
end
end
function cnt=bracketsCheck(v)
    t =  regexp(v, '[(\w*)]','split');
    y=1;cnt=0;
    while y<length(t)+1
        tt = isempty(t{y});
        if tt==0
            cnt=cnt+1;
        end
        y=y+1;
    end
end
function setMsxOptions(obj,varargin)
solver=obj.getMsxSolver;
areaunits=obj.getMsxAreaUnits;
rateunits=obj.getMsxRateUnits;
rtol=obj.getMsxRtol;
atol=obj.getMsxAtol;
timestep=obj.getMsxTimeStep;
coupling=obj.getMsxCoupling;
compiler=obj.getMsxCompiler;

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
                
[info,tline] = readAllFile(obj.MsxFile);
fid2 = fopen([obj.pathfile,obj.MsxTempFile],'w');
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
            if strcmp(a{u},'[OPTIONS]')
                fprintf(fid2,'[OPTIONS]');
                sect=1;
                break;     
            elseif strcmp(a{u},'[SPECIES]')
                fprintf(fid2,'[SPECIES]');
                sect=0;
                break;
            end
            % section [OPTIONS]
            if (sect==1) 
                fprintf(fid2,['AREA_UNITS',blanks(5),'%s\n'],areaunits);
                fprintf(fid2,['RATE_UNITS',blanks(5),'%s\n'],rateunits);
                fprintf(fid2,['SOLVER',blanks(5),'%s\n'],solver);
                fprintf(fid2,['COUPLING',blanks(5),'%s\n'],coupling);
                fprintf(fid2,['COMPILER',blanks(5),'%s\n'],compiler);
                fprintf(fid2,['TIMESTEP',blanks(5),'%d\n'],timestep);
                fprintf(fid2,['RTOL',blanks(5),'%s\n'],rtol);
                fprintf(fid2,['ATOL',blanks(5),'%s\n\n'],atol);
                sect=2;
                break;
            elseif sect==0
                fprintf(fid2,tline{t});
                break;
            end
            u=u+1;
        end
    end
    if sect~=2
        fprintf(fid2,'\n');
    end
end
obj.MsxUnload;
obj.msx(obj.MsxTempFile,1);
end