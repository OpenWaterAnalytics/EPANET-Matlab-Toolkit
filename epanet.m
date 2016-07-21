classdef epanet <handle
    %EPANET-Matlab Toolkit v2.1: A Matlab Class for EPANET and EPANET-MSX
    %libraries
    %
    %
    %   How to run: 
    %   d=epanet('networks/Net1_Rossman2000.inp');
    %   
    %   To select a different DLL version:  
    %   d=epanet('networks/Net1_Rossman2000.inp','epanet2');
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
        ControlLevelValues; %The control level values
        ControlLinkIndex; % Set of control types in links
        ControlNodeIndex; % Set of control types in nodes
        ControlRules; %Retrieves the parameters of all control statements
        ControlRulesCount;  % Number of controls
        Controls;
        ControlSettings; % Settings for the controls
        ControlTypes; % Set of control types
        ControlTypesIndex; %Index of the control types
        CurveCount; % Number of curves
        CurvesInfo;
        EnergyEfficiencyUnits; % Units for efficiency
        EnergyUnits; %Units for energy
        Errcode; %Code for the EPANET error message
        HeadCurveIndex;
        InputFile;  % Name of the input file
        Iterations;
        LibEPANET; %EPANET library dll
        LibEPANETpath; %EPANET library dll path
        LinkBulkReactionCoeff; % Bulk reaction coefficient of each link
        LinkCount; % Number of links
        LinkDiameter; % Diameter of each link
        LinkFlowUnits; %Units of flow
        LinkFrictionFactorUnits; %Units for friction factor
        LinkIndex; % Index of links
        LinkInitialSetting; %Initial settings of links
        LinkInitialStatus; %Initial status of links
        LinkLength; % Length of links
        LinkLengthsUnits; % Units of length
        LinkLengthUnits; % Units of length
        LinkMinorLossCoeff; %Minor loss coefficient of links
        LinkMinorLossCoeffUnits; %Minor loss coefficient units
        LinkNameID; % Name ID of links
        LinkPipeCount;% Number of pipes
        LinkPipeDiameterUnits; % Units for pipe diameters
        LinkPipeIndex; % Index of pipe links
        LinkPipeNameID; % Name ID of pipe links
        LinkPipeRoughnessCoeffUnits; %Pipe roughness coefficient units
        LinkPumpCount;%Number of pumps
        LinkPumpIndex; % Index of pumps
        LinkPumpNameID; % Name ID of pumps
        LinkPumpPatternIndex;
        LinkPumpPatternNameID;
        LinkPumpPowerUnits; %Units of power
        LinkPumpType;
        LinkPumpTypeCode;
        LinkRoughnessCoeff; % Roughness coefficient of links
        LinkType; %ID of link type
        LinkTypeIndex; %Index of link type
        LinkValveCount;% Number of valves
        LinkValveIndex; % Index of valves
        LinkValveNameID; % ID name of valves
        LinkVelocityUnits; % Units for velocity
        LinkWallReactionCoeff; %Wall reaction coefficient of links
        NodeBaseDemands;            % Base demands of nodes
        NodeCoordinates;            % Coordinates for each node (long/lat & intermediate pipe coordinates)
        NodeCount;                  % Number of nodes
        NodeDemandPatternIndex;     % Index of demand patterns
        NodeDemandPatternNameID;
        NodeElevations;             % Elevation of nodes
        NodeElevationUnits;         % Units for elevation
        NodeEmitterCoeff;           % Emmitter Coefficient of nodes
        NodeEmitterCoefficientUnits;% Units for emitter coefficient
        NodeHeadUnits; % Nodal head units
        NodeIndex;                  % Index of nodes
        NodeInitialQuality;         %Initial quality of nodes
        NodeJunctionCount;          % Number of junctions
        NodeJunctionIndex; % Index of node junctions
        NodeJunctionNameID; %Name ID of node junctions
        NodeNameID; %Name ID of all nodes
        NodeDemandCategoriesNumber;
        NodePatternIndex;
        NodePressureUnits; % PUnits for Pressure
        NodeReservoirCount;         % Number of reservoirs
        NodeReservoirIndex; %Index of reservoirs
        NodeReservoirNameID;%Name ID of reservoirs
        NodesConnectingLinksID; %Name IDs of nodes which connect links
        NodesConnectingLinksIndex; %Indices of nodes which connect links
        NodeSourcePatternIndex; %Index of pattern for node sources
        NodeSourceQuality; %Quality of node sources
        NodeSourceTypeIndex; %Index of source type
        NodeTankBulkReactionCoeff; %Bulk reaction coefficients in tanks
        NodeTankCount;              % Number of tanks
        NodeTankDiameter; %Tank Diameters
        NodeTankDiameterUnits;      % Units for tank diameters
        NodeTankIndex; %Indices of Tanks
        NodeTankInitialLevel; %Initial water level in tanks
        NodeTankInitialWaterVolume; %Initial water volume in tanks
        NodeTankMaximumWaterLevel; % Maximum water level in tanks
        NodeTankMaxVolume;
        NodeTankMinimumFraction; % Fraction of the total tank volume devoted to the inlet/outlet compartment
        NodeTankMinimumWaterLevel; % Minimum water level
        NodeTankMinimumWaterVolume; % Minimum water volume
        NodeTankMixingModelCode; % Code of mixing model (MIXED:0, 2COMP:1, FIFO:2, LIFO:3)
        NodeTankMixingModelType; % Type of mixing model (MIXED, 2COMP, FIFO, or LIFO)
        NodeTankMixZoneVolume;  % Mixing zone volume
        NodeTankNameID; % Name ID of Tanks
        NodeTankReservoirCount;     % Number of tanks and reservoirs
        NodeTankVolumeCurveIndex; %Index of curve for tank volumes
        NodeTankVolumeUnits; %Units for volume
        NodeType; %ID of node type
        NodeTypeIndex; %Index of nodetype
        OptionsAccuracyValue; %Convergence value (0.001 is default)
        OptionsEmitterExponent; %Exponent of pressure at an emmiter node (0.5 is default)
        OptionsHeadloss; %Headloss formula (Hazen-Williams, Darcy-Weisbach or Chezy-Manning)
        OptionsHydraulics; %Save or Use hydraulic soltion. *** Not implemented ***
        OptionsMaxTrials; % Maximum number of trials (40 is default)
        OptionsPattern; % *** Not implemented ***
        OptionsPatternDemandMultiplier; %Multiply demand values (1 is default)
        OptionsQualityTolerance; %Tolerance for water quality (0.01 is default)
        OptionsSpecificGravity; %*** Not implemented ***
        OptionsUnbalanced; %*** Not implemented ***
        OptionsUnbalancedContinueN; % *** Not implemented ***
        OptionsViscosity; %*** Not implemented ***
        Pattern; % get all patterns
        PatternAverageValue;
        PatternCount;  % Number of patterns
        PatternDemandsUnits; %Units for demands
        PatternIndex; %Indices of the patterns
        PatternLengths; %Length of the patterns
        PatternNameID; %ID of the patterns
        QualityChemName;
        QualityChemUnits;
        QualityCode; % Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        QualityReactionCoeffBulkUnits; % Bulk reaction coefficient units
        QualityReactionCoeffWallUnits; % Wall reaction coefficient units
        QualitySourceMassInjectionUnits; % Units for source mass injection
        QualityTraceNodeIndex; %Index of trace node (0 if QualityCode<3)
        QualityType; % Water quality analysis type (None/Chemical/Age/Trace)
        QualityUnits; %Units for quality concentration.
        QualityWaterAgeUnits; %Units for water age
        RelativeError;
        TimeHaltFlag;
        TimeHTime;
        TimeHydraulicStep; %Hydraulic time step
        TimeNextEvent; %find the lesser of the hydraulic time step length, or the time to next fill/empty
        TimePatternStart; %Pattern start time
        TimePatternStep; %Pattern Step
        TimeQualityStep; %Quality Step
        TimeReportingPeriods; % Reporting periods
        TimeReportingStart; %Start time for reporting
        TimeReportingStep; %Reporting time step
        TimeRuleControlStep; % Time step for evaluating rule-based controls
        TimeSimulationDuration; %Simulation duration
        TimeStarting_Time;
        TimeStatisticsIndex; %Index of time series post-processing type ('NONE':0,'AVERAGE':1,'MINIMUM':2,'MAXIMUM':3, 'RANGE':4)
        TimeStatisticsType; %Type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
        Units_SI_Metric;
        Units_US_Customary;
        Version; % EPANET version
        
        % Parameters used with EPANET MSX
        MSXLibEPANET;           
        MSXLibEPANETPath;
        MSXConstantsNameID;
        MSXConstantsValue;
        MSXConstantsCount;
        MSXConstantsIndex;
        MSXParametersCount;
        MSXPatternsCount;
        MSXSpeciesCount;
        MSXLinkInitqualValue;
        MSXNodeInitqualValue;
        MSXFile;
        MSXTempFile;
        MSXParametersNameID;
        MSXParametersIndex;
        MSXParametersPipesValue;
        MSXParametersTanksValue;
        MSXPatternsNameID;
        MSXPatternsIndex;
        MSXPatternsLengths;
        MSXPattern;
        MSXEquationsPipes;
        MSXSources;
        MSXSourceLevel;
        MSXSourceNodeNameID;
        MSXSourcePatternNameID;
        MSXSourcePatternIndex;
        MSXSourceType;
        MSXSourceTypeCode;
        MSXSpeciesATOL;
        MSXSpeciesIndex;
        MSXSpeciesNameID;
        MSXSpeciesRTOL;
        MSXSpeciesType;
        MSXSpeciesUnits;
        MSXEquationsTanks;
        MSXEquationsTerms; 
        
        % Parameters used when the Binary mode is used 
        Bin;
        BinControlLinksID;
        BinControlNodesID;
        BinControlRulesCount;
        BinControlsInfo;
        BinCountInitialQualitylines;
        BinCountPatternlines;
        BinCountReactionlines;
        BinCountStatuslines;
        BinCurveAllLines;
        BinCurveCount; % Number of curves
        BinCurveNameID; %ID name of curves
        BinCurveTypes; % Type of curves 
        BinCurveXvalue; %X-value of curves
        BinCurveYvalue; %Y-value of curves
        BinLinkBulkReactionCoeff;
        BinLinkCount;
        BinLinkDiameters;
        BinLinkFlowUnits;
        BinLinkFromNode;
        BinLinkGlobalBulkReactionCoeff;
        BinLinkGlobalWallReactionCoeff;
        BinLinkInitialStatus;
        BinLinkInitialStatusNameID;
        BinLinkLengths;
        BinLinkNameID;
        BinLinkPipeCount;
        BinLinkPipeDiameters;
        BinLinkPipeIndex;
        BinLinkPipeLengths;
        BinLinkPipeMinorLoss;
        BinLinkPipeNameID;
        BinLinkPipeRoughness;
        BinLinkPipeStatus;
        BinLinkPumpCount;
        BinLinkPumpCurveNameID;
        BinLinkPumpIndex;
        BinLinkPumpNameID;
        BinLinkPumpNameIDPower;
        BinLinkPumpPatterns;
        BinLinkPumpPower;
        BinLinkPumpStatus;  
        BinLinkPumpStatusNameID;   
        BinLinkRoughnessCoeff;
        BinLinkSettings;
        BinLinkToNode;
        BinLinkType;
        BinLinkValveCount;
        BinLinkValveDiameters;
        BinLinkValveIndex;
        BinLinkValveMinorLoss;
        BinLinkValveNameID;
        BinLinkValveSetting;
        BinLinkValveStatus;
        BinLinkValveStatusNameID;
        BinLinkValveType;
        BinLinkWallReactionCoeff;
        BinNodeBaseDemands;
        BinNodeCoordinates;
        BinNodeCount;  
        BinNodeDemandPatternNameID;
        BinNodeElevations;
        BinNodeInitialQuality;
        BinNodeJunctionCount; 
        BinNodeJunctionElevation;
        BinNodeJunctionIndex;
        BinNodeJunctionNameID;
        BinNodeJunctionsBaseDemands;
        BinNodeJunctionsBaseDemandsID;
        BinNodeNameID;
        BinNodePressureUnits;
        BinNodeResDemandPatternNameID;
        BinNodeReservoirCount; 
        BinNodeReservoirElevation;
        BinNodeReservoirIndex;
        BinNodeReservoirNameID;
        BinNodeSourcePatternIndex;
        BinNodeSourcePatternNameID;
        BinNodeSourceQuality;
        BinNodeSourceType;
        BinNodeSourceTypeCode;
        BinNodeSourceTypeIndex;
        BinNodeTankCount;
        BinNodeTankDiameter;
        BinNodeTankElevation;
        BinNodeTankIndex;
        BinNodeTankInitLevel;
        BinNodeTankMaxLevel;
        BinNodeTankMinimumFraction;
        BinNodeTankMinLevel;
        BinNodeTankMinVol;
        BinNodeTankMixID;
        BinNodeTankMixModel;
        BinNodeTankNameID;
        BinNodeTankReservoirCount;   
        BinNodeType;
        BinOptionsAccuracyValue;
        BinOptionsDiffusivity;
        BinOptionsEmitterExponent;
        BinOptionsHeadloss;
        BinOptionsMaxTrials;
        BinOptionsPattern;
        BinOptionsPatternDemandMultiplier;
        BinOptionsQualityTolerance;
        BinOptionsSpecificGravity;
        BinOptionsUnbalanced;
        BinOptionsViscosity;
        BinPatternCount;
        BinPatternLengths;
        BinPatternMatrix;
        BinPatternNameID;
        BinPatternValue; 
        BinQualityCode;
        BinQualityTraceNodeID;
        BinQualityTraceNodeIndex;
        BinQualityType;% Water quality analysis code (None:0/Chemical:1/Age:2/Trace:3)
        BinQualityUnits;
        BinRulesControlLinksID;
        BinRulesControlNodesID;
        BinRulesControlsInfo;
        BinRulesCount;
        BinUnits_SI_Metric;
        BinTempfile;
        BinTimeHydraulicStep;
        BinTimePatternStart;
        BinTimePatternStep;
        BinTimeQualityStep;
        BinTimeReportingStart;
        BinTimeReportingStep;
        BinTimeSimulationDuration;
        BinTimeStatistics;
        BinTimeStatisticsIndex;
        BinUnits;
        BinUnits_US_Customary;
    end
    properties (Constant = true)
        classversion='2.1a';
        
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
        
        MSXTYPEAREAUNITS={'FT2','M2','CM2'}; % sets the units used to express pipe wall surface area
        MSXTYPERATEUNITS={'SEC','MIN','HR','DAY'}; % is the units in which all reaction rate terms are expressed
        MSXTYPESOLVER={'EUL','RK5','ROS2'}; % is the choice of numerical integration method used to solve the reaction system
        MSXTYPECOUPLING={'FULL','NONE'}; % is the choice of numerical integration method used to solve the reaction system
        MSXTYPECOMPILER={'NONE','VC','GC'}; % is the choice of numerical integration method used to solve the reaction system
    end
    methods
        function obj = epanet(varargin)
            %Constructor of the EPANET Class
            try unloadlibrary('epanet2');catch e; end
            try unloadlibrary('epanetmsx');catch e; end
            % DLLs
            if strcmp(computer('arch'),'win64')% if no DLL is given, select one automatically
                if exist('64bit')==7 % name is a folder.
                    obj.LibEPANETpath = [pwd,'\64bit\'];
                else
                    warning('Folder "64bit" does not exit.');return;
                end
            elseif strcmp(computer('arch'),'win32')
                if exist('32bit')==7
                    obj.LibEPANETpath = [pwd,'\32bit\'];
                else
                    warning('Folder "32bit" does not exit.');return;
                end
            end
            warning on;   
            obj.InputFile=varargin{1}; % Get name of INP file
            % Bin functions
            if nargin==2
                if strcmp(upper(varargin{2}),'BIN')
                    obj.BinTempfile=[obj.InputFile(1:end-4),'_temp.inp'];
                    copyfile(obj.InputFile,obj.BinTempfile);
                    value=obj.getBinCurvesInfo;
                    if ~isempty(value.BinCurveNameID), obj.remAddBinCurvesID(obj.BinTempfile);end
                    obj.InputFile=obj.BinTempfile;
                    obj.Bin=0;
                    obj = BinUpdateClass(obj);
                    obj.saveBinInpFile;
                    return;
                end
            end
            obj.Bin=1;
            if nargin==2
                if ~isempty(find(obj.InputFile==' '))
                    warning(['File "', obj.InputFile, '" is not a valid']);return;
                end
                obj.LibEPANET=varargin{2}; % Get DLL LibEPANET (e.g. epanet20012x86 for 32-bit)
                obj.LibEPANETpath = [pwd,'\'];
                warning off;
                try  loadlibrary([obj.LibEPANETpath,obj.LibEPANET],[obj.LibEPANETpath,obj.LibEPANET,'.h']); 
                catch e
                   warning on; 
                   obj.Errcode=-1;
                   warning(['File "', obj.LibEPANET, '" is not a valid win application.']);return;
                end
                warning on;
                
            elseif nargin==1
                obj.LibEPANET = 'epanet2';
                [~,inp]=fileparts(obj.InputFile);
                if ~isempty(find(inp==' '))
                    warning(['File "', obj.InputFile, '" is not a valid.']);return;
                end
            end
            if ~exist(obj.InputFile,'file')
                warning(['File "', obj.InputFile, '" does not exist in folder.']);return;
            end
            if strcmp(computer('arch'),'win64') || strcmp(computer('arch'),'win32')
                %Load EPANET Library
                warning off;
                ENLoadLibrary(obj.LibEPANETpath,obj.LibEPANET);
                warning on;
                %Open the file
                obj.Errcode=ENopen(obj.InputFile,[obj.InputFile(1:end-4),'.txt'],[obj.InputFile(1:end-4),'.bin'],obj.LibEPANET);
                if obj.Errcode~=0
                    warning('Could not open the file, please check INP file.');return;
                end
                %Save the temporary input file
                obj.BinTempfile=[obj.InputFile(1:end-4),'_temp.inp'];
                obj.saveInputFile(obj.BinTempfile,1); %create a new INP file (Working Copy) using the SAVE command of EPANET
                obj.closeNetwork;  %ENclose; %Close input file
                %Load temporary file
                obj.Errcode=ENopen(obj.BinTempfile,[obj.InputFile(1:end-4),'_temp.txt'], [obj.InputFile(1:end-4),'_temp.bin'],obj.LibEPANET);
                if obj.Errcode~=0
                    warning('Could not open the file, please check INP file.');return;
                end
            end
            
            % Get some link data
            [obj.LinkDiameter,obj.LinkLength,obj.LinkRoughnessCoeff,obj.LinkMinorLossCoeff,obj.LinkInitialStatus,...
            obj.LinkInitialSetting,obj.LinkBulkReactionCoeff,obj.LinkWallReactionCoeff,obj.NodesConnectingLinksIndex,...
            obj.LinkTypeIndex] = obj.getLinksInfo;
            % Get some node data
            [obj.NodeElevations,obj.NodeDemandPatternIndex,obj.NodeEmitterCoeff,obj.NodeInitialQuality,...
            obj.NodeSourceQuality,obj.NodeSourcePatternIndex,obj.NodeSourceTypeIndex,obj.NodeTypeIndex] = obj.getNodesInfo; 
            %Get all the countable network parameters
            obj.NodeCount = obj.getNodeCount;
            obj.NodeTankReservoirCount = obj.getNodeTankReservoirCount;
            obj.LinkCount = obj.getLinkCount;
            obj.PatternCount = obj.getPatternCount;
            obj.CurveCount = obj.getCurveCount;
            obj.ControlRulesCount = obj.getControlRulesCount;
            obj.NodeJunctionCount = obj.NodeCount-obj.NodeTankReservoirCount; %obj.getNodeJunctionCount;
            % Get type of the parameters
            obj.LinkType=obj.TYPELINK(obj.LinkTypeIndex+1); 
            obj.NodeType=obj.TYPENODE(obj.NodeTypeIndex+1);
            %Get all the countable network parameters
            obj.LinkPipeCount = sum(strcmp(obj.LinkType,'PIPE'))+sum(strcmp(obj.LinkType,'CVPIPE')); %obj.getLinkPipeCount;
            obj.LinkPumpCount = sum(strcmp(obj.LinkType,'PUMP')); %obj.getLinkPumpCount;
            obj.NodeReservoirCount = sum(strcmp(obj.NodeType,'RESERVOIR')); %obj.getNodeReservoirCount;
            obj.NodeTankCount = obj.NodeCount-obj.NodeJunctionCount-obj.NodeReservoirCount; %obj.getNodeTankCount;
            obj.LinkValveCount = obj.LinkCount-obj.LinkPipeCount-obj.LinkPumpCount;%obj.getLinkValveCount;
            %Get all the controls
            obj.Controls = obj.getControls;
            %Get the flow units
            obj.LinkFlowUnits = obj.getFlowUnits;
            %Get all the link data
            obj.LinkNameID = obj.getLinkNameID;
            obj.LinkIndex = 1:obj.LinkCount;   %obj.getLinkIndex;
            obj.LinkPipeIndex = 1:obj.LinkPipeCount; %find(strcmp(obj.LinkType,'PIPE'));
            obj.LinkPumpIndex = obj.LinkPipeCount+1:obj.LinkPipeCount+obj.LinkPumpCount; %find(strcmp(obj.LinkType,'PUMP'));
            obj.LinkValveIndex = find(obj.LinkTypeIndex>2);
            obj.LinkPipeNameID = obj.LinkNameID(obj.LinkPipeIndex);
            obj.LinkPumpNameID = obj.LinkNameID(obj.LinkPumpIndex);
            obj.LinkValveNameID = obj.LinkNameID(obj.LinkValveIndex);
            %Get all the node data
            obj.NodeNameID = obj.getNodeNameID;
            tmp(:,1)=obj.NodeNameID(obj.NodesConnectingLinksIndex(:,1)');
            tmp(:,2) = obj.NodeNameID(obj.NodesConnectingLinksIndex(:,2)');
            obj.NodesConnectingLinksID=tmp;
            obj.NodeIndex = 1:obj.NodeCount; %obj.getNodeIndex;
            obj.NodeReservoirIndex = find(obj.NodeTypeIndex==1); %find(strcmp(obj.NodeType,'RESERVOIR'));
            obj.NodeTankIndex = find(obj.NodeTypeIndex==2); %find(strcmp(obj.NodeType,'TANK'));
            obj.NodeJunctionIndex = 1:obj.NodeJunctionCount; %find(strcmp(obj.NodeType,'JUNCTION'));
            obj.NodeReservoirNameID=obj.NodeNameID(obj.NodeReservoirIndex);
            obj.NodeTankNameID=obj.NodeNameID(obj.NodeTankIndex);
            obj.NodeJunctionNameID=obj.NodeNameID(obj.NodeJunctionIndex);
            obj.NodePatternIndex=obj.getNodePatternIndex;
            obj.NodeBaseDemands = obj.getNodeBaseDemands;
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
            obj.PatternIndex = 1:obj.PatternCount; %obj.getPatternIndex;
            obj.PatternLengths = obj.getPatternLengths;
            obj.Pattern = obj.getPattern;
            %Get quality types
            obj.QualityCode = obj.getQualityCode;
            obj.QualityTraceNodeIndex = obj.getQualityTraceNodeIndex;
%             obj.QualityType = obj.getQualityType;
            n = obj.getQualityInfo;
            obj.QualityChemUnits = n.QualityChemUnits;
            obj.QualityChemName= n.QualityChemUnits;
            
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
            try%New version dev2.1.dll LibEPANET
                obj.TimeStarting_Time = obj.getTimeStarting_Time;
                obj.TimeHTime = obj.getTimeHTime;
                obj.TimeHaltFlag = obj.getTimeHaltFlag;
                obj.TimeNextEvent = obj.getTimeNextEvent;
                obj.NodeTankMaxVolume = obj.getNodeTankMaxVolume;
                obj.NodeBaseDemands = obj.getNodeBaseDemands;
                obj.NodeDemandCategoriesNumber = obj.getNodeDemandCategoriesNumber;
                obj.PatternAverageValue = obj.getPatternAverageValue;
                n = obj.getStatistic;
                obj.RelativeError = n.RelativeError;
                obj.Iterations = n.Iterations;
                obj.NodeDemandPatternNameID = obj.getNodeDemandPatternNameID;
                obj.NodeDemandPatternIndex = obj.getNodeDemandPatternIndex;
                obj.HeadCurveIndex = obj.getHeadCurveIndex;
                obj.LinkPumpPatternNameID = obj.getLinkPumpPatternNameID;
                obj.LinkPumpPatternIndex = obj.getLinkPumpPatternIndex;
                obj.LinkPumpTypeCode = obj.getLinkPumpTypeCode;
                obj.LinkPumpType = obj.getLinkPumpType;
                obj.CurvesInfo = obj.getCurvesInfo; % New version dev2.1
            catch e
            end
            %Get data from raw file (for information which cannot be
            %accessed by the epanet library)
            value=obj.getNodeCoordinates;
            %Get coordinates
            obj.NodeCoordinates{1} = value{1};
            obj.NodeCoordinates{2} = value{2};
            obj.NodeCoordinates{3} = value{3};
            obj.NodeCoordinates{4} = value{4};
            
            %     US Customary - SI metric
            if find(strcmp(obj.LinkFlowUnits, obj.TYPEUNITS))<6
                obj.Units_US_Customary=1;
            else
                obj.Units_SI_Metric=1;
            end

            if obj.Units_US_Customary==1;
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
        function Errcode = loadEPANETFile(obj,varargin)
           [Errcode] = ENopen(varargin{1},[varargin{1}(1:end-4),'.rpt'],[varargin{1}(1:end-4),'.bin'],obj.LibEPANET); 
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
                    [obj.Errcode, obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)] = ENgetcontrol(i,obj.LibEPANET);
                    obj.ControlTypes(i)={obj.TYPECONTROL(obj.ControlTypesIndex(i)+1)};
                    value{i}={obj.ControlTypes{i},obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)};
                end
            else
                value=-1;
            end
        end
        function value = getNodeCount(obj)
            % Retrieves the number of nodes
            [obj.Errcode, value] = ENgetcount(0,obj.LibEPANET);
        end
        function value = getNodeTankReservoirCount(obj)
            % Retrieves the number of tanks
            [obj.Errcode, value] = ENgetcount(1,obj.LibEPANET);
        end
        function value = getLinkCount(obj)
            % Retrieves the number of links
            [obj.Errcode, value] = ENgetcount(2,obj.LibEPANET);
        end
        function value = getPatternCount(obj)
            % Retrieves the number of patterns
            [obj.Errcode, value] = ENgetcount(3,obj.LibEPANET);
        end
        function value = getCurveCount(obj)
            % Retrieves the number of curves
            [obj.Errcode, value] = ENgetcount(4,obj.LibEPANET);
        end
        function value = getControlRulesCount(obj)
            % Retrieves the number of controls
            [obj.Errcode, value] = ENgetcount(5,obj.LibEPANET);
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
        function errmssg = getError(obj,Errcode)
            %Retrieves the text of the message associated with a particular error or warning code.
            [errmssg , Errcode] = ENgeterror(Errcode,obj.LibEPANET);
        end
        function value = getFlowUnits(obj)
            %Retrieves flow units used to express all flow rates.
            [obj.Errcode, flowunitsindex] = ENgetflowunits(obj.LibEPANET);
            value=obj.TYPEUNITS(flowunitsindex+1);
        end
        function value = getLinkNameID(obj,varargin)
            % Retrieves the ID label(s) of all links, or the IDs of an index set of links
            if isempty(varargin)
                value{obj.getLinkCount}=[];
                for i=1:obj.getLinkCount
                    [obj.Errcode, value{i}]=ENgetlinkid(i,obj.LibEPANET);
                end
            else
                k=1;
                if ~isempty(varargin{1})
                    value{length(varargin{1})}=[];
                    for i=varargin{1}
                        [obj.Errcode, value{k}]=ENgetlinkid(i,obj.LibEPANET);
                        if obj.Errcode==204, error(obj.getError(obj.Errcode)), return; end   
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
                    [obj.Errcode, value(k)] = ENgetlinkindex(varargin{1}{j},obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetlinkindex(varargin{1},obj.LibEPANET);
            end
        end
        function value = getLinkPipeIndex(obj)
            %Retrieves the pipe indices
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PIPE'));
%             if isempty(value), value=-1; end
        end
        function value = getLinkPumpIndex(obj)
            %Retrieves the pump indices
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PUMP'));
%             if isempty(value), value=-1; end
        end
        function value = getLinkValveIndex(obj)
            %Retrieves the valve indices
            value = obj.getLinkPipeCount+obj.getLinkPumpCount+1:obj.getLinkCount;
%             if isempty(value), value=-1; end
        end
        function value = getLinkNodesIndex(obj)
            %Retrieves the indexes of the from/to nodes of all links.
            value(obj.getLinkCount,1:2)=[nan nan];
            for i=1:obj.getLinkCount
                [obj.Errcode,linkFromNode,linkToNode] = ENgetlinknodes(i,obj.LibEPANET);
                value(i,:)= [linkFromNode,linkToNode];
            end
        end
        function value = getNodesConnectingLinksID(obj)
            %Retrieves the id of the from/to nodes of all links.
            obj.NodesConnectingLinksIndex=obj.getLinkNodesIndex;
            if obj.getLinkCount
                value(:,1)=obj.getNodeNameID(obj.NodesConnectingLinksIndex(:,1)');
                value(:,2) = obj.getNodeNameID(obj.NodesConnectingLinksIndex(:,2)');
            else
                value=-1;
            end
        end
        function value = getLinkType(obj, varargin)
            %Retrieves the link-type code for all links.
            indices = getLinkIndices(obj,varargin);
            value=obj.TYPELINK(obj.getLinkTypeIndex(indices)+1);
        end
        function value = getLinkTypeIndex(obj, varargin)
            %Retrieves the link-type code for all links.
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode,value(j)] = ENgetlinktype(i,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end  
                j=j+1;
            end
        end
        function [LinkDiameter,LinkLength,LinkRoughnessCoeff,LinkMinorLossCoeff,LinkInitialStatus,...
                LinkInitialSetting,LinkBulkReactionCoeff,LinkWallReactionCoeff,NodesConnectingLinksIndex,...
                LinkTypeIndex] = getLinksInfo(obj)
            for i=1:obj.getLinkCount
                [~, LinkDiameter(i)] = ENgetlinkvalue(i,0,obj.LibEPANET);
                [~, LinkLength(i)] = ENgetlinkvalue(i,1,obj.LibEPANET);
                [~, LinkRoughnessCoeff(i)] = ENgetlinkvalue(i,2,obj.LibEPANET);
                [~, LinkMinorLossCoeff(i)] = ENgetlinkvalue(i,3,obj.LibEPANET);
                [~, LinkInitialStatus(i)] = ENgetlinkvalue(i,4,obj.LibEPANET);
                [~, LinkInitialSetting(i)] = ENgetlinkvalue(i,5,obj.LibEPANET);
                [~, LinkBulkReactionCoeff(i)] = ENgetlinkvalue(i,6,obj.LibEPANET);
                [~, LinkWallReactionCoeff(i)] = ENgetlinkvalue(i,7,obj.LibEPANET);
                [~,NodesConnectingLinksIndex(i,1),NodesConnectingLinksIndex(i,2)] = ENgetlinknodes(i,obj.LibEPANET);
                [~,LinkTypeIndex(i)] = ENgetlinktype(i,obj.LibEPANET);
            end
        end
        function value = getLinkDiameter(obj, varargin)
            indices = getLinkIndices(obj,varargin);j=1;
            %Retrieves the value of all link diameters
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,0,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkLength(obj, varargin)
            %Retrieves the value of all link lengths
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,1,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end 
                j=j+1;
            end
        end
        function value = getLinkRoughnessCoeff(obj, varargin)
            %Retrieves the value of all link roughness
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,2,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkMinorLossCoeff(obj, varargin)
            %Retrieves the value of all link minor loss coefficients
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,3,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkInitialStatus(obj, varargin)
            %Retrieves the value of all link initial status
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,4,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkInitialSetting(obj, varargin)
            %Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,5,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkBulkReactionCoeff(obj, varargin)
            %Retrieves the value of all link bulk reaction coefficients
            value=[];
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,6,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkWallReactionCoeff(obj, varargin)
            %Retrieves the value of all link wall reaction coefficients
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,7,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkFlows(obj, varargin)
            %Retrieves the value of all computed link flow rates
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,8,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkVelocity(obj, varargin)
            %Retrieves the value of all computed link velocities
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,9,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkHeadloss(obj, varargin)
            %Retrieves the value of all computed link headloss
            value=[];
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,10,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkStatus(obj, varargin)
            %Retrieves the value of all computed link status (0 = closed, 1 = open)
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,11,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkSettings(obj, varargin)
            %Retrieves the value of all computed link roughness for pipes or actual speed for pumps or actual setting for valves
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,12,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkPumpEnergy(obj, varargin)
            %Retrieves the value of all computed energy in kwatts
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,13,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkQuality(obj, varargin)
            %New version dev2.1
            indices = getLinkIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetlinkvalue(i,14,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getLinkPumpPatternIndex(obj)
            %New version dev2.1
            value=zeros(1,obj.getLinkPumpCount);v=1;
            for i=obj.getLinkPumpIndex
                [obj.Errcode, value(v)] = ENgetlinkvalue(i,15,obj.LibEPANET);
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
                    [obj.Errcode, value{i}]=ENgetnodeid(i,obj.LibEPANET);
                end
            else
                k=1;
                value{length(varargin{1})}=[];
                for i=varargin{1}
                    [obj.Errcode, value{k}]=ENgetnodeid(i,obj.LibEPANET);
                    if obj.Errcode==203, value=NaN; return; end   
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
                    [obj.Errcode, value(k)] = ENgetnodeindex(varargin{1}{j},obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetnodeindex(varargin{1},obj.LibEPANET);
            end
        end
        function value = getNodeReservoirIndex(obj)
            %Retrieves the indices of reservoirs
            if obj.getNodeReservoirCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'RESERVOIR'));
            else
                value=0;
            end
        end
        function value = getNodeJunctionIndex(obj)
            %Retrieves the indices of junctions
            if obj.getNodeJunctionCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'JUNCTION'));
            else
                value=0;
            end
        end
        function value = getNodeType(obj, varargin)
            %Retrieves the node-type code for all nodes
            indices = getNodeIndices(obj,varargin);
            value=obj.TYPENODE(obj.getNodeTypeIndex(indices)+1);
        end
        function value = getNodeTypeIndex(obj, varargin)
            %Retrieves the node-type code for all nodes
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode,value(j)] = ENgetnodetype(i,obj.LibEPANET);  
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end  
                j=j+1;
            end
        end
        function [NodeElevations,NodeDemandPatternIndex,NodeEmitterCoeff,NodeInitialQuality,...
                NodeSourceQuality,NodeSourcePatternIndex,NodeSourceTypeIndex,NodeTypeIndex] = getNodesInfo(obj)
            for i=1:obj.getNodeCount
                [~, NodeElevations(i)] = ENgetnodevalue(i,0,obj.LibEPANET);
                [~, NodeDemandPatternIndex(i)] = ENgetnodevalue(i,2,obj.LibEPANET);
                [~, NodeEmitterCoeff(i)] = ENgetnodevalue(i,3,obj.LibEPANET);
                [~, NodeInitialQuality(i)] = ENgetnodevalue(i,4,obj.LibEPANET);
                [~, NodeSourceQuality(i)] = ENgetnodevalue(i,5,obj.LibEPANET);
                [~, NodeSourcePatternIndex(i)] = ENgetnodevalue(i,6,obj.LibEPANET);
                [~, NodeSourceTypeIndex(i)] = ENgetnodevalue(i,7,obj.LibEPANET);
                [~, NodeTypeIndex(i)] = ENgetnodetype(i,obj.LibEPANET);
            end
        end
        function value = getNodeElevations(obj, varargin)
            %Retrieves the value of all node elevations
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,0,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeBaseDemands(obj, varargin)
            %New version dev2.1
            chckfunctions=libfunctions(obj.LibEPANET);
            if sum(strcmp(chckfunctions,'ENgetbasedemand'))
                numdemands = obj.getNodeDemandCategoriesNumber;
                val=zeros(max(numdemands),obj.getNodeCount);
                for i=obj.getNodeIndex
                    v=1;
                    for u=1:numdemands(i)
                        [obj.Errcode, val(v,i)] = ENgetbasedemand(i,u,obj.LibEPANET);v=v+1;
                    end
                end
                for i=1:size(val,1)
                    value{i} = val(i,:);
                end
            else%version epanet20012
                %Retrieves the value of all node base demands
                indices = getNodeIndices(obj,varargin);j=1;
                for i=indices
                    [obj.Errcode, value(j)] = ENgetnodevalue(i,1,obj.LibEPANET); 
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                    j=j+1;
                end   
            end
        end
        function value = getNodeDemandCategoriesNumber(obj, varargin)
            %New version dev2.1
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnumdemands(i,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeDemandPatternIndex(obj)
            %New version dev2.1
            numdemands = obj.getNodeDemandCategoriesNumber;
            val=zeros(max(numdemands),obj.getNodeCount);
            for i=obj.getNodeIndex
                v=1;
                for u=1:numdemands(i)
                    [obj.Errcode, val(v,i)] = ENgetdemandpattern(i,u,obj.LibEPANET);v=v+1;
                end             
            end
            for i=1:size(val,1)
                value{i} = val(i,:);
            end
        end
        function value = getNodeDemandPatternNameID(obj, varargin)
            %New version dev2.1
            value={};
            v = obj.getNodeDemandPatternIndex;
            m = obj.getPatternNameID;
            if ~isempty(varargin)
                numdemands = obj.getNodeDemandCategoriesNumber(varargin{1});
            else
                numdemands = obj.getNodeDemandCategoriesNumber;
            end
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                for u=1:numdemands(i)
                    if v{u}(i)~=0 
                        val{u,j}= char(m(v{u}(i))); 
                    else
                        val{u,j}= '';
                    end
                    j=j+1;
                end
                if numdemands(i)==0
                    val{1,j}= []; j=j+1;
                end
            end
            for i=1:size(val,1)
                value{i} = val(i,:);
            end
        end
        function value = getStatistic(obj)
            %New version dev2.1
            %Input: none
            %Output: *iter = # of iterations to reach solution
            %*relerr = convergence error in solution
            %returns error code
            [obj.Errcode, value.Iterations] = ENgetstatistic(0,obj.LibEPANET);
            [obj.Errcode, value.RelativeError] = ENgetstatistic(1,obj.LibEPANET);
        end    
        function value = getNodePatternIndex(obj, varargin)
            %Retrieves the value of all node demand pattern indices
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,2,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeEmitterCoeff(obj, varargin)
            %Retrieves the value of all node emmitter coefficients
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,3,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeInitialQuality(obj, varargin)
            %Retrieves the value of all node initial quality
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,4,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeSourceQuality(obj, varargin)
            %Retrieves the value of all nodes source quality
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,5,obj.LibEPANET); 
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeSourcePatternIndex(obj, varargin)
            %Retrieves the value of all node source pattern index
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,6,obj.LibEPANET); 
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeSourceTypeIndex(obj, varargin)
            %Retrieves the value of all node source index
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,7,obj.LibEPANET); 
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeSourceType(obj, varargin)
            %Retrieves the value of all node source type
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, temp] = ENgetnodevalue(i,7,obj.LibEPANET);
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                if ~isnan(temp)
                    value(j)=obj.TYPESOURCE(temp+1);
                else
                    value{j}=temp;
                end
                j=j+1;
            end
        end
        function value = getNodeTankInitialLevel(obj, varargin)
            %Retrieves the value of all tank initial water levels
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,8,obj.LibEPANET); 
                if obj.Errcode==251, value(j)=NaN; end
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeActualDemand(obj, varargin)
            %Retrieves the computed value of all actual demands
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,9,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeActualDemandSensingNodes(obj, varargin)
            %Retrieves the computed demand values at some sensing nodes
            value=zeros(1,length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.Errcode, value(v)] = ENgetnodevalue(i,9,obj.LibEPANET);v=v+1;
            end
        end
        function value = getNodeHydaulicHead(obj, varargin)
            %Retrieves the computed values of all hydraulic heads
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,10,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodePressure(obj, varargin)
            %Retrieves the computed values of all node pressures
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,11,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeActualQuality(obj, varargin)
            %Retrieves the computed values of the actual quality for all nodes
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,12,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeMassFlowRate(obj, varargin)
            %Retrieves the computed mass flow rates per minute of chemical sources
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,13,obj.LibEPANET); 
%                 if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeActualQualitySensingNodes(obj,varargin)
            %Retrieves the computed quality values at some sensing nodes
            value=zeros(1,length(varargin{1}));v=1;
            for i=varargin{1}
                [obj.Errcode, value(v)] = ENgetnodevalue(i,12,obj.LibEPANET);v=v+1;
            end
        end
        function value = getNodeTankInitialWaterVolume(obj, varargin)
            %Retrieves the tank initial volume
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,14,obj.LibEPANET); 
                if obj.Errcode==251, value(j)=NaN; end
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMixiningModel(obj)
            %Retrieves the tank mixing mode (mix1, mix2, fifo, lifo)
            obj.NodeTankInitialWaterVolume=nan(1,obj.getNodeCount);
            obj.NodeTankMixingModelCode=nan(1,obj.getNodeCount);
            obj.NodeTankMixingModelType={};
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.Errcode, obj.NodeTankMixingModelCode(i)] = ENgetnodevalue(i, 15,obj.LibEPANET);
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
        function value = getNodeTankMixZoneVolume(obj, varargin)
            %Retrieves the tank mixing zone volume
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,16,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankDiameter(obj, varargin)
            %Retrieves the tank diameters
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,17,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMinimumWaterVolume(obj, varargin)
            %Retrieves the tank minimum volume
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,18,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankVolumeCurveIndex(obj, varargin)
            %Retrieves the tank volume curve index
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,19,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMinimumWaterLevel(obj, varargin)
            %Retrieves the tank minimum water level
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,20,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMaximumWaterLevel(obj, varargin)
            %Retrieves the tank maximum water level
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,21,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMinimumFraction(obj, varargin)
            %Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,22,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankBulkReactionCoeff(obj, varargin)
            %Retrieves the tank bulk rate coefficient
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,23,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankVolume(obj, varargin)
            %New version dev2.1
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,24,obj.LibEPANET); 
                if obj.Errcode==251, value(j)=NaN; end
                if obj.Errcode==203, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
            end
        end
        function value = getNodeTankMaxVolume(obj, varargin)
            %New version dev2.1
            indices = getNodeIndices(obj,varargin);j=1;
            for i=indices
                [obj.Errcode, value(j)] = ENgetnodevalue(i,25,obj.LibEPANET); 
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                j=j+1;
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
            [obj.Errcode, value] = ENgetoption(0,obj.LibEPANET);
        end
        function value = getOptionsAccuracyValue(obj)
            % Retrieve the analysis convergence criterion (0.001)
            [obj.Errcode, value] = ENgetoption(1,obj.LibEPANET);
        end
        function value = getOptionsQualityTolerance(obj)
            % Retrieve the water quality analysis tolerance
            [obj.Errcode, value] = ENgetoption(2,obj.LibEPANET);
        end
        function value = getOptionsEmitterExponent(obj)
            % Retrieve power exponent for the emmitters (0.5)
            [obj.Errcode, value] = ENgetoption(3,obj.LibEPANET);
        end
        function value = getOptionsPatternDemandMultiplier(obj)
            % Retrieve the demand multiplier (x1)
            [obj.Errcode, value] = ENgetoption(4,obj.LibEPANET);
        end
        function value = getPatternNameID(obj,varargin)
            %Retrieves the ID label of all or some time patterns indices
            value={};
            if obj.getPatternCount
                if isempty(varargin) 
                    value{obj.getPatternCount}=[];
                    for i=1:obj.getPatternCount
                        [obj.Errcode, value{i}]=ENgetpatternid(i,obj.LibEPANET);
                    end
                else
                    k=1;
                    for i=varargin{1}
                        [obj.Errcode, value{k}]=ENgetpatternid(i,obj.LibEPANET);
    %                     if obj.Errcode==205, error(obj.getError(obj.Errcode)), return; end   
                        k=k+1;
                    end
                end
            end
        end
        function value = getCurveNameID(obj,varargin)
            %Retrieves ID of a curve with specific index
            %New version dev2.1
            if obj.getCurveCount
                if isempty(varargin) 
                    value{obj.getCurveCount}=[];
                    for i=1:obj.getCurveCount
                        [obj.Errcode, value{i}]=ENgetcurveid(i,obj.LibEPANET);
                    end
                else
                    k=1;
                    for i=varargin{1}
                        [obj.Errcode, value{k}]=ENgetcurveid(i,obj.LibEPANET);
                        if obj.Errcode==206, error('Input Error 205: function call  refers to undefined curve.'), return; end   
                        k=k+1;
                    end
                end
            end
        end
        function value = getCurveLengths(obj,varargin)
            %Retrieves number of points in a curve 
            %New version dev2.1
            if isempty(varargin)
                tmpCurves=1:obj.getCurveCount;
                for i=tmpCurves
                    [obj.Errcode, value(i)]=ENgetcurvelen(i,obj.LibEPANET);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = ENgetcurvelen(obj.getCurveIndex(varargin{1}{j}),obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetcurvelen(obj.getCurveIndex(varargin{1}),obj.LibEPANET);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.Errcode, value(k)]=ENgetcurvelen(i,obj.LibEPANET);
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
                    [obj.Errcode, value(k)] = ENgetcurveindex(varargin{1}{j},obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetcurveindex(varargin{1},obj.LibEPANET);
            end
        end
        function setCurve(obj,index,curveVector)
            %Sets x,y values for a specific curve
            %New version dev2.1
            nfactors=size(curveVector,1);%x = number of points in curve
            [obj.Errcode] = ENsetcurve(index, curveVector(:,1), curveVector(:,2), nfactors,obj.LibEPANET);
        end
        function setCurveValue(obj,index, curvePnt, value)
            %Retrieves x,y point for a specific point number and curve
            %New version dev2.1
            x=value(1); y=value(2);
            [obj.Errcode] = ENsetcurvevalue(index, curvePnt, x, y, obj.LibEPANET);
        end
        function value = getPatternIndex(obj,varargin)
            %Retrieves the index of all or some time patterns IDs
            if isempty(varargin)
                value=1:obj.getPatternCount;
            elseif isa(varargin{1},'cell')
                k=1;
                value{length(varargin{1})}=[];
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = ENgetpatternindex(varargin{1}{j},obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetpatternindex(varargin{1},obj.LibEPANET);
            end
        end
        function value = getPatternLengths(obj,varargin)
            %Retrieves the number of time periods in all or some patterns
            if isempty(varargin)
                tmpPatterns=1:obj.getPatternCount;
                for i=tmpPatterns
                    [obj.Errcode, value(i)]=ENgetpatternlen(i,obj.LibEPANET);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = ENgetpatternlen(obj.getPatternIndex(varargin{1}{j}),obj.LibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = ENgetpatternlen(obj.getPatternIndex(varargin{1}),obj.LibEPANET);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.Errcode, value(k)]=ENgetpatternlen(i,obj.LibEPANET);
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
                    [obj.Errcode, value(i,j)] = ENgetpatternvalue(i, j,obj.LibEPANET);
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
            [obj.Errcode, value] = ENgetpatternvalue(patternIndex, patternStep,obj.LibEPANET);
        end
        function value = getQualityType(obj,varargin)
            %Retrieves the type of water quality analysis type
%             [obj.Errcode, obj.QualityCode,obj.QualityTraceNodeIndex] = ENgetqualinfo(obj.LibEPANET); % bug
%             value=obj.TYPEQUALITY(obj.QualityCode+1);
            if nargin>1
                obj.saveInputFile(obj.BinTempfile,1);
            else
                obj.saveInputFile(obj.BinTempfile);
            end
            value = obj.getBinQualType;
%             value = {obj.getBinOptionsInfo.BinQualityType};
        end 
        function value = getQualityInfo(obj)
            [obj.Errcode, ~,value.QualityChemName,value.QualityChemUnits,~] = ENgetqualinfo(obj.LibEPANET);
        end
        function value = getQualityCode(obj)
            %Retrieves the code of water quality analysis type
            [obj.Errcode, value,obj.QualityTraceNodeIndex] = ENgetqualtype(obj.LibEPANET);
        end
        function value = getQualityTraceNodeIndex(obj)
            %Retrieves the trace node index of water quality analysis type
            [obj.Errcode, obj.QualityCode,value] = ENgetqualtype(obj.LibEPANET);
        end
        function value = getTimeSimulationDuration(obj)
            %Retrieves the value of simulation duration
            [obj.Errcode, value] = ENgettimeparam(0,obj.LibEPANET);
        end
        function value = getTimeHydraulicStep(obj)
            %Retrieves the value of the hydraulic time step
            [obj.Errcode, value] = ENgettimeparam(1,obj.LibEPANET);
        end
        function value = getTimeQualityStep(obj)
            %Retrieves the value of the water quality time step
            [obj.Errcode, value] = ENgettimeparam(2,obj.LibEPANET);
        end
        function value = getTimePatternStep(obj)
            %Retrieves the value of the pattern time step
            [obj.Errcode, value] = ENgettimeparam(3,obj.LibEPANET);
        end
        function value = getTimePatternStart(obj)
            %Retrieves the value of pattern start time
            [obj.Errcode, value] = ENgettimeparam(4,obj.LibEPANET);
        end
        function value = getTimeReportingStep(obj)
            %Retrieves the value of the reporting time step
            [obj.Errcode, value] = ENgettimeparam(5,obj.LibEPANET);
        end
        function value = getTimeReportingStart(obj)
            %Retrieves the value of the reporting start time
            [obj.Errcode, value] = ENgettimeparam(6,obj.LibEPANET);
        end
        function value = getTimeRuleControlStep(obj)
            %Retrieves the time step for evaluating rule-based controls
            [obj.Errcode, value] = ENgettimeparam(7,obj.LibEPANET);
        end
        function value = getTimeStatisticsType(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.Errcode, obj.TimeStatisticsIndex] = ENgettimeparam(8,obj.LibEPANET);
            value=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
        end
        function value = getTimeStatisticsIndex(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.Errcode, value] = ENgettimeparam(8,obj.LibEPANET);
        end
        function value = getTimeReportingPeriods(obj)
            %Retrieves the number of reporting periods saved to the binary
            [obj.Errcode, value] = ENgettimeparam(9,obj.LibEPANET);
        end
        %%%%% New version dev2.1 %%%%%
        function value = getTimeStarting_Time(obj)
            %Retrieves the number of start time
            [obj.Errcode, value] = ENgettimeparam(10,obj.LibEPANET);
        end
        function value = getTimeHTime(obj)
            %Retrieves the number of htime
            [obj.Errcode, value] = ENgettimeparam(11,obj.LibEPANET);
        end
        function value = getTimeHaltFlag(obj)
            %Retrieves the number of  halt flag
            [obj.Errcode, value] = ENgettimeparam(12,obj.LibEPANET);
        end
        function value = getTimeNextEvent(obj)
            %Retrieves the number of next event 
            [obj.Errcode, value] = ENgettimeparam(13,obj.LibEPANET);
        end
        function value = getCurvesInfo(obj)
            %New version dev2.1
            %Input:   curveIndex = curve index
            %Output:  *nValues = number of points on curve
            %         *xValues = values for x
            %         *yValues = values for y
            %Returns: error code
            %Purpose: retrieves end nodes of a specific link
            for i=1:obj.getCurveCount
                [obj.Errcode, value.CurveNameID{i}, value.CurveNvalue{i}, value.CurveXvalue{i}, value.CurveYvalue{i}] = ENgetcurve(i,obj.LibEPANET);
            end
        end
        function value = getConnectivityMatrix(obj, varargin)
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
                [obj.Errcode] = ENaddcurve(varargin{1},obj.LibEPANET);
                valueIndex = getCurveIndex(obj,varargin{1});
            elseif nargin==3
                [obj.Errcode] = ENaddcurve(varargin{1},obj.LibEPANET);
                valueIndex = getCurveIndex(obj,varargin{1});
                setCurve(obj,valueIndex,varargin{2});
            end
        end
        function value = getCurveValue(obj,index,varargin)
            %New version dev2.1
            %Retrieves x,y point for a specific point number and curve
            if nargin<3
                tmplen=obj.getCurveLengths;
                value=zeros(tmplen(index),2);
                for i=1:tmplen(index)
                    [obj.Errcode, value(i,1), value(i,2)] = ENgetcurvevalue(index, i,obj.LibEPANET);
                end
            else
                [obj.Errcode, x, y] = ENgetcurvevalue(index, varargin{1},obj.LibEPANET);
                %if obj.Errcode, error(obj.getError(obj.Errcode)), return; end
                value = [x y];
            end
        end
        function value = getHeadCurveIndex(obj)
            %New version dev2.1
            %Retrieves index of a head curve for specific link index
            v=1;
            for i=obj.getLinkPumpIndex
                [obj.Errcode, value(v)] = ENgetheadcurveindex(i,obj.LibEPANET);
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
                    [obj.Errcode, value(v)] = ENgetpumptype(i,obj.LibEPANET);
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
        function value = getPatternAverageValue(obj)
            %New version dev2.1
            for i=obj.getPatternIndex
                [obj.Errcode, value(i)] = ENgetaveragepatternvalue(i,obj.LibEPANET);
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
        function value = getLibFunctions(obj)
            value = libfunctions(obj.LibEPANET)'';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = getVersion(obj)
            % Retrieve the current EPANET LibEPANET
            [obj.Errcode, value] = ENgetversion(obj.LibEPANET);
        end
        function value = getComputedHydraulicTimeSeries(obj,varargin)
            % Compute hydraulic simulation and retrieve all time-series
            obj.openHydraulicAnalysis;
            obj.initializeHydraulicAnalysis
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep;
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
            clear initlinkmatrix initnodematrix;
            k=1;tstep=1;
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
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeQualityStep;
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
            clear initnodematrix;
            k=1;t=1;tleft=1;
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
            [obj.Errcode] = ENsolveH(obj.LibEPANET);
        end
        function solveCompleteQuality(obj)
            [obj.Errcode] = ENsolveQ(obj.LibEPANET);
        end
        function valueIndex = addPattern(obj,varargin)
            valueIndex=-1;
            if nargin==2
                [obj.Errcode] = ENaddpattern(varargin{1},obj.LibEPANET);
                valueIndex = getPatternIndex(obj,varargin{1});
            elseif nargin==3
                [obj.Errcode] = ENaddpattern(varargin{1},obj.LibEPANET);
                valueIndex = getPatternIndex(obj,varargin{1});
                setPattern(obj,valueIndex,varargin{2});
            end
        end
        function setControl(obj,controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel)
            % Example: d.setControl(1,1,13,1,11,150)
            if controlRuleIndex<=obj.getControlRulesCount
                [obj.Errcode] = ENsetcontrol(controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel,obj.LibEPANET);
            else
                disp('New rules cannot be added in this LibEPANET')
            end
        end
        function setLinkDiameter(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i,0, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end            
        end
        function setLinkLength(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 1, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkRoughnessCoeff(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 2, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkMinorLossCoeff(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 3, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkInitialStatus(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices% Cannot set status for a check valve
                [obj.Errcode] = ENsetlinkvalue(i, 4, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkInitialSetting(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 5, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkBulkReactionCoeff(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 6, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkWallReactionCoeff(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 7, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkStatus(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 11, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setLinkSettings(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getLinkIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetlinkvalue(i, 12, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeElevations(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, 0, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeBaseDemands(obj, value, varargin)
            if nargin==3, 
                indices = value; value=varargin{1};j=1;
                for i=indices
                    [obj.Errcode] = ENsetnodevalue(i, 1, value(j),obj.LibEPANET); j=j+1;
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                end
                return;
            end
            %New version dev2.1
            chckfunctions=libfunctions(obj.LibEPANET);
            if sum(strcmp(chckfunctions,'ENsetbasedemand'))
                NodeNumDemandC=obj.getNodeDemandCategoriesNumber;
                for i=1:obj.getNodeJunctionCount
                    for u=1:NodeNumDemandC(i)
                        [obj.Errcode] = ENsetbasedemand(i, NodeNumDemandC(u), value{NodeNumDemandC(u)}(i),obj.LibEPANET);
                    end
                end
            else %version epanet20012
                if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
                j=1;
                for i=indices
                    [obj.Errcode] = ENsetnodevalue(i, 1, value(j),obj.LibEPANET); j=j+1;
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                end
            end
        end
        function setNodeCoordinates(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            if ~isempty(varargin)
                for i=indices
                    [obj.Errcode] = ENsetcoord(i,value(1),value(2),obj.LibEPANET);
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                end
            else
                for i=1:length(value)
                    x=value{1}(i);
                    y=value{2}(i);
                    [obj.Errcode] = ENsetcoord(i,x,y,obj.LibEPANET);
                end
            end
        end
        function setNodeDemandPatternIndex(obj, value, varargin)
            if nargin==3, 
                indices = value; value=varargin{1};j=1;
                for i=indices
                    [obj.Errcode] = ENsetnodevalue(i, 2, value(j),obj.LibEPANET); j=j+1;
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                end
                return;
            end
            if iscell(value)
                value=value{1};
            end
            for i=1:length(value)
                [obj.Errcode] = ENsetnodevalue(i, 2, value(i),obj.LibEPANET);
            end
        end
        function setNodeEmitterCoeff(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, 3, value(j),obj.LibEPANET); j=j+1;
%                 if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeInitialQuality(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, 4, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeTankInitialLevel(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 8, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 8, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankMixingModelType(obj, value, varargin)
            if nargin==3
                code=strfind(strcmpi(varargin{1},obj.TYPEMIXMODEL),1)-1;
                [obj.Errcode] = ENsetnodevalue(value, 15, code,obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                code=strfind(strcmpi(value(i),obj.TYPEMIXMODEL),1)-1;
                [obj.Errcode] = ENsetnodevalue(i, 15, code,obj.LibEPANET);
            end
        end
        function setNodeTankDiameter(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 17, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 17, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankMinimumWaterLevel(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 20, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 20, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankMinimumWaterVolume(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 18, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 18, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankMaximumWaterLevel(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 21, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 21, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankMinimumFraction(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 22, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 22, value(i),obj.LibEPANET);
            end
        end
        function setNodeTankBulkReactionCoeff(obj, value, varargin)
            if nargin==3
                [obj.Errcode] = ENsetnodevalue(value, 23, varargin{1},obj.LibEPANET);
                if obj.Errcode, error(obj.getError(obj.Errcode)); end 
                return;
            end
            for i=obj.getNodeTankIndex
                [obj.Errcode] = ENsetnodevalue(i, 23, value(i),obj.LibEPANET);
            end
        end
        function setNodeSourceQuality(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, 5, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeSourcePatternIndex(obj, value, varargin)
            if nargin==3, indices = value; value=varargin{1}; else  indices = getNodeIndices(obj,varargin); end
            j=1;
            for i=indices
                [obj.Errcode] = ENsetnodevalue(i, 6, value(j),obj.LibEPANET); j=j+1;
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
            end
        end
        function setNodeSourceType(obj, index, value)
            value=find(strcmpi(obj.TYPESOURCE,value)==1)-1;
            [obj.Errcode] = ENsetnodevalue(index, 7, value,obj.LibEPANET);
        end
        function setOptionsMaxTrials(obj,value)
            [obj.Errcode] = ENsetoption(0,value,obj.LibEPANET);
        end
        function setOptionsAccuracyValue(obj,value)
            [obj.Errcode] = ENsetoption(1,value,obj.LibEPANET);
        end
        function setOptionsQualityTolerance(obj,value)
            [obj.Errcode] = ENsetoption(2,value,obj.LibEPANET);
        end
        function setOptionsEmitterExponent(obj,value)
            [obj.Errcode] = ENsetoption(3,value,obj.LibEPANET);
        end
        function setOptionsPatternDemandMultiplier(obj,value)
            [obj.Errcode] = ENsetoption(4,value,obj.LibEPANET);
        end
        function setTimeSimulationDuration(obj,value)
            [obj.Errcode] = ENsettimeparam(0,value,obj.LibEPANET);
        end
        function setTimeHydraulicStep(obj,value)
            % Hstep = min(Pstep,Hstep)
            % Hstep = min(Rstep,Hstep)
            % Hstep = min(Qstep,Hstep)
            [obj.Errcode] = ENsettimeparam(1,value,obj.LibEPANET);
        end
        function setTimeQualityStep(obj,value)
            % Qstep = min(Qstep,Hstep)
            [obj.Errcode] = ENsettimeparam(2,value,obj.LibEPANET);
        end
        function setTimePatternStep(obj,value)
            [obj.Errcode] = ENsettimeparam(3,value,obj.LibEPANET);
        end
        function setTimePatternStart(obj,value)
            [obj.Errcode] = ENsettimeparam(4,value,obj.LibEPANET);
        end
        function setTimeReportingStep(obj,value)
            [obj.Errcode] = ENsettimeparam(5,value,obj.LibEPANET);
        end
        function setTimeReportingStart(obj,value)
            [obj.Errcode] = ENsettimeparam(6,value,obj.LibEPANET);
        end
        function setTimeRuleControlStep(obj,value)
            [obj.Errcode] = ENsettimeparam(7,value,obj.LibEPANET);
        end
        function setTimeStatisticsType(obj,value)
            %'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
            tmpindex=find(strcmpi(obj.TYPESTATS,value)==1)-1;
            [obj.Errcode] = ENsettimeparam(8,tmpindex,obj.LibEPANET);
        end
        function setTimeHTime(obj,value)
            % epanet20100
            [obj.Errcode] = ENsettimeparam(11,value,obj.LibEPANET);
        end
        function setTimeHaltFlag(obj,value)
            % epanet20100
            [obj.Errcode] = ENsettimeparam(12,value,obj.LibEPANET);
        end
%         function value = setTimeReportingPeriods(obj)
%             [obj.Errcode, value] = ENgettimeparam(9,obj.LibEPANET);
%         function value = setTimeStarting_Time(obj)
%             [obj.Errcode, value] = ENgettimeparam(10,obj.LibEPANET);
%         function value = setTimeNextEvent(obj)
%             [obj.Errcode, value] = ENgettimeparam(13,obj.LibEPANET);
%         end
        function setPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.Errcode] = ENsetpattern(index, patternVector, nfactors,obj.LibEPANET);
        end
        function setPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.Errcode] = ENsetpattern(i, patternMatrix(i,:), nfactors,obj.LibEPANET);
            end
        end
        function setPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.Errcode] = ENsetpatternvalue(index, patternTimeStep, patternFactor,obj.LibEPANET);
        end
        function setQualityType(obj,varargin)
            qualcode=0;chemname='';chemunits='';tracenode='';
            if find(strcmpi(varargin,'none')==1)
                [obj.Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.LibEPANET);
            elseif find(strcmpi(varargin,'age')==1)
                qualcode=2;
                [obj.Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.LibEPANET);
            elseif find(strcmpi(varargin,'chem')==1)
                qualcode=1;
                chemname=varargin{1};
                if nargin<3
                    chemunits='mg/L';
                else
                    chemunits=varargin{2};
                end
                [obj.Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.LibEPANET);
            elseif find(strcmpi(varargin,'trace')==1)
                qualcode=3;
                tracenode=varargin{2};
                [obj.Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.LibEPANET);
            else
                qualcode=1;
                chemname=varargin{1};
                if nargin<3
                    chemunits='mg/L';
                else
                    chemunits=varargin{2};
                end
                [obj.Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,obj.LibEPANET);                
            end
        end
        function setReportFormatReset(obj)
            [obj.Errcode]=ENresetreport(obj.LibEPANET);
        end
        function setReportStatus(obj,value)
            %'yes','no','full'
            statuslevel=find(strcmpi(obj.TYPEREPORT,value)==1)-1;
            [obj.Errcode] = ENsetstatusreport(statuslevel,obj.LibEPANET);
        end
        function setReport(obj,value)
            [obj.Errcode] = ENsetreport(value,obj.LibEPANET);
        end
        function closeNetwork(obj)
            [obj.Errcode] = ENclose(obj.LibEPANET);
        end
        function closeHydraulicAnalysis(obj)
            [obj.Errcode] = ENcloseH(obj.LibEPANET);
        end
        function closeQualityAnalysis(obj)
            [obj.Errcode] = ENcloseQ(obj.LibEPANET);
        end
        function saveHydraulicFile(obj,hydname)
            [obj.Errcode]=ENsavehydfile(hydname,obj.LibEPANET);
        end
        function useHydraulicFile(obj,hydname)
            [obj.Errcode]=ENusehydfile(hydname,obj.LibEPANET);
        end
        function initializeHydraulicAnalysis(obj)
            [obj.Errcode] = ENinitH(1,obj.LibEPANET);
        end
        function initializeQualityAnalysis(obj)
            [obj.Errcode] = ENinitQ(1,obj.LibEPANET);
        end
        function tstep = nextHydraulicAnalysisStep(obj)
            [obj.Errcode, tstep] = ENnextH(obj.LibEPANET);
        end
        function tstep = nextQualityAnalysisStep(obj)
            [obj.Errcode, tstep] = ENnextQ(obj.LibEPANET);
        end
        function openHydraulicAnalysis(obj)
            [obj.Errcode] = ENopenH(obj.LibEPANET);
        end
        function openQualityAnalysis(obj)
            [obj.Errcode] = ENopenQ(obj.LibEPANET);
        end
        function tstep = runHydraulicAnalysis(obj)
            [obj.Errcode, tstep] = ENrunH(obj.LibEPANET);
        end
        function tstep = runQualityAnalysis(obj)
            [obj.Errcode, tstep] = ENrunQ(obj.LibEPANET);
        end
        function saveHydraulicsOutputReportingFile(obj)
            [obj.Errcode] = ENsaveH(obj.LibEPANET);
        end
        function tleft=stepQualityAnalysisTimeLeft(obj)
            [obj.Errcode, tleft] = ENstepQ(obj.LibEPANET);
        end
        function Errcode = saveInputFile(obj,inpname,varargin)
            if strcmp(inpname,obj.BinTempfile) && nargin<3 %&& ~isempty(varargin)
                [addSectionCoordinates,addSectionRules] = obj.getBinCoordRuleSections;
                [obj.Errcode] = ENsaveinpfile(inpname,obj.LibEPANET);
                % Open epanet input file
                [~,info] = obj.readInpFile;
                endSectionIndex=find(~cellfun(@isempty,regexp(info,'END','match')));
                coordSectionIndex=find(~cellfun(@isempty,regexp(info,'COORDINATES','match')));
%                 verticesSectionIndex=find(~cellfun(@isempty,regexp(info,'VERTICES','match')));
                info(endSectionIndex)='';
                f1=fopen([obj.BinTempfile],'w');
                if ~isempty(coordSectionIndex)
                    fprintf(f1, '%s\n', info{1:coordSectionIndex-1});
                else
                    fprintf(f1, '%s\n', info{:});
                end
                if ~isempty(addSectionRules)
                    fprintf(f1, '%s\n', addSectionRules{:});
                end
                if ~isempty(addSectionCoordinates) % && isempty(coordSectionIndex)
                    fprintf(f1, '%s\n', addSectionCoordinates{:});
                else
                    fprintf(f1, '[END]\n');
                end
                fclose(f1);return;
            end
            [Errcode] = ENsaveinpfile(inpname,obj.LibEPANET);
            % The code below is because of a bug in EPANET 2.00.12
            % When saving using ENsaveinpfile, it does not save the type of the curves.
            obj.remAddBinCurvesID(inpname);
        end
        function writeLineInReportFile(obj, line)
            [obj.Errcode] = ENwriteline (line,obj.LibEPANET);
        end
        function writeReport(obj)
            %Writes a formatted text report on simulation results to the Report file
            [obj.Errcode]=ENreport(obj.LibEPANET);
        end
        function unload(obj)
            ENclose(obj.LibEPANET);
            ENMatlabCleanup(obj.LibEPANET);
%             if exist([obj.BinTempfile(1:end-4),'.bin'])==2
%                 delete([obj.BinTempfile(1:end-4),'.bin']);
%             end
%             delete(obj.BinTempfile);
%             if exist([obj.BinTempfile(1:end-4),'.txt'])==2
%                 delete([obj.BinTempfile(1:end-4),'.txt']);
%             end
%             [p,f]=fileparts(obj.InputFile);
%             if exist([p,'/',f,'.txt'])==2
%                 delete([p,'/',f,'.txt']);
%             end
%             if exist(obj.MSXTempFile)==2
%                 delete(obj.MSXTempFile);
%             end
            disp('EPANET Class is unloaded')
        end
        function loadMSXFile(obj,msxname,varargin)
            if isempty(varargin)
                MSXMatlabSetup(obj,msxname);
            else
                MSXMatlabSetup(obj,msxname,varargin);
            end
        end
        function runMSXexe(obj)
            [~,mm]=system(['cmd /c for %A in ("',pwd,'") do @echo %~sA']);
            mmPwd=regexp(mm,'\s','split');
            inpfile=obj.BinTempfile;
            rptfile=[inpfile(1:length(inpfile)-4),'.txt'];
            if strcmp(computer('arch'),'win64')
                    folder='64bit';
                r = sprintf('%s\\%s\\epanetmsx.exe %s %s %s',mmPwd{1},folder,inpfile,obj.MSXTempFile,rptfile);
            elseif strcmp(computer('arch'),'win32')
                    folder='32bit';
                r = sprintf('%s\\%s\\epanetmsx.exe %s %s %s',mmPwd{1},folder,inpfile,obj.MSXTempFile,rptfile);
            end
            system(r);
        end
        function value = getMSXEquationsTerms(obj)
            [value,~,~] = getEquations(obj.MSXFile);
        end
        function value = getMSXEquationsPipes(obj)
            [~,value,~] = getEquations(obj.MSXFile);
        end
        function value = getMSXEquationsTanks(obj)
            [~,~,value] = getEquations(obj.MSXFile);
        end
        function value = getMSXTimeStep(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.timestep;
        end
        function value = getMSXSolver(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.solver;
        end
        function value = getMSXAreaUnits(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.areaunits;
        end        
        function value = getMSXRateUnits(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.rateunits;
        end 
        function value = getMSXRtol(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.rtol;
        end 
        function value = getMSXAtol(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.atol;
        end
        function value = getMSXCoupling(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.coupling;
        end
        function value = getMSXCompiler(obj)
            [value] = getMSXOptions(obj.MSXFile);
            value = value.compiler;
        end
        function value = getMSXSpeciesCount(obj)
            % Species, Constants, Parameters, Patterns
            [obj.Errcode, value] = MSXgetcount(3,obj.MSXLibEPANET);
        end
        function value = getMSXConstantsCount(obj)
            [obj.Errcode, value] = MSXgetcount(6,obj.MSXLibEPANET);
        end
        function value = getMSXParametersCount(obj)
            [obj.Errcode, value] = MSXgetcount(5,obj.MSXLibEPANET);
        end
        function value = getMSXPatternsCount(obj)
            [obj.Errcode, value] = MSXgetcount(7,obj.MSXLibEPANET);
        end
        function value = getMSXSpeciesNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMSXSpeciesCount
                    [obj.Errcode, len] = MSXgetIDlen(3,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{i}]=MSXgetID(3,i,len,obj.MSXLibEPANET);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(3,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{k}]=MSXgetID(3,i,len,obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXSpeciesType(obj)
            if obj.getMSXSpeciesCount
                for i=1:obj.getMSXSpeciesCount
                    [obj.Errcode,value{i},~,~,~] = MSXgetspecies(i,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXSpeciesUnits(obj)
            if obj.getMSXSpeciesCount
                for i=1:obj.getMSXSpeciesCount
                    [obj.Errcode,~,value{i},~,~] = MSXgetspecies(i,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXSpeciesATOL(obj)
            if obj.getMSXSpeciesCount
                for i=1:obj.getMSXSpeciesCount
                    [obj.Errcode,~,~,value,~] = MSXgetspecies(i,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXSpeciesRTOL(obj)
            if obj.getMSXSpeciesCount
                for i=1:obj.getMSXSpeciesCount
                    [obj.Errcode,~,~,~,value] = MSXgetspecies(i,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXSpeciesIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMSXSpeciesCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = MSXgetindex(3,varargin{1}{j},obj.MSXLibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET,3,varargin{1});
            end
        end
        function value = getMSXConstantsNameID(obj)
            if obj.getMSXConstantsCount
                for i=1:obj.getMSXConstantsCount
                    [obj.Errcode, len] = MSXgetIDlen(6,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{i}] = MSXgetID(6,i,len,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXConstantsValue(obj)
            if obj.getMSXConstantsCount
                for i=1:obj.getMSXConstantsCount
                    [obj.Errcode, value(i)] = MSXgetconstant(i,obj.MSXLibEPANET);
                end
            else
                value=-1;
            end
        end
        function value = getMSXConstantsIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMSXConstantsCount;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = MSXgetindex(6,varargin{1}{j},obj.MSXLibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET,6,varargin{1});
            end
        end
        function value = getMSXParametersNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMSXParametersCount
                    [obj.Errcode, len] = MSXgetIDlen(5,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{i}]=MSXgetID(5,i,len,obj.MSXLibEPANET);
                end
                if ~obj.getMSXParametersCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(5,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{k}]=MSXgetID(5,i,len,obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXParametersIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMSXParametersCount;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = MSXgetindex(5,varargin{1}{j},obj.MSXLibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = MSXgetindex(5,varargin{1},obj.MSXLibEPANET);
            end
        end
        function value = getMSXParametersTanksValue(obj)
            value={};
            if ~obj.getMSXParametersCount, value=0;return;end
            if ~length(obj.NodeTankIndex), value=0;return;end
            for i=1:length(obj.NodeTankIndex)
                for j=1:obj.MSXParametersCount%
                    [obj.Errcode, value{obj.NodeTankIndex(i)}(j)] = MSXgetparameter(0,obj.NodeTankIndex(i),j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXParametersPipesValue(obj)
            if ~obj.getMSXParametersCount
                value=0;return;
            end
            for i=1:obj.getLinkPipeCount
                for j=1:obj.getMSXParametersCount
                    [obj.Errcode, value{i}(j)] = MSXgetparameter(1,i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXPatternsNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMSXPatternsCount
                    [obj.Errcode, len] = MSXgetIDlen(7,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{i}]=MSXgetID(7,i,len,obj.MSXLibEPANET);
                end
                if ~obj.getMSXPatternsCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.Errcode, len] = MSXgetIDlen(7,i,obj.MSXLibEPANET);
                    [obj.Errcode, value{k}]=MSXgetID(7,i,len,obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXPatternsIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getMSXPatternsCount;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET,7,varargin{1});
                    if obj.Errcode
                        value{k}=0;
                    end
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = MSXgetindex(obj.MSXLibEPANET,7,varargin{1});
                if obj.Errcode
                    value=0;
                end
            end
        end
        function value = getMSXPatternsLengths(obj,varargin)
            if isempty(varargin)
                if obj.getMSXPatternsCount
                    for i=obj.getMSXPatternsIndex
                        [obj.Errcode, value(i)]=MSXgetpatternlen(i,obj.MSXLibEPANET);
                    end
                else
                    value=-1;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.Errcode, value(k)] = MSXgetpatternlen(obj.getMSXPatternsIndex(varargin{1}{j}),obj.MSXLibEPANET);
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.Errcode, value] = MSXgetpatternlen(obj.getMSXPatternsIndex(varargin{1}),obj.MSXLibEPANET);
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.Errcode, value(k)]=MSXgetpatternlen(i,obj.MSXLibEPANET);
                    k=k+1;
                end
            end
        end
        function value = getMSXNodeInitqualValue(obj)
            if obj.getMSXSpeciesCount==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, value{i}(j)] = MSXgetinitqual(0,i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXLinkInitqualValue(obj)
            if obj.getMSXSpeciesCount==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.getLinkCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, value{i}(j)] = MSXgetinitqual(1,i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXSources(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode, obj.MSXSourceType{i}{j},obj.MSXSourceLevel{i}(j),obj.MSXSourcePatternIndex{i}(j)] = MSXgetsource(i,j,obj.MSXLibEPANET);
                    obj.MSXSourceTypeCode{i}(j)=find(strcmp(obj.TYPESOURCEMSX,obj.MSXSourceType{i}{j}))-2;
                end
            end
            SnodeID=obj.getMSXSourceNodeNameID;
           % value={obj.MSXSourceType,obj.MSXSourceTypeCode,obj.MSXSourceLevel,obj.MSXSourcePatternIndex,SnodeID};
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
                    [obj.Errcode, value{i}{j},~,~] = MSXgetsource(i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXSourceLevel(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode,~,value{i}(j),~] = MSXgetsource(i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXSourcePatternIndex(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMSXSpeciesCount
                    [obj.Errcode,~,~,value{i}(j)] = MSXgetsource(i,j,obj.MSXLibEPANET);
                end
            end
        end
        function value = getMSXPattern(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getMSXPatternsLengths);
            value=nan(obj.getMSXPatternsCount,tmpmaxlen);
            for i=1:obj.getMSXPatternsCount
                tmplength=obj.getMSXPatternsLengths(i);
                for j=1:tmplength
                    [obj.Errcode, value(i,j)] = MSXgetpatternvalue(i,j,obj.MSXLibEPANET);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                
            end
        end
        function value = getMSXPatternValue(obj,patternIndex, patternStep) %Mass flow rate per minute of a chemical source
            [obj.Errcode, value] = MSXgetpatternvalue(patternIndex, patternStep, obj.MSXLibEPANET);
        end
        function value = getMSXSpeciesConcentration(obj, type, index, species)
            [obj.Errcode, value] = MSXgetqual(type, index, species, obj.MSXLibEPANET);
        end
        function value = getMSXComputedQualityNode(obj,varargin)
            if obj.getMSXSpeciesCount==0
                value=0;
                return;
            end
            if ~isempty(varargin)
                if length(varargin)==1
                    ss=varargin{1};%index node
                    uu=1:obj.getMSXSpeciesCount;
                elseif length(varargin)==2
                    ss=varargin{1};%index node
                    uu=varargin{2};%index species
                end
            else
                ss=1:obj.getNodeCount;%index node
                uu=1:obj.getMSXSpeciesCount;
            end
            % Obtain a hydraulic solution
            obj.solveMSXCompleteHydraulics(obj.MSXLibEPANET);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.initializeMSXQualityAnalysis(0);
            % Retrieve species concentration at node
            k=2; tleft=1;t=0;i=1;nl=1;
            value.Time(k,:)=0;
            timeSmle=obj.getTimeSimulationDuration;%bug at time
            while(tleft>0 && obj.Errcode==0 && timeSmle~=t)
                [t, tleft]=obj.stepMSXQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<3600 || t==3600
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMSXNodeInitqualValue{ss}(j);
                        end
                    else
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMSXSpeciesConcentration(0, ss, j);%node code0
                        end
                    end
                else
                    if t<3600 || t==3600
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMSXNodeInitqualValue{(nl)}(j);
                            end
                        end
                    else
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMSXSpeciesConcentration(0, (nl), j);%node code0
                            end
                        end
                    end
                    nl=nl+1;
                end
                value.Time(k,:)=t;
                k=k+1;
            end
        end
        function value = getMSXComputedQualityLink(obj,varargin)
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
            % Obtain a hydraulic solution
            obj.solveMSXCompleteHydraulics(obj.MSXLibEPANET);
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.initializeMSXQualityAnalysis(0);
            
            % Retrieve species concentration at node
            k=2;tleft=1;i=1;
            value.Time(k,:)=0;
            while(tleft>0 && obj.Errcode==0)
                [t, tleft]=obj.stepMSXQualityAnalysisTimeLeft;
                if ~isempty(varargin)
                    if t<3600 || t==3600
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMSXLinkInitqualValue{ss}(j);
                        end
                    else
                        for j=uu
                            value.Quality{j,i}(k,:)=obj.getMSXSpeciesConcentration(1, ss, j);
                        end
                    end
                else
                    if t<3600 || t==3600
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMSXLinkInitqualValue{(nl)}(j);
                            end
                        end
                    else
                        for nl=ss
                            for j=uu
                                value.Quality{nl}{j}(k)=obj.getMSXSpeciesConcentration(1, (nl), j);%link code1
                            end
                        end
                    end
                    nl=nl+1;
                end
                value.Time(k,:)=t;
                k=k+1;
            end
        end
        function plotMSXSpeciesNodeConcentration(obj,varargin)
            s=obj.getMSXComputedQualityNode(varargin{1},varargin{2});
            nodesID=obj.getNodeNameID;
            SpeciesNameID=obj.getMSXSpeciesNameID;
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
        function plotMSXSpeciesLinkConcentration(obj,varargin)
            s=obj.getMSXComputedQualityLink(varargin{1},varargin{2});
            linksID=obj.getLinkNameID;
            SpeciesNameID=obj.getMSXSpeciesNameID;
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
        function value = getMSXError(obj,Errcode)
            [obj.Errcode, value] = MSXgeterror(Errcode,obj.MSXLibEPANET);
        end
        function solveMSXCompleteHydraulics(obj,varargin)
            [obj.Errcode] = MSXsolveH(obj.MSXLibEPANET);
        end
        function solveMSXCompleteQuality(obj,varargin)
            [obj.Errcode] = MSXsolveQ(obj.MSXLibEPANET);
        end
        function writeMSXReport(obj,varargin)
            [obj.Errcode]=MSXreport(obj.MSXLibEPANET);
        end
        function index = addMSXPattern(obj,varargin)
            index=-1;
            if nargin==2
                [obj.Errcode] = MSXaddpattern(varargin{1},obj.MSXLibEPANET);
                [obj.Errcode, index] = MSXgetindex(obj.MSXLibEPANET,7,varargin{1});
            elseif nargin==3
                [obj.Errcode] = MSXaddpattern(varargin{1},obj.MSXLibEPANET);
                [obj.Errcode, index] = MSXgetindex(obj.MSXLibEPANET,7,varargin{1});
                setMSXPattern(obj,index,varargin{2});
            end
        end
        function setMSXSources(obj, node, species, type, level, pat)
            MSXsetsource(node, species, type, level, pat, obj.MSXLibEPANET);
        end
        function setMSXConstantsValue(obj, value)
            for i=1:length(value)
                [obj.Errcode] = MSXsetconstant(i, value(i), obj.MSXLibEPANET);
            end
        end
        function setMSXParametersTanksValue(obj, NodeTankIndex, paramindex, value)
            if ~sum(NodeTankIndex==obj.NodeTankIndex)
                fprintf('>> Invalid Tank Index <<\n');obj.NodeTankIndex
                return;
            end
            [obj.Errcode] = MSXsetparameter(0, NodeTankIndex, paramindex, value, obj.MSXLibEPANET);
        end
        function setMSXParametersPipesValue(obj, pipeIndex, value)
            for i=1:length(value)
                [obj.Errcode] = MSXsetparameter(1, pipeIndex, i, value(i), obj.MSXLibEPANET);
            end
        end
        function setMSXNodeInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.Errcode] = MSXsetinitqual(0, i, j, value{i}(j), obj.MSXLibEPANET);
                end
            end
        end
        function setMSXLinkInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.Errcode] = MSXsetinitqual(1, i, j, value{i}(j), obj.MSXLibEPANET);
                end
            end
        end
        function setMSXPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.Errcode] = MSXsetpattern(index, patternVector, nfactors, obj.MSXLibEPANET);
        end
        function setMSXPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.Errcode] = MSXsetpattern(i, patternMatrix(i,:), nfactors, obj.MSXLibEPANET);
            end
        end
        function setMSXPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.Errcode] = MSXsetpatternvalue(index, patternTimeStep, patternFactor, obj.MSXLibEPANET);
        end
        function setMSXTimeStep(obj,timestep)
            setMSXOptions(obj,'timestep',timestep);
        end
        function setMSXAreaUnitsFT2(obj)
            setMSXOptions(obj,'areaunits','FT2');
        end
        function setMSXAreaUnitsM2(obj)
            setMSXOptions(obj,'areaunits','M2');
        end
        function setMSXAreaUnitsCM2(obj)
            setMSXOptions(obj,'areaunits','CM2');
        end
        function setMSXRateUnitsSEC(obj)
            setMSXOptions(obj,'rateunits','SEC');
        end
        function setMSXRateUnitsMIN(obj)
            setMSXOptions(obj,'rateunits','MIN');
        end
        function setMSXRateUnitsHR(obj)
            setMSXOptions(obj,'rateunits','HR');
        end        
        function setMSXRateUnitsDAY(obj)
            setMSXOptions(obj,'rateunits','DAY');
        end
        function setMSXSolverEUL(obj)
            setMSXOptions(obj,'solver','EUL');
        end
        function setMSXSolverRK5(obj)
            setMSXOptions(obj,'solver','RK5');
        end        
        function setMSXSolverROS2(obj)
            setMSXOptions(obj,'solver','ROS2');
        end
        function setMSXCouplingFULL(obj)
            setMSXOptions(obj,'coupling','FULL');
        end
        function setMSXCouplingNONE(obj)
            setMSXOptions(obj,'coupling','NONE');
        end
        function setMSXCompilerNONE(obj)
            setMSXOptions(obj,'compiler','NONE');
        end
        function setMSXCompilerVC(obj)
            setMSXOptions(obj,'compiler','VC');
        end
        function setMSXCompilerGC(obj)
            setMSXOptions(obj,'compiler','GC');
        end
        function setMSXAtol(obj,atol)
            setMSXOptions(obj,'atol',atol);
        end
        function setMSXRtol(obj,rtol)
            setMSXOptions(obj,'rtol',rtol);
        end        
        function saveMSXQualityFile(obj,outfname)
            [obj.Errcode]=MSXsaveoutfile(outfname,obj.MSXLibEPANET);
        end
        function useMSXHydraulicFile(obj,hydname)
            [obj.Errcode]=MSXusehydfile(hydname,obj.MSXLibEPANET);
        end
        function initializeMSXQualityAnalysis(obj,flag)
            [obj.Errcode] = MSXinit(flag,obj.MSXLibEPANET);
        end
        function [t, tleft]= stepMSXQualityAnalysisTimeLeft(obj)
            [obj.Errcode, t, tleft] = MSXstep(obj.MSXLibEPANET);
        end
        function saveMSXFile(obj,msxname)
            [obj.Errcode] = MSXsavemsxfile(msxname,obj.MSXLibEPANET);
        end
        function msx = writeMSXFile(obj,msx)
            % Input Arguments
            % msx={};
            % msx.msxFile = 'networks/MSXFileName.msx';
            % % section Title
            % msx.titleDescription{1} = 'Example: MSXFile.';
            % % section Options
            % msx.options{1}='FT2'; %AREA_UNITS FT2/M2/CM2
            % msx.options{2}='DAY'; %TIME_UNITS SEC/MIN/HR/DAY
            % msx.options{3}='EUL'; %SOLVER EUL/RK5/ROS2
            % msx.options{4}='NONE'; %COUPLING FULL/NONE
            % msx.options{5}='NONE'; %COMPILER NONE/VC/GC
            % msx.options{6}=3600; %TIMESTEP in seconds
            % msx.options{7}=0.01;  %ATOL value
            % msx.options{8}=0.001;  %RTOL value
            % % section Species
            % % <type> <specieID> <units> (<atol> <rtol>)
            % msx.species{1}={'BULK'}; %type BULK/WALL
            % msx.species{2}={'CL2'}; %specieID
            % msx.species{3}={'MG'}; %units UG/MG
            % msx.species{4}={0.01}; %atol
            % msx.species{5}={0.001}; %rtol
            % 
            % % section Coefficients 
            % % CONSTANT name value % PARAMETER name value
            % msx.coefficients{1}={'PARAMETER','PARAMETER'}; 
            % msx.coefficients{2}={'Kb','Kw'}; 
            % msx.coefficients{3}={0.3,1}; 
            % 
            % % section Terms
            % % <termID> <expression>
            % msx.terms{1}={'Kf'}; % termID
            % msx.terms{2}={'1.5826e-4 * RE^0.88 / D'}; % expression
            % 
            % % section Pipes
            % % EQUIL <specieID> <expression>
            % % RATE <specieID> <expression>
            % % FORMULA <specieID> <expression>
            % msx.pipes{1} ={'RATE'}; %type
            % msx.pipes{2} ={'CL2'}; %specieID
            % msx.pipes{3} ={'-Kb*CL2 - (4/D)*Kw*Kf/(Kw+Kf)*CL2'}; %expression
            % 
            % % section Tanks
            % % EQUIL <specieID> <expression>
            % % RATE <specieID> <expression>
            % % FORMULA <specieID> <expression>
            % msx.tanks{1} ={'RATE'}; %type
            % msx.tanks{2} ={'CL2'}; %specieID
            % msx.tanks{3} ={'-Kb*CL2'}; %expression

            % % section Sources
            % % <type> <nodeID> <specieID> <strength> (<patternID>)
            % msx.sources{1}={''}; %CONC/MASS/FLOW/SETPOINT
            % msx.sources{2}={''}; %nodeID
            % msx.sources{3}={''}; %specieID
            % msx.sources{4}={''}; %strength
            % msx.sources{5}={''}; %patternID
            % 
            % % section Quality Global
            % % GLOBAL <specieID> <value>
            % msx.global{1} = {''};
            % msx.global{2} = {''};%specieID
            % msx.global{3} = {''};%value
            % % others
            % % NODE <nodeID> <bulkSpecieID> <value>
            % % LINK <linkID> <wallSpecieID> <value>
            % msx.quality{1} = {''}; %NODE/LINK
            % msx.quality{2} = {''}; %ID
            % msx.quality{3} = {''}; %bulkSpecieID/wallSpecieID
            % msx.quality{4} = {''}; %value
            % 
            % % section Parameters
            % % PIPE <pipeID> <paramID> <value>
            % % TANK <tankID> <paramID> <value>
            % msx.parameters{1} = {''};
            % msx.parameters{2} = {''};
            % msx.parameters{3} = {''};
            % msx.parameters{4} = {''};
            % 
            % % section Patterns
            % % <patternID> <multiplier> <multiplier> 
            % msx.patterns{1} = {''}; %patternID
            % msx.patterns{2} = {''}; %multiplier
            space=5;
            f = fopen(msx.msxFile,'w');
            fprintf(f,'[TITLE]\n');
            if isfield(msx,'titleDescription')
                for i=1:length(msx.titleDescription)
                    fprintf(f,msx.titleDescription{i});
                end
            end

            fprintf(f,'\n\n[OPTIONS]\n');
            options = {'AREA_UNITS', 'RATE_UNITS', 'SOLVER', 'COUPLING', 'COMPILER',...
                'TIMESTEP', 'ATOL', 'RTOL'};
            spaces=blanks(space);
            if isfield(msx,'options')
                for i=1:length(msx.options)
                    fprintf(f,num2str(options{i}));
                    fprintf(f,spaces);
                    fprintf(f,num2str(msx.options{i}));
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[SPECIES]\n');
            if isfield(msx,'species')
                for u=1:length(msx.species{1})
                    for i=1:length(msx.species)
                        if isnumeric(msx.species{i}{u})
                            fprintf(f,num2str(msx.species{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.species{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end
            fprintf(f,'\n[COEFFICIENTS]\n');
            if isfield(msx,'coefficients')
                for u=1:length(msx.coefficients{1})
                    for i=1:length(msx.coefficients)
                        if isnumeric(msx.coefficients{i}{u})
                            fprintf(f,num2str(msx.coefficients{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.coefficients{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[TERMS]\n');
            if isfield(msx,'terms')
                for u=1:length(msx.terms{1})
                    for i=1:length(msx.terms)
                        if isnumeric(msx.terms{i}{u})
                            fprintf(f,num2str(msx.terms{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.terms{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[PIPES]\n');
            if isfield(msx,'pipes')
                for u=1:length(msx.pipes{1})
                    for i=1:length(msx.pipes)
                        if isnumeric(msx.pipes{i}{u})
                            fprintf(f,num2str(msx.pipes{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.pipes{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[TANKS]\n');
            if isfield(msx,'tanks')
                for u=1:length(msx.tanks{1})
                    for i=1:length(msx.tanks)
                        if isnumeric(msx.tanks{i}{u})
                            fprintf(f,num2str(msx.tanks{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.tanks{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[SOURCES]\n');
            if isfield(msx,'sources')
                for u=1:length(msx.sources{1})
                    for i=1:length(msx.sources)
                        if isnumeric(msx.sources{i}{u})
                            fprintf(f,num2str(msx.sources{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.sources{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'\n[QUALITY]\n');
            if isfield(msx,'global')
                for u=1:length(msx.global{1})
                    for i=1:length(msx.global)
                        if isnumeric(msx.global{i}{u})
                            fprintf(f,num2str(msx.global{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.global{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end
            if isfield(msx,'quality')
                for u=1:length(msx.quality{1})
                    for i=1:length(msx.quality)
                        if isnumeric(msx.quality{i}{u})
                            fprintf(f,num2str(msx.quality{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.quality{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'[PARAMETERS]\n');
            if isfield(msx,'parameters')
                for u=1:length(msx.parameters{1})
                    for i=1:length(msx.parameters)
                        if isnumeric(msx.parameters{i}{u})
                            fprintf(f,num2str(msx.parameters{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.parameters{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'[PATTERNS]\n');
            if isfield(msx,'patterns')
                for u=1:length(msx.patterns{1})
                    for i=1:length(msx.patterns)
                        if isnumeric(msx.patterns{i}{u})
                            fprintf(f,num2str(msx.patterns{i}{u}));
                        else
                            fprintf(f,num2str(char(msx.patterns{i}{u})));
                        end
                        fprintf(f,spaces);
                    end
                    fprintf(f,'\n');
                end
            end

            fprintf(f,'[REPORT]\n');
            fprintf(f,'NODES ALL\n');
            fprintf(f,'LINKS ALL\n');
            fclose(f)
        end
        function unloadMSX(obj)
            MSXclose(obj);
            MSXMatlabCleanup(obj);
        end
        function obj = BinUpdateClass(obj)
            sect=0;i=1;t=1;q=1;
            typecode=0;x=1;b=1;d=1;
            if obj.Bin
                obj.saveInputFile([obj.BinTempfile]);
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
                        obj.BinCountStatuslines=0;
                        continue;
                        % [PATTERNS]
                    elseif strcmpi(tok(1:5),'[PATT')
                        sect=10;
                        obj.BinPatternLengths=[];
                        obj.BinPatternNameID={};
                        obj.BinPatternValue={};
                        obj.BinCountPatternlines=0;d=1;h=1;
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
                        obj.BinCountInitialQualitylines=0;
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
                        obj.BinCountReactionlines=0;                        
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
                    obj.BinCountInitialQualitylines=d;
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
                    obj.BinCountReactionlines=d;
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
        function [Errcode]=setBinNodeInitialQuality(obj,varargin)
            parameter=varargin{1};
            zz=obj.BinNodeCount-obj.BinCountInitialQualitylines+1;
            sections={'[QUALITY]','[SOURCES]'};
            [Errcode]=setBinParam2(obj,parameter,sections,zz);   
        end
        function [Errcode]=setBinLinkReactionCoeff(obj,varargin)
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
           fid = fopen([obj.BinTempfile],'w');start=0;
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
               Errcode=closeOpenNetwork(obj);
           end
        end
        function [Errcode]=setBinQualityChem(obj,varargin)
            sections={'[OPTIONS]','[REPORT]'};
            indexParameter=1;
            parameter='Quality            	chem   mg/L';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [Errcode]=setBinQualityNone(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	None';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [Errcode]=setBinQualityAge(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1;
            parameter='Quality            	Age';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [Errcode]=setBinQualityTrace(obj,varargin)
            sections={'[OPTIONS]','[COORDINATES]'};
            indexParameter=1; 
            if ~sum(strcmp(varargin{1},obj.BinNodeNameID))
                warning('Invalid property found.');
                return
            end
            parameter=['Quality            	Trace ',varargin{1}];
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);          
        end
        function [Errcode]=setBinTimeSimulationDuration(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=1;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeHydraulicStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeQualityStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimePatternStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimePatternStart(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=5;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeReportingStep(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=6;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeReportingStart(obj,varargin)
            parameter=varargin{1};
            sections={'[TIMES]','[REPORT]'};
            indexParameter=7;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeStatisticsNone(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='None';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeStatisticsAverage(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='AVERAGE';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeStatisticsMinimum(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='MINIMUM';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeStatisticsMaximum(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='MAXIMUM';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinTimeStatisticsRange(obj,varargin)
            sections={'[TIMES]','[REPORT]'};
            indexParameter=8;
            parameter='RANGE';
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankInitLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankMinLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankMaxLevel(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=5;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankDiameter(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=6;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeTankMinVol(obj,varargin)
            parameter=varargin{1};
            sections={'[TANKS]','[PIPES]'};
            indexParameter=7;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinLinkGlobalWallReactionCoeff(obj,varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]','[MIXING]'};
            [Errcode]=setBinParam(obj,3,parameter,sections);        
        end
        function [Errcode]=setBinLinkGlobalBulkReactionCoeff(obj,varargin)
            parameter=varargin{1};
            sections={'[REACTIONS]','[MIXING]'};
            [Errcode]=setBinParam(obj,1,parameter,sections);        
        end
        function BinClose(obj)
            fclose all;
            if exist([obj.BinTempfile(1:end-4),'.bin'])==2
                delete([obj.BinTempfile(1:end-4),'.bin'])
            end
            delete(obj.BinTempfile)
            delete([obj.BinTempfile(1:end-4),'.txt'])
            if exist([obj.BinTempfile(1:end-4),'.msx'])==2
                delete([obj.BinTempfile(1:end-4),'.msx'])
            end
        end
        function [Errcode]=setBinLinkValvesParameters(obj,varargin)
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
                            warning('Invalid property found.');Errcode=-1;
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
                            warning('Invalid argument found.');Errcode=-1;
                            return;
                        end   
                        zz=abs(obj.BinLinkPumpCount+obj.BinLinkValveCount-obj.BinCountStatuslines);
                        sections={'[STATUS]','[PATTERNS]','valve'};
                        setBinParam2(obj,Status,sections,zz);   
                    otherwise
                        warning('Invalid property found.');Errcode=-1;
                        return;
                end
            end            
            [tlines]=regexp( fileread(obj.BinTempfile), '\n', 'split');
            fid = fopen([obj.BinTempfile],'w');
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
                Errcode=closeOpenNetwork(obj);
            end      
        end
        function [Errcode]=setBinNodeResDemandPatternNameID(obj,varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]','[TANKS]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeReservoirElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[RESERVOIRS]','[TANKS]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);        
        end
        function [Errcode]=setBinNodeJunctionElevation(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=2;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinNodeJunctionsBaseDemands(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=3;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinNodeDemandPatternNameID(obj,varargin)
            parameter=varargin{1};
            sections={'[JUNCTIONS]','[RESERVOIRS]','[DEMANDS]','[STATUS]','[EMITTERS]'};
            indexParameter=4;
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinPattern(obj,varargin)
            idpattern=varargin{1};
            %ind=obj.getPatternIndex(idpattern);
            values=varargin{2}; %obj.getPatternValue{ind};
            sections={'[PATTERNS]','[CURVES]'};
            v=obj.getBinPatternsInfo;
            patterns=v.BinPatternNameID;
            if sum(strcmp(idpattern,patterns))
                Errcode=setBinParam(obj,idpattern,values,sections);
            else
                warning('Invalid argument found.');Errcode=-1;
                return
            end
        end
        function [Errcode]=addBinPattern(obj,varargin)
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
                Errcode=setBinParam2(obj,values,sections,zz,newidpattern);
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
            end
        end
        function [Errcode]=setBinNodeSourceQuality(obj,varargin)
            sections={'[SOURCES]','[MIXING]'};
            values=varargin{1};
            [Errcode]=setBinParam(obj,11,values,sections);
        end
        function saveBinInpFile(obj)
            [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
            f = fopen([obj.BinTempfile],'w');
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
        function [Errcode]=addBinJunction(obj,varargin)
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
            Errcode=addLinkWarnings(obj,typecode,varargin{7},toNode); 
            if Errcode==-1, return; end
            [Errcode]=addNode(obj,0,newID,X,Y,newElevation,newBaseDemand,newDemandPattern);
            if Errcode~=200 && obj.Bin==1
                return;
            end
            if typecode==1
                newPipeID=varargin{7}; 
                newLength=varargin{9};
                newDiameter=varargin{10};
                newRoughness=varargin{11};
                [Errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{7}; 
                newCurveIDofPump=varargin{9}; 
                newCurveXvalue=varargin{10}; 
                newCurveYvalue=varargin{11}; 
                newCurveType=varargin{12};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [Errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{7}; 
                newValveDiameter=varargin{9}; 
                newValveSetting=varargin{10}; 
                [Errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [Errcode]=addBinReservoir(obj,varargin)
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
            Errcode=addLinkWarnings(obj,typecode,varargin{5},toNode); 
            if Errcode==-1, return; end            
            [Errcode]=addNode(obj,1,newID,X,Y,newElevation);
            if Errcode
                return;
            end
            if strcmp(upper(varargin{end}),'PIPE')
                newPipeID=varargin{5}; 
                newLength=varargin{7};
                newDiameter=varargin{8};
                newRoughness=varargin{9};
                [Errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{5}; 
                newCurveIDofPump=varargin{7}; 
                newCurveXvalue=varargin{8}; 
                newCurveYvalue=varargin{9}; 
                newCurveType=varargin{10};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [Errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{5}; 
                newValveDiameter=varargin{7}; 
                newValveSetting=varargin{8}; 
                [Errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [Errcode]=addBinTank(obj,varargin)
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
            Errcode=addLinkWarnings(obj,typecode,varargin{11},toNode); 
            if Errcode==-1, return; end            
            [Errcode]=addNode(obj,2,newID,X,Y,MaxLevel,Diameter,Initlevel,newElevation,initqual,MinLevel,MinVol);
            if Errcode
                return;
            end
            if strcmp(upper(varargin{end}),'PIPE')
                newPipeID=varargin{11}; 
                newLength=varargin{13};
                newDiameter=varargin{14};
                newRoughness=varargin{15};
                [Errcode]=addBinPipe(obj,newPipeID,fromNode,toNode,newLength,newDiameter,newRoughness);
            elseif typecode==2
                newPumpID=varargin{11}; 
                newCurveIDofPump=varargin{13}; 
                newCurveXvalue=varargin{14}; 
                newCurveYvalue=varargin{15}; 
                newCurveType=varargin{16};  % PUMP, EFFICIENCY, VOLUME, HEADLOSS  
                [Errcode]=addBinPump(obj,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            else 
                newValveID=varargin{11}; 
                newValveDiameter=varargin{13}; 
                newValveSetting=varargin{14}; 
                [Errcode]=addLink(obj,typecode,newValveID,fromNode,toNode,newValveDiameter,newValveSetting);
            end
        end
        function [Errcode]=addBinPipe(obj,newLink,fromNode,toNode,newLength,newDiameter,newRoughness)
            [Errcode]=addLink(obj,1,newLink,fromNode,toNode,newLength,newDiameter,newRoughness);
        end
        function [Errcode]=addBinPump(obj,newPumpID,fromNode,toNode,varargin)
            Errcode=-1;
            if length(varargin)==4
                newCurveIDofPump=varargin{1};
                newCurveXvalue=varargin{2};
                newCurveYvalue=varargin{3};
                newCurveType=varargin{4};
                if strcmp(upper(newCurveType),'PUMP')
                    Errcode=addBinCurvePump(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Head
                elseif strcmp(upper(newCurveType),'EFFICIENCY')
                    Errcode=addBinCurveEfficiency(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Efficiency
                elseif strcmp(upper(newCurveType),'VOLUME')
                    Errcode=addBinCurveVolume(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Heigh-Volume
                elseif strcmp(upper(newCurveType),'HEADLOSS')
                    Errcode=addBinCurveHeadloss(obj,newCurveIDofPump,newCurveXvalue,newCurveYvalue);%Flow-Headloss
                end
                [Errcode]=addLink(obj,2,newPumpID,fromNode,toNode,newCurveIDofPump,newCurveXvalue,newCurveYvalue,newCurveType);
            elseif length(varargin)==1
                power=varargin{1};
                [Errcode]=addLink(obj,2,newPumpID,fromNode,toNode,power);
            end
        end
        function [Errcode]=addBinValvePRV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,3,newLink,fromNode,toNode,diameter,setting); % Pressure Reducing Valve
        end
        function [Errcode]=addBinValvePSV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,4,newLink,fromNode,toNode,diameter,setting); % Pressure Sustaining Valve
        end
        function [Errcode]=addBinValvePBV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,5,newLink,fromNode,toNode,diameter,setting); % Pressure Breaker Valve
        end
        function [Errcode]=addBinValveFCV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,6,newLink,fromNode,toNode,diameter,setting); % Flow Control Valve
        end
        function [Errcode]=addBinValveTCV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,7,newLink,fromNode,toNode,diameter,setting); % Throttle Control Valve
        end
        function [Errcode]=addBinValveGPV(obj,newLink,fromNode,toNode,diameter,setting)
            [Errcode]=addLink(obj,8,newLink,fromNode,toNode,diameter,setting); % General Purpose Valve
        end
        function [Errcode]=addBinCurvePump(obj,newCurveID,varargin)
            CurveX=varargin{1};
            CurveY=varargin{2};
            [Errcode]=addCurve(obj,newCurveID,CurveX,CurveY,0);  %ID Flow-OptionsHeadloss
        end
        function [Errcode]=addBinCurveEfficiency(obj,newCurveID,CurveX,CurveY)
            [Errcode]=addCurve(obj,newCurveID,CurveX,CurveY,1);  %ID Flow-Efficiency
        end
        function [Errcode]=addBinCurveVolume(obj,newCurveID,CurveX,CurveY)
            [Errcode]=addCurve(obj,newCurveID,CurveX,CurveY,2);  %ID Heigh-Volume
        end
        function [Errcode]=addBinCurveHeadloss(obj,newCurveID,CurveX,CurveY)
            [Errcode]=addCurve(obj,newCurveID,CurveX,CurveY,3);  %ID Flow-OptionsHeadloss
        end
        function [Errcode]=addBinControl(obj,x,status,y_t_c,param,z,varargin)
            if nargin==6
                [Errcode]=addNewControl(obj,x,status,y_t_c,param,z);
            elseif nargin==5
                [Errcode]=addNewControl(obj,x,status,y_t_c,param);
            else
                [Errcode]=addNewControl(obj,x,status,y_t_c);
            end
        end
        function [Errcode]=removeBinNodeID(obj,NodeID)
            [Errcode]=rmNode(obj,NodeID);
        end
        function [Errcode]=removeBinCurveID(obj,CurveID)
            [Errcode]=rmCurveID(obj,CurveID);
        end
        function [Errcode]=removeBinLinkID(obj,LinkID)
            [Errcode]=rmLink(obj,LinkID);
        end
        function [Errcode]=removeBinControlLinkID(obj,ID)
            [Errcode]=rmControl(obj,1,ID);
        end
        function [Errcode]=removeBinRulesControlLinkID(obj,ID)
            [Errcode]=rmRulesControl(obj,1,ID);
        end
        function [Errcode]=removeBinRulesControlNodeID(obj,ID)
            [Errcode]=rmRulesControl(obj,1,ID);
        end
        function [Errcode]=removeBinControlNodeID(obj,ID)
            [Errcode]=rmControl(obj,0,ID);
        end
        function [Errcode]=setBinFlowUnitsGPM(obj)
            [Errcode]=Options(obj,'GPM'); %gallons per minute
        end
        function [Errcode]=setBinFlowUnitsLPS(obj)
            [Errcode]=Options(obj,'LPS'); %liters per second
        end
        function [Errcode]=setBinFlowUnitsMGD(obj)
            [Errcode]=Options(obj,'MGD'); %million gallons per day
        end
        function [Errcode]=setBinFlowUnitsIMGD(obj)
            [Errcode]=Options(obj,'IMGD'); %Imperial mgd
        end
        function [Errcode]=setBinFlowUnitsCFS(obj)
            [Errcode]=Options(obj,'CFS'); %cubic feet per second
        end
        function [Errcode]=setBinFlowUnitsAFD(obj)
            [Errcode]=Options(obj,'AFD'); %acre-feet per day
        end
        function [Errcode]=setBinFlowUnitsLPM(obj)
            [Errcode]=Options(obj,'LPM'); %liters per minute
        end
        function [Errcode]=setBinFlowUnitsMLD(obj)
            [Errcode]=Options(obj,'MLD'); %million liters per day
        end
        function [Errcode]=setBinFlowUnitsCMH(obj)
            [Errcode]=Options(obj,'CMH'); %cubic meters per hour
        end
        function [Errcode]=setBinFlowUnitsCMD(obj)
            [Errcode]=Options(obj,'CMD'); %cubic meters per day
        end
        function [Errcode]=setBinHeadlossHW(obj)
            [Errcode]=Options(obj,'','H-W');  %Hazen-Wiliams
        end
        function [Errcode]=setBinHeadlossDW(obj)
            [Errcode]=Options(obj,'','D-W');  %Darcy-Weisbach
        end
        function [Errcode]=setBinHeadlossCM(obj)
            [Errcode]=Options(obj,'','C-M');  %Chezy-Manning
        end        
        function [Errcode]=setBinLinkPipeLengths(obj,varargin)
            parameter=varargin{1};
            indexParameter=4;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections);
        end
        function [Errcode]=setBinLinkPipeDiameters(obj,varargin)
            parameter=varargin{1};
            indexParameter=5;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinLinkPipeRoughness(obj,varargin)
            parameter=varargin{1};
            indexParameter=6;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinLinkPipeMinorLoss(obj,varargin)
            parameter=varargin{1};
            indexParameter=7;
            sections={'[PIPES]', '[PUMPS]'};
            [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinLinkPipeStatus(obj,varargin)
           indexParameter=8;
           if sum(strcmp(lower(varargin{1}),'closed')+strcmp(lower(varargin{1}),'open')+strcmp(lower(varargin{1}),'cv'))==obj.BinLinkPipeCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
           end
           sections={'[PIPES]', '[PUMPS]'};
           [Errcode]=setBinParam(obj,indexParameter,parameter,sections); 
        end
        function [Errcode]=setBinLinkPumpStatus(obj,varargin)
            if sum(strcmp(lower(varargin{1}),'closed')+strcmp(lower(varargin{1}),'open'))==obj.BinLinkPumpCount
                parameter=varargin{1};
            else
                warning('Invalid argument found.');Errcode=-1;
                return;
            end   
            zz=abs(obj.BinLinkPumpCount+1+obj.BinLinkValveCount-obj.BinCountStatuslines);
            sections={'[STATUS]','[PATTERNS]','pump'};
            [Errcode]=setBinParam2(obj,parameter,sections,zz);   
        end
        function [Errcode]=setBinLinkPipesParameters(obj,varargin)
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
                        if sum(strcmp(lower(varargin{2*i}),'closed')+strcmp(lower(varargin{2*i}),'open'))==obj.BinLinkPipeCount
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
            fid = fopen([obj.BinTempfile], 'w');
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
                Errcode=closeOpenNetwork(obj);
            end
        end
        function [Errcode]=setBinNodeJunctionsParameters(obj,varargin)
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
            fid = fopen([obj.BinTempfile], 'w');
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
                Errcode=closeOpenNetwork(obj);
            end
        end
        function [Errcode]=setBinNodeTankParameters(obj,varargin)
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
%                         v=obj.getBinNodesInfo;
                        mixm=upper(MixModel(find(cellfun('isempty',(varargin{2*i}))==0)));
                        for u=1:length(mixm)
                            if ~sum(strcmp(mixm(u),{'MIX1', 'FIFO','LIFO','MIXED','2COMP'}))
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
            fid = fopen([obj.BinTempfile], 'w');
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
               Errcode=closeOpenNetwork(obj);
           end
        end
        function [Errcode]=setBinNodeReservoirParameters(obj,varargin)
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
            fid = fopen([obj.BinTempfile], 'w');
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
                Errcode=closeOpenNetwork(obj);
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
            [fid,binfile] = runEPANETexe(obj);
            if fid~=-1
                data = fread(fid,'int32');
                fclose(fid);
                value.BinNumberReportingPeriods = data(end-2);
                clear data;
                fid1 = fopen(binfile, 'r');

                % Seek to the 10th byte ('J'), read 5
                fseek(fid1, 0, 'bof');
                value.Binmagicnumber=fread(fid1, 1, 'uint32');
                value.BinLibEPANET=fread(fid1, 1, 'uint32');
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
            if fid==-1
                fprintf('"Run was unsuccessful."\n');
            else
                fprintf('"Run was successful."\n');
            end
        end
        function [info,tline,allines] = readInpFile(obj,varargin)
            if ~sum(strcmp(who,'varargin'))
                [info,tline,allines] = readAllFile(obj.BinTempfile);
            else
                [info,tline,allines] = readAllFile(obj.BinTempfile);
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
                        value.BinCountInitialQualitylines=0;
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
                    value.BinCountInitialQualitylines=d;
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
                    if strcmpi(tok(1:5),'[OPTI')
                        sect=3;
                    end
                end
                if sect==0
                    continue;
                elseif sect==3
                    if strcmp(regexp(tline,'\w*QUALITY*\w','match'),'QUALITY')
                        tl=regexp(tline,'\s*','split');mm=2;
                        if isempty(tl{1})
                            mm=3;
                        end
                        valueQualityType=tl(mm);
                    end
                end
            end
        end
        function [valueCoord,valueRule] = getBinCoordRuleSections(obj)
            % Open epanet input file
            [info]=regexp( fileread(obj.InputFile), '\n', 'split'); 
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
                    if strcmpi(tok(1:5),'[COOR')
                        sect=1;
                        valueCoord{d}=tline;
                        continue;
                        % [RULES] section
                    elseif strcmpi(tok(1:5),'[RULE')
                        sect=2;
                        valueRule{dRule}=tline;
                        continue;
                    elseif strcmpi(tok(1:5),'[ENER')
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
        function value = getNodeCoordinates(obj, varargin)
            if isempty(varargin)
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
                        [~,index] = ENgetlinkindex(char(A{1}),obj.LibEPANET);
                        nvert(index) = nvert(index) + 1;
                        vertx{index}(nvert(index)) = A{2};
                        verty{index}(nvert(index)) = A{3};
                    end
                end
            end
            try
                indices = getNodeIndices(obj,varargin);j=1;
                for i=indices
                    [obj.Errcode,vx(j),vy(j)]=ENgetcoord(i,obj.LibEPANET);
                    if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                    j=j+1;
                end
            catch e
            end
            if isempty(varargin)
                value{1} = vx;
                value{2} = vy;
                value{3} = vertx;
                value{4} = verty;
            else
                if obj.Errcode, error(obj.getError(obj.Errcode)), return; end   
                value = [vx' vy'];
            end
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
                        value.BinCountStatuslines=0;
                        continue;                        
                    elseif strcmpi(tok(1:5),'[REAC')
                        sect=5;d=1;
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
                    value.BinCountStatuslines=d;
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
function [Errcode] = ENwriteline (line,LibEPANET)
    [Errcode]=calllib(LibEPANET,'ENwriteline',line);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode] = ENaddpattern(patid,LibEPANET)
    Errcode=calllib(LibEPANET,'ENaddpattern',patid);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode] = ENclose(LibEPANET)
    [Errcode]=calllib(LibEPANET,'ENclose');
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode] = ENcloseH(LibEPANET)
    [Errcode]=calllib(LibEPANET,'ENcloseH');
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode, value] = ENgetnodevalue(index, paramcode,LibEPANET)
    value=single(0);
    %p=libpointer('singlePtr',value);
    index=int32(index);
    paramcode=int32(paramcode);
    [Errcode, value]=calllib(LibEPANET,'ENgetnodevalue',index, paramcode,value);
    if Errcode==240
        value=NaN;
    end
end
function [Errcode, value] = ENgetbasedemand(index,numdemands,LibEPANET)
    %epanet20100
    [Errcode,value]=calllib(LibEPANET,'ENgetbasedemand',index,numdemands,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode, value] = ENgetnumdemands(index,LibEPANET)
    %epanet20100
    [Errcode,value]=calllib(LibEPANET,'ENgetnumdemands',index,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode, value] = ENgetdemandpattern(index,numdemands,LibEPANET)
    %epanet20100
    [Errcode,value]=calllib(LibEPANET,'ENgetdemandpattern',index,numdemands,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode, value] = ENgetstatistic(code,LibEPANET)
    %epanet20100
    [Errcode,value]=calllib(LibEPANET,'ENgetstatistic',code,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode] = ENcloseQ(LibEPANET)
    [Errcode]=calllib(LibEPANET,'ENcloseQ');
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end

function [Errcode, ctype,lindex,setting,nindex,level] = ENgetcontrol(cindex,LibEPANET)
    [Errcode, ctype,lindex,setting,nindex,level]=calllib(LibEPANET,'ENgetcontrol',cindex,0,0,0,0,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [Errcode, count] = ENgetcount(countcode,LibEPANET)
    [Errcode,count]=calllib(LibEPANET,'ENgetcount',countcode,0);
    if Errcode
        ENgeterror(Errcode,LibEPANET);
    end
end
function [errmsg, e] = ENgeterror(Errcode,LibEPANET)
    if Errcode
        errmsg = char(32*ones(1,79));
        [e,errmsg] = calllib(LibEPANET,'ENgeterror',Errcode,errmsg,79);
    else
        e=0;
        errmsg='';
    end
end
function [Errcode,flowunitsindex] = ENgetflowunits(LibEPANET)
[Errcode, flowunitsindex]=calllib(LibEPANET,'ENgetflowunits',0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,id] = ENgetlinkid(index,LibEPANET)
id=char(32*ones(1,31));
[Errcode,id]=calllib(LibEPANET,'ENgetlinkid',index,id);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,index] = ENgetlinkindex(id,LibEPANET)
[Errcode,~,index]=calllib(LibEPANET,'ENgetlinkindex',id,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,from,to] = ENgetlinknodes(index,LibEPANET)
[Errcode,from,to]=calllib(LibEPANET,'ENgetlinknodes',index,0,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, type] = ENgetlinktype(index,LibEPANET)
[Errcode,type]=calllib(LibEPANET,'ENgetlinktype',index,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetlinkvalue(index, paramcode,LibEPANET)
[Errcode,value]=calllib(LibEPANET,'ENgetlinkvalue',index, paramcode, 0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,id] = ENgetnodeid(index,LibEPANET)
id=char(32*ones(1,31));
[Errcode,id]=calllib(LibEPANET,'ENgetnodeid',index,id);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,index] = ENgetnodeindex(id,LibEPANET)
[Errcode, ~, index]=calllib(LibEPANET,'ENgetnodeindex',id,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, type] = ENgetnodetype(index,LibEPANET)
[Errcode,type]=calllib(LibEPANET,'ENgetnodetype',index,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetoption(optioncode,LibEPANET)
[Errcode,value]=calllib(LibEPANET,'ENgetoption',optioncode,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, id] = ENgetpatternid(index,LibEPANET)
id=char(32*ones(1,31));
[Errcode,id]=calllib(LibEPANET,'ENgetpatternid',index,id);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, id] = ENgetcurveid(index,LibEPANET)
%New version dev2.1
id=char(32*ones(1,31));
[Errcode,id]=calllib(LibEPANET,'ENgetcurveid',index,id);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, index] = ENgetpatternindex(id,LibEPANET)
[Errcode,~, index]=calllib(LibEPANET,'ENgetpatternindex',id,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, len] = ENgetpatternlen(index,LibEPANET)
[Errcode,len]=calllib(LibEPANET,'ENgetpatternlen',index,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetpatternvalue(index, period,LibEPANET)
[Errcode,value]=calllib(LibEPANET,'ENgetpatternvalue',index, period, 0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,qualcode,tracenode] = ENgetqualtype(LibEPANET)
[Errcode,qualcode,tracenode]=calllib(LibEPANET,'ENgetqualtype',0,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, timevalue] = ENgettimeparam(paramcode,LibEPANET)
[Errcode,timevalue]=calllib(LibEPANET,'ENgettimeparam',paramcode,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, LibEPANET] = ENgetversion(LibEPANET)
[Errcode,LibEPANET]=calllib(LibEPANET,'ENgetversion',0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENinitH(flag,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENinitH',flag);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENinitQ(saveflag,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENinitQ',saveflag);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
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
function ENLoadLibrary(LibEPANETpath,LibEPANET)
if ~libisloaded(LibEPANET)
    loadlibrary([LibEPANETpath,LibEPANET],[LibEPANETpath,LibEPANET,'.h'])
end
if libisloaded(LibEPANET)
    LibEPANETString = 'EPANET loaded sucessfuly.';
    disp(LibEPANETString);
else
    warning('There was an error loading the EPANET library (DLL).')
end
end
function [Errcode, tstep] = ENnextH(LibEPANET)
[Errcode,tstep]=calllib(LibEPANET,'ENnextH',int32(0));
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, tstep] = ENnextQ(LibEPANET)
[Errcode,tstep]=calllib(LibEPANET,'ENnextQ',int32(0));
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
tstep = double(tstep);
end
function [Errcode] = ENopen(inpname,repname,binname,LibEPANET) %DE
    Errcode=calllib(LibEPANET,'ENopen',inpname,repname,binname);
    if Errcode
       [~,errmsg] = calllib(LibEPANET,'ENgeterror',Errcode,char(32*ones(1,79)),79);
       warning(errmsg);
    end
end
function [Errcode] = ENopenH(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENopenH');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENopenQ(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENopenQ');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENreport(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENreport');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENresetreport(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENresetreport');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, t] = ENrunH(LibEPANET)
[Errcode,t]=calllib(LibEPANET,'ENrunH',int32(0));
t = double(t);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, t] = ENrunQ(LibEPANET)
t=int32(0);
[Errcode,t]=calllib(LibEPANET,'ENrunQ',t);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsaveH(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsaveH');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsavehydfile(fname,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsavehydfile',fname);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsaveinpfile(inpname,LibEPANET)
Errcode=calllib(LibEPANET,'ENsaveinpfile',inpname);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetcontrol(cindex,ctype,lindex,setting,nindex,level,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetcontrol',cindex,ctype,lindex,setting,nindex,level);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetlinkvalue(index, paramcode, value,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetlinkvalue',index, paramcode, value);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetnodevalue(index, paramcode, value,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetnodevalue',index, paramcode, value);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetoption(optioncode,value,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetoption',optioncode,value);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetpattern(index, factors, nfactors,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetpattern',index,factors,nfactors);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetpatternvalue(index, period, value,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetpatternvalue',index, period, value);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetqualtype',qualcode,chemname,chemunits,tracenode);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetreport(command,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetreport',command);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetstatusreport(statuslevel,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsetstatusreport',statuslevel);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsettimeparam(paramcode, timevalue,LibEPANET)
paramcode=int32(paramcode);
timevalue=int32(timevalue);
[Errcode]=calllib(LibEPANET,'ENsettimeparam',paramcode,timevalue);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsolveH(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsolveH');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsolveQ(LibEPANET)
[Errcode]=calllib(LibEPANET,'ENsolveQ');
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, tleft] = ENstepQ(LibEPANET)
tleft=int32(0);
[Errcode,tleft]=calllib(LibEPANET,'ENstepQ',tleft);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
tleft=double(tleft);
end
function [Errcode] = ENusehydfile(hydfname,LibEPANET)
[Errcode]=calllib(LibEPANET,'ENusehydfile',hydfname);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetcurve(index, x, y, nfactors,LibEPANET)
% New version dev2.1
[Errcode]=calllib(LibEPANET,'ENsetcurve',index,x,y,nfactors);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, x, y] = ENgetcurvevalue(index, period,LibEPANET)
% New version dev2.1
[Errcode,x, y]=calllib(LibEPANET,'ENgetcurvevalue',index, period, 0, 0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, x, y] = ENsetcurvevalue(index, pnt, x, y, LibEPANET)
% New version dev2.1
% index  = curve index
% pnt    = curve's point number
% x      = curve x value
% y      = curve y value                            
% sets x,y point for a specific point and curve 
[Errcode]=calllib(LibEPANET,'ENsetcurvevalue',index, pnt, x, y);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, index] = ENgetcurveindex(id,LibEPANET)
% New version dev2.1
[Errcode,~, index]=calllib(LibEPANET,'ENgetcurveindex',id,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENaddcurve(cid,LibEPANET)
% New version dev2.1
Errcode=calllib(LibEPANET,'ENaddcurve',cid);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, ids, nvalue, xvalue, yvalue] = ENgetcurve(value,LibEPANET)
[~,~,nvalue,~,~]=calllib(LibEPANET,'ENgetcurve',value,char(32*ones(1,31)),0,0,0);
[Errcode,ids,~, xvalue, yvalue]=calllib(LibEPANET,'ENgetcurve',value,char(32*ones(1,31)),0,zeros(1,nvalue),zeros(1,nvalue));
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, len] = ENgetcurvelen(index,LibEPANET)
% New version dev2.1
[Errcode,len]=calllib(LibEPANET,'ENgetcurvelen',index,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetheadcurveindex(pumpindex,LibEPANET)
% New version dev2.1
[Errcode,value]=calllib(LibEPANET,'ENgetheadcurveindex',pumpindex,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetpumptype(pumpindex,LibEPANET)
% New version dev2.1
[Errcode,value]=calllib(LibEPANET,'ENgetpumptype',pumpindex,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, value] = ENgetaveragepatternvalue(index,LibEPANET) 
% return  average pattern value
% New version dev2.1
[Errcode,value]=calllib(LibEPANET,'ENgetaveragepatternvalue',index,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode, x, y] = ENgetcoord(index,LibEPANET)
% New version dev2.1
[Errcode,x,y]=calllib(LibEPANET,'ENgetcoord',index,0,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetcoord(index,x,y,LibEPANET)
% New version dev2.1
[Errcode]=calllib(LibEPANET,'ENsetcoord',index,x,y);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode] = ENsetbasedemand(index, demandIdx, value, LibEPANET)
% New version dev2.1
[Errcode]=calllib(LibEPANET,'ENsetbasedemand',index, demandIdx, value);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [Errcode,qualcode,chemname,chemunits,tracenode] = ENgetqualinfo(LibEPANET)
chm=char(32*ones(1,31));
[Errcode,qualcode,chemname,chemunits,tracenode]=calllib(LibEPANET,'ENgetqualinfo',0,chm,chm,0);
if Errcode
    ENgeterror(Errcode,LibEPANET);
end
end
function [obj] = MSXMatlabSetup(obj,msxname,varargin)
if ~isempty(varargin)
    if varargin{1}{1}~=1
        if nargin==3 
            obj.MSXLibEPANETPath=char(varargin{1});
            obj.MSXLibEPANETPath=[fileparts(obj.MSXLibEPANETPath),'\'];
            if isempty(varargin{1})
                obj.MSXLibEPANETPath='';
            end
        end
    end
else
    obj.MSXLibEPANETPath=obj.LibEPANETpath;
end
obj.MSXLibEPANET='epanetmsx'; % Get DLL LibEPANET (e.g. epanet20012x86 for 32-bit)
if ~libisloaded(obj.MSXLibEPANET)
    loadlibrary([obj.MSXLibEPANETPath,obj.MSXLibEPANET],[obj.MSXLibEPANETPath,[obj.MSXLibEPANET,'.h']]);
end

obj.MSXFile = char(msxname);
%Save the temporary msx file
if ~isempty(varargin)
    if varargin{1}{1}==1
        if ~iscell(varargin)
            obj.MSXTempFile=obj.MSXFile;
        end
    end
else
    obj.MSXTempFile=[obj.MSXFile(1:end-4),'_temp.msx'];
    copyfile(obj.MSXFile,obj.MSXTempFile);
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
function [Errcode] = MSXopen(obj,varargin)
[Errcode] = calllib(obj.MSXLibEPANET,'MSXopen',obj.MSXTempFile);
if Errcode
    MSXerror(Errcode,obj.MSXLibEPANET);
end
if (Errcode == 520)
    disp('current MSX project will be closed and the new project will be opened');
    [Errcode] = MSXclose(obj);
    if Errcode
        MSXerror(Errcode,obj.MSXLibEPANET);
    else
        [Errcode] = calllib(obj.MSXLibEPANET,'MSXopen',obj.MSXTempFile);
        if Errcode
            MSXerror(Errcode,obj.MSXLibEPANET);
        end
    end
end
end
function [Errcode] = MSXclose(obj)
[Errcode] = calllib(obj.MSXLibEPANET,'MSXclose');
if Errcode
    MSXerror(Errcode,obj.MSXLibEPANET);
end
end
function [e] = MSXerror(Errcode,MSXLibEPANET)
len=80;
errstring=char(32*ones(1,len+1));
[e,errstring] = calllib(MSXLibEPANET,'MSXgeterror',Errcode,errstring,len);
disp(errstring);
end
function [Errcode, count] = MSXgetcount(code,MSXLibEPANET)
count=0;
[Errcode,count] = calllib(MSXLibEPANET,'MSXgetcount',code,count);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, index] = MSXgetindex(MSXLibEPANET,varargin)
index =0;
if ~isnumeric(varargin{1})
    varargin{1}=varargin{2};
    varargin{2}=varargin{3};
end
[Errcode,~,index]=calllib(MSXLibEPANET,'MSXgetindex',varargin{1},varargin{2},index);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, id] = MSXgetID(type,index,len,MSXLibEPANET)
id=char(32*ones(1,len+1));
[Errcode,id]=calllib(MSXLibEPANET,'MSXgetID',type,index,id,len);
id=id(1:len);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, len] = MSXgetIDlen(type,index,MSXLibEPANET)
len=0;
[Errcode,len]=calllib(MSXLibEPANET,'MSXgetIDlen',type,index,len);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, type, units, atol, rtol] = MSXgetspecies(index,MSXLibEPANET)
type=0; rtol=0; atol=0;
units=char(32*ones(1,16));
[Errcode,type,units,atol,rtol]=calllib(MSXLibEPANET,'MSXgetspecies',index,type,units,atol,rtol);
switch type
    case 0
        type='BULK';   % for a bulk water species
    case 1
        type='WALL';   % for a pipe wall surface species
end
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetconstant(index,MSXLibEPANET)
value=0;
[Errcode,value]=calllib(MSXLibEPANET,'MSXgetconstant',index,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetparameter(type,index,param,MSXLibEPANET)
value=0;
[Errcode,value]=calllib(MSXLibEPANET,'MSXgetparameter',type,index,param,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, patlen] = MSXgetpatternlen(patindex,MSXLibEPANET)
patlen=0;
[Errcode,patlen]=calllib(MSXLibEPANET,'MSXgetpatternlen',patindex,patlen);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetpatternvalue(patindex,period,MSXLibEPANET)
value=0;
[Errcode,value]=calllib(MSXLibEPANET,'MSXgetpatternvalue',patindex,period,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, value] = MSXgetinitqual(obj,index,species,MSXLibEPANET)
value=0;
[Errcode,value]=calllib(MSXLibEPANET,'MSXgetinitqual',obj,index,species,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, type, level, pat] = MSXgetsource(node,species,MSXLibEPANET)
type=0;
level=0;
pat=0;
[Errcode,type,level,pat]=calllib(MSXLibEPANET,'MSXgetsource',node,species,type,level,pat);
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
    MSXerror(Errcode,MSXLibEPANET);
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
function [Errcode] = MSXsaveoutfile(outfname,MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET,'MSXsaveoutfile',outfname);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsavemsxfile(msxname,MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET,'MSXsavemsxfile',msxname);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetconstant(index, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetconstant',index,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetparameter(type, index, param, value, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetparameter',type,index,param,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetinitqual(type,index,species,value,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetinitqual',type,index,species,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetpattern(index, factors, nfactors, MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetpattern',index,factors,nfactors);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetpatternvalue(pat, period, value,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetpatternvalue',pat,period,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsolveQ(MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsolveQ');
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsolveH(MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsolveH');
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXaddpattern(patid,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXaddpattern',patid);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXusehydfile(hydfname,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXusehydfile',hydfname);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode, t, tleft] = MSXstep(MSXLibEPANET)
t=int32(0);
tleft=int32(0);
[Errcode,t,tleft]=calllib(MSXLibEPANET,'MSXstep',t,tleft);
t = double(t);
tleft = double(tleft);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXinit(flag,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXinit',flag);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXreport(MSXLibEPANET)
[Errcode] = calllib(MSXLibEPANET,'MSXreport');
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [e, errmsg] = MSXgeterror(Errcode,MSXLibEPANET)
errmsg = char(32*ones(1,80));
[e,errmsg] = calllib(MSXLibEPANET,'MSXgeterror',Errcode,errmsg,80);
if e
    MSXerror(e);
end
end
function [Errcode, value] = MSXgetqual(type, index, species,MSXLibEPANET)
value=0;
[Errcode,value]=calllib(MSXLibEPANET,'MSXgetqual',type,index,species,value);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
end
end
function [Errcode] = MSXsetsource(node,species,type,level,pat,MSXLibEPANET)
[Errcode]=calllib(MSXLibEPANET,'MSXsetsource',node,species,type,level,pat);
if Errcode
    MSXerror(Errcode,MSXLibEPANET);
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
fclose(fid);
end
function value=getMSXOptions(msxname)

if isempty(msxname)
    warning('Please load MSX File.');
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
fclose(fid);
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
    chckfunctions=libfunctions(obj.LibEPANET);
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
for i=find(~cellfun(@isempty,legendString))
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
function [Errcode]=setBinParam(obj,indexParameter,parameter,sections,varargin)
    ok=0;Errcode=0;
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
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
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
    fid = fopen([obj.BinTempfile], 'w');
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
        Errcode=closeOpenNetwork(obj);
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
function [Errcode]=setBinParam2(obj,parameter,sections,zz,varargin)
    Errcode=0;
    if strcmp(sections{1},'[STATUS]')
        value =obj.getBinLinksInfo;
        if strcmp(sections{3},'pump')
            nameID=value.BinLinkPumpStatusNameID;
            cntlv=obj.BinLinkPumpCount;
        elseif strcmp(sections{3},'valve')
            nameID=value.BinLinkValveStatusNameID;
            cntlv=obj.BinLinkValveCount;
            if strcmp(upper(parameter),'NONE'), Errcode=-1;return;end
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
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
    fid = fopen([obj.BinTempfile], 'w');
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
        Errcode=closeOpenNetwork(obj);
    end
end
function value = getBinParam(obj,sections,varargin)  
    warning off;
    [tlines]=regexp( fileread([obj.BinTempfile]), '\n', 'split');
    if strcmp(sections{1},'[SOURCES]')
        value.BinNodeSourcePatternIndex = nan(1,obj.BinNodeCount);
        value.BinNodeSourceQuality = nan(1,obj.BinNodeCount);
        value.BinNodeSourceTypeCode = nan(1,obj.BinNodeCount);
        value.BinNodeSourceType = cell(1,obj.BinNodeCount);
        value.BinNodeSourcePatternNameID = cell(1,obj.BinNodeCount);
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
                    value.BinCountPatternlines=d;
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
function [Errcode]=addCurve(obj,newCurveID,varargin)
v=obj.getBinCurvesInfo;Errcode=0;
CurveX=varargin{1};
CurveY=varargin{2};
typecode=varargin{3};
% PUMP 0 EFFICIENCY 1 VOLUME 2 HEADLOSS 3
for i=1:length(CurveX)
    if i+1<length(CurveX)+1
        if CurveX(i)>=CurveX(i+1)
            if strfind([0 1 3],typecode)
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
i=1;exists=0;
while i<length(v.BinCurveNameID)+1
    exists(i) = strcmp(newCurveID,v.BinCurveNameID{i});
    i=i+1;
end
if sum(exists)>0
    s = sprintf('Curve "%s" already exists.',newCurveID);
    warning(s);Errcode=-1;
    return
end
sect=0;
% Open and read inpname
% Read all file and save in variable info
[~,info] = obj.readInpFile;
% write
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
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
function Errcode=addNode(obj,typecode,varargin)
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
        for i=1:nodes.BinNodeCount
            t = strcmpi(nodes.BinNodeNameID(i),l);
            if sum(t)==0
                s = sprintf('Node %s disconnected.',char(nodes.BinNodeNameID(i)));
                warning(s);Errcode=-1;
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
                Errcode=-1;
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
        warning('There is no such object in the network.');Errcode=-1;
        return;
    end
    for i=1:length(nodes.BinNodeNameID)
        exists(i) = strcmp(newID,char(nodes.BinNodeNameID(i)));
    end
    if (sum(exists)==1)
        s = sprintf('Node "%s" already exists.',newID);
        warning(s);Errcode=-1;
        return;
    end
    % Get type of node
    A = [0 1 2];
    code=strfind(A,typecode);
    if length(code)==0
        warning('There is no such typecode(0-2)!');Errcode=-1;
        return;
    end
    % check section in inpname, [JUNCTIONS], [RESERVOIRS], [TANKS]
    stank_check=1;
    sreservoir_check=1;
    sjunction_check=1;
    % Open and read inpname
    % Read all file and save in variable info
    [~,info,~] = obj.readInpFile;
    fid2 = fopen([obj.BinTempfile],'w');
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
        Errcode=closeOpenNetwork(obj);
    end
end
function Errcode=addLink(obj,typecode,newLink,fromNode,toNode,varargin)
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
        curveID='';
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
[Errcode]=addLinkWarnings(obj,typecode,newLink,toNode);
crvs = obj.getBinCurvesInfo;
% Open and read inpname
% Read all file and save in variable info
[~,info] = obj.readInpFile;
fid2 = fopen([obj.BinTempfile],'w');
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
                    fprintf(fid2,'\n%s%s%s%s%s%s%s%s%s',newLink,sps,fromNode,sps,...
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
    Errcode=closeOpenNetwork(obj);
end
end
function [Errcode] = rmNode(obj,NodeID)
% Remove node from the network.
% Check if id new already exists
nodes = obj.getBinNodesInfo;Errcode=0;
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
    Errcode=-1;
    return;
end
if sum(strcmp(NodeID,nodes.BinNodeReservoirNameID)) || sum(strcmp(NodeID,nodes.BinNodeTankNameID))
    if (length(char(nodes.BinNodeReservoirNameID))+length(char(nodes.BinNodeTankNameID))-1)==0;
        warning('This tank/reservoir has not removed.');
        Errcode=-1;
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
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function Errcode=rmRulesControl(obj,type,id)
% Remove control from the network.
exists=0;Errcode=0;
exists1=0;
rulescontrols = obj.getBinRulesControlsInfo;
if type==1
    if length(rulescontrols.BinRulesControlLinksID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);Errcode=-1;
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
        warning(s);Errcode=-1;
        return;
    end
    for i=1:length(rulescontrols.BinRulesControlNodesID)
        if type==0
            exists1(i) = strcmp( rulescontrols.BinRulesControlNodesID{i}{length(rulescontrols.BinRulesControlNodesID{1})},char(id));
        else
            warning('Type is NODE(0) or LINK(1)');Errcode=-1;
            return;
        end
    end
end
[~,info] = obj.readInpFile;
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function Errcode=rmControl(obj,type,id)
% Remove control from the network.
exists=0;Errcode=0;
exists1=0;
controls = obj.getBinControlsInfo;
if type==1
    if length(controls.BinControlLinksID)==0
        s = sprintf('There is no such object in the network.');
        warning(s);Errcode=-1;
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
        warning(s);Errcode=-1;
        return;
    end
    for i=1:length(controls.BinControlNodesID)
        if type==0
            exists1(i) = strcmp(controls.BinControlNodesID{i},char(id));
        else
            warning('Type is NODE(0) or LINK(1)');Errcode=-1;
            return;
        end
    end
end
[~,info] = obj.readInpFile;
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function [Errcode] = rmLink(obj,LinkID)
% Remove link from the network.
% Check if id new already exists
links = obj.getBinLinksInfo;Errcode=0;
if length(links.BinLinkNameID)==0
    s = sprintf('There is no such object in the network.');
    warning(s);Errcode=-1;
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
    warning(s);Errcode=-1;
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
fid2 = fopen([obj.BinTempfile],'w');

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
            Errcode=0;
        else
            Errcode=-1;
        end
    end
end
if obj.Bin==1
    Errcode=closeOpenNetwork(obj);
end
end
function [Errcode]=addNewControl(obj,x,status,y_t_c,param,z,varargin)
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
        warning(s);Errcode=-1;
        return;
    end
end
% Check if id new already exists
links = obj.getBinLinksInfo; Errcode=0;
if length(char(links.BinLinkNameID))==0
    s = sprintf('There is no such object in the network.');
    warning(s);Errcode=-1;
    return;
end
i=1;
while i<length(char(links.BinLinkNameID))+1
    exists(i) = strcmp(x,char(links.BinLinkNameID(i)));
    i=i+1;
end
if (sum(exists)~=1)
    s = sprintf('There is no such object in the network.');
    warning(s);Errcode=-1;
    return;
end
type_n='[CONTROLS]';
[~,info] = obj.readInpFile;
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function [Errcode]=rmCurveID(obj,CurveID,varargin)
% Check if id new already exists
valueC=obj.getBinCurvesInfo;Errcode=0;
if length(valueC.BinCurveNameID)==0
    s = sprintf('There is no such object in the network.');
    warning(s);Errcode=-1;
    return;
end
i=1;
while i<length(valueC.BinCurveNameID)+1
    exists(i) = strcmp(CurveID,char(valueC.BinCurveNameID(i)));
    i=i+1;
end
if (sum(exists)==0)
    s = sprintf('There is no such object in the network.');
    warning(s);Errcode=-1;
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
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function [Errcode]=Options(obj,newFlowUnits,headloss,varargin)
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
newUnits_US_Customary=0;
newUnits_SI_Metric=0;
switch newFlowUnits
    case 'CFS'
        newUnits_US_Customary=1;
    case 'GPM'
        newUnits_US_Customary=1;
    case 'MGD'
        newUnits_US_Customary=1;
    case 'IMGD'
        newUnits_US_Customary=1;
    case 'AFD'
        newUnits_US_Customary=1;
    case 'LPS'
        newUnits_SI_Metric=1;
    case 'LPM'
        newUnits_SI_Metric=1;
    case 'MLD'
        newUnits_SI_Metric=1;
    case 'CMH'
        newUnits_SI_Metric=1;
    case 'CMD'
        newUnits_SI_Metric=1;
end
if newUnits_US_Customary==value.BinUnits_US_Customary
    changes=0; newUnits_US_Customary=0;
    newUnits_SI_Metric=0; % feet to feet
elseif newUnits_SI_Metric==value.BinUnits_SI_Metric
    changes=1; newUnits_US_Customary=0;
    newUnits_SI_Metric=0; % meter to meter
elseif value.BinUnits_US_Customary==1 && newUnits_US_Customary==0
    changes=1; % feet to meter or cubic feet to cubic meter
elseif value.BinUnits_US_Customary==0 && newUnits_US_Customary==1
    changes=2; % meter to feet or cubic meter to cubic feet
end
Units=newUnits_US_Customary+newUnits_SI_Metric;
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

[info] = readAllFile(obj.BinTempfile);
fid2 = fopen([obj.BinTempfile],'w');
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
    Errcode=closeOpenNetwork(obj);
end
end
function Errcode=closeOpenNetwork(obj)
    obj.closeNetwork;  %Close input file 
    Errcode=ENopen([obj.BinTempfile],[obj.BinTempfile(1:end-4),'.txt'],[obj.BinTempfile(1:end-4),'.bin'],obj.LibEPANET);
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
function [fid,binfile] = runEPANETexe(obj)
    [~,mm]=system(['cmd /c for %A in ("',pwd,'") do @echo %~sA']);
    mmPwd=regexp(mm,'\s','split');
    pp=[mmPwd{1},'/'];
    inpfile=[pp,obj.BinTempfile];
    rptfile=[inpfile(1:length(inpfile)-4),'.txt'];
    binfile=[inpfile(1:length(inpfile)-4),'.bin'];
    if exist(binfile)==2, fclose all; delete(binfile); end
    if strcmp(computer('arch'),'win64')
            folder='64bit';
        r = sprintf('%s\\%s\\epanet2d.exe %s %s %s',mmPwd{1},folder,inpfile,rptfile,binfile);
    elseif strcmp(computer('arch'),'win32')
            folder='32bit';
        r = sprintf('%s\\%s\\epanet2d.exe %s %s %s',mmPwd{1},folder,inpfile,rptfile,binfile);
    end
    system(r);
    fid = fopen(binfile,'r');
end
function value = getBinComputedTimeSeries(obj,indParam,varargin)
    [fid,binfile] = runEPANETexe(obj);
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
        clear data;
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
    if fid==-1
        fprintf('"Run was unsuccessful."\n');
    else
        fprintf('"Run was successful."\n');
    end
end
function Errcode=addLinkWarnings(obj,typecode,newLink,toNode)
% Check if id new already exists
Nodes = obj.getBinNodesInfo;
Errcode=0;
if length(Nodes.BinNodeNameID)==0
    Errcode=-1;
    return;
end
existsTo=0;
for i=1:length(Nodes.BinNodeNameID) 
    existsTo(i) = strcmp(toNode,char(Nodes.BinNodeNameID{i}));
end
if sum(existsTo)~=1
    s = sprintf('There is no node "%s" in the network.',toNode);
    warning(s);Errcode=-1;
    return;
end
A = [0 1 2 3 4 5 6 7 8];
code = strfind(A,typecode);
if length(code)==0
    warning('There is no such typecode(0-8)');
    Errcode=-1;
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
        Errcode=-1;
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
        Errcode=-1;
        warning(s);
        return;
    end
end

if typecode==2
    crvs = obj.getBinCurvesInfo;
    if isempty(char(crvs.BinCurveNameID))
        s = sprintf('No head curve supplied for pump %s.',newLink);
        Errcode=-1;
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
function setMSXOptions(obj,varargin)
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
                
[info,tline] = readAllFile(obj.MSXFile);
fid2 = fopen([obj.MSXTempFile],'w');
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
obj.unloadMSX;
obj.loadMSXFile(obj.MSXTempFile,1);
end
function indices = getLinkIndices(obj, varargin)
    indices = getIndices(obj.getLinkCount, varargin{1});
end
function indices = getNodeIndices(obj, varargin)
    indices = getIndices(obj.getNodeCount, varargin{1});
end
function indices = getIndices(cnt, varargin)
    if isempty(varargin{1})
        indices=1:cnt;
    else
        indices=varargin{1}{1};
    end 
end