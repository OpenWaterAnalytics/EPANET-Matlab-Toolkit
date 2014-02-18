classdef epanet <handle
    %epanet EPANET-Matlab Class: A Matlab Class for EPANET and EPANET-MSX
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
    %   upadated version of the software
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
    %   Licensed under the EUPL, version 1.1 or - as soon they will be
    %   approved by the European Commission - subsequent versions of the
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
        %EPANET-MSX
        MsxConstantsNameID; 
        MsxConstantsValue;
        MsxConstantsCount;
        MsxConstantsIndex;
        MsxParametersCount;
        MsxPatternsCount;
        MsxSpeciesCount;
        MsxLinkInitqualValue;
        MsxNodeInitqualValue;
        MSXFile;
        MSXPathFile;
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
        MsxSpeciesATOL;
        MsxSpeciesIndex;
        MsxSpeciesNameID;
        MsxSpeciesRTOL;
        MsxSpeciesType;
        MsxSpeciesUnits;
        MsxEquationsTanks;
        MsxEquationsTerms;
%         MsxComputedQualityNode;
%         MsxComputedQualityLink;
        %EPANET
        NodeCoordinates; % Coordinates for each node (long/lat & intermediate pipe coordinates)
        NodeJunctionsCount; %Number of junctions
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
        LinkVelocityUnits; % Units for velocity
        ControlLevelValues; %The control level values
        ControlLinkIndex; % Set of control types in links
        ControlNodeIndex; % Set of control types in nodes
        ControlSettings; % Settings for the controls
        ControlTypes; % Set of control types
        ControlTypesIndex; %Index of the control types
        ControlRules; %Retrieves the parameters of all control statements
        ControlRulesCount;  % Number of controls
        TimeRuleControlStep; % Time step for evaluating rule-based controls
        CurveCount;    % Number of curves
        CurveNameID; %ID name of curves
        CurveTypes; % Type of curves
        CurveXvalue; %X-value of curves
        CurveYvalue; %Y-value of curves
        PatternCount;  % Number of patterns
        PatternDemandsUnits; %Units for demands
        PatternID; %ID of the patterns
        PatternIndex; %Indices of the patterns
        PatternLengths; %Length of the patterns
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
        errcode; %Code for the EPANET error message
        pathfile;   % The path of the input file
        pathfileMsx;
        inputfile;  % Name of the input file
        version; %EPANET library version
    end
    properties (Constant = true)
        TYPECONTROL={'LOWLEVEL','HIGHLEVEL', 'TIMER', 'TIMEOFDAY'}; % Constants for control: 'LOWLEVEL','HILEVEL', 'TIMER', 'TIMEOFDAY'
        TYPECURVE={'PUMP','EFFICIENCY','VOLUME','HEADLOSS'}; % Constants for pumps: 'PUMP','EFFICIENCY','VOLUME','HEADLOSS'
        TYPELINK={'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV', 'VALVE'}; % Constants for valves: 'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV', 'VALVE'
        TYPEMIXMODEL={'MIX1','MIX2', 'FIFO','LIFO'}; % Constants for mixing models: 'MIX1','MIX2', 'FIFO','LIFO'
        TYPENODE={'JUNCTION','RESERVOIR', 'TANK'}; % Contants for nodes: 'JUNCTION','RESERVOIR', 'TANK'
        TYPEPUMP={'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'}; % Constants for pumps: 'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'
        TYPEQUALITY={'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'}; % Constants for quality: 'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'
        TYPEREPORT={'YES','NO','FULL'}; % Constants for report: 'YES','NO','FULL'
        TYPESOURCE={'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'}; % Constants for sources: 'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'
        TYPESTATS={'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'}; % Constants for statistics: 'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
        TYPEUNITS={'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'}; % Constants for units: 'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'
    end
    methods
        function obj   = epanet(pathfile,varargin)
            %   The epanet object constructor. Example call:
            %   n=epanet('Net1.inp');
            
            path(path,genpath(pwd)); % Set current folder and subfolders in PATH
            if nargin==2 % Get name of INP file
                inpfile=varargin{1};
            elseif nargin==1
                inpfile=pathfile;
            end
            obj.inputfile=inpfile;
            %Load EPANET Library
            [obj.errcode]=obj.epanetLoadLibrary;
            %Open the file
            obj.LoadInpFile([pwd,'\NETWORKS\',inpfile], '', '');
            %Set path of temporary file
            obj.pathfile=[pwd,'\RESULTS\temp.inp'];
            %Save the temporary input file
            obj.saveInputFile(obj.pathfile);
            %Close input file
            ENclose;
            %Load temporary file
            obj.LoadInpFile(obj.pathfile,[pwd,'\RESULTS\temp.txt'], [pwd,'\RESULTS\temp.out']);
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
            obj.NodeJunctionsCount = obj.getNodeJunctionsCount;
            obj.LinkPipeCount = obj.getLinkPipeCount;
            obj.LinkPumpCount = obj.getLinkPumpCount;
            obj.LinkValveCount = obj.getLinkValveCount; 
            %Get all the controls
            obj.getControls;
            %Get the flow units
            obj.LinkFlowUnits = obj.getFlowUnits;
            %Get all the link data
            obj.LinkNameID = obj.getLinkNameID;
            obj.LinkIndex = obj.getLinkIndex;
            obj.LinkPipeIndex = find(strcmp(obj.LinkType,'PIPE'));
            obj.LinkPumpIndex = find(strcmp(obj.LinkType,'PUMP'));
            obj.LinkValveIndex = find(strcmp(obj.LinkType,'VALVE'));
            obj.LinkTypeIndex = obj.getLinkTypeIndex;
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
            %Get all the node data
            obj.NodesConnectingLinksIndex = obj.getLinkNodesIndex;
            obj.NodesConnectingLinksID = obj.getNodesConnectingLinksID;
            obj.NodeNameID = obj.getNodeNameID;
            obj.NodeIndex = obj.getNodeIndex;
            obj.NodeTypeIndex = obj.getNodeTypeIndex;
            obj.NodeElevations = obj.getNodeElevations;
            obj.NodeBaseDemands = obj.getNodeBaseDemands;
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
            obj.PatternID = obj.getPatternID;
            obj.PatternIndex = obj.getPatternIndex;
            obj.PatternLengths = obj.getPatternLengths;
            %Get quality types
            obj.QualityCode = obj.getQualityCode;
            obj.QualityTraceNodeIndex = obj.getQualityTraceNodeIndex;
            obj.QualityType = obj.getQualityType;
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
            %Get current EPANET version
            obj.version = obj.getVersion;
            %Get data from raw file (for information which cannot be
            %accessed by the epanet library)
            value=obj.getInputFileInfo;
            %Get curves
            obj.CurveNameID=value.CurveID;
            obj.CurveXvalue=value.CurveX;
            obj.CurveYvalue=value.CurveY;
            obj.CurveTypes =value.typeCurve;
            %Get coordinates
            obj.NodeCoordinates{1} = value.vx;
            obj.NodeCoordinates{2} = value.vy;
            obj.NodeCoordinates{3} = value.vertx;
            obj.NodeCoordinates{4} = value.verty;
            %Get Units
            obj.LinkFlowUnits=value.LinkFlowUnits;
            obj.QualityUnits=value.QualityUnits;
            obj.OptionsHeadloss=value.OptionsHeadloss;
            obj.PatternDemandsUnits=value.PatternDemandsUnits;
            obj.LinkPipeDiameterUnits=value.LinkPipeDiameterUnits;
            obj.NodeTankDiameterUnits=value.NodeTankDiameterUnits;
            obj.EnergyEfficiencyUnits=value.EnergyEfficiencyUnits;
            obj.NodeElevationUnits=value.NodeElevationUnits;
            obj.NodeEmitterCoefficientUnits=value.NodeEmitterCoefficientUnits;
            obj.EnergyUnits=value.EnergyUnits;
            obj.LinkFrictionFactorUnits=value.LinkFrictionFactorUnits;
            obj.NodeHeadUnits=value.NodeHeadUnits;
            obj.LinkLengthUnits=value.LinkLengthUnits;
            obj.LinkMinorLossCoeffUnits=value.LinkMinorLossCoeffUnits;
            obj.LinkPumpPowerUnits=value.LinkPumpPowerUnits;
            obj.NodePressureUnits=value.NodePressureUnits;
            obj.QualityReactionCoeffBulkUnits=value.QualityReactionCoeffBulkUnits;
            obj.QualityReactionCoeffWallUnits=value.QualityReactionCoeffWallUnits;
            obj.LinkPipeRoughnessCoeffUnits=value.LinkPipeRoughnessCoeffUnits;
            obj.QualitySourceMassInjectionUnits=value.QualitySourceMassInjectionUnits;
            obj.LinkVelocityUnits=value.LinkVelocityUnits;
            obj.NodeTankVolumeUnits=value.NodeTankVolumeUnits;
            obj.QualityWaterAgeUnits=value.QualityWaterAgeUnits;
        end % End of epanet class constructor
        function errcode = LoadInpFile(obj,inpname,repname,binname,varargin)
            [errcode] = ENopen(inpname,repname,binname);
        end
        function errcode = epanetLoadLibrary(obj)
            [errcode] = ENLoadLibrary;
        end
        function plot(obj,varargin)
           ENplot(obj,varargin{:});
        end
        function value = getControls(obj)
            %Retrieves the parameters of all control statements
            if obj.getControlRulesCount
                cnt=obj.getControlRulesCount;
                obj.ControlTypes{cnt}=[];
                obj.ControlTypesIndex(cnt)=NaN;
                obj.ControlLinkIndex(cnt)=NaN;
                obj.ControlSettings(cnt)=NaN;
                obj.ControlNodeIndex(cnt)=NaN;
                obj.ControlLevelValues(cnt)=NaN;
                for i=1:obj.getControlRulesCount
                    [obj.errcode, obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)] = ENgetcontrol(i);
                    obj.ControlTypes(i)={obj.TYPECONTROL(obj.ControlTypesIndex(i)+1)};
                    value{i}={obj.ControlTypes{i},obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)};
                end
            else
                value=-1;
            end
            %value=obj.ControlRules;
        end
        function value = getNodeCount(obj)
            % Retrieves the number of nodes
            [obj.errcode, value] = ENgetcount(0);
        end
        function value = getNodeTankReservoirCount(obj)
            % Retrieves the number of tanks
            [obj.errcode, value] = ENgetcount(1);
        end
        function value = getLinkCount(obj)
            % Retrieves the number of links
            [obj.errcode, value] = ENgetcount(2);
        end
        function value = getPatternCount(obj)
            % Retrieves the number of patterns
            [obj.errcode, value] = ENgetcount(3);
        end
        function value = getCurveCount(obj)
            % Retrieves the number of curves
            [obj.errcode, value] = ENgetcount(4);
        end
        function value = getControlRulesCount(obj)
            % Retrieves the number of controls
            [obj.errcode, value] = ENgetcount(5);
        end
        function value = getNodeTankCount(obj)
            value = sum(strcmp(obj.getNodeType,'TANK'));
        end
        function value = getNodeReservoirCount(obj)
            value = sum(strcmp(obj.getNodeType,'RESERVOIR'));
        end
        function value = getNodeJunctionsCount(obj)
            value = sum(strcmp(obj.getNodeType,'JUNCTION'));
        end
        function value = getLinkPipeCount(obj)
            value = sum(strcmp(obj.getLinkType,'PIPE'))+sum(strcmp(obj.getLinkType,'CVPIPE'));
        end
        function value = getLinkPumpCount(obj)
            value = sum(strcmp(obj.getLinkType,'PUMP'));
        end
        function value = getLinkValveCount(obj)
            value = obj.getLinkCount - (obj.getLinkPipeCount + obj.getLinkPumpCount);
        end      
        function value = getError(obj,errcode)
            %Retrieves the text of the message associated with a particular error or warning code.
            [obj.errcode, value] = ENgeterror(errcode);
        end
        function value = getFlowUnits(obj)
            %Retrieves flow units used to express all flow rates.
            [obj.errcode, flowunitsindex] = ENgetflowunits();
            obj.LinkFlowUnits=obj.TYPEUNITS(flowunitsindex+1);
            value=obj.LinkFlowUnits;
        end
        function value = getLinkNameID(obj,varargin)
            % Retrieves the ID label(s) of all links, or the IDs of an index set of links
            if isempty(varargin)
                value{obj.getLinkCount}=[];
                for i=1:obj.getLinkCount
                    [obj.errcode, value{i}]=ENgetlinkid(i);
                end
            else
                k=1;
                value{length(varargin{1})}=[];
                for i=varargin{1}
                    [obj.errcode, value{k}]=ENgetlinkid(i);
                    k=k+1;
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
                    [obj.errcode, value(k)] = ENgetlinkindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetlinkindex(varargin{1});
            end
        end
        function value = getLinkPipeIndex(obj)
            %Retrieves the pipe index
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PIPE'));
            if isempty(value), value=-1; end
        end        
        function value = getLinkPumpIndex(obj)
            %Retrieves the pipe index
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'PUMP'));
            if isempty(value), value=-1; end
        end
        function value = getLinkValveIndex(obj)
            %Retrieves the pipe index
            tmpLinkTypes=obj.getLinkType;
            value = find(strcmp(tmpLinkTypes,'VALVE'));
            if isempty(value), value=-1; end
        end
        function value = getLinkNodesIndex(obj)
            %Retrieves the indexes of the from/to nodes of all links.
            value(obj.getLinkCount,1:2)=[nan nan];
            for i=1:obj.getLinkCount
                [obj.errcode,linkFromNode,linkToNode] = ENgetlinknodes(i);
                value(i,:)= [linkFromNode,linkToNode];
            end
        end
        function value = getNodesConnectingLinksID(obj)
            %Retrieves the id of the from/to nodes of all links.
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
                [obj.errcode,value(i)] = ENgetlinktype(i);
                if value(i)>2
                    value(i)=9; %Valve
                elseif value(i)==1
                    value(i)=1; %cvpipe pipe
                end
            end           
        end
        function value = getLinkDiameter(obj)
            %Retrieves the value of all link diameters
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,0);
            end
        end
        function value = getLinkLength(obj)
            %Retrieves the value of all link lengths
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,1);
            end
        end
        function value = getLinkRoughnessCoeff(obj)
            %Retrieves the value of all link roughness
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,2);
            end
        end
        function value = getLinkMinorLossCoeff(obj)
            %Retrieves the value of all link minor loss coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,3);
            end
        end
        function value = getLinkInitialStatus(obj)
            %Retrieves the value of all link initial status
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,4);
            end
        end
        function value = getLinkInitialSetting(obj)
            %Retrieves the value of all link roughness for pipes or initial speed for pumps or initial setting for valves
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,5);
            end
        end
        function value = getLinkBulkReactionCoeff(obj)
            %Retrieves the value of all link bulk reaction coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,6);
            end
        end
        function value = getLinkWallReactionCoeff(obj)
            %Retrieves the value of all link wall reaction coefficients
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,7);
            end
        end
        function value = getLinkFlows(obj)
            %Retrieves the value of all computed link flow rates
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,8);
            end
        end
        function value = getLinkVelocity(obj)
            %Retrieves the value of all computed link velocities
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,9);
            end
        end
        function value = getLinkHeadloss(obj)
            %Retrieves the value of all computed link headloss
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,10);
            end
        end
        function value = getLinkStatus(obj)
            %Retrieves the value of all computed link status (0 = closed, 1 = open)
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,11);
            end
        end
        function value = getLinkSettings(obj)
            %Retrieves the value of all computed link roughness for pipes or actual speed for pumps or actual setting for valves
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,12);
            end
        end
        function value = getLinkPumpEnergy(obj) 
            %Retrieves the value of all computed energy in kwatts
            value=zeros(1,obj.getLinkCount);
            for i=1:obj.getLinkCount
                [obj.errcode, value(i)] = ENgetlinkvalue(i,13);
            end
        end
        function value = getNodeNameID(obj,varargin)
            %Retrieves the ID label of all nodes or some nodes with a specified index.
            if isempty(varargin)
                value{obj.getNodeCount}=[];
                for i=1:obj.getNodeCount
                    [obj.errcode, value{i}]=ENgetnodeid(i);
                end
            else
                k=1;
                value{length(varargin{1})}=[];
                for i=varargin{1}
                    [obj.errcode, value{k}]=ENgetnodeid(i);
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
                    [obj.errcode, value(k)] = ENgetnodeindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetnodeindex(varargin{1});
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
            if obj.getNodeJunctionsCount
                tmpNodeTypes=obj.getNodeType;
                value = find(strcmp(tmpNodeTypes,'JUNCTION'));
            else
                value=0;
            end
        end
        function value = getNodeType(obj)
            %Retrieves the node-type code for all nodes
            for i=1:obj.getNodeCount
                [obj.errcode,obj.NodeTypeIndex(i)] = ENgetnodetype(i);
                value(i)=obj.TYPENODE(obj.NodeTypeIndex(i)+1);
            end
        end
        function value = getNodeTypeIndex(obj)
            %Retrieves the node-type code for all nodes
            for i=1:obj.getNodeCount
                [obj.errcode,value(i)] = ENgetnodetype(i);
            end
        end
        function value = getNodeElevations(obj)
            %Retrieves the value of all node elevations
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,0);
            end
        end
        function value = getNodeBaseDemands(obj)
            %Retrieves the value of all node base demands
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,1);
            end
        end
        function value = getNodeDemandPatternIndex(obj)
            %Retrieves the value of all node demand pattern indices
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,2);
            end
        end
        function value = getNodeEmitterCoeff(obj)
            %Retrieves the value of all node emmitter coefficients
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,3);
            end
        end
        function value = getNodeInitialQuality(obj)
            %Retrieves the value of all node initial quality
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,4);
            end
        end
        function value = getNodeSourceQuality(obj)
            %Retrieves the value of all nodes source quality
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,5);
            end
        end
        function value = getNodeSourcePatternIndex(obj)
            %Retrieves the value of all node source pattern index
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,6);
            end
        end
        function value = getNodeSourceType(obj)
            %Retrieves the value of all node source type
            value=cell(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, temp] = ENgetnodevalue(i,7);
                if ~isnan(temp)
                    value(i)=obj.TYPESOURCE(temp+1);
                end
            end
        end
        function value = getNodeTankInitialLevel(obj)
            %Retrieves the value of all tank initial water levels
            value=nan(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,8);
            end
        end
        function value = getNodeActualDemand(obj)
            %Retrieves the computed value of all actual demands
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,9);
            end
        end
        function value = getNodeActualDemandSensingNodes(obj,varargin)
            %Retrieves the computed demand values at some sensing nodes
            value=zeros(1,obj.getNodeCount);
            for i=1:length(varargin{1})
                [obj.errcode, value(varargin{1}(i))] = ENgetnodevalue(varargin{1}(i),9);
            end
        end
        function value = getNodeHydaulicHead(obj)
            %Retrieves the computed values of all hydraulic heads
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,10);
            end
        end
        function value = getNodePressure(obj)
            %Retrieves the computed values of all node pressures
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,11);
            end
        end
        function value = getNodeActualQuality(obj)
            %Retrieves the computed values of the actual quality for all nodes
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,12);
            end
        end
        function value = getNodeMassFlowRate(obj)
            %Retrieves the computed mass flow rates per minute of chemical sources
            value=zeros(1,obj.getNodeCount);
            for i=1:obj.getNodeCount
                [obj.errcode, value(i)] = ENgetnodevalue(i,13);
            end
        end
        function value = getNodeActualQualitySensingNodes(obj,varargin)
            %Retrieves the computed quality values at some sensing nodes
            value=zeros(1,obj.getNodeCount);
            for i=1:length(varargin{1})
                [obj.errcode, value(varargin{1}(i))] = ENgetnodevalue(varargin{1}(i),12);
            end
        end
        function value = getNodeTankInitialWaterVolume(obj)
            %Retrieves the tank initial volume
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i,14);
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
                    [obj.errcode, obj.NodeTankMixingModelCode(i)] = ENgetnodevalue(i, 15);
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
                    [obj.errcode, value(i)] = ENgetnodevalue(i,16);
                end
            end
        end
        function value = getNodeTankDiameter(obj)
            %Retrieves the tank diameters
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 17);
                end
            end
        end
        function value = getNodeTankMinimumWaterVolume(obj)
            %Retrieves the tank minimum volume
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 18);
                end
            end
        end 
        function value = getNodeTankVolumeCurveIndex(obj)
            %Retrieves the tank volume curve index
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 19);
                end
            end
        end
        function value = getNodeTankMinimumWaterLevel(obj)
            %Retrieves the tank minimum water level
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 20);
                end
            end
        end
        function value = getNodeTankMaximumWaterLevel(obj)
            %Retrieves the tank maximum water level
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 21);
                end
            end
        end
        function value = getNodeTankMinimumFraction(obj)
            %Retrieves the tank Fraction of total volume occupied by the inlet/outlet zone in a 2-compartment tank
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 22);
                end
            end
        end
        function value = getNodeTankBulkReactionCoeff(obj)
            %Retrieves the tank bulk rate coefficient
            value=zeros(1,obj.getNodeCount);
            if obj.getNodeTankCount
                for i=obj.getNodeTankIndex
                    [obj.errcode, value(i)] = ENgetnodevalue(i, 23);
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
            [obj.errcode, value] = ENgetoption(0);
        end
        function value = getOptionsAccuracyValue(obj)
            % Retrieve the analysis convergence criterion (0.001)
            [obj.errcode, value] = ENgetoption(1);
        end
        function value = getOptionsQualityTolerance(obj)
            % Retrieve the water quality analysis tolerance
            [obj.errcode, value] = ENgetoption(2);
        end
        function value = getOptionsEmitterExponent(obj)
            % Retrieve power exponent for the emmitters (0.5)
            [obj.errcode, value] = ENgetoption(3);
        end
        function value = getOptionsPatternDemandMultiplier(obj)
            % Retrieve the demand multiplier (x1)
            [obj.errcode, value] = ENgetoption(4);
        end
        function value = getPatternID(obj,varargin)
            %Retrieves the ID label of all or some time patterns indices
            if isempty(varargin)
                value{obj.getPatternCount}=[];
                for i=1:obj.getPatternCount
                    [obj.errcode, value{i}]=ENgetpatternid(i);
                end
            else
                k=1;
                value{length(varargin{1})}=[];
                for i=varargin{1}
                    [obj.errcode, value{k}]=ENgetpatternid(i);
                    k=k+1;
                end
            end
        end
        function value = getPatternIndex(obj,varargin)
            %Retrieves the index of all or some time patterns IDs
            if isempty(varargin)
                value=1:obj.getPatternCount;
            elseif isa(varargin{1},'cell')
                k=1;
                value{length(varargin{1})}=[];
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetpatternindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetpatternindex(varargin{1});
            end
        end
        function value = getPatternLengths(obj,varargin)
            %Retrieves the number of time periods in all or some patterns
            if isempty(varargin)
                tmpPatterns=1:obj.getPatternCount;
                for i=tmpPatterns
                    [obj.errcode, value(i)]=ENgetpatternlen(i);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = ENgetpatternlen(obj.getPatternIndex(varargin{1}{j}));
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = ENgetpatternlen(obj.getPatternIndex(varargin{1}));
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errcode, value(k)]=ENgetpatternlen(i);
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
                    [obj.errcode, value(i,j)] = ENgetpatternvalue(i, j);
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
            [obj.errcode, value] = ENgetpatternvalue(patternIndex, patternStep);
        end
        function value = getQualityType(obj)
            %Retrieves the type of water quality analysis type
            [obj.errcode, obj.QualityCode,obj.QualityTraceNodeIndex] = ENgetqualtype();
            value=obj.TYPEQUALITY(obj.QualityCode+1);
        end
        function value = getQualityCode(obj)
            %Retrieves the code of water quality analysis type
            [obj.errcode, value,obj.QualityTraceNodeIndex] = ENgetqualtype();
        end
        function value = getQualityTraceNodeIndex(obj)
            %Retrieves the trace node index of water quality analysis type
            [obj.errcode, obj.QualityCode,value] = ENgetqualtype();
        end
        function value = getTimeSimulationDuration(obj)
            %Retrieves the value of simulation duration
            [obj.errcode, value] = ENgettimeparam(0);
        end
        function value = getTimeHydraulicStep(obj)
            %Retrieves the value of the hydraulic time step
            [obj.errcode, value] = ENgettimeparam(1);
        end
        function value = getTimeQualityStep(obj)
            %Retrieves the value of the water quality time step
            [obj.errcode, value] = ENgettimeparam(2);
        end
        function value = getTimePatternStep(obj)
            %Retrieves the value of the pattern time step
            [obj.errcode, value] = ENgettimeparam(3);
        end
        function value = getTimePatternStart(obj)
            %Retrieves the value of pattern start time
            [obj.errcode, value] = ENgettimeparam(4);
        end
        function value = getTimeReportingStep(obj)
            %Retrieves the value of the reporting time step
            [obj.errcode, value] = ENgettimeparam(5);
        end
        function value = getTimeReportingStart(obj)
            %Retrieves the value of the reporting start time
            [obj.errcode, value] = ENgettimeparam(6);
        end
        function value = getTimeRuleControlStep(obj)
            %Retrieves the time step for evaluating rule-based controls
            [obj.errcode, value] = ENgettimeparam(7);
        end
        function value = getTimeStatisticsType(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.errcode, obj.TimeStatisticsIndex] = ENgettimeparam(8);
            value=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
        end
        function value = getTimeStatisticsIndex(obj)
            %Retrieves the type of time series post-processing ('NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE')
            [obj.errcode, value] = ENgettimeparam(8);
        end
        function value = getTimeReportingPeriods(obj)
            %Retrieves the number of reporting periods saved to the binary
            %output file
            [obj.errcode, value] = ENgettimeparam(9);
        end
        function value = getVersion(obj)
            % Retrieve the current EPANET version
            [obj.errcode, value] = ENgetversion();
        end
        function value = getComputedHydraulicTimeSeries(obj)
            % Compute hydraulic simulation and retrieve all time-series
            obj.openHydraulicAnalysis;
            obj.initializeHydraulicAnalysis
            tstep=1;
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep+1;
            initnodematrix=zeros(totalsteps, obj.getNodeCount);
            initlinkmatrix=zeros(totalsteps, obj.getLinkCount);
            value.Time=zeros(totalsteps,1);
            value.Pressure=initnodematrix;
            value.Demand=initnodematrix;
            value.Head=initnodematrix;
            value.Flow=initlinkmatrix;
            value.Velocity=initlinkmatrix;
            value.HeadLoss=initlinkmatrix;
            value.Status=initlinkmatrix;
            value.Setting=initlinkmatrix;
            value.Energy=initlinkmatrix;
            k=1;
            while (tstep>0)
                t=obj.runHydraulicAnalysis;
                value.Time(k,:)=t;
                value.Pressure(k,:)=obj.getNodePressure;
                value.Demand(k,:)=obj.getNodeActualDemand;
                value.Head(k,:)=obj.getNodeHydaulicHead;
                value.Flow(k,:)=obj.getLinkFlows;
                value.Velocity(k,:)=obj.getLinkVelocity;
                value.HeadLoss(k,:)=obj.getLinkHeadloss;
                value.Status(k,:)=obj.getLinkStatus;
                value.Setting(k,:)=obj.getLinkSettings;
                value.Energy(k,:)=obj.getLinkPumpEnergy;
                tstep = obj.nextHydraulicAnalysisStep;
                k=k+1;
            end
            obj.closeHydraulicAnalysis;
        end
        function value = getComputedQualityTimeSeries(obj,varargin)
            % Compute Quality simulation and retrieve all or some time-series
            obj.solveCompleteHydraulics
            obj.solveCompleteQuality
            obj.openQualityAnalysis
            obj.initializeQualityAnalysis
            tleft=1;
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeQualityStep;
            initnodematrix=zeros(totalsteps, obj.getNodeCount);
            if size(varargin,2)==0
                varargin={'time', 'quality', 'mass'};
            end
            if find(strcmpi(varargin,'time'))
                value.Time=zeros(totalsteps,1);
            end
            if find(strcmpi(varargin,'quality'))
                value.Quality=initnodematrix;
            end
            if find(strcmpi(varargin,'mass'))
                value.MassFlowRate=initnodematrix;
            end
            if find(strcmpi(varargin,'demand'))
                value.Demand=initnodematrix;
            end
            if find(strcmpi(varargin,'qualitySensingNodes'))
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
                    value.Quality(k,:)=obj.getNodeActualQualitySensingNodes(varargin{2});
                end
                tleft = obj.stepQualityAnalysisTimeLeft;
                k=k+1;
                if t==obj.getTimeSimulationDuration
                    t=obj.getTimeSimulationDuration+1;
                end
            end
            obj.closeQualityAnalysis;
        end
        function solveCompleteHydraulics(obj)
            [obj.errcode] = ENsolveH();
        end
        function solveCompleteQuality(obj)
            [obj.errcode] = ENsolveQ();
        end
        function valueIndex = addPattern(obj,varargin)
            valueIndex=-1;
            if nargin==2
                [obj.errcode] = ENaddpattern(varargin{1});
                valueIndex = getPatternIndex(obj,varargin{1});
            elseif nargin==3
                [obj.errcode] = ENaddpattern(varargin{1});
                valueIndex = getPatternIndex(obj,varargin{1});
                setPattern(obj,valueIndex,varargin{2});
            end
        end
        function setControl(obj,controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel)
            % Example: d.setControl(1,1,13,1,11,150) controlRuleIndex must
            % exist
            if controlRuleIndex<=obj.getControlRulesCount
                [obj.errcode] = ENsetcontrol(controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel);
                obj.ControlTypes={};
                for i=1:obj.getControlRulesCount
                    [obj.errcode, obj.ControlTypesIndex(i),obj.ControlLinkIndex(i),obj.ControlSettings(i),obj.ControlNodeIndex(i),obj.ControlLevelValues(i)] = ENgetcontrol(i);
                    obj.ControlTypes(i)=obj.TYPECONTROL(obj.ControlTypesIndex(i)+1);
                end
                obj.ControlRules={obj.ControlTypes,obj.ControlTypesIndex,obj.ControlLinkIndex,obj.ControlSettings,obj.ControlNodeIndex,obj.ControlLevelValues};
            else
                disp('New rules cannot be added in this version')
            end
        end
        function setLinkDiameter(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 0, value(i));
            end
        end
        function setLinkLength(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 1, value(i));
            end
        end
        function setLinkRoughnessCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 2, value(i));
            end
        end
        function setLinkMinorLossCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 3, value(i));
            end
        end
        function setLinkInitialStatus(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 4, value(i));
            end
        end
        function setLinkInitialSetting(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 5, value(i));
            end
        end
        function setLinkBulkReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 6, value(i));
            end
        end
        function setLinkWallReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 7, value(i));
            end
        end
        function setLinkStatus(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 11, value(i));
            end
        end
        function setLinkSettings(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetlinkvalue(i, 12, value(i));
            end
        end
        function setNodeElevations(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 0, value(i));
            end
            %             [obj.errcode, obj.NodeElevations(index)] =
            %             ENgetnodevalue(index, 0);
        end
        function setNodeBaseDemands(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 1, value(i));
            end
            %[obj.errcode, obj.NodeBaseDemands(index)] =
            %ENgetnodevalue(index, 1);
        end
        function setNodeDemandPatternIndex(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 2, value(i));
            end
        end
        function setNodeEmitterCoeff(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 3, value(i));
            end
        end
        function setNodeInitialQuality(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 4, value(i));
            end
        end
        function setNodeTankLevelInitial(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 8, value(i));
            end
        end
        function setNodeTankMixingModelType(obj, value)
            for i=obj.getNodeTankIndex
                code=strfind(strcmpi(value(i),obj.TYPEMIXMODEL),1)-1;
                [obj.errcode] = ENsetnodevalue(i, 15, code);
            end
        end
        function setNodeTankDiameter(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 17, value(i));
            end
        end
        function setNodeTankMinimumWaterLevel(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 20, value(i));
            end
        end
        function setNodeTankMinimumWaterVolume(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 18, value(i));
            end
        end 
        function setNodeTankMaximumWaterLevel(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 21, value(i));
            end
        end
        function setNodeTankMinimumFraction(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 22, value(i));
            end
        end
        function setNodeTankBulkReactionCoeff(obj, value)
            for i=obj.getNodeTankIndex
                [obj.errcode] = ENsetnodevalue(i, 23, value(i));
            end
        end
        function setNodeSourceQuality(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 5, value(i));
            end
        end
        function setNodeSourcePatternIndex(obj, value)
            for i=1:length(value)
                [obj.errcode] = ENsetnodevalue(i, 6, value(i));
            end
        end
        function setNodeSourceType(obj, index, value)
            value=find(strcmpi(obj.TYPESOURCE,value)==1)-1;
            [obj.errcode] = ENsetnodevalue(index, 7, value);
        end
        function setOptionsMaxTrials(obj,value)
            [obj.errcode] = ENsetoption(0,value);
            [obj.errcode, obj.OptionsMaxTrials] = ENgetoption(0);
        end
        function setOptionsAccuracyValue(obj,value)
            [obj.errcode] = ENsetoption(1,value);
            [obj.errcode, obj.OptionsAccuracyValue] = ENgetoption(1);
        end
        function setOptionsQualityTolerance(obj,value)
            [obj.errcode] = ENsetoption(2,value);
            [obj.errcode, obj.OptionsQualityTolerance] = ENgetoption(2);
        end
        function setOptionsEmitterExponent(obj,value)
            [obj.errcode] = ENsetoption(3,value);
            [obj.errcode, obj.OptionsEmitterExponent] = ENgetoption(3);
        end
        function setOptionsPatternDemandMultiplier(obj,value)
            [obj.errcode] = ENsetoption(4,value);
            [obj.errcode, obj.OptionsPatternDemandMultiplier] = ENgetoption(4);
        end
        function setTimeSimulationDuration(obj,value)
            [obj.errcode] = ENsettimeparam(0,value);
            [obj.errcode, obj.TimeSimulationDuration] = ENgettimeparam(0);
        end
        function setTimeHydraulicStep(obj,value)
            [obj.errcode] = ENsettimeparam(1,value);
            [obj.errcode, obj.TimeHydraulicStep] = ENgettimeparam(1);
        end
        function setTimeQualityStep(obj,value)
            [obj.errcode] = ENsettimeparam(2,value);
            [obj.errcode, obj.TimeQualityStep] = ENgettimeparam(2);
        end
        function setTimePatternStep(obj,value)
            [obj.errcode] = ENsettimeparam(3,value);
            [obj.errcode, obj.TimePatternStep] = ENgettimeparam(3);
        end
        function setTimePatternStart(obj,value)
            [obj.errcode] = ENsettimeparam(4,value);
            [obj.errcode, obj.TimePatternStart] = ENgettimeparam(4);
        end
        function setTimeReportingStep(obj,value)
            [obj.errcode] = ENsettimeparam(5,value);
            [obj.errcode, obj.TimeReportingStep] = ENgettimeparam(5);
        end
        function setTimeReportingStart(obj,value)
            [obj.errcode] = ENsettimeparam(6,value);
            [obj.errcode, obj.TimeReportingStart] = ENgettimeparam(6);
        end
        function setTimeStatisticsType(obj,value)
            %'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
            tmpindex=find(strcmpi(obj.TYPESTATS,value)==1)-1;
            [obj.errcode] = ENsettimeparam(8,tmpindex);
            [obj.errcode, obj.TimeStatisticsIndex] = ENgettimeparam(8);
        end
        function setTimeRuleControlStep(obj,value)
            [obj.errcode] = ENsettimeparam(7,value);
            [obj.errcode, obj.TimeRuleControlStep] = ENgettimeparam(7);
        end
        function setPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errcode] = ENsetpattern(index, patternVector, nfactors);
        end
        function setPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errcode] = ENsetpattern(i, patternMatrix(i,:), nfactors);
            end
        end
        function setPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.errcode] = ENsetpatternvalue(index, patternTimeStep, patternFactor);
        end
        function setQualityType(obj,varargin)
            qualcode=0;
            chemname='';
            chemunits='';
            tracenode='';
            if find(strcmpi(varargin,'none')==1)
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'age')==1)
                qualcode=2;
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'chem')==1)
                qualcode=1;
                chemname=varargin{1};
                chemunits=varargin{2};
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'trace')==1)
                qualcode=3;
                tracenode=varargin{2};
                [obj.errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            end
        end
        function setReportFormatReset(obj)
            [obj.errcode]=ENresetreport();
        end
        function setReportStatus(obj,value) 
            %'yes','no','full'
            statuslevel=find(strcmpi(obj.TYPEREPORT,value)==1)-1;
            [obj.errcode] = ENsetstatusreport(statuslevel);
        end
        function setReport(obj,value)
            [obj.errcode] = ENsetreport(value);
        end
        function closeHydraulicAnalysis(obj)
            [obj.errcode] = ENcloseH();
        end
        function closeQualityAnalysis(obj)
            [obj.errcode] = ENcloseQ();
        end
        function saveHydraulicFile(obj,hydname)
            [obj.errcode]=ENsavehydfile(hydname);
        end
        function useHydraulicFile(obj,hydname)
            [obj.errcode]=ENusehydfile(hydname);
        end
        function initializeHydraulicAnalysis(obj)
            [obj.errcode] = ENinitH(1);
        end
        function initializeQualityAnalysis(obj)
            [obj.errcode] = ENinitQ(1);
        end
        function tstep = nextHydraulicAnalysisStep(obj)
            [obj.errcode, tstep] = ENnextH();
        end
        function tstep = nextQualityAnalysisStep(obj)
            [obj.errcode, tstep] = ENnextQ();
        end
        function openHydraulicAnalysis(obj)
            [obj.errcode] = ENopenH();
        end
        function openQualityAnalysis(obj)
            [obj.errcode] = ENopenQ();
        end
        function tstep = runHydraulicAnalysis(obj)
            [obj.errcode, tstep] = ENrunH();
        end
        function tstep = runQualityAnalysis(obj)
            [obj.errcode, tstep] = ENrunQ();
        end
        function saveHydraulicsOutputReportingFile(obj)
            [obj.errcode] = ENsaveH();
        end
        function tleft=stepQualityAnalysisTimeLeft(obj)
            [obj.errcode, tleft] = ENstepQ();
        end
        function saveInputFile(obj,inpname)
            [obj.errcode] = ENsaveinpfile(inpname);
            % The code below is because of a bug in EPANET 2.00.12
            % When saving using ENsaveinpfile, it does not save the type of the curves.
            value=obj.getCurveInfo;
            if ~isempty(value.CurvesID)
                obj.remAddCurvesID(value.CurvesID,value.Clines);
            end
        end
        function writeLineInReportFile(obj, line)
            [obj.errcode] = ENwriteline (line);
        end
        function writeReport(obj)
            %Writes a formatted text report on simulation results to the Report file
            [obj.errcode]=ENreport();
        end
        function unload(varargin)
            ENclose;
            ENMatlabCleanup;
        end
        function setFlowUnitsGPM(obj,newinpname)
            Options(obj,newinpname,'GPM') %gallons per minute
        end
        function setFlowUnitsLPS(obj,newinpname)
            Options(obj,newinpname,'LPS') %liters per second
        end
        function setFlowUnitsMGD(obj,newinpname)
            Options(obj,newinpname,'MGD') %million gallons per day
        end
        function setFlowUnitsIMGD(obj,newinpname)
            Options(obj,newinpname,'IMGD') %Imperial mgd
        end
        function setFlowUnitsCFS(obj,newinpname)
            Options(obj,newinpname,'CFS') %cubic feet per second
        end
        function setFlowUnitsAFD(obj,newinpname)
            Options(obj,newinpname,'AFD') %acre-feet per day
        end
        function setFlowUnitsLPM(obj,newinpname)
            Options(obj,newinpname,'LPM') %liters per minute
        end
        function setFlowUnitsMLD(obj,newinpname)
            Options(obj,newinpname,'MLD') %million liters per day
        end
        function setFlowUnitsCMH(obj,newinpname)
            Options(obj,newinpname,'CMH') %cubic meters per hour
        end
        function setFlowUnitsCMD(obj,newinpname)
            Options(obj,newinpname,'CMD') %cubic meters per day
        end
        function setHeadlossHW(obj,newinpname)
            Options(obj,newinpname,'','H-W')  %Hazen-Wiliams
        end
        function setHeadlossDW(obj,newinpname)
            Options(obj,newinpname,'','D-W')  %Darcy-Weisbach
        end
        function setHeadlossCM(obj,newinpname)
            Options(obj,newinpname,'','C-M')  %Chezy-Manning
        end
        function msx(obj,msxname)
            [obj] = MSXMatlabSetup(obj,msxname);
        end
        function value = getMsxEquationsTerms(obj)
            [value,~,~] = getEquations(obj.MSXFile);
        end
        function value = getMsxEquationsPipes(obj)
            [~,value,~] = getEquations(obj.MSXFile);
        end
        function value = getMsxEquationsTanks(obj)
            [~,~,value] = getEquations(obj.MSXFile);
        end
        function value = getMsxTimeStep(obj)
            [value] = MsxTimeStep(obj.MSXFile);
        end
        function value = getMsxSpeciesCount(obj)
            % Species, Constants, Parameters, Patterns
            [obj.errcode, value] = MSXgetcount(3);
        end
        function value = getMsxConstantsCount(obj)
            [obj.errcode, value] = MSXgetcount(6);
        end
        function value = getMsxParametersCount(obj)
            [obj.errcode, value] = MSXgetcount(5);
        end
        function value = getMsxPatternsCount(obj)
            [obj.errcode, value] = MSXgetcount(7);
        end
        function value = getMsxSpeciesNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode, len] = MSXgetIDlen(3,i);
                    [obj.errcode, value{i}]=MSXgetID(3,i,len);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(3,i);
                    [obj.errcode, value{k}]=MSXgetID(3,i,len);
                    k=k+1;
                end
            end
        end  
        function value = getMsxSpeciesType(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,value{i},~,~,~] = MSXgetspecies(i);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesUnits(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,value{i},~,~] = MSXgetspecies(i);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesATOL(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,value,~] = MSXgetspecies(i);
                end
            else
                value=-1;
            end
        end
        function value = getMsxSpeciesRTOL(obj)
            if obj.getMsxSpeciesCount
                for i=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,~,value] = MSXgetspecies(i);
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
                    [obj.errcode, value(k)] = MSXgetindex(3,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(3,varargin{1});
            end
        end
        function value = getMsxConstantsNameID(obj)
            if obj.getMsxConstantsCount
                for i=1:obj.getMsxConstantsCount
                    [obj.errcode, len] = MSXgetIDlen(6,i);
                    [obj.errcode, value{i}] = MSXgetID(6,i,len);
                end
            else
                value=-1;
            end
        end
        function value = getMsxConstantsValue(obj)
            if obj.getMsxConstantsCount
                for i=1:obj.getMsxConstantsCount
                    [obj.errcode, value(i)] = MSXgetconstant(i);
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
                    [obj.errcode, value(k)] = MSXgetindex(6,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(6,varargin{1});
            end
        end
        function value = getMsxParametersNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxParametersCount
                    [obj.errcode, len] = MSXgetIDlen(5,i);
                    [obj.errcode, value{i}]=MSXgetID(5,i,len);
                end
                if ~obj.getMsxParametersCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(5,i);
                    [obj.errcode, value{k}]=MSXgetID(5,i,len);
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
                    [obj.errcode, value(k)] = MSXgetindex(5,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetindex(5,varargin{1});
            end
        end
        function value = getMsxParametersTanksValue(obj)
            value={};
            if ~obj.getMsxParametersCount, value=0;return;end
            if ~length(obj.NodeTankIndex), value=0;return;end
            for i=1:length(obj.getNodeTankCount)
                for j=1:obj.MsxParametersCount
                    [obj.errcode, value{obj.NodeTankIndex(i)}(j)] = MSXgetparameter(0,obj.NodeTankIndex(i),j);
                end
            end
        end
        function value = getMsxParametersPipesValue(obj)
            if ~obj.getMsxParametersCount
                value=0;return;
            end
            for i=1:obj.getLinkPipeCount
                for j=1:obj.getMsxParametersCount
                    [obj.errcode, value{i}(j)] = MSXgetparameter(1,i,j);
                end
            end
        end
        function value = getMsxPatternsNameID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getMsxPatternsCount
                    [obj.errcode, len] = MSXgetIDlen(7,i);
                    [obj.errcode, value{i}]=MSXgetID(7,i,len);
                end
                if ~obj.getMsxPatternsCount
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errcode, len] = MSXgetIDlen(7,i);
                    [obj.errcode, value{k}]=MSXgetID(7,i,len);
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
                    [obj.errcode, len] = MSXgetIDlen(7,j);
                    [obj.errcode, value{k}] = MSXgetID(7, obj.MsxPatternsIndex,len);
                    if obj.errcode
                        value{k}=0;
                    end
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, obj.MsxPatternsIndex] = MSXgetindex(obj,7,varargin{1});
                [obj.errcode, len] = MSXgetIDlen(7,obj.MsxPatternsIndex);
                [obj.errcode, value] = MSXgetID(7, obj.MsxPatternsIndex,len);
                if obj.errcode
                    value=0;
                end
            end
        end
        function value = getMsxPatternsLengths(obj,varargin)
            if isempty(varargin)
                if obj.getMsxPatternsCount
                    for i=obj.getMsxPatternsIndex
                        [obj.errcode, value(i)]=MSXgetpatternlen(i);
                    end
                else
                    value=-1;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errcode, value(k)] = MSXgetpatternlen(obj.getMsxPatternsIndex(varargin{1}{j}));
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errcode, value] = MSXgetpatternlen(obj.getMsxPatternsIndex(varargin{1}));
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errcode, value(k)]=MSXgetpatternlen(i);
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
                    [obj.errcode, value{i}(j)] = MSXgetinitqual(0,i,j);
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
                    [obj.errcode, value{i}(j)] = MSXgetinitqual(1,i,j);
                end
            end
        end
        function value = getMsxSources(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, obj.MsxSourceType{i}{j},obj.MsxSourceLevel{i}(j),obj.MsxSourcePatternIndex{i}(j)] = MSXgetsource(i,j);
                end
            end
            value={obj.MsxSourceType,obj.MsxSourceLevel,obj.MsxSourcePatternIndex,obj.getMsxSourceNodeNameID};
        end
        function value = getMsxSourceNodeNameID(obj)
             value = obj.getNodeNameID;
        end
        function value = getMsxSourceType(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode, value{i}{j},~,~] = MSXgetsource(i,j);
                end
            end
        end
        function value = getMsxSourceLevel(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,value{i}(j),~] = MSXgetsource(i,j);
                end
            end
        end
        function value = getMsxSourcePatternIndex(obj)
            for i=1:obj.getNodeCount
                for j=1:obj.getMsxSpeciesCount
                    [obj.errcode,~,~,value{i}(j)] = MSXgetsource(i,j);
                end
            end
        end
        function value = getMsxPattern(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getMsxPatternsLengths);
            value=nan(obj.getMsxPatternsCount,tmpmaxlen);
            for i=1:obj.getMsxPatternsCount
                tmplength=obj.getMsxPatternsLengths(i);
                for j=1:tmplength
                    [obj.errcode, value(i,j)] = MSXgetpatternvalue(i, j);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                
            end
        end
        function value = getMsxPatternValue(obj,patternIndex, patternStep) %Mass flow rate per minute of a chemical source
            [obj.errcode, value] = MSXgetpatternvalue(patternIndex, patternStep);
        end
        function value = getMsxSpeciesConcentration(obj, type, index, species)
            [obj.errcode, value] = MSXgetqual(type, index, species);
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
            end
            % Obtain a hydraulic solution
            obj.MsxSolveCompleteHydraulics();
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.MsxInitializeQualityAnalysis(0);
            % Retrieve species concentration at node
            k=1; tleft=1;t=0;i=1;
            timeSmle=obj.getTimeSimulationDuration;%bug at time
            while(tleft>0 && obj.errcode==0 && timeSmle~=t)
                [t, tleft]=obj.MsxStepQualityAnalysisTimeLeft();
                for j=uu
                    value.Quality{j,i}(k,:)=obj.getMsxSpeciesConcentration(0, ss, j);%node code0
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
            end
            % Obtain a hydraulic solution
            obj.MsxSolveCompleteHydraulics();
            % Run a step-wise water quality analysis without saving
            % RESULTS to file
            obj.MsxInitializeQualityAnalysis(0);

            % Retrieve species concentration at node
            k=1;tleft=1;i=1;
            while(tleft>0 && obj.errcode==0)
                [t, tleft]=obj.MsxStepQualityAnalysisTimeLeft();
                for j=uu
                    value.Quality{j,i}(k,:)=obj.getMsxSpeciesConcentration(1, ss, j); 
                end
                value.Time(k,:)=t;
                k=k+1;
            end
        end
        function MsxPlotConcentrationSpeciesOfNodes(obj,varargin)
            s=obj.getMsxComputedQualityNode(varargin{1},varargin{2});
            nodesID=obj.getNodeNameID;
            SpeciesNameID=obj.getMsxSpeciesNameID;
%             SpCnt=obj.getMsxSpeciesCount;
%             NodCnt=obj.getNodeCount;
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
%             SpCnt=obj.getMsxSpeciesCount;
%             LinkCnt=obj.getLinkCount;
%             for l=1:LinkCnt
            for l=varargin{1}
                linkID=linksID(l);
                figure('Name',['LINK ',char(linkID)]);
%                 for i=1:SpCnt
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
            [obj.errcode, value] = MSXgeterror(errcode);
        end
        function MsxSolveCompleteHydraulics(obj)
            [obj.errcode] = MSXsolveH();
        end
        function MsxSolveCompleteQuality(obj)
            [obj.errcode] = MSXsolveQ();
        end
        function MsxWriteReport(obj)
            [obj.errcode]=MSXreport();
        end
        function index = MsxAddPattern(obj,varargin)
            index=-1;
            if nargin==2
                [obj.errcode] = MSXaddpattern(varargin{1});
                [obj.errcode, index] = MSXgetindex(obj,7,varargin{1});
            elseif nargin==3
                [obj.errcode] = MSXaddpattern(varargin{1});
                [obj.errcode, index] = MSXgetindex(obj,7,varargin{1});
                setMsxPattern(obj,index,varargin{2});
            end
        end
        function setMsxSources(obj, node, species, type, level, pat)
            MSXsetsource(node, species, type, level, pat);
        end
        function setMsxConstantsValue(obj, value)
            for i=1:length(value)
                [obj.errcode] = MSXsetconstant(i, value(i));
            end
        end
        function setMsxParametersTanksValue(obj, NodeTankIndex, value)
            if ~sum(NodeTankIndex==obj.NodeTankIndex)
                fprintf('>> Invalid Tank Index <<\n');obj.NodeTankIndex
                return;
            end
            for i=1:length(value)
                [obj.errcode] = MSXsetparameter(0, NodeTankIndex, i, value(i));
            end
        end
        function setMsxParametersPipesValue(obj, pipeIndex, value)
            for i=1:length(value)
                [obj.errcode] = MSXsetparameter(1, pipeIndex, i, value(i));
            end
        end
        function setMsxNodeInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.errcode] = MSXsetinitqual(0, i, j, value{i}(j));
                end
            end
        end
        function setMsxLinkInitqualValue(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.errcode] = MSXsetinitqual(1, i, j, value{i}(j));
                end
            end
        end
        function setMsxPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errcode] = MSXsetpattern(index, patternVector, nfactors);
        end
        function setMsxPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errcode] = MSXsetpattern(i, patternMatrix(i,:), nfactors);
            end
        end
        function setMsxPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.errcode] = MSXsetpatternvalue(index, patternTimeStep, patternFactor);
        end
        function MsxSaveQualityFile(obj,outfname)
            [obj.errcode]=MSXsaveoutfile(outfname);
        end
        function MsxUseHydraulicFile(obj,hydname)
            [obj.errcode]=MSXusehydfile(hydname);
        end
        function MsxInitializeQualityAnalysis(obj,flag)
            [obj.errcode] = MSXinit(flag);
        end
        function [t, tleft]=MsxStepQualityAnalysisTimeLeft(obj)
            [obj.errcode, t, tleft] = MSXstep();
        end
        function MsxSaveFile(obj,msxname)
            [obj.errcode] = MSXsavemsxfile(msxname);
        end
        function MsxUnload(varargin)
            MSXclose;
            MSXMatlabCleanup;
        end
        function NodeCoordinates = getCoordinates(obj)
            [vx,vy,vertx,verty]  = getNodeCoord(obj);
            NodeCoordinates{1} = vx;
            NodeCoordinates{2} = vy;
            NodeCoordinates{3} = vertx;
            NodeCoordinates{4} = verty;
        end
        function value = getInputFileInfo(obj)
            value=InputFileInfo(obj);
        end
        function value = getCurveInfo(obj)
            [value.CurvesID,value.CurveX,value.CurveY,value.Clines,value.typeCurve]=CurveInfo(obj);
        end
        function value = getLinksInfo(obj)
            value=LinksInfo(obj);
        end
        function value = getNodesInfo(obj)
            value=NodesInfo(obj);
        end
        function value = getPumpInfo(obj)
            value=PumpInfo(obj);
        end
        function value = getControlsInfo(obj)
            value=ControlsInfo(obj);
        end
        function value = getFlowUnitsHeadlossFormula(obj)
            value=FlowUnitsHeadlossFormula(obj);
        end
        function remAddCurvesID(obj,CurveID,Clines)
            remAddCurve(obj,CurveID,Clines);
        end
    end
end
function [errcode] = ENwriteline (line)
[errcode]=calllib('epanet2','ENwriteline',line);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENaddpattern(patid)
errcode=calllib('epanet2','ENaddpattern',patid);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENclose()
[errcode]=calllib('epanet2','ENclose');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENcloseH()
[errcode]=calllib('epanet2','ENcloseH');
if errcode
    ENerror(errcode);
end
end
function [errcode, value] = ENgetnodevalue(index, paramcode)
[errcode, value]=calllib('epanet2','ENgetnodevalue',index, paramcode, 0);
if errcode==240
    value=NaN;
end
end
function [errcode] = ENcloseQ()
[errcode]=calllib('epanet2','ENcloseQ');
if errcode
    ENerror(errcode);
end
end
function [e] = ENerror(errcode)
errstring=char(32*ones(1,80));
len=80;
[e,errstring] = calllib('epanet2','ENgeterror',errcode,errstring,len);
disp(errstring);
end
function [errcode, ctype,lindex,setting,nindex,level] = ENgetcontrol(cindex)
[errcode, ctype,lindex,setting,nindex,level]=calllib('epanet2','ENgetcontrol',cindex,0,0,0,0,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, count] = ENgetcount(countcode)
[errcode,count]=calllib('epanet2','ENgetcount',countcode,0);
if errcode
    ENerror(errcode);
end
end
function [e, errmsg] = ENgeterror(errcode)
errmsg = char(32*ones(1,80));
[e,errmsg] = calllib('epanet2','ENgeterror',errcode,errmsg,80);
if e
    ENerror(e);
end
end
function [errcode,flowunitsindex] = ENgetflowunits()
[errcode, flowunitsindex]=calllib('epanet2','ENgetflowunits',0);
if errcode
    ENerror(errcode);
end
end
function [errcode,id] = ENgetlinkid(index)
id=char(32*ones(1,17));
[errcode,id]=calllib('epanet2','ENgetlinkid',index,id);
if errcode
    ENerror(errcode);
end
end
function [errcode,index] = ENgetlinkindex(id)
[errcode,~,index]=calllib('epanet2','ENgetlinkindex',id,0);
if errcode
    ENerror(errcode);
end
end
function [errcode,from,to] = ENgetlinknodes(index)
[errcode,from,to]=calllib('epanet2','ENgetlinknodes',index,0,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, type] = ENgetlinktype(index)
[errcode,type]=calllib('epanet2','ENgetlinktype',index,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, value] = ENgetlinkvalue(index, paramcode)
[errcode,value]=calllib('epanet2','ENgetlinkvalue',index, paramcode, 0);
if errcode
    ENerror(errcode);
end
end
function [errcode,id] = ENgetnodeid(index)
id=char(32*ones(1,17));
[errcode,id]=calllib('epanet2','ENgetnodeid',index,id);
if errcode
    ENerror(errcode);
end
end
function [errcode,index] = ENgetnodeindex(id)
[errcode, ~, index]=calllib('epanet2','ENgetnodeindex',id,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, type] = ENgetnodetype(index)
[errcode,type]=calllib('epanet2','ENgetnodetype',index,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, value] = ENgetoption(optioncode)
[errcode,value]=calllib('epanet2','ENgetoption',optioncode,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, id] = ENgetpatternid(index)
id=char(32*ones(1,31));
[errcode,id]=calllib('epanet2','ENgetpatternid',index,id);
if errcode
    ENerror(errcode);
end
end
function [errcode, index] = ENgetpatternindex(id)
[errcode,~, index]=calllib('epanet2','ENgetpatternindex',id,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, len] = ENgetpatternlen(index)
[errcode,len]=calllib('epanet2','ENgetpatternlen',index,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, value] = ENgetpatternvalue(index, period)
[errcode,value]=calllib('epanet2','ENgetpatternvalue',index, period, 0);
if errcode
    ENerror(errcode);
end
end
function [errcode,qualcode,tracenode] = ENgetqualtype()
[errcode,qualcode,tracenode]=calllib('epanet2','ENgetqualtype',0,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, timevalue] = ENgettimeparam(paramcode)
[errcode,timevalue]=calllib('epanet2','ENgettimeparam',paramcode,0);
if errcode
    ENerror(errcode);
end
end
function [errcode, version] = ENgetversion()
[errcode,version]=calllib('epanet2','ENgetversion',0);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENinitH(flag)
[errcode]=calllib('epanet2','ENinitH',flag);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENinitQ(saveflag)
[errcode]=calllib('epanet2','ENinitQ',saveflag);
if errcode
    ENerror(errcode);
end
end
function ENMatlabCleanup()
%     if nargin == 1
% ENDLLNAME='epanet2';
%     end;
% Load library
if libisloaded('epanet2')
    unloadlibrary('epanet2');
else
    errstring =['Library ', 'epanet2', '.dll was not loaded'];
    disp(errstring);
end;
end
function [errcode] = ENLoadLibrary()
% try
%     unloadlibrary('epanet2')
% catch err
% end
errcode=0;
if ~libisloaded('epanet2')
    loadlibrary('epanet2','epanet2.h')
     pause(0.1);  
end
end
function [errcode, tstep] = ENnextH()
[errcode,tstep]=calllib('epanet2','ENnextH',0);
if errcode
    ENerror(errcode);
end
end
function [errcode, tstep] = ENnextQ()
[errcode,tstep]=calllib('epanet2','ENnextQ',0);
if errcode
    ENerror(errcode);
end
tstep = double(tstep);
end
function [errcode] = ENopen(inpname,repname,binname,varargin)
errcode=calllib('epanet2','ENopen',inpname,repname,binname);
% while errcode~=0
%    try
        errcode=calllib('epanet2','ENopen',inpname,repname,binname);
        if errcode==302
            unloadlibrary('epanet2');
            ENLoadLibrary;
        end
      pause(0.5);     
    % catch err
 %   end
%end
end
function [errcode] = ENopenH()
[errcode]=calllib('epanet2','ENopenH');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENopenQ()
[errcode]=calllib('epanet2','ENopenQ');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENreport()
[errcode]=calllib('epanet2','ENreport');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENresetreport()
[errcode]=calllib('epanet2','ENresetreport');
if errcode
    ENerror(errcode);
end
end
function [errcode, t] = ENrunH()
[errcode,t]=calllib('epanet2','ENrunH',0);
t = double(t);
if errcode
    ENerror(errcode);
end
end
function [errcode, t] = ENrunQ()
[errcode,t]=calllib('epanet2','ENrunQ',0);
if errcode
    ENerror(errcode);
end
t = double(t);
end
function [errcode] = ENsaveH()
[errcode]=calllib('epanet2','ENsaveH');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsavehydfile(fname)
[errcode]=calllib('epanet2','ENsavehydfile',fname);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsaveinpfile(inpname)
errcode=calllib('epanet2','ENsaveinpfile',inpname);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetcontrol(cindex,ctype,lindex,setting,nindex,level)
[errcode]=calllib('epanet2','ENsetcontrol',cindex,ctype,lindex,setting,nindex,level);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetlinkvalue(index, paramcode, value)
[errcode]=calllib('epanet2','ENsetlinkvalue',index, paramcode, value);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetnodevalue(index, paramcode, value)
[errcode]=calllib('epanet2','ENsetnodevalue',index, paramcode, value);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetoption(optioncode,value)
[errcode]=calllib('epanet2','ENsetoption',optioncode,value);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetpattern(index, factors, nfactors)
[errcode]=calllib('epanet2','ENsetpattern',index,factors,nfactors);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetpatternvalue(index, period, value)
[errcode]=calllib('epanet2','ENsetpatternvalue',index, period, value);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode)
[errcode]=calllib('epanet2','ENsetqualtype',qualcode,chemname,chemunits,tracenode);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetreport(command)
[errcode]=calllib('epanet2','ENsetreport',command);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsetstatusreport(statuslevel)
[errcode]=calllib('epanet2','ENsetstatusreport',statuslevel);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsettimeparam(paramcode, timevalue)
paramcode=int32(paramcode);
timevalue=int32(timevalue);
[errcode]=calllib('epanet2','ENsettimeparam',paramcode,timevalue);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsolveH()
[errcode]=calllib('epanet2','ENsolveH');
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENsolveQ()
[errcode]=calllib('epanet2','ENsolveQ');
if errcode
    ENerror(errcode);
end
end
function [errcode, tleft] = ENstepQ()
tleft=int32(0);
[errcode,tleft]=calllib('epanet2','ENstepQ',tleft);
if errcode
    ENerror(errcode);
end
end
function [errcode] = ENusehydfile(hydfname)
[errcode]=calllib('epanet2','ENusehydfile',hydfname);
if errcode
    ENerror(errcode);
end
end
function ENplot(obj,varargin)

% Initiality
highlightnode=0;
highlightlink=0;
highlightnodeindex=[];
highlightlinkindex=[];
Node=char('no');
Link=char('no');
fontsize=10;
selectColorNode={''};
selectColorLink={''};

for i=1:(nargin/2)
    argument =lower(varargin{2*(i-1)+1});
    switch argument
        case 'nodes' % Nodes
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            Node=varargin{2*i};
        case 'links' % Links
            if ~strcmp(lower(varargin{2*i}),'yes') && ~strcmp(lower(varargin{2*i}),'no')
                warning('Invalid argument.');
                return
            end
            Link=varargin{2*i};
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
        otherwise
            warning('Invalid property found.');
            return
    end
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
cla
% Get node names and x, y coordiantes
%     NodeCoordinates = obj.getCoordinates; nodes = obj.getNodesInfo;
value=obj.getInputFileInfo;

NodeCoordinates{1} = value.vx;
NodeCoordinates{2} = value.vy;
NodeCoordinates{3} = value.vertx;
NodeCoordinates{4} = value.verty;
if isa(highlightnode,'cell')
    for i=1:length(highlightnode)
        n = strcmp(value.NodesAll,highlightnode{i});
        if sum(n)==0
            warning('Undefined node with id "%s" in function call therefore the index is zero.', char(highlightnode{i}));
        else
            highlightnodeindex(i) = strfind(n,1);
        end
    end
end

if isa(highlightlink,'cell')
    for i=1:length(highlightlink)
        n = strcmp(value.LinksAll,highlightlink{i});
        if sum(n)==0
            warning('Undefined link with id "%s" in function call therefore the index is zero.', char(highlightlink{i}));
        else
            highlightlinkindex(i) = strfind(n,1);
        end
    end
end

for i=1:value.LinkCount
    FromNode=strfind(strcmp(value.FromNode{i},value.NodesAll),1);
    ToNode=strfind(strcmp(value.ToNode{i},value.NodesAll),1);
    
    if FromNode
        x1 = double(NodeCoordinates{1}(FromNode));
        y1 = double(NodeCoordinates{2}(FromNode));
    end
    if ToNode
        x2 = double(NodeCoordinates{1}(ToNode));
        y2 = double(NodeCoordinates{2}(ToNode));
    end
    
    hh=strfind(highlightlinkindex,i);
    
    %         h(:,4)=line([x1,x2],[y1,y2],'LineWidth',1);
    h(:,4)=line([x1 NodeCoordinates{3}{i} x2],[y1 NodeCoordinates{4}{i} y2],'LineWidth',1);
    
    legendString{4} = char('Pipes');
    % Plot Pumps
    if sum(strfind(value.LinkPumpIndex,i))
        colornode = 'm';
        if length(hh) && isempty(selectColorLink)
            colornode = 'r';
        end
        h(:,5)=plot((x1+x2)/2,(y1+y2)/2,'mv','LineWidth',2,'MarkerEdgeColor','m',...
            'MarkerFaceColor','m',...
            'MarkerSize',5);
        plot((x1+x2)/2,(y1+y2)/2,'mv','LineWidth',2,'MarkerEdgeColor',colornode,...
            'MarkerFaceColor',colornode,...
            'MarkerSize',5);
        
        legendString{5} = char('Pumps');
    end
    
    % Plot Valves
    if sum(strfind(value.LinkValveIndex,i))
        colornode = 'k';
        if length(hh) && isempty(selectColorLink)
            colornode = 'r';
        end
        h(:,6)=plot((x1+x2)/2,(y1+y2)/2,'k*','LineWidth',2,'MarkerEdgeColor',colornode,...
            'MarkerFaceColor',colornode,'MarkerSize',7);
        legendString{6} = char('Valves');
    end
    
    % Show Link id
    if (strcmp(lower(Link),'yes') && ~length(hh))
        text((x1+x2)/2,(y1+y2)/2,value.LinksAll(i),'Fontsize',fontsize);
    end
    
    if length(hh) && isempty(selectColorLink)
        line([x1,x2],[y1,y2],'LineWidth',2,'Color','r');
        text((x1+x2)/2,(y1+y2)/2,value.LinksAll(i),'Fontsize',fontsize);
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
                line([x1 NodeCoordinates{3}{i} x2],[y1 NodeCoordinates{4}{i} y2],'LineWidth',2,'Color',nm{1}{1});
            else
                line([x1 NodeCoordinates{3}{i} x2],[y1 NodeCoordinates{4}{i} y2],'LineWidth',2,'Color',nm{1});
            end
        else
            line([x1 NodeCoordinates{3}{i} x2],[y1 NodeCoordinates{4}{i} y2],'LineWidth',2,'Color',char(selectColorLink(hh)));
        end
        text((x1+x2)/2,(y1+y2)/2,value.LinksAll(i),'Fontsize',fontsize);
    end
    
    hold on
end

% Coordinates for node FROM
for i=1:value.NodeCount
    [x] = double(NodeCoordinates{1}(i));
    [y] = double(NodeCoordinates{2}(i));
    
    hh=strfind(highlightnodeindex,i);
    h(:,1)=plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerSize',5);
    legendString{1}= char('Junctions');
    
    % Plot Reservoirs
    if sum(strfind(value.NodeReservoirIndex,i))
        colornode = 'g';
        if length(hh) && isempty(selectColorNode)
            colornode = 'r';
        end
        h(:,2)=plot(x,y,'s','LineWidth',2,'MarkerEdgeColor','g',...
            'MarkerFaceColor','g',...
            'MarkerSize',13);
        plot(x,y,'s','LineWidth',2,'MarkerEdgeColor', colornode,...
            'MarkerFaceColor', colornode,...
            'MarkerSize',13);
        legendString{2} = char('Reservoirs');
    end
    % Plot Tanks
    if sum(strfind(value.NodeTankIndex,i))
        colornode='c';
        if length(hh) && isempty(selectColorNode)
            colornode='r';
        elseif length(hh) && ~isempty(selectColorNode)
            colornode= 'b';
        end
        h(:,3)=plot(x,y,'p','LineWidth',2,'MarkerEdgeColor','c',...
            'MarkerFaceColor','c',...
            'MarkerSize',16);
        
        plot(x,y,'p','LineWidth',2,'MarkerEdgeColor',colornode,...
            'MarkerFaceColor',colornode,...
            'MarkerSize',16);
        
        legendString{3} = char('Tanks');
    end
    
    % Show Node id
    if (strcmp(lower(Node),'yes') && ~length(hh))
        text(x,y,value.NodesAll(i),'Fontsize',fontsize);%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
    end
    
    if length(hh) && isempty(selectColorNode)
        plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','r',...
            'MarkerFaceColor','r',...
            'MarkerSize',10)
        text(x,y,value.NodesAll(i),'Fontsize',fontsize)%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
    elseif length(hh) && ~isempty(selectColorNode)
        try 
            tt=length(selectColorNode{hh});
        catch err
            tt=2;
        end
       if tt>1
            if length(selectColorNode(hh))==1
                nm{1}=selectColorNode(hh);
            else
                nm=selectColorNode(hh);
            end
            if iscell(nm{1}) 
                plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',nm{1}{1},'MarkerFaceColor',nm{1}{1},'MarkerSize',10)
            else
                plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',nm{1},'MarkerFaceColor',nm{1},'MarkerSize',10)
            end
       else
        plot(x, y,'o','LineWidth',2,'MarkerEdgeColor',char(selectColorNode(hh)),'MarkerFaceColor',char(selectColorNode(hh)),...
            'MarkerSize',10)
       end
        text(x,y,value.NodesAll(i),'Fontsize',fontsize)%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
    end
    hold on
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
[xmax,~]=max(NodeCoordinates{1});
[xmin,~]=min(NodeCoordinates{1});
[ymax,~]=max(NodeCoordinates{2});
[ymin,~]=min(NodeCoordinates{2});

%     xmax=yxmax(1); ymax=yxmax(2); xmin=yxmin(1); ymin=yxmin(2);
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
end

function [vx,vy,vertx,verty] = getNodeCoord(obj)
% Initialize
nodes=obj.getNodesInfo;
links=obj.getLinksInfo;

vx = NaN(nodes.NodeCount,1);
vy = NaN(nodes.NodeCount,1);
vertx = cell(links.LinkCount,1);
verty = cell(links.LinkCount,1);
nvert = zeros(links.LinkCount,1);
% Open epanet input file
[fid,message]=fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end

sect = 0;%i=1;
% Read each line from input file.
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    
    % Get first token in the line
    tok = strtok(tline);
    
    % Skip blank Clines and comments
    if isempty(tok), continue, end
    if (tok(1) == ';'), continue, end
    
    % Check if at start of a new COOR or VERT section
    if (tok(1) == '[')
        % [COORDINATES] section
        if strcmpi(tok(1:5),'[COOR')
            sect = 1;
            continue;
            % [VERTICES] section
        elseif strcmpi(tok(1:5),'[VERT')
            sect = 2;
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
        
        % Coordinates
    elseif sect == 1
        A = textscan(tline,'%s %f %f');
        % get the node index
        a=strcmp(A{1},nodes.NodesAll);
        index=strfind(a,1);
        if length(index)==0
            return;
        end
        vx(index) = A{2};
        vy(index) = A{3};
        % Vertices
    elseif sect == 2
        A = textscan(tline,'%s %f %f');
        [errcode,index] = ENgetlinkindex(char(A{1}));
        if errcode ~=0
            return;
        end
        nvert(index) = nvert(index) + 1;
        vertx{index}(nvert(index)) = A{2};
        verty{index}(nvert(index)) = A{3};
    end
end
end
function [obj] = MSXMatlabSetup(obj,msxname)
if ~libisloaded('epanetmsx')
    loadlibrary('epanetmsx','epanetmsx.h');
end
obj.MSXPathFile = [pwd,'\NETWORKS\',msxname];
obj.MSXFile = msxname;
%Open the file
[obj.errcode] = MSXopen(obj.MSXPathFile);
%Set path of temporary file
obj.pathfileMsx=[pwd,'\RESULTS\temp.msx'];
if ~strcmpi(msxname,'temp.msx')
    dir_struct = dir(strcat(obj.pathfileMsx));
    [~,sorted_index] = sortrows({dir_struct.name}');
    if length(sorted_index), delete(obj.pathfileMsx); end
end
%Save the temporary input file
obj.MsxSaveFile(obj.pathfileMsx);
%Close msx file
MSXclose;
[obj.errcode] = MSXopen(obj.pathfileMsx);
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
function [errcode] = MSXopen(pathfile)

[errcode] = calllib('epanetmsx','MSXopen',pathfile);
if errcode
    MSXerror(errcode);
end
if (errcode == 520)
    disp('current MSX project will be closed and the new project will be opened');
    [errcode] = MSXclose();
    if errcode
        MSXerror(errcode);
    else
        [errcode] = calllib('epanetmsx','MSXopen',pathfile);
        if errcode
            MSXerror(errcode);
        end
    end
end
end
function [errcode] = MSXclose()
[errcode] = calllib('epanetmsx','MSXclose');
if errcode
    MSXerror(errcode);
end
end
function [e] = MSXerror(errcode)
e=0; len=80;
errstring=char(32*ones(1,len+1));
[e,errstring] = calllib('epanetmsx','MSXgeterror',errcode,errstring,len);
disp(errstring);
end
function [errcode, count] = MSXgetcount(code)
count=0;
[errcode,count] = calllib('epanetmsx','MSXgetcount',code,count);
if errcode
    MSXerror(errcode);
end
end
function [errcode, index] = MSXgetindex(varargin)
index =0;
if ~isnumeric(varargin{1})
    varargin{1}=varargin{2};
    varargin{2}=varargin{3};
end
[errcode,id,index]=calllib('epanetmsx','MSXgetindex',varargin{1},varargin{2},index);
if errcode
    MSXerror(errcode);
end
end
function [errcode, id] = MSXgetID(type,index,len)
id=char(32*ones(1,len+1));
[errcode,id]=calllib('epanetmsx','MSXgetID',type,index,id,len);
id=id(1:len);
if errcode
    MSXerror(errcode);
end
end
function [errcode, len] = MSXgetIDlen(type,index)
len=0;
[errcode,len]=calllib('epanetmsx','MSXgetIDlen',type,index,len);
if errcode
    MSXerror(errcode);
end
end
function [errcode, type, units, atol, rtol] = MSXgetspecies(index)
type=0; rtol=0; atol=0;
units=char(32*ones(1,16));
[errcode,type,units,atol,rtol]=calllib('epanetmsx','MSXgetspecies',index,type,units,atol,rtol);
switch type
    case 0
        type='BULK';   % for a bulk water species
    case 1
        type='WALL';   % for a pipe wall surface species
end
if errcode
    MSXerror(errcode);
end
end
function [errcode, value] = MSXgetconstant(index)
value=0;
[errcode,value]=calllib('epanetmsx','MSXgetconstant',index,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode, value] = MSXgetparameter(type,index,param)
value=0;
[errcode,value]=calllib('epanetmsx','MSXgetparameter',type,index,param,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode, patlen] = MSXgetpatternlen(patindex)
patlen=0;
[errcode,patlen]=calllib('epanetmsx','MSXgetpatternlen',patindex,patlen);
if errcode
    MSXerror(errcode);
end
end
function [errcode, value] = MSXgetpatternvalue(patindex,period)
value=0;
[errcode,value]=calllib('epanetmsx','MSXgetpatternvalue',patindex,period,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode, value] = MSXgetinitqual(obj,index,species)
value=0;
[errcode,value]=calllib('epanetmsx','MSXgetinitqual',obj,index,species,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode, type, level, pat] = MSXgetsource(node,species)
type=0;
level=0;
pat=0;
[errcode,type,level,pat]=calllib('epanetmsx','MSXgetsource',node,species,type,level,pat);
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
    MSXerror(errcode);
end
end
function MSXMatlabCleanup()
% Unload library
if libisloaded('epanetmsx')
    unloadlibrary('epanetmsx');
else
    errstring =['Library ', 'epanetmsx', '.dll was not loaded'];
    disp(errstring);
end;
end
function [errcode] = MSXsaveoutfile(outfname)
[errcode] = calllib('epanetmsx','MSXsaveoutfile',outfname);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsavemsxfile(msxname)
[errcode] = calllib('epanetmsx','MSXsavemsxfile',msxname);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetconstant(index, value)
[errcode]=calllib('epanetmsx','MSXsetconstant',index,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetparameter(type, index, param, value)
[errcode]=calllib('epanetmsx','MSXsetparameter',type,index,param,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetinitqual(type,index,species,value)
[errcode]=calllib('epanetmsx','MSXsetinitqual',type,index,species,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetpattern(index, factors, nfactors)
[errcode]=calllib('epanetmsx','MSXsetpattern',index,factors,nfactors);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetpatternvalue(pat, period, value)
[errcode]=calllib('epanetmsx','MSXsetpatternvalue',pat,period,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsolveQ()
[errcode]=calllib('epanetmsx','MSXsolveQ');
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsolveH()
[errcode]=calllib('epanetmsx','MSXsolveH');
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXaddpattern(patid)
[errcode]=calllib('epanetmsx','MSXaddpattern',patid);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXusehydfile(hydfname)
[errcode]=calllib('epanetmsx','MSXusehydfile',hydfname);
if errcode
    MSXerror(errcode);
end
end
function [errcode, t, tleft] = MSXstep()
t=int32(0);
tleft=int32(0);
[errcode,t,tleft]=calllib('epanetmsx','MSXstep',t,tleft);
t = double(t);
tleft = double(tleft);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXinit(flag)
[errcode]=calllib('epanetmsx','MSXinit',flag);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXreport()
[errcode] = calllib('epanetmsx','MSXreport');
if errcode
    MSXerror(errcode);
end
end
function [e, errmsg] = MSXgeterror(errcode)
errmsg = char(32*ones(1,80));
[e,errmsg] = calllib('epanetmsx','MSXgeterror',errcode,errmsg,80);
if e
    MSXerror(e);
end
end
function [errcode, value] = MSXgetqual(type, index, species)
value=0;
[errcode,value]=calllib('epanetmsx','MSXgetqual',type,index,species,value);
if errcode
    MSXerror(errcode);
end
end
function [errcode] = MSXsetsource(node,species,type,level,pat)
[errcode]=calllib('epanetmsx','MSXsetsource',node,species,type,level,pat);
if errcode
    MSXerror(errcode);
end
end
function value=InputFileInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.pathfile,'rt');
    if fid < 0
        disp(message)
        return
    end
    %Nodes
    value.JunctionsID={};
    value.ReservoirsID={};
    value.TanksID={};
    value.NodeReservoirIndex=0;
    value.NodeTankIndex=0;
    value.NodeJunctionIndex=0;
    value.NodeCount=0;
    value.sectTanks=0;
    value.sectJunctions=0;
    value.sectReservoirs=0;
    k=1;p=1;r=1;

    %Links
    value.PipesID={};
    value.PumpsID={};
    value.ValvesID={};
    value.LinkPipeIndex=0;
    value.LinkPumpIndex=0;
    value.LinkValveIndex=0;
    value.FromNode={};
    value.ToNode={};
    value.LinkCount=0;
    value.sectPipes=0;
    value.sectPumps=0;
    value.sectValves=0;
    sect=0; i=1;t=1;q=1;

    %Curves
    value.typeCurve=[];
    typecode=0;
    value.CurveID={};
    value.CurveX={};
    value.CurveY={};
    value.Clines={}; sectCurve=0; x=1; b=1;

    %Controls
    value.controlsInfo={};
    value.linksID={};
    value.nodesID={};
    value.sectControls=0; d=1;

    %Options
    value.LinkFlowUnits={};
    value.OptionsHeadloss={};
    value.sectOptions=0;
    value.SImetric=0;
    value.UScustomary=0;
    value.QualityUnits={}; g=1; f=1;

    while 1
        tline = fgetl(fid);
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
                value.sectJunctions=value.sectJunctions+1;
                continue;
                % [RESERVOIRS] section
            elseif strcmpi(tok(1:5),'[RESE')
                sect=2;
                value.sectReservoirs=value.sectReservoirs+1;
                continue;
                % [TANKS] section
            elseif strcmpi(tok(1:5),'[TANK')
                sect=3;
                value.sectTanks=value.sectTanks+1;
                continue;
                % [PIPES] section
            elseif strcmpi(tok(1:5),'[PIPE')
                sect=4;
                value.sectPipes=value.sectPipes+1;
                continue;
                % [PUMPS] section
            elseif strcmpi(tok(1:5),'[PUMP')
                sect=5;
                value.sectPumps=value.sectPumps+1;
                continue;
                % [VALVES] section
            elseif strcmpi(tok(1:5),'[VALV')
                sect=6;
                value.sectValves=value.sectValves+1;
                continue;
                % [CURVES] section
            elseif strcmpi(tok(1:5),'[CURV')
                sect=7;
                sectCurve= sectCurve+1;
                continue;
                % [CONTROLS] section
            elseif strcmpi(tok(1:5),'[CONT')
                sect=8;
                value.sectControls=value.sectControls+1;
                continue;
                % [OPTIONS] section
            elseif strcmpi(tok(1:5),'[OPTI')
                sect=9;
                value.sectOptions=value.sectOptions+1;
                value.LinksAll=[value.PipesID value.PumpsID value.ValvesID];
                value.NodesAll=[value.JunctionsID value.ReservoirsID value.TanksID];
                value.vx = NaN(value.NodeCount,1);
                value.vy = NaN(value.NodeCount,1);
                value.vertx = cell(value.LinkCount,1);
                value.verty = cell(value.LinkCount,1);
                nvert = zeros(value.LinkCount,1);
                continue;
                % [COORDINATES] section
            elseif strcmpi(tok(1:5),'[COOR')
                sect=10;
                continue;
                % [VERTICES] section
            elseif strcmpi(tok(1:5),'[VERT')
                sect=11;
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
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.JunctionsID{k}=atline{1};
            value.NodeJunctionIndex(k)=k;
            k=k+1;
            value.NodeCount=value.NodeCount+1;
        elseif sect==2
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.ReservoirsID{r}=atline{1};
            value.NodeReservoirIndex(r)=k;
            k=k+1;
            r=r+1;
            value.NodeCount=value.NodeCount+1;
        elseif sect==3
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.TanksID{p}=atline{1};
            value.NodeTankIndex(p)=k;
            k=k+1;
            p=p+1;
            value.NodeCount=value.NodeCount+1;

            % Links
        elseif sect==4
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.PipesID{t}=atline{1};
            value.LinkPipeIndex(t)=t;
            value.FromNode{t}=atline{2};
            value.ToNode{t}=atline{3};
            t=t+1;
            value.LinkCount= value.LinkCount+1;
        elseif sect==5
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.PumpsID{q}=atline{1};
            value.LinkPumpIndex(q)=t;
            value.FromNode{t}=atline{2};
            value.ToNode{t}=atline{3};
            t=t+1;
            q=q+1;
            value.LinkCount= value.LinkCount+1;
        elseif sect==6
            a = regexp(tline, '\s*','split');uu=1;
            for tt=1:length(a)
                if isempty(a{tt})
                    %skip
                elseif (a{tt}==';')
                    %skip
                else
                    atline{uu}=a{tt}; uu=uu+1;
                end
            end
            value.ValvesID{i}=atline{1};
            value.LinkValveIndex(i)=t;
            value.FromNode{t}=atline{2};
            value.ToNode{t}=atline{3};
            t=t+1;
            i=i+1;
            value.LinkCount= value.LinkCount+1;
            % Curves
        elseif sect==7
            ee=regexp(tline,'\w*EFFICIENCY*\w','match');
            nn=regexp(tline,'\w*VOLUME*\w','match');
            kk=regexp(tline,'\w*HEADLOSS*\w','match');

            if strcmp(ee,'EFFICIENCY'), typecode=1;   % EFFICIENCY
                value.Clines{b}=tline;b=b+1;continue;
            elseif strcmp(nn,'VOLUME'), typecode=2;   % VOLUME
                value.Clines{b}=tline;b=b+1;continue;
            elseif strcmp(kk,'HEADLOSS'), typecode=3; % HEADLOSS
                value.Clines{b}=tline;b=b+1;continue;
            elseif (~length(strcmp(nn,'VOLUME')) || ~length(strcmp(ee,'EFFICIENCY')) || ~length(strcmp(kk,'HEADLOSS'))) &&  (tok(1)==';'), typecode=0; % HEADLOSS
                value.Clines{b}=tline;b=b+1;continue;
            else
                value.typeCurve(x)=typecode;
            end
            a = textscan(tline,'%s %f %f');
            value.CurveID(x)=a{1};
            value.CurveX(x)=a(2);
            value.CurveY(x)=a(3);
            value.Clines{b}=tline;
            x=x+1;b=b+1;
            % Controls
        elseif sect==8
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
            value.controlsInfo{d}=atline;
            value.linksID{d}=atline{2};
            t = regexp(tline, '\w*TIME\w*','match');
            if length(t)==0
                value.nodesID{d}=atline{6};
            end
            d=d+1;
            % Options
        elseif sect==9
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
            if strcmp(upper(atline{1}),'UNITS')
                value.LinkFlowUnits{f}=atline{2};
            end
            if strcmp(upper(atline{1}),'HEADLOSS')
                value.OptionsHeadloss{g}=atline{2};
                g=g+1;
            end
            if strcmp(upper(atline{1}),'PRESSURE')
                value.NodePressureUnits=atline{2};
            end
            if strcmp(upper(atline{1}),'QUALITY')
                if ~strcmp(upper(atline{2}),'NONE')
                    if length(atline)==3
                        value.QualityUnits=atline{3};
                    end
                else
                    value.QualityUnits='NONE';
                end
            end
            f=f+1;
            % Coordinates
        elseif sect==10
            A = textscan(tline,'%s %f %f');
            % get the node index
            a=strcmp(A{1},value.NodesAll);
            index=strfind(a,1);
            if length(index)==0
                return;
            end
            value.vx(index) = A{2};
            value.vy(index) = A{3};
            % Vertices
        elseif sect==11
            A = textscan(tline,'%s %f %f');
            [errcode,index] = ENgetlinkindex(char(A{1}));
            if errcode ~=0
                return;
            end
            nvert(index) = nvert(index) + 1;
            value.vertx{index}(nvert(index)) = A{2};
            value.verty{index}(nvert(index)) = A{3};
        end
    end

    %     US Customary - SI metric
    switch char(value.LinkFlowUnits)
        case 'CFS'
            value.UScustomary=1;
        case 'GPM'
            value.UScustomary=1;
        case 'MGD'
            value.UScustomary=1;
        case 'IMGD'
            value.UScustomary=1;
        case 'AFD'
            value.UScustomary=1;
        case 'LPS'
            value.SImetric=1;
        case 'LPM'
            value.SImetric=1;
        case 'MLD'
            value.SImetric=1;
        case 'CMH'
            value.SImetric=1;
        case 'CMD'
            value.SImetric=1;
    end

    if value.UScustomary==1;
        value.PatternDemandsUnits=value.LinkFlowUnits;
        value.LinkPipeDiameterUnits='inches';
        value.NodeTankDiameterUnits='feet';
        value.EnergyEfficiencyUnits='percent';
        value.NodeElevationUnits='feet';
        value.NodeEmitterCoefficientUnits='flow units @ 1 psi drop';
        value.EnergyUnits='kwatt-hours';
        value.LinkFrictionFactorUnits='unitless';
        value.NodeHeadUnits='feet';
        value.LinkLengthUnits='feet';
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
        value.PatternDemandsUnits=value.LinkFlowUnits;
        value.LinkPipeDiameterUnits='millimeters';
        value.NodeTankDiameterUnits='meters';
        value.EnergyEfficiencyUnits='percent';
        value.NodeElevationUnits='meters';
        value.NodeEmitterCoefficientUnits='flow units @ 1 meter drop';
        value.EnergyUnits='kwatt-hours';
        value.LinkFrictionFactorUnits='unitless';
        value.NodeHeadUnits='meters';
        value.LinkLengthUnits='meters';
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
function [curvesID,X,Y,Clines,typeCurve,sectCurve]=CurveInfo(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end

typeCurve=[];
typecode=0;
curvesID={};
X={};
Y={};
Clines={};
sect=0; i=1; u=1; sectCurve=0;
% Read each line from msx file.
while 1
    tline = fgetl(fid);
    a=regexp(tline,'\s*','split');
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
        %             typecode=0;
        if (tok(1) == ';'), continue, end
    end
    
    if (tok(1) == '[')
        % [CURVES] section
        if strcmpi(tok(1:5),'[CURV')
            sect = 1;
            sectCurve= sectCurve+1;
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
        
        if strcmp(ee,'EFFICIENCY'), typecode=1;   % EFFICIENCY
            Clines{u}=tline;u=u+1;continue;
        elseif strcmp(nn,'VOLUME'), typecode=2;   % VOLUME
            Clines{u}=tline;u=u+1;continue;
        elseif strcmp(kk,'HEADLOSS'), typecode=3; % HEADLOSS
            Clines{u}=tline;u=u+1;continue;
        elseif (~length(strcmp(nn,'VOLUME')) || ~length(strcmp(ee,'EFFICIENCY')) || ~length(strcmp(kk,'HEADLOSS'))) &&  (tok(1)==';'), typecode=0; % HEADLOSS
            Clines{u}=tline;u=u+1;continue;
        else
            typeCurve(i)=typecode;
        end
        a = textscan(tline,'%s %f %f');
        curvesID(i)=a{1};
        X(i)=a(2);
        Y(i)=a(3);
        Clines{u}=tline;
        i=i+1;u=u+1;
    end
end
end
function value=LinksInfo(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end
value.PipesID={};
value.PumpsID={};
value.ValvesID={};
value.LinkPipeIndex=0;
value.LinkPumpIndex=0;
value.LinkValveIndex=0;
value.FromNode={};
value.ToNode={};
value.LinkCount=0;
value.sectPipes=0;
value.sectPumps=0;
value.sectValves=0;
sect=0; i=1;t=1;q=1;
while 1
    tline = fgetl(fid);
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
            value.sectPipes=value.sectPipes+1;
            continue;
            % [PUMPS] section
        elseif strcmpi(tok(1:5),'[PUMP')
            sect=2;
            value.sectPumps=value.sectPumps+1;
            continue;
            % [VALVES] section
        elseif strcmpi(tok(1:5),'[VALV')
            sect=3;
            value.sectValves=value.sectValves+1;
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
        
        % Nodes
    elseif sect == 1
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.PipesID{t}=atline{1};
        value.LinkPipeIndex(t)=t;
        value.FromNode{t}=atline{2};
        value.ToNode{t}=atline{3};
        t=t+1;
    elseif sect == 2
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.PumpsID{q}=atline{1};
        value.LinkPumpIndex(q)=t;
        value.FromNode{t}=atline{2};
        value.ToNode{t}=atline{3};
        t=t+1;
        q=q+1;
    elseif sect == 3
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.ValvesID{i}=atline{1};
        value.LinkValveIndex(i)=t;
        value.FromNode{t}=atline{2};
        value.ToNode{t}=atline{3};
        t=t+1;
        i=i+1;
    end
    value.LinkCount= value.LinkCount+1;
end
value.LinksAll=[value.PipesID value.PumpsID value.ValvesID];
end
function value=NodesInfo(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end

value.JunctionsID={};
value.ReservoirsID={};
value.TanksID={};
value.NodeReservoirIndex=0;
value.NodeTankIndex=0;
value.NodeJunctionIndex=0;
value.NodeCount=0;

value.sectTanks=0;
value.sectJunctions=0;
value.sectReservoirs=0;
sect=0; t=1;i=1;q=1;
while 1
    tline = fgetl(fid);
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
            value.sectJunctions=value.sectJunctions+1;
            continue;
            % [RESERVOIRS] section
        elseif strcmpi(tok(1:5),'[RESE')
            sect=2;
            value.sectReservoirs=value.sectReservoirs+1;
            continue;
            % [TANKS] section
        elseif strcmpi(tok(1:5),'[TANK')
            sect=3;
            value.sectTanks=value.sectTanks+1;
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
        
        % Nodes
    elseif sect == 1
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.JunctionsID{i}=atline{1};
        value.NodeJunctionIndex(i)=i;
        i=i+1;
        
    elseif sect == 2
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.ReservoirsID{t}=atline{1};
        value.NodeReservoirIndex(t)=i;
        i=i+1;
        t=t+1;
        
    elseif sect == 3
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.TanksID{q}=atline{1};
        value.NodeTankIndex(q)=i;
        i=i+1;
        q=q+1;
    end
    value.NodeCount=value.NodeCount+1;
end

value.NodesAll=[value.JunctionsID value.ReservoirsID value.TanksID];
end
function value=ControlsInfo(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end
value.controlsInfo={};
value.linksID={};
value.nodesID={};
value.sectControls=0;

sect=0;i=1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    % Get first token in the line
    tok = strtok(tline);
    % Skip blank Clines and comments
    if isempty(tok), continue, end
    if (tok(1) == ';'), continue, end
    if (tok(1) == '[')
        % [CONTROLS] section
        if strcmpi(tok(1:5),'[CONT')
            sect=1;
            value.sectControls=value.sectControls+1;
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
        
        % Controls
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
        value.controlsInfo{i}=atline;
        value.linksID{i}=atline{2};
        t = regexp(tline, '\w*TIME\w*','match');
        if length(t)==0
            value.nodesID{i}=atline{6};
        end
        i=i+1;
    end
    
end
end
function value=FlowUnitsHeadlossFormula(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end
value.LinkFlowUnits={};
value.OptionsHeadloss={};
value.sectOptions=0;
value.SImetric=0;
value.UScustomary=0;
value.QualityUnits={};

sect=0;i=1;t=1;
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
            value.sectOptions=value.sectOptions+1;
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
        if strcmp(upper(atline{1}),'UNITS')
            value.LinkFlowUnits{i}=atline{2};
        end
        if strcmp(upper(atline{1}),'HEADLOSS')
            value.OptionsHeadloss{t}=atline{2};
            t=t+1;
        end
        if strcmp(upper(atline{1}),'PRESSURE')
            value.NodePressureUnits=atline{2};
        end
        if strcmp(upper(atline{1}),'QUALITY')
            if ~strcmp(upper(atline{2}),'NONE')
                if length(atline)==3
                    value.QualityUnits=atline{3};
                end
            else
                value.QualityUnits='NONE';
            end
        end
        i=i+1;
    end
    
end

%     US Customary - SI metric
switch char(value.LinkFlowUnits)
    case 'CFS'
        value.UScustomary=1;
    case 'GPM'
        value.UScustomary=1;
    case 'MGD'
        value.UScustomary=1;
    case 'IMGD'
        value.UScustomary=1;
    case 'AFD'
        value.UScustomary=1;
    case 'LPS'
        value.SImetric=1;
    case 'LPM'
        value.SImetric=1;
    case 'MLD'
        value.SImetric=1;
    case 'CMH'
        value.SImetric=1;
    case 'CMD'
        value.SImetric=1;
end

if value.UScustomary==1;
    value.PatternDemandsUnits=value.LinkFlowUnits;
    value.LinkPipeDiameterUnits='inches';
    value.NodeTankDiameterUnits='feet';
    value.EnergyEfficiencyUnits='percent';
    value.NodeElevationUnits='feet';
    value.NodeEmitterCoefficientUnits='flow units @ 1 psi drop';
    value.EnergyUnits='kwatt-hours';
    value.LinkFrictionFactorUnits='unitless';
    value.NodeHeadUnits='feet';
    value.LinkLengthUnits='feet';
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
    value.PatternDemandsUnits=value.LinkFlowUnits;
    value.LinkPipeDiameterUnits='millimeters';
    value.NodeTankDiameterUnits='meters';
    value.EnergyEfficiencyUnits='percent';
    value.NodeElevationUnits='meters';
    value.NodeEmitterCoefficientUnits='flow units @ 1 meter drop';
    value.EnergyUnits='kwatt-hours';
    value.LinkFrictionFactorUnits='unitless';
    value.NodeHeadUnits='meters';
    value.LinkLengthUnits='meters';
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
function value=PumpInfo(obj)
% Open epanet input file
[fid,message] = fopen(obj.pathfile,'rt');
if fid < 0
    disp(message)
    return
end

value.PumpsID={};
value.fromNodePumps={};
value.toNodePumps={};
value.PumpCurveID={};
sect=0; i=1; value.sectPump=0;
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
        % [PUMPS] section
        if strcmpi(tok(1:5),'[PUMP')
            sect=1;
            value.sectPump=value.sectPump+1;
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
        
        % Pumps
    elseif sect == 1
        a = regexp(tline, '\s*','split');uu=1;
        for tt=1:length(a)
            if isempty(a{tt})
                %skip
            elseif (a{tt}==';')
                %skip
            else
                atline{uu}=a{tt}; uu=uu+1;
            end
        end
        value.PumpsID{i}=atline{1};
        value.fromNodePumps{i}=atline{2};
        value.toNodePumps{i}=atline{3};
        value.PumpCurveID{i}=atline{5};
        i=i+1;
    end
end
end
function remAddCurve(obj,CurvesID,Cline,varargin)

% Open and read inpname
% Read all file and save in variable info
info = readAllFile(obj);
fid2 = fopen(obj.pathfile,'w');

sect=0;
t=1;
nn=0; sps = {'                  '};
while t<length(info)+1
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
        % skip
    elseif isempty(c)
        % skip
    else
        u=1;
        while u < length(a)+1
            
            if strcmp(a{u},'[CURVES]')
                fprintf(fid2,'[CURVES]');
                nn=0;
                sect=1;
                fprintf(fid2,'\n');break;
            elseif sect==0
                fprintf(fid2,'%s%s',c,sps{:});
                fprintf(fid2,'\n');break;
            end
            
            % section [CURVES]
            if (sect==1) && (nn==0)
                while ~strcmp(a{u},'[CONTROLS]')
                    t=t+1;
                    c = info{t};
                    a = regexp(c, '\s*','split');
                end
                for i=1:length(Cline)
                    fprintf(fid2,'%s\n',char(Cline{i}),sps{:});
                end
                fprintf(fid2,'%s%s',a{u},sps{:});nn=1;sect=0;
                fprintf(fid2,'\n');break;
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps{:});
                end
            end
            u=u+1;
        end
    end
    t=t+1;
    %         fprintf(fid2,'\n');
end
fclose all;
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
function info = readAllFile(obj)
fid = fopen(obj.pathfile,'r+');
t=1;
while ~feof(fid)
    tline=fgetl(fid);
    info{t} = tline;
    t = t+1;
end
end
function Options(obj,newinpname,newFlowUnits,headloss,varargin)
% Notes: Flow units codes are as follows: CFS	cubic feet per second
% GPM	gallons per minute MGD	million gallons per day IMGD
% Imperial mgd AFD	acre-feet per day LPS	liters per second LPM
% liters per minute MLD	million liters per day CMH	cubic meters per
% hour CMD	cubic meters per day
value=obj.getFlowUnitsHeadlossFormula;
previousFlowUnits=value.LinkFlowUnits;
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
if newUScustomary==value.UScustomary || newSImetric==value.SImetric
    changes=0; newUScustomary=0;
    newSImetric=0;% feet to feet, meter to meter
elseif value.UScustomary==1 && newUScustomary==0
    changes=1; % feet to meter or cubic feet to cubic meter
elseif value.UScustomary==0 && newUScustomary==1
    changes=2; % meter to feet or cubic meter to cubic feet
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Units=newUScustomary+newSImetric;
variables=who;nheadl=0;
if ~sum(strcmp('headloss',variables))
    headloss=value.OptionsHeadloss;
    nheadl=1;
end
nodes = obj.getNodesInfo;
links = obj.getLinksInfo;
controls=obj.getControlsInfo;
curves=obj.getCurveInfo;

% Open and read inpname
% Read all file and save in variable info
info = readAllFile(obj);
fid2 = fopen(obj.pathfile,'w');

sections=[0 0 0 0 0 0 0 0 0 0];

nn=0;pp=1;sps={'                  '};
for t = 1:length(info)
    c = info{t};
    a = regexp(c, '\s*','split');
    if isempty(a)
        % skip
    elseif isempty(c)
        % skip
    else
        u=1;
        while u < length(a)+1
            if strcmp(a{u},'[JUNCTIONS]') && Units
                fprintf(fid2,'[JUNCTIONS]');
                sections=[1 0 0 0 0 0 0 0 0 0];
                break;
            elseif strcmp(a{u},'[RESERVOIRS]') && Units
                fprintf(fid2,'[RESERVOIRS]');
                nn=0;pp=1;
                sections=[0 1 0 0 0 0 0 0 0 0];
                break;
            elseif strcmp(a{u},'[TANKS]') && Units
                fprintf(fid2,'[TANKS]');
                nn=0;pp=1;
                sections=[0 0 1 0 0 0 0 0 0 0];
                break;
            elseif strcmp(a{u},'[PIPES]') && Units
                fprintf(fid2,'[PIPES]');
                nn=0;pp=1;
                sections=[0 0 0 1 0 0 0 0 0 0];
                break;
            elseif strcmp(a{u},'[PUMPS]') && Units
                fprintf(fid2,'[PUMPS]');
                nn=0;pp=1;
                sections=[0 0 0 0 1 0 0 0 0 0];
                break;
            elseif strcmp(a{u},'[VALVES]') && Units
                fprintf(fid2,'[VALVES]');
                nn=0;pp=1;
                sections=[0 0 0 0 0 1 0 0 0 0];
                break;
            elseif strcmp(a{u},'[DEMANDS]') && ((Units || ~changes) && nheadl)
                fprintf(fid2,'[DEMANDS]');
                nn=0;pp=1;
                sections=[0 0 0 0 0 0 1 0 0 0];
                break;
            elseif strcmp(a{u},'[CURVES]') && ((Units || ~changes) && nheadl)
                fprintf(fid2,'[CURVES]');
                nn=0;pp=1;ww=1;
                sections=[0 0 0 0 0 0 0 1 0 0];
                break;
            elseif strcmp(a{u},'[CONTROLS]') && Units
                fprintf(fid2,'[CONTROLS]');
                nn=0;pp=1;
                sections=[0 0 0 0 0 0 0 0 1 0];
                break;
            elseif strcmp(a{u},'[OPTIONS]')
                fprintf(fid2,'[OPTIONS]');
                sections=[0 0 0 0 0 0 0 0 0 1];nn=0;
                break;
            end
            
            if length(strfind(sections,1))==0
                sect=0;
            else
                sect=strfind(sections,1);
            end
            
            % section [JUNCTIONS]
            if (sect==1) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(nodes.JunctionsID))+1
                        if strcmp(a{mm},nodes.JunctionsID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                            end
                        end
                    else
                        nn=1;
                    end
                end
                break;
                % section [RESERVOIRS]
            elseif (sect==2) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(nodes.ReservoirsID))+1
                        if strcmp(a{mm},nodes.ReservoirsID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                            end
                        end
                    else
                        nn=1;
                    end
                end
                break;
                % section [TANKS]
            elseif (sect==3) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(nodes.TanksID))+1
                        if strcmp(a{mm},nodes.TanksID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+4})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+5})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+6})*0.02831685),sps{:});
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+4})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+5})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+6})*35.3147),sps{:});
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                % section [PIPES]
            elseif (sect==4) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(links.PipesID))+1
                        if strcmp(a{mm},links.PipesID{pp})
                            pp=pp+1;
                            for mm=mm:(mm+2)
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            end
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*25.4),sps{:});
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.03937007874),sps{:});
                            end
                            for mm=(mm+3):(mm+5)
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
            % section [PUMPS]
            elseif (sect==5) && (nn==0)
                mm=1; 
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(links.PumpsID))+1
                        if strcmp(a{mm},links.PumpsID{pp})
                            pp=pp+1;
                            for mm=mm:length(a(1:end-1))
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            end
                            power=regexp(c,'\w*POWER*\w','match');
                            if strcmpi(upper(power),'POWER')
                                if changes==1 
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.745699882507324),sps{:});
                                elseif changes==2
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})/0.745699882507324),sps{:});
                                end
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                
            % section [VALVES]
            elseif (sect==6) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(links.ValvesID))+1
                        if strcmp(a{mm},links.ValvesID{pp})
                            pp=pp+1;
                            for mm=mm:(mm+2)
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            end
                            prv=regexp(c,'\w*PRV*\w','match');
                            psv=regexp(c,'\w*PSV*\w','match');
                            pbv=regexp(c,'\w*PBV*\w','match');
                            fcv=regexp(c,'\w*FCV*\w','match');
                            if changes==1
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*25.4),sps{:});
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps{:});
                                if strcmpi(upper(prv),'PRV') || strcmpi(upper(psv),'PRV') || strcmpi(upper(pbv),'PRV')
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})*0.3048),sps{:});
                                elseif strcmpi(upper(fcv),'FCV')
                                    setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                end
                            elseif changes==2
                                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.03937007874),sps{:});
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps{:});
                                if strcmpi(upper(prv),'PRV') || strcmpi(upper(psv),'PRV') || strcmpi(upper(pbv),'PRV')
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+3})/0.3048),sps{:});
                                elseif strcmpi(upper(fcv),'FCV')
                                    setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                end
                            end
                            for mm=(mm+4):length(a)
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                
                % section [DEMANDS]
            elseif (sect==7) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(char(nodes.JunctionsID))+1
                        if strcmp(a{mm},nodes.JunctionsID{pp})
                            pp=pp+1;
                            fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                            setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            fprintf(fid2,'%s%s',char(a{mm+2}),sps{:});
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                % section [CURVES]
            elseif (sect==8) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(curves.Clines)+1 && ~isempty(char(a)) % PUMP % EFFICIENCY % VOLUME
                        if ww<length(curves.CurvesID)+1
                            if curves.typeCurve(ww)==0
                                if strcmp(a{1},';PUMP:')
                                    fprintf(fid2,c);pp=pp+1;break;
                                end
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                if changes==1
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps{:});
                                elseif changes==2
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*3.281),sps{:});
                                else
                                    fprintf(fid2,'%s%s',char(a{mm+2}),sps{:});
                                end
                            elseif curves.typeCurve(ww)==1
                                ee=regexp(c,'\w*EFFICIENCY*\w','match');
                                if length(strcmp(ee,'EFFICIENCY'))
                                    fprintf(fid2,c);pp=pp+1;break;
                                end
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                                fprintf(fid2,'%s%s',char(a{mm+2}),sps{:});
                            elseif curves.typeCurve(ww)==2
                                gg=regexp(c,'\w*VOLUME*\w','match');
                                if length(strcmp(gg,'VOLUME'))
                                    fprintf(fid2,c);pp=pp+1;break;
                                end
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                if changes==1
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.831685e-02),sps{:});
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})*0.3048),sps{:});
                                else
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+2})),sps{:});
                                end
                            elseif curves.typeCurve(ww)==3 % HEADLOSS
                                kk=regexp(c,'\w*HEADLOSS*\w','match');
                                if length(strcmp(kk,'HEADLOSS'))
                                    fprintf(fid2,c);pp=pp+1;break;
                                end
                                fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                if changes==1
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                                elseif changes==2
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                                else
                                    fprintf(fid2,'%s%s',char(a{mm+1}),sps{:});
                                end
                                mm=mm+1;
                                setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
                            end
                            ww=ww+1;
                        end
                        pp=pp+1;
                    else
                        if ~(ww<length(curves.CurvesID)+1), nn=1; end
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                
                % section [CONTROLS]
            elseif (sect==9) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if pp<length(controls.controlsInfo)+1
                        pp=pp+1;
                        if length(a)>4
                            if strcmp(upper(a{mm+6}),'BELOW') || strcmp(upper(a{mm+6}),'ABOVE')
                                for mm=mm:(mm+6)
                                    fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                end
                                if changes==1
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
                                elseif changes==2
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:});
                                end
                            else
                                for mm=mm:length(a)
                                    fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                end
                            end
                        end
                    else
                        nn=1;
                        fprintf(fid2,'%s%s',char(a{1}),sps{:});
                    end
                end
                break;
                
                % section [OPTIONS]
            elseif (sect==10) && (nn==0)
                mm=1;
                if mm < length(a)+1
                    if isempty(a{mm})
                        % skip
                        mm=mm+1;
                    end
                    if strcmp(upper(a{mm}),'UNITS') && (Units || ~changes)
                        fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                        if nheadl
                            fprintf(fid2,'%s%s',char(newFlowUnits),sps{:});nn=1;
                        else
                            fprintf(fid2,'%s%s',char(previousFlowUnits),sps{:});
                        end
                    elseif strcmp(upper(a{mm}),'HEADLOSS')
                        fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                        fprintf(fid2,'%s%s',char(headloss),sps{:});
                        nn=1;
                    else
                        fprintf(fid2,c);
                    end
                end
                break;
                
                
            elseif isempty(a{u}) && nn==0
            else
                if isempty(a{u}) && nn==1
                else
                    fprintf(fid2,'%s%s',a{u},sps{:});
                end
            end
            u=u+1;
        end
    end
    fprintf(fid2,'\n');
end
fclose all;
copyfile([pwd,'\RESULTS\','temp.inp'],[pwd,'\NETWORKS\',newinpname]);
% obj=epanet('temp2.inp');
% obj.saveInputFile([pwd,'\NETWORKS\',newinpname])
end
function setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
if strcmp(previousFlowUnits,'GPM')
    switch newFlowUnits %(GPM)
        case 'CFS' 
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00222816399286988),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00144),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00119905),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.004419191),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0630902),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.785412),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.005450993),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2271247),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*5.450993),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'CFS')
    switch newFlowUnits
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*448.8312),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.6463169),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.5381711),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.983471),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*28.31685),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1899.011),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.446576),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*101.9406),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2446.576),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'MGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.547229),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*694.4445),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.8326738),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.068883),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*43.81264),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2628.758),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.785412),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*157.7255),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3785.412),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'IMGD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.858145),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*833.9936),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.200951),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.685577),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*52.61681),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3157.008),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.546092),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*189.4205),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4546.092),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'AFD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.5041667),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*226.2857),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3258514),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.271328),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*14.27641),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*856.5846),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.233482),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*51.39508),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1233.482),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'LPS')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.03531466),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*15.85032),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.02282446),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01900533),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.07004562),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*60),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0864),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.6),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*86.4),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'LPM')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0005885777),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.264172),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0003804078),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0003167556),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0011674272),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01666667),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00144),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.06),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.44),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'MLD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.4087345),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*183.4528),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.264172),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2199692),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.8107132),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*11.57407),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*694.4445),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*41.66667),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1000),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'CMH')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.009809635),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.402868),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.006340129),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00527926),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01945712),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.2777778),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*16.66667),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.024),sps{:});
        case 'CMD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*24),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
elseif strcmp(previousFlowUnits,'CMD')
    switch newFlowUnits
        case 'CFS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0004087345),sps{:});
        case 'GPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.1834528),sps{:});
        case 'MGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.000264172),sps{:});
        case 'IMGD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0002199692),sps{:});
        case 'AFD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0008107132),sps{:});
        case 'LPS'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.01157407),sps{:});
        case 'LPM'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.6944444),sps{:});
        case 'MLD'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.001),sps{:});
        case 'CMH'
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.04166667),sps{:});
        otherwise
            fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});
    end
end
end
function value=MsxTimeStep(msxname)
% Open epanet input file
[fid,message] = fopen(msxname,'rt');
if fid < 0
    disp(message)
    return
end

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
            value=str2num(atline{2});
        end
    end
    
end
end
