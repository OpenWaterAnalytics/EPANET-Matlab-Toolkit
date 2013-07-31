%{
 Copyright 2013 KIOS Research Center for Intelligent Systems and Networks, University of Cyprus (www.kios.org.cy)

 Licensed under the EUPL, Version 1.1 or ï¿½ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence");
 You may not use this work except in compliance with the Licence.
 You may obtain a copy of the Licence at:

 http://ec.europa.eu/idabc/eupl

 Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the Licence for the specific language governing permissions and limitations under the Licence.
%}

classdef Epanet <handle 
    % The EPANET class v.1.04
    properties %(Attributes)
        InputFile;
        PathFile;
        CountNodes;
        CountTanksReservoirs;
        CountLinks;
        CountPatterns;
        CountCurves;
        CountControls;
        CountReservoirs;
        CountJunctions;
        CountTanks;
        CountPipes;
        CountPumps;
        CountValves;
        ControlType;
        ControlLink;
        ControlSeetting;
        ControlNodeIndex;
        ControlTypeIndex;
        ControlLevel;
        ControlAll;
        CoordinatesXY;
        CurvesNameID;
        CurvesX;
        CurvesY;
        flowUnits;
        UnitsType;
        LinkNameID;
        LinkIndex;
        LinkTypeIndex;
        LinkPipeNameID;
        LinkPipeIndex;
        LinkPumpNameID;
        LinkPumpIndex;
        LinkTypeID;
        LinkDiameter;
        LinkRoughness;
        LinkMinorLossCoef;
        LinkInitialStatus;
        LinkInitialSetting;
        LinkBulkReactionCoeff;
        LinkWallReactionCoeff;
        LinkLength;
        NodeNameID;
        NodeIndex;
        NodeTypeIndex;
        NodeJunctionNameID;
        NodeJunctionIndex;
        NodeReservoirNameID;
        NodeReservoirIndex;
        NodeTypeID;
        NodeElevations;
        NodeBaseDemands;
        NodeDemandPatternIndex;
        NodeEmitterCoeff;
        NodeInitialQuality;
        NodeSourceQuality;
        NodeSourcePatternIndex;
        NodeSourceTypeIndex;
        NodesConnectingLinksIndex;
        NodesConnectingLinksID;
        TankNameID;
        TankIndex;
        TankInitialLevel;
        TankInitialWaterVolume;
        TankMixingModelCode;
        TankMixZoneVolume;
        TankDiameter;
        TankMinimumWaterVolume;
        TankVolumeCurveIndex;
        TankMinimumWaterLevel;
        TankMaximumWaterLevel;
        TankMinimumFraction;
        TankBulkReactionCoeff;
        TankMixingModelType;
        OptionsTrials;
        OptionsAccuracy;
        OptionsTolerance;
        OptionsEmmiterExponent;
        OptionsDemandMult;
        PatternID;
        PatternIndex;
        PatternLengths;
        QualityCode;
        QualityType;
        QualityTraceNodeIndex;
        TimeSimulationDuration;
        TimeHydraulicStep;
        TimeQualityStep;
        TimePatternStep;
        TimePatternStart;
        TimeReportingStart;
        TimeReportingStep;
        TimeRuleControlStep;
        TimeStatisticsIndex;
        TimeStatisticsType;
        TimeReportingPeriods;
        Version;
        ValveNameID;
        ValveIndex;
        errorCode;
        TYPECONTROL={'LOWLEVEL','HILEVEL', 'TIMER', 'TIMEOFDAY'};
        TYPEUNITS={'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'};
        TYPELINK={'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV', 'VALVE'};
        TYPENODE={'JUNCTION','RESERVOIR', 'TANK'};
        TYPEMIXMODEL={'MIX1','MIX2', 'FIFO','LIFO'};
        TYPEQUALITY={'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'};
        TYPESTATS={'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'};
        TYPESOURCE={'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'};
        TYPEREPORT={'YES','NO','FULL'};
        TYPEPUMP={'CONSTANT_HORSEPOWER', 'POWER_FUNCTION', 'CUSTOM'};

        
        %%%%%%%%%%%%%%%%%% EPANET - MSX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MsxFile;
        MsxPathFile;
        CountSpeciesMsx;
        CountConstantsMsx;
        CountParametersMsx;
        CountPatternsMsx;
       
        TermsFormula;
        PipesFormula;
        TanksFormula;
        
        ConstantNameIDMsx;
        ConstantValueMsx;
        
        ParameterNameIDMsx;
        ParametersIndexMsx;
        TankParameterValueMsx;
        PipeParameterValueMsx;

        PatternIDMsx;
        PatternIndexMsx;
        PatternLengthsMsx;
                
        SpeciesNameIDMsx;
        SpeciesIndexMsx;
        SpeciesTypeMsx;
        SpeciesUnitsMsx;
        SpeciesAtolMsx;
        SpeciesRtolMsx;

        NodeInitqualValueMsx;
        LinkInitqualValueMsx;
        
        SourceTypeMsx;
        SourceLevelMsx;
        SourceNodeIDMsx;
        SourcePatternIndexMsx; 
        SourcePatternIDMsx;
%         SourceSpeciesNameIDMsx;
        SourceAllMsx;
        
        %% Units %%
        FlowUnits;
        ConcentrationUnits;
        Headloss;
        DemandsUnits;
        DiameterPipesUnits;
        DiameterTanksUnits;
        EfficiencyUnits;
        ElevationUnits;
        EmitterCoeffUnits;
        EnergyUnits;
        FrictionFactorUnits;
        HeadUnits;
        LengthUnits;
        MinorLossCoeffUnits;
        PowerUnits;
        PressureUnits;
        ReactionCoeffBulkUnits;
        ReactionCoeffWallUnits;
        RoughnessCoeffUnits;
        SourceMassInjectionUnits;
        VelocityUnits;
        VolumeUnits;
        WaterAgeUnits;
    end
    
    methods
        function obj =  Epanet(pathfile,varargin)
            path(path,genpath(pwd));

            if nargin==2
                inpfile=varargin{1};
%                 obj.PathFile=pathfile;
            elseif nargin==1
                inpfile=pathfile;
%                 obj.PathFile=which(pathfile);
            end
                
            %%%%%%%%%%%%%%%%%
            if ~strcmp(inpfile,'tempInpFile.inp')
                copyfile(inpfile,[pwd,'\Results\tempInpFile.inp']);
            end
            tmp=1;  
            save([pwd,'\Results\','tmpInp'],'tmp','-mat')
            %%%%%%%%%%%%%%%%%
            
            %ENMatlabSetup  
            [obj.errorCode]=ENMatlabSetup('epanet2','epanet2.h');
            
            obj.PathFile=which('tempInpFile.inp');
            
            %ENopen  
            [obj.errorCode] = ENopen(obj.PathFile, strcat('tempInpFile.inp','.rpt'), strcat('tempInpFile.inp','.out'));
            obj.InputFile=inpfile;
            
            % ENsaveinpfile
            obj.saveInputFile('tempInpFile.inp'); % OK %
            movefile('tempInpFile.inp',[pwd,'\Results']);
        
            %ENgetcount  
            [obj.errorCode, obj.CountNodes] = ENgetcount(0);
            [obj.errorCode, obj.CountTanksReservoirs] = ENgetcount(1);
            [obj.errorCode, obj.CountLinks] = ENgetcount(2);
            [obj.errorCode, obj.CountPatterns] = ENgetcount(3);
            [obj.errorCode, obj.CountCurves] = ENgetcount(4);
            [obj.errorCode, obj.CountControls] = ENgetcount(5);
            tmpNodeTypes=obj.getNodeType;
            tmpLinkTypes=obj.getLinkType;
            
            obj.NodeReservoirIndex = find(strcmp(tmpNodeTypes,'RESERVOIR'));
            obj.TankIndex = find(strcmp(tmpNodeTypes,'TANK'));
            obj.NodeJunctionIndex = find(strcmp(tmpNodeTypes,'JUNCTION'));
            obj.LinkPipeIndex = find(strcmp(tmpLinkTypes,'PIPE'));
            obj.LinkPumpIndex = find(strcmp(tmpLinkTypes,'PUMP'));
            obj.ValveIndex = find(strcmp(tmpLinkTypes,'VALVE'));
            
            obj.CountReservoirs=sum(strcmp(tmpNodeTypes,'RESERVOIR'));
            obj.CountTanks=sum(strcmp(tmpNodeTypes,'TANK'));
            obj.CountJunctions=sum(strcmp(tmpNodeTypes,'JUNCTION'));
            obj.CountPipes=sum(strcmp(tmpLinkTypes,'PIPE'))+sum(strcmp(tmpLinkTypes,'CVPIPE'));
            obj.CountPumps=sum(strcmp(tmpLinkTypes,'PUMP'));
            obj.CountValves=obj.CountLinks - (obj.CountPipes + obj.CountPumps);
         
            %ENgetcontrol
            if obj.CountControls
                obj.ControlType{obj.CountControls}=[];
                obj.ControlTypeIndex(obj.CountControls)=NaN;
                obj.ControlLink(obj.CountControls)=NaN;
                obj.ControlSeetting(obj.CountControls)=NaN;
                obj.ControlNodeIndex(obj.CountControls)=NaN;
                obj.ControlLevel(obj.CountControls)=NaN;
            end
            %tmpControlType={'LOWLEVEL','HILEVEL', 'TIMER', 'TIMEOFDAY'};
            for i=1:obj.CountControls
                [obj.errorCode, obj.ControlTypeIndex(i),obj.ControlLink(i),obj.ControlSeetting(i),obj.ControlNodeIndex(i),obj.ControlLevel(i)] = ENgetcontrol(i);
                obj.ControlType(i)=obj.TYPECONTROL(obj.ControlTypeIndex(i)+1);
            end
            obj.ControlAll={obj.ControlType,obj.ControlTypeIndex,obj.ControlLink,obj.ControlSeetting,obj.ControlNodeIndex,obj.ControlLevel};
            
            %ENgetflowunits
            [obj.errorCode, obj.flowUnits] = ENgetflowunits();
            %tmpUnits={'CFS', 'GPM', 'MGD', 'IMGD', 'AFD', 'LPS', 'LPM', 'MLD', 'CMH', 'CMD'};
            obj.UnitsType=obj.TYPEUNITS(obj.flowUnits+1);
            
            
            %tmpLinkType={'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'};
            obj.LinkTypeID={};
            obj.NodesConnectingLinksID={};
            for i=1:obj.CountLinks
                %ENgetlinkid
                [obj.errorCode, obj.LinkNameID{i}] = ENgetlinkid(i);
                %ENgetlinkindex
                [obj.errorCode, obj.LinkIndex(i)] = ENgetlinkindex(obj.LinkNameID{i});
                %ENgetlinknodes
                [obj.errorCode,linkFromNode,linkToNode] = ENgetlinknodes(i);
                obj.NodesConnectingLinksIndex(i,:)= [linkFromNode,linkToNode];
                obj.NodesConnectingLinksID(i,1) = obj.getLinkID(linkFromNode);
                obj.NodesConnectingLinksID(i,2) = obj.getLinkID(linkToNode);
                %ENgetlinktype
                [obj.errorCode,obj.LinkTypeIndex(i)] = ENgetlinktype(i);
                obj.LinkTypeID(i)=obj.TYPELINK(obj.LinkTypeIndex(i)+1);
                %ENgetlinkvalue
                [obj.errorCode, obj.LinkDiameter(i)] = ENgetlinkvalue(i,0);
                [obj.errorCode, obj.LinkLength(i)] = ENgetlinkvalue(i,1);
                [obj.errorCode, obj.LinkRoughness(i)] = ENgetlinkvalue(i,2);
                [obj.errorCode, obj.LinkMinorLossCoef(i)] = ENgetlinkvalue(i,3);
                [obj.errorCode, obj.LinkInitialStatus(i)] = ENgetlinkvalue(i,4);
                [obj.errorCode, obj.LinkInitialSetting(i)] = ENgetlinkvalue(i,5);
                [obj.errorCode, obj.LinkBulkReactionCoeff(i)] = ENgetlinkvalue(i,6);
                [obj.errorCode, obj.LinkWallReactionCoeff(i)] = ENgetlinkvalue(i,7);
            end
            
            %tmpNodeType={'JUNCTION','RESERVOIR', 'TANK'};
            %tmpSourceType={'CONCEN','MASS', 'SETPOINT', 'FLOWPACED'};
            obj.NodeTypeID={};
            for i=1:obj.CountNodes
                %ENgetnodeid
                [obj.errorCode, obj.NodeNameID{i}] = ENgetnodeid(i);
                %ENgetnodeindex
                [obj.errorCode, obj.NodeIndex(i)] = ENgetnodeindex(obj.NodeNameID{i});
                %ENgetnodetype
                [obj.errorCode, obj.NodeTypeIndex(i)] = ENgetnodetype(i);
                obj.NodeTypeID(i)=obj.TYPENODE(obj.NodeTypeIndex(i)+1);
                %ENgetnodevalue
                [obj.errorCode, obj.NodeElevations(i)] = ENgetnodevalue(i, 0);
                [obj.errorCode, obj.NodeBaseDemands(i)] = ENgetnodevalue(i, 1);
                [obj.errorCode, obj.NodeDemandPatternIndex(i)] = ENgetnodevalue(i, 2);
                [obj.errorCode, obj.NodeEmitterCoeff(i)] = ENgetnodevalue(i, 3);
                [obj.errorCode, obj.NodeInitialQuality(i)] = ENgetnodevalue(i, 4);
                [obj.errorCode, obj.NodeSourceQuality(i)] = ENgetnodevalue(i, 5);
                [obj.errorCode, obj.NodeSourcePatternIndex(i)] = ENgetnodevalue(i, 6);
                [obj.errorCode, obj.NodeSourceTypeIndex(i)] = ENgetnodevalue(i, 7);
            end
                        
            obj.NodeReservoirNameID=obj.NodeNameID(obj.NodeReservoirIndex);
            obj.TankNameID=obj.NodeNameID(obj.TankIndex);
            obj.NodeJunctionNameID=obj.NodeNameID(obj.NodeJunctionIndex);
            obj.LinkPipeNameID=obj.LinkNameID(obj.LinkPipeIndex);            
            obj.LinkPumpNameID=obj.LinkNameID(obj.LinkPumpIndex);            
            obj.ValveNameID=obj.LinkNameID(obj.ValveIndex);   
            
            
            %ENgetnodevalue (for Tanks)
            obj.TankInitialLevel=nan(1,obj.CountNodes);
            obj.TankInitialWaterVolume=nan(1,obj.CountNodes);
            obj.TankMixingModelCode=nan(1,obj.CountNodes);
            obj.TankMixZoneVolume=nan(1,obj.CountNodes);
            obj.TankDiameter=nan(1,obj.CountNodes);
            obj.TankMinimumWaterVolume=nan(1,obj.CountNodes);
            obj.TankVolumeCurveIndex=nan(1,obj.CountNodes);
            obj.TankMinimumWaterLevel=nan(1,obj.CountNodes);
            obj.TankMaximumWaterLevel=nan(1,obj.CountNodes);
            obj.TankMinimumFraction=nan(1,obj.CountNodes);
            obj.TankBulkReactionCoeff=nan(1,obj.CountNodes);
            tmpTanks=find(strcmpi(obj.NodeTypeID,'TANK')==1);
            %tmpMixModel={'MIX1','MIX2', 'FIFO','LIFO'};
            obj.TankMixingModelType={};
            for i=tmpTanks
                [obj.errorCode, obj.TankInitialLevel(i)] = ENgetnodevalue(i, 8);
                [obj.errorCode, obj.TankInitialWaterVolume(i)] = ENgetnodevalue(i, 14);
                [obj.errorCode, obj.TankMixingModelCode(i)] = ENgetnodevalue(i, 15);
                obj.TankMixingModelType(i)=obj.TYPEMIXMODEL(obj.TankMixingModelCode(i)+1);
                [obj.errorCode, obj.TankMixZoneVolume(i)] = ENgetnodevalue(i, 16);
                [obj.errorCode, obj.TankDiameter(i)] = ENgetnodevalue(i, 17);
                [obj.errorCode, obj.TankMinimumWaterVolume(i)] = ENgetnodevalue(i, 18);
                [obj.errorCode, obj.TankVolumeCurveIndex(i)] = ENgetnodevalue(i, 19);
                [obj.errorCode, obj.TankMinimumWaterLevel(i)] = ENgetnodevalue(i, 20);
                [obj.errorCode, obj.TankMaximumWaterLevel(i)] = ENgetnodevalue(i, 21);
                [obj.errorCode, obj.TankMinimumFraction(i)] = ENgetnodevalue(i, 22);
                [obj.errorCode, obj.TankBulkReactionCoeff(i)] = ENgetnodevalue(i, 23);
            end
            
            %ENgetoption
            [obj.errorCode, obj.OptionsTrials] = ENgetoption(0);
            [obj.errorCode, obj.OptionsAccuracy] = ENgetoption(1);
            [obj.errorCode, obj.OptionsTolerance] = ENgetoption(2);
            [obj.errorCode, obj.OptionsEmmiterExponent] = ENgetoption(3);
            [obj.errorCode, obj.OptionsDemandMult] = ENgetoption(4);
            
            
            obj.PatternID={};
            for i=1:obj.getCountPatterns
                %ENgetpatternid
                [obj.errorCode, obj.PatternID{i}] = ENgetpatternid(i);
                %ENgetpatterindex
                [obj.errorCode, obj.PatternIndex(i)] = ENgetpatternindex(obj.PatternID{i});
                %Engetpatternlen
                [obj.errorCode, obj.PatternLengths(i)] = ENgetpatternlen(i);
            end
            
            %ENgetqualtype
            %tmpQuality={'NONE', 'CHEM', 'AGE', 'TRACE', 'MULTIS'};
            [obj.errorCode, obj.QualityCode,obj.QualityTraceNodeIndex] = ENgetqualtype();
            obj.QualityType=obj.TYPEQUALITY(obj.QualityCode+1);
            
            %ENgettimeparam
            [obj.errorCode, obj.TimeSimulationDuration] = ENgettimeparam(0);
            [obj.errorCode, obj.TimeHydraulicStep] = ENgettimeparam(1);
            [obj.errorCode, obj.TimeQualityStep] = ENgettimeparam(2);
            [obj.errorCode, obj.TimePatternStep] = ENgettimeparam(3);
            [obj.errorCode, obj.TimePatternStart] = ENgettimeparam(4);
            [obj.errorCode, obj.TimeReportingStep] = ENgettimeparam(5);
            [obj.errorCode, obj.TimeReportingStart] = ENgettimeparam(6);
            [obj.errorCode, obj.TimeRuleControlStep] = ENgettimeparam(7);
            [obj.errorCode, obj.TimeStatisticsIndex] = ENgettimeparam(8);
            %tmpStats={'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'};
            obj.TimeStatisticsType=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
            [obj.errorCode, obj.TimeReportingPeriods] = ENgettimeparam(9);
            
            %ENgetVersion
            [obj.errorCode, obj.Version] = ENgetversion();
            
            
            %Curves
            value=getCurveInfo(obj);
            obj.CurvesNameID=value.CurvesID;
            obj.CurvesX=value.CurveX;
            obj.CurvesY=value.CurveY;
            
            %Coordinates
            [obj.CoordinatesXY{1},obj.CoordinatesXY{2},obj.CoordinatesXY{3},obj.CoordinatesXY{4}]=getNodeCoord(obj);      
            
            %Units
            value=obj.getFlowUnitsHeadlossFormula;
            obj.FlowUnits=value.FlowUnits;
            obj.ConcentrationUnits=value.ConcentrationUnits;
            obj.Headloss=value.Headloss;
            obj.DemandsUnits=value.DemandsUnits;
            obj.DiameterPipesUnits=value.DiameterPipesUnits;
            obj.DiameterTanksUnits=value.DiameterTanksUnits;
            obj.EfficiencyUnits=value.EfficiencyUnits;
            obj.ElevationUnits=value.ElevationUnits;
            obj.EmitterCoeffUnits=value.EmitterCoeffUnits;
            obj.EnergyUnits=value.EnergyUnits;
            obj.FrictionFactorUnits=value.FrictionFactorUnits;
            obj.HeadUnits=value.HeadUnits;
            obj.LengthUnits=value.LengthUnits;
            obj.MinorLossCoeffUnits=value.MinorLossCoeffUnits;
            obj.PowerUnits=value.PowerUnits;
            obj.PressureUnits=value.PressureUnits;
            obj.ReactionCoeffBulkUnits=value.ReactionCoeffBulkUnits;
            obj.ReactionCoeffWallUnits=value.ReactionCoeffWallUnits;
            obj.RoughnessCoeffUnits=value.RoughnessCoeffUnits;
            obj.SourceMassInjectionUnits=value.SourceMassInjectionUnits;
            obj.VelocityUnits=value.VelocityUnits;
            obj.VolumeUnits=value.VolumeUnits;
            obj.WaterAgeUnits=value.WaterAgeUnits;
        end       
        %CONTROLS: EPANET cannot add new controls
       
        %ENsolveH
        function solveCompleteHydraulics(obj)
            [obj.errorCode] = ENsolveH();
        end
        
        %ENsolveQ
        function solveCompleteQuality(obj)
            [obj.errorCode] = ENsolveQ();
        end
        
        %%%%%%%%%%%%%%%%% ADD FUNCTIONS %%%%%%%%%%%%%%%%%
        
        %ENaddpattern
        function valueIndex = addPattern(obj,varargin)
            valueIndex=-1;
            if nargin==2
                [obj.errorCode] = ENaddpattern(varargin{1});
                valueIndex = getPatternIndex(obj,varargin{1});
            elseif nargin==3
                [obj.errorCode] = ENaddpattern(varargin{1});
                valueIndex = getPatternIndex(obj,varargin{1});
                setPattern(obj,valueIndex,varargin{2});
            end
        end
        
        
        %%%%%%%%%%%%%%%%% SET FUNCTIONS %%%%%%%%%%%%%%%%%
        
        
        %ENsetcontrol
        function setControl(obj,controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel)
            % Example: d.setControl(1,1,13,1,11,150)
            % controlRuleIndex must exist
            if controlRuleIndex<=obj.getCountControls
                [obj.errorCode] = ENsetcontrol(controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel);
                obj.ControlType={};
                for i=1:obj.getCountControls
                    [obj.errorCode, obj.ControlTypeIndex(i),obj.ControlLink(i),obj.ControlSeetting(i),obj.ControlNodeIndex(i),obj.ControlLevel(i)] = ENgetcontrol(i);
                    obj.ControlType(i)=obj.TYPECONTROL(obj.ControlTypeIndex(i)+1);
                end
                obj.ControlAll={obj.ControlType,obj.ControlTypeIndex,obj.ControlLink,obj.ControlSeetting,obj.ControlNodeIndex,obj.ControlLevel};
            else
                disp('New rules cannot be added in this version')
            end
        end
        
        %setLinkID --- not exist
        
        %ENsetlinkvalue
        function setLinkDiameter(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 0, value(i));
            end
        end
%         function setLinkDiameter(obj,index, value)
%             [obj.errorCode] = ENsetlinkvalue(index, 0, value);
%         end
        function setLinkLength(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 1, value(i));
            end
        end
        function setLinkRoughness(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 2, value(i));
            end
        end
        function setLinkMinorLossCoeff(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 3, value(i));
            end
        end
        function setLinkInitialStatus(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 4, value(i));
            end
        end
        function setLinkInitialSettings(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 5, value(i));
            end
        end
        function setLinkBulkReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 6, value(i));
            end
        end
        function setLinkWallReactionCoeff(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 7, value(i));
            end
        end
        function setLinkStatus(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 11, value(i));
            end
        end
        function setLinkSettings(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetlinkvalue(i, 12, value(i));
            end
        end
        
        
        %ENsetnodevalue
        function setNodeElevation(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 0, value(i));
            end
%             [obj.errorCode, obj.NodeElevations(index)] = ENgetnodevalue(index, 0);
        end
        function setNodeBaseDemand(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 1, value(i));
            end
            %[obj.errorCode, obj.NodeBaseDemands(index)] = ENgetnodevalue(index, 1);
        end
        function setNodeDemandPatternIndex(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 2, value(i));
            end
        end
        function setNodeEmitterCoeff(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 3, value(i));
            end
        end
        function setNodeInitialQuality(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 4, value(i));
            end
        end
        function setTankLevelInitial(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 8, value(i));
            end
        end
        
        %ENsetoption
        function setOptionTrial(obj,value)
            [obj.errorCode] = ENsetoption(0,value);
            [obj.errorCode, obj.OptionsTrials] = ENgetoption(0);
        end
        function setOptionAccuracy(obj,value)
            [obj.errorCode] = ENsetoption(1,value);
            [obj.errorCode, obj.OptionsAccuracy] = ENgetoption(1);
        end
        function setOptionTolerance(obj,value)
            [obj.errorCode] = ENsetoption(2,value);
            [obj.errorCode, obj.OptionsTolerance] = ENgetoption(2);
        end
        function setOptionEmitterExponent(obj,value)
            [obj.errorCode] = ENsetoption(3,value);
            [obj.errorCode, obj.OptionsEmmiterExponent] = ENgetoption(3);
        end
        function setOptionDemandMult(obj,value)
            [obj.errorCode] = ENsetoption(4,value);
            [obj.errorCode, obj.OptionsDemandMult] = ENgetoption(4);
        end
        
        %ENsettimeparam
        function setTimeSimulationDuration(obj,value)
            [obj.errorCode] = ENsettimeparam(0,value);
            [obj.errorCode, obj.TimeSimulationDuration] = ENgettimeparam(0);
            
        end
        function setTimeHydraulicStep(obj,value)
            [obj.errorCode] = ENsettimeparam(1,value);
            [obj.errorCode, obj.TimeHydraulicStep] = ENgettimeparam(1);
        end
        function setTimeQualityStep(obj,value)
            [obj.errorCode] = ENsettimeparam(2,value);
            [obj.errorCode, obj.TimeQualityStep] = ENgettimeparam(2);
        end
        function setTimePatternStep(obj,value)
            [obj.errorCode] = ENsettimeparam(3,value);
            [obj.errorCode, obj.TimePatternStep] = ENgettimeparam(3);
        end
        function setTimePatternStart(obj,value)
            [obj.errorCode] = ENsettimeparam(4,value);
            [obj.errorCode, obj.TimePatternStart] = ENgettimeparam(4);
        end
        function setTimeReportingStep(obj,value)
            [obj.errorCode] = ENsettimeparam(5,value);
            [obj.errorCode, obj.TimeReportingStep] = ENgettimeparam(5);
        end
        function setTimeReportingStart(obj,value)
            [obj.errorCode] = ENsettimeparam(6,value);
            [obj.errorCode, obj.TimeReportingStart] = ENgettimeparam(6);
        end
        
        function setTimeStatistics(obj,value)
            %'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'
            tmpindex=find(strcmpi(obj.TYPESTATS,value)==1)-1;
            [obj.errorCode] = ENsettimeparam(8,tmpindex);
            [obj.errorCode, obj.TimeStatisticsIndex] = ENgettimeparam(8);
        end
        
        %ENsetpattern
        function setPattern(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errorCode] = ENsetpattern(index, patternVector, nfactors);
        end
        function setPatternMatrix(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errorCode] = ENsetpattern(i, patternMatrix(i,:), nfactors);
            end
        end
        %ENsetpatternvalue
        function setPatternValue(obj,index, patternTimeStep, patternFactor)
            [obj.errorCode] = ENsetpatternvalue(index, patternTimeStep, patternFactor);
        end
        
        
        function setQualityType(obj,varargin)
            qualcode=0;
            chemname='';
            chemunits='';
            tracenode='';
            if find(strcmpi(varargin,'none')==1)
                [obj.errorCode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'age')==1)
                qualcode=2;
                [obj.errorCode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'chem')==1)
                qualcode=1;
                chemname=varargin{1};
                chemunits=varargin{2};
                [obj.errorCode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            elseif find(strcmpi(varargin,'trace')==1)
                qualcode=3;
                tracenode=varargin{2};
                [obj.errorCode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode);
            end
        end
        
        %ENresetreport
        function setReportFormatReset(obj)
            [obj.errorCode]=ENresetreport();
        end
        
        %ENsetstatusreport
        function setReportStatus(obj,value) %'yes','no','full'
            statuslevel=find(strcmpi(obj.TYPEREPORT,value)==1)-1;
            [obj.errorCode] = ENsetstatusreport(statuslevel);
        end
        
        %ENsetreport
        function setReport(obj,value)
            [obj.errorCode] = ENsetreport(value);
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%% GET FUNCTIONS %%%%%%%%%%%%%%%%%
        
        function value = getControls(obj)
            for i=1:obj.getCountControls
                [obj.errorCode, obj.ControlTypeIndex(i),obj.ControlLink(i),obj.ControlSeetting(i),obj.ControlNodeIndex(i),obj.ControlLevel(i)] = ENgetcontrol(i);
                obj.ControlType(i)=obj.TYPECONTROL(obj.ControlTypeIndex(i)+1);
            end
            obj.ControlAll={obj.ControlType,obj.ControlTypeIndex,obj.ControlLink,obj.ControlSeetting,obj.ControlNodeIndex,obj.ControlLevel};
            value=obj.ControlAll;
        end
        
        %ENgetcount
        function value  =  getCountNodes(obj)
            % Nodes, Tanks/Reservoirs, Links, Patterns, Curves, Controls
            [obj.errorCode, value] = ENgetcount(0);
        end
        function value  =  getCountTanksReservoirs(obj)
            [obj.errorCode, value] = ENgetcount(1);
        end
        function value  =  getCountLinks(obj)
            [obj.errorCode, value] = ENgetcount(2);
        end
        function value  =  getCountPatterns(obj)
            [obj.errorCode, value] = ENgetcount(3);
        end
        function value  =  getCountCurves(obj)
            [obj.errorCode, value] = ENgetcount(4);
        end
        function value  =  getCountControls(obj)
            [obj.errorCode, value] = ENgetcount(5);
        end
        
        
        %ENgeterror
        function value = getError(obj,errcode)
            [obj.errorCode, value] = ENgeterror(errcode);
        end
        
        %ENgetflowunits
        function value = getFlowUnits(obj)
            [obj.errorCode, obj.flowUnits] = ENgetflowunits();
            obj.UnitsType=obj.TYPEUNITS(obj.flowUnits+1);
            value=obj.UnitsType;
        end
        
        %ENgetlinkid
        function value = getLinkID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountLinks
                    [obj.errorCode, value{i}]=ENgetlinkid(i);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, value{k}]=ENgetlinkid(i);
                    k=k+1;
                end
            end
        end
        
        %ENgetlinkindex
        function value = getLinkIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountLinks;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = ENgetlinkindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, value] = ENgetlinkindex(varargin{1});
            end
        end
        
        %ENgetlinknodes
        function value = getLinkNodes(obj)
            for i=1:obj.getCountNodes
                [obj.errorCode,linkFromNode,linkToNode] = ENgetlinknodes(i);
                value(i,:)= [linkFromNode,linkToNode];
            end
        end
        
        %ENgetlinknodes
        function value = getLinkNodesIndex(obj)
            for i=1:obj.getCountLinks
                [obj.errorCode,linkFromNode,linkToNode] = ENgetlinknodes(i);
                value(i,:)= [linkFromNode,linkToNode];
            end
        end
        
        %ENgetlinktype
        function value = getLinkType(obj)
            for i=1:obj.CountLinks
                [obj.errorCode,obj.LinkTypeIndex(i)] = ENgetlinktype(i);
                if obj.LinkTypeIndex(i)>2
                    obj.LinkTypeIndex(i)=9; %Valve  
                elseif obj.LinkTypeIndex(i)==1
                    obj.LinkTypeIndex(i)=1; %cvpipe pipe
                end
                value(i)=obj.TYPELINK(obj.LinkTypeIndex(i)+1);
            end
        end
        
        %ENgetlinkvalue
        function value = getLinkDiameter(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,0);
            end
        end
        function value = getLinkLength(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,1);
            end
        end
        function value = getLinkRoughness(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,2);
            end
        end
        function value = getLinkMinorLossCoeff(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,3);
            end
        end
        function value = getLinkInitialStatus(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,4);
            end
        end
        function value = getLinkInitialSettings(obj)
            %Roughness for pipes,initial speed for pumps,initial setting
            %for valves
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,5);
            end
        end
        function value = getLinkBulkReactionCoeff(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,6);
            end
        end
        function value = getLinkWallReactionCoeff(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,7);
            end
        end
        function value = getLinkFlows(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,8);
            end
        end
        function value = getLinkVelocity(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,9);
            end
        end
        function value = getLinkHeadloss(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,10);
            end
        end
        function value = getLinkStatus(obj)
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,11);
            end
        end
        function value = getLinkSettings(obj) %Roughness for pipes, actual speed for pumps, actual setting for valves
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,12);
            end
        end
        function value = getLinkEnergy(obj) %in kwatts
            value=zeros(1,obj.CountLinks);
            for i=1:obj.CountLinks
                [obj.errorCode, value(i)] = ENgetlinkvalue(i,13);
            end
        end
        
        
        %ENgetnodeid
        function value = getNodeID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountNodes
                    [obj.errorCode, value{i}]=ENgetnodeid(i);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, value{k}]=ENgetnodeid(i);
                    k=k+1;
                end
            end
        end
        
        %ENgetnodeindex
        function value = getNodeIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountNodes;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = ENgetnodeindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, value] = ENgetnodeindex(varargin{1});
            end
        end
        
        %ENgetnodetype
        function value = getNodeType(obj)
            for i=1:obj.getCountNodes
                [obj.errorCode,obj.NodeTypeIndex(i)] = ENgetnodetype(i);
                value(i)=obj.TYPENODE(obj.NodeTypeIndex(i)+1);
            end
        end
        
        %ENgetnodevalue
        function value = getNodeElevation(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,0);
            end
        end
        function value = getNodeBaseDemand(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,1);
            end
        end
        function value = getNodeDemandPatternIndex(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,2);
            end
        end
        function value = getNodeEmitterCoeff(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,3);
            end
        end
        function value = getNodeInitialQuality(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,4);
            end
        end
        
        
        
        function value = getTankLevelInitial(obj)
            value=nan(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,8);
            end
        end
        
        function value = getNodeActualDemand(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,9);
            end
        end
        function value = getNodeActualDemandSensingNodes(obj,varargin)
            value=zeros(1,obj.getCountNodes);
            for i=1:length(varargin{1})
                [obj.errorCode, value(varargin{1}(i))] = ENgetnodevalue(varargin{1}(i),12);
            end
        end
        function value = getNodeHydaulicHead(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,10);
            end
        end
        function value = getNodePressure(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,11);
            end
        end
        function value = getNodeActualQuality(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,12);
            end
        end
        function value = getNodeMassFlowRate(obj) %Mass flow rate per minute of a chemical source
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.CountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,13);
            end
        end
        function value = getNodeActualQualitySensingNodes(obj,varargin)
            value=zeros(1,obj.getCountNodes);
            for i=1:length(varargin{1})
                [obj.errorCode, value(varargin{1}(i))] = ENgetnodevalue(varargin{1}(i),12);
            end
        end
        
        %ENgetoption
        function value = getOptionTrial(obj)
            [obj.errorCode, value] = ENgetoption(0);
        end
        function value = getOptionAccuracy(obj)
            [obj.errorCode, value] = ENgetoption(1);
        end
        function value = getOptionTolerance(obj)
            [obj.errorCode, value] = ENgetoption(2);
        end
        function value = getOptionEmitterExponent(obj)
            [obj.errorCode, value] = ENgetoption(3);
        end
        function value = getOptionDemandMult(obj)
            [obj.errorCode, value] = ENgetoption(4);
        end
        
        
        
        %ENgetpatternid
        function value = getPatternID(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountPatterns
                    [obj.errorCode, value{i}]=ENgetpatternid(i);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, value{k}]=ENgetpatternid(i);
                    k=k+1;
                end
            end
        end
        
        %ENgetpatternindex
        function value = getPatternIndex(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountPatterns;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = ENgetpatternindex(varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, value] = ENgetpatternindex(varargin{1});
            end
        end
        
        %ENgetpatternlen
        function value = getPatternLength(obj,varargin)
            if isempty(varargin)
                tmpPatterns=1:obj.getCountPatterns;
                for i=tmpPatterns
                    [obj.errorCode, value(i)]=ENgetpatternlen(i);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = ENgetpatternlen(obj.getPatternIndex(varargin{1}{j}));
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, value] = ENgetpatternlen(obj.getPatternIndex(varargin{1}));
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errorCode, value(k)]=ENgetpatternlen(i);
                    k=k+1;
                end
            end
        end
        
        %ENgetpatternvalue
        function value = getPattern(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getPatternLength);
            value=nan(obj.getCountPatterns,tmpmaxlen);
            for i=1:obj.getCountPatterns
                tmplength=obj.getPatternLength(i);
                for j=1:tmplength
                    [obj.errorCode, value(i,j)] = ENgetpatternvalue(i, j);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                    
            end
        end
        
        %ENgetpatternvalue
        function value = getPatternValue(obj,patternIndex, patternStep) %Mass flow rate per minute of a chemical source
            [obj.errorCode, value] = ENgetpatternvalue(patternIndex, patternStep);
        end
        
        
        %ENgetqualtype
        function value = getQualityType(obj)
            [obj.errorCode, obj.QualityCode,obj.QualityTraceNodeIndex] = ENgetqualtype();
            value=obj.TYPEQUALITY(obj.QualityCode+1);
        end
        
        
        %ENgettimeparam
        function value = getTimeSimulationDuration(obj)
            [obj.errorCode, value] = ENgettimeparam(0);
        end
        function value = getTimeHydraulicStep(obj)
            [obj.errorCode, value] = ENgettimeparam(1);
        end
        function value = getTimeQualityStep(obj)
            [obj.errorCode, value] = ENgettimeparam(2);
        end
        function value = getTimePatternStep(obj)
            [obj.errorCode, value] = ENgettimeparam(3);
        end
        function value = getTimePatternStart(obj)
            [obj.errorCode, value] = ENgettimeparam(4);
        end
        function value = getTimeReportingStep(obj)
            [obj.errorCode, value] = ENgettimeparam(5);
        end
        function value = getTimeReportingStart(obj)
            [obj.errorCode, value] = ENgettimeparam(6);
        end
        function value = getTimeStatistics(obj)
            [obj.errorCode, obj.TimeStatisticsIndex] = ENgettimeparam(8);
            %tmpStats={'NONE','AVERAGE','MINIMUM','MAXIMUM', 'RANGE'};
            value=obj.TYPESTATS(obj.TimeStatisticsIndex+1);
        end
        function value = getTimeReportingPeriods(obj)
            [obj.errorCode, value] = ENgettimeparam(9);
        end
        
        
        
        %ENreport
        function getReport(obj)
            [obj.errorCode]=ENreport();
        end
        
        function value=getVersion(obj)
            [obj.errorCode, value] = ENgetversion();
        end
        
         function value=getComputedHydraulicTimeSeries(obj)
            obj.openHydraulicAnalysis;
            obj.initializeHydraulicAnalysis
            tstep=1; 
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeHydraulicStep+1;
            initnodematrix=zeros(totalsteps, obj.getCountNodes);
            initlinkmatrix=zeros(totalsteps, obj.getCountLinks);
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
                value.Energy(k,:)=obj.getLinkEnergy;   
                tstep = obj.nextHydraulicAnalysisStep;
                k=k+1;

            end
            obj.closeHydraulicAnalysis;
        end       
        
        
        
        function value=getComputedQualityTimeSeries(obj,varargin)
            obj.openQualityAnalysis
            obj.initializeQualityAnalysis
            tleft=1; 
            totalsteps=obj.getTimeSimulationDuration/obj.getTimeQualityStep;
            initnodematrix=zeros(totalsteps, obj.getCountNodes);
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
            k=1;
            while (tleft>0)
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
                %tstep=obj.nextQualityAnalysisStep;
                k=k+1;
            end
            obj.closeQualityAnalysis;
        end
        
        
        %%%%%%%%%%%%%%%%% OPERATIONS %%%%%%%%%%%%%%%%%%%
        
        %ENclose & ENMatlabCleanup
        function unload(varargin)
            ENclose;
            ENMatlabCleanup;
        end
        
        
        %ENcloseH
        function closeHydraulicAnalysis(obj)
            [obj.errorCode] = ENcloseH();
        end
        
        %ENcloseQ
        function closeQualityAnalysis(obj)
            [obj.errorCode] = ENcloseQ();
        end
        
        
        
        %ENsavehydfile
        function saveHydraulicFile(obj,hydname)
            [obj.errorCode]=ENsavehydfile(hydname);
        end
        
        
        %ENusehydfile
        function useHydraulicFile(obj,hydname)
            [obj.errorCode]=ENusehydfile(hydname);
        end
        
        
        %ENinitH
        function initializeHydraulicAnalysis(obj)
            [obj.errorCode] = ENinitH(1);
        end
        
        %ENinitQ
        function initializeQualityAnalysis(obj)
            [obj.errorCode] = ENinitQ(1);
        end
        
        %ENnextH
        function tstep = nextHydraulicAnalysisStep(obj)
            [obj.errorCode, tstep] = ENnextH();
        end
        
        %ENnextQ
        function tstep = nextQualityAnalysisStep(obj)
            [obj.errorCode, tstep] = ENnextQ();
        end
        
        %ENopenH
        function openHydraulicAnalysis(obj)
            [obj.errorCode] = ENopenH();
        end
        
        %ENopenQ
        function openQualityAnalysis(obj)
            [obj.errorCode] = ENopenQ();
        end
        
        %ENrunH
        function tstep = runHydraulicAnalysis(obj)
            [obj.errorCode, tstep] = ENrunH();
        end
        
        %ENrunQ
        function tstep = runQualityAnalysis(obj)
            [obj.errorCode, tstep] = ENrunQ();
        end
        
        
        %ENsaveH
        function saveHydraulicsOutputReportingFile(obj)
            [obj.errorCode] = ENsaveH();
        end
        
        
        %ENstepQ
        function tleft=stepQualityAnalysisTimeLeft(obj)
            [obj.errorCode, tleft] = ENstepQ();
        end
        
        %ENsaveinpfile
        function saveInputFile(obj,inpname)
            [obj.errorCode] = ENsaveinpfile(inpname);
        end
        
        %ENwriteline
        function writeLineInReportFile(obj, line)
            [obj.errorCode] = ENwriteline (line);
        end
        

        function LoadMSX(obj,msxname)
           [obj] = MSXMatlabSetup(obj,msxname);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = getNodeSourceQuality(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,5);
            end
        end
        function setNodeSourceQuality(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 5, value(i));
            end
        end
        function value = getNodeSourcePatternIndex(obj)
            value=zeros(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, value(i)] = ENgetnodevalue(i,6);
            end
        end
        function setNodeSourcePatternIndex(obj, value)
            for i=1:length(value)
                [obj.errorCode] = ENsetnodevalue(i, 6, value(i));
            end
        end
        function value = getNodeSourceType(obj)
            %value=zeros(1,obj.getCountNodes);
            value=cell(1,obj.getCountNodes);
            for i=1:obj.getCountNodes
                [obj.errorCode, temp] = ENgetnodevalue(i,7);
                if ~isnan(temp)
                    value(i)=obj.TYPESOURCE(temp+1);
                end
            end
        end
        function setNodeSourceType(obj, index, value)
            value=find(strcmpi(obj.TYPESOURCE,value)==1)-1;
            [obj.errorCode] = ENsetnodevalue(index, 7, value);
        end
        function value = getTimeRuleControlStep(obj)
            [obj.errorCode, value] = ENgettimeparam(7);
        end
        function setTimeRuleControlStep(obj,value)
            [obj.errorCode] = ENsettimeparam(7,value);
            [obj.errorCode, obj.TimeRuleControlStep] = ENgettimeparam(7);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value=getPumpType(obj)
            tmpLinkTypes=obj.getLinkType;
            indices=find(strcmp(tmpLinkTypes,'PUMP')==1);
            value=cell(1,length(tmpLinkTypes));
            for index=indices
                [obj.errorCode, v] = ENgetpumptype(index);
                value(index)=obj.TYPEPUMP(v+1);
            end
        end
        
        %%%% New Functions
        % Get
        function value=getCurveInfo(obj)
         [value.CurvesID,value.CurveX,value.CurveY]=CurveInfo(obj);
        end
        function value=getLinksInfo(obj)
            value=LinksInfo(obj);
        end
        function value=getNodesInfo(obj)
            value=NodesInfo(obj);
        end
        function value=getPumpInfo(obj)
            value=PumpInfo(obj);
        end
        function value=getControlsInfo(obj)
            value=ControlsInfo(obj);
        end
        function value=getFlowUnitsHeadlossFormula(obj)
            value=FlowUnitsHeadlossFormula(obj);
        end
        
        % Add
        function addNewCurve(obj,newCurveID,CurveX,CurveY)
        	addCurve(obj,newCurveID,CurveX,CurveY);  
        end
        function addPipe(obj,newLink,fromNode,toNode)
            addLink(obj,1,newLink,fromNode,toNode);
        end 
        function addPump(obj,newLink,fromNode,toNode,curveID)
            addLink(obj,2,newLink,fromNode,toNode,curveID);
        end
        function addValvePRV(obj,newLink,fromNode,toNode)
            addLink(obj,3,newLink,fromNode,toNode); % Pressure Reducing Valve 
        end
        function addValvePSV(obj,newLink,fromNode,toNode)
            addLink(obj,4,newLink,fromNode,toNode); % Pressure Sustaining Valve 
        end
        function addValvePBV(obj,newLink,fromNode,toNode)
            addLink(obj,5,newLink,fromNode,toNode); % Pressure Breaker Valve 
        end
        function addValveFCV(obj,newLink,fromNode,toNode)
            addLink(obj,6,newLink,fromNode,toNode); % Flow Control Valve
        end
        function addValveTCV(obj,newLink,fromNode,toNode)
            addLink(obj,7,newLink,fromNode,toNode); % Throttle Control Valve 
        end
        function addValveGPV(obj,newLink,fromNode,toNode)
            addLink(obj,8,newLink,fromNode,toNode); % General Purpose Valve 
        end
        
        
        function addJunction(obj,newID,X,Y)
            addNode(obj,0,newID,X,Y)
        end 
        function addReservoir(obj,newID,X,Y)
            addNode(obj,1,newID,X,Y)
        end
        function addTank(obj,newID,X,Y)
            addNode(obj,2,newID,X,Y)
        end
        function addControl(obj,x,status,y_t_c,param,z,varargin)
            if nargin==6
                addNewControl(obj,x,status,y_t_c,param,z)
            elseif nargin==5
                addNewControl(obj,x,status,y_t_c,param)
            else
                addNewControl(obj,x,status,y_t_c)
            end
        end
        
        % Remove
        function removeCurveID(obj,CurveID)
        	rmCurveID(obj,CurveID);
        end
        function removeControlLinkID(obj,ID)
        	rmControl(obj,1,ID);
        end
        function removeControlNodeID(obj,ID)
        	rmControl(obj,0,ID);
        end
        function warn1=removeLinkID(obj,LinkID)
            warn1=rmLink(obj,LinkID);
        end
        function removeNodeID(obj,NodeID)
            rmNode(obj,NodeID);
        end
        % Set Flow Units
        function setFlowUnitsGPM(obj)
            Options(obj,'GPM') %gallons per minute
        end
        function setFlowUnitsLPS(obj)
            Options(obj,'LPS') %liters per second
        end        
        function setFlowUnitsMGD(obj)
            Options(obj,'MGD') %million gallons per day
        end  
        function setFlowUnitsIMGD(obj)
            Options(obj,'IMGD') %Imperial mgd
        end
        function setFlowUnitsCFS(obj)
            Options(obj,'CFS') %cubic feet per second
        end
        function setFlowUnitsAFD(obj)
            Options(obj,'AFD') %acre-feet per day
        end
        function setFlowUnitsLPM(obj)
            Options(obj,'LPM') %liters per minute
        end
        function setFlowUnitsMLD(obj)
            Options(obj,'MLD') %million liters per day
        end    
        function setFlowUnitsCMH(obj)
            Options(obj,'CMH') %cubic meters per hour
        end   
        function setFlowUnitsCMD(obj)
            Options(obj,'CMD') %cubic meters per day
        end           	
        % Set HeadLoss Formula
        function setHeadlossHW(obj)
            Options(obj,'','H-W')  %Hazen-Wiliams
        end        
        function setHeadlossDW(obj)
            Options(obj,'','D-W')  %Darcy-Weisbach
        end  
        function setHeadlossCM(obj)
            Options(obj,'','C-M')  %Chezy-Manning
        end
        
%       function value=getCurveID(obj)
%             tmpLinkTypes=obj.getLinkType;
%             indices=find(strcmp(tmpLinkTypes,'PUMP')==1);
%             value=cell(1,length(tmpLinkTypes));
%             for index=indices
%                 [obj.errorCode, value{index}] = ENgetheadcurve(index);
%             end
%        end
        
%         %% FOR FUTURE VERSIONS
%         function value=getCurves(obj)
%             tmpLinkTypes=obj.getLinkType;
%             indices=find(strcmp(tmpLinkTypes,'PUMP')==1);
%             value=cell(1,length(tmpLinkTypes));
%             for index=indices
%                 [errcode, nValues, xValues, yValues] = ENgetcurve(index);
%                 value(index)=[nValues, xValues, yValues];
%             end
%         end        
       
        % ENplot
        function plot(obj,varargin)    
            ENplot(obj,varargin{:});
        end
        
        function CoordinatesXY = getCoordinates(inpname)
            [vx,vy,vertx,verty]  = getNodeCoord(inpname);
            CoordinatesXY{1} = vx;
            CoordinatesXY{2} = vy;
            CoordinatesXY{3} = vertx;
            CoordinatesXY{4} = verty;
        end
        
        %%%%%%%%%%%%%%%%%%% EPANET - MSX %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% SOLVE FUNCTIONS %%%%%%%%%%%%%%%%%
        %MSXsolveH
        function solveCompleteHydraulicsMsx(obj)
            [obj.errorCode] = MSXsolveH();
        end
        
        %MSXsolveQ
        function solveCompleteQualityMsx(obj)
            [obj.errorCode] = MSXsolveQ();
        end
        
        %%%%%%%%%%%%%%%%% ADD FUNCTIONS %%%%%%%%%%%%%%%%%
        
        %MSXaddpattern
        function valueIndex = addPatternMsx(obj,varargin)
            valueIndex=-1;
            if nargin==2
                [obj.errorCode] = MSXaddpattern(varargin{1});
                [obj.errorCode, valueIndex] = MSXgetindex(obj,7,varargin{1}); 
            elseif nargin==3
                [obj.errorCode] = MSXaddpattern(varargin{1});
                [obj.errorCode, valueIndex] = MSXgetindex(obj,7,varargin{1}); 
                setPatternMsx(obj,valueIndex,varargin{2});
            end
        end
        
        %%%%%%%%%%%%%%%%% SET FUNCTIONS %%%%%%%%%%%%%%%%%
        %MSXsetsource
        function setSourceMsx(obj, node, species, type, level, pat)
            MSXsetsource(node, species, type, level, pat);
        end

        %MSXsetconstant  	
        function setConstantValueMsx(obj, value)
            for i=1:length(value)
                [obj.errorCode] = MSXsetconstant(i, value(i));
            end
        end
        
        %MSXsetparameter
        function setParameterTankValueMsx(obj, tankIndex, value)
            if ~sum(tankIndex==obj.TankIndex)
                fprintf('>> Invalid Tank Index <<\n');obj.TankIndex
                return;
            end
            for i=1:length(value)
                [obj.errorCode] = MSXsetparameter(0, tankIndex, i, value(i));
            end
        end
        
        function setParameterPipeValueMsx(obj, pipeIndex, value)
            for i=1:length(value)
                [obj.errorCode] = MSXsetparameter(1, pipeIndex, i, value(i));
            end
        end
        
        %MSXsetinitqual
        function setInitqualNodeValueMsx(obj, value)
            for i=1:length(value)
                for j=1:size(value{1},2)
                    [obj.errorCode] = MSXsetinitqual(0, i, j, value{i}(j));
                end
            end
        end
        
        function setInitqualLinkValueMsx(obj, value)
            for i=1:size(value,1)
                for j=1:size(value{1},2)
                    [obj.errorCode] = MSXsetinitqual(1, i, j, value{i}(j));
                end
            end
        end
        
        %MSXsetpattern
        function setPatternMsx(obj,index,patternVector)
            nfactors=length(patternVector);
            [obj.errorCode] = MSXsetpattern(index, patternVector, nfactors);
        end
        
        function setPatternMatrixMsx(obj,patternMatrix)
            nfactors=size(patternMatrix,2);
            for i=1:size(patternMatrix,1)
                [obj.errorCode] = MSXsetpattern(i, patternMatrix(i,:), nfactors);
            end
        end
        
        %MSXsetpatternvalue
        function setPatternValueMsx(obj,index, patternTimeStep, patternFactor)
            [obj.errorCode] = MSXsetpatternvalue(index, patternTimeStep, patternFactor);
        end
     
        %%%%%%%%%%%%%%%%% GET FUNCTIONS %%%%%%%%%%%%%%%%%
        
        %MSXgetsource
        function value = getSourcesMsx(obj)
            for i=1:obj.CountNodes
                for j=1:obj.CountSpeciesMsx 
                   [obj.errorCode, obj.SourceTypeMsx{i}{j},obj.SourceLevelMsx{i}(j),obj.SourcePatternIndexMsx{i}(j)] = MSXgetsource(i,j);
                end
            end
            obj.SourceAllMsx={obj.SourceTypeMsx,obj.SourceLevelMsx,obj.SourcePatternIndexMsx,obj.SourceNodeIDMsx};
            value=obj.SourceAllMsx;
        end
        
        %MSXgetinitqual
        function value = getInitqualNodeValueMsx(obj)
            if obj.getCountSpeciesMsx==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.CountNodes
                for j=1:obj.getCountSpeciesMsx
                   [obj.errorCode, value{i}(j)] = MSXgetinitqual(0,i,j);   
                end
            end
        end
        
        function value = getInitqualLinkValueMsx(obj)
            if obj.getCountSpeciesMsx==0
                value{1}(1)=0;
                return;
            end
            for i=1:obj.CountLinks
                for j=1:obj.getCountSpeciesMsx
                   [obj.errorCode, value{i}(j)] = MSXgetinitqual(1,i,j);   
                end
            end
        end
            
        % MSXgetparameter
        function value = getParameterTankValueMsx(obj)
            value={};
            if ~obj.getCountParametersMsx
                value=0;return;
            end
            for i=1:length(obj.CountTanks)
                for j=1:obj.CountParametersMsx
                   [obj.errorCode, value{obj.TankIndex(i)}(j)] = MSXgetparameter(0,obj.TankIndex(i),j);   
                end
            end
        end
        
        function value = getParameterPipeValueMsx(obj)
            if ~obj.getCountParametersMsx
                value=0;return;
            end
            for i=1:obj.CountPipes
                for j=1:obj.CountParametersMsx
                   [obj.errorCode, value{i}(j)] = MSXgetparameter(1,i,j);   
                end
            end
        end
            
        %MSXgetconstant
        function value = getConstantValueMsx(obj)
            if ~obj.getCountConstantsMsx
                value=0;
            end
            for i=1:obj.getCountConstantsMsx
                [obj.errorCode, len] = MSXgetIDlen(6,i);
                [obj.errorCode, obj.ConstantNameIDMsx{i}] = MSXgetID(6,i,len);
                [obj.errorCode, value(i)] = MSXgetconstant(i);
            end
        end
            
        %MSXgetcount
        function value  =  getCountSpeciesMsx(obj)
            % Species, Constants, Parameters, Patterns 
            [obj.errorCode, value] = MSXgetcount(3);
        end
        function value  =  getCountConstantsMsx(obj)
            [obj.errorCode, value] = MSXgetcount(6);
        end
        function value  =  getCountParametersMsx(obj)
            [obj.errorCode, value] = MSXgetcount(5);
        end
        function value  =  getCountPatternsMsx(obj)
            [obj.errorCode, value] = MSXgetcount(7);
        end        
        
        %MSXgeterror
        function value = getErrorMsx(obj,errcode)
            [obj.errorCode, value] = MSXgeterror(errcode);
        end
               
        %Species ID
        function value = getSpeciesIDMsx(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountSpeciesMsx
                    [obj.errorCode, len] = MSXgetIDlen(3,i);
                    [obj.errorCode, value{i}]=MSXgetID(3,i,len);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, len] = MSXgetIDlen(3,i);
                    [obj.errorCode, value{k}]=MSXgetID(3,i,len);
                    k=k+1;
                end
            end
        end
        
        %Constants ID
        function value = getConstantsIDMsx(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountConstantsMsx
                    [obj.errorCode, len] = MSXgetIDlen(6,i);
                    [obj.errorCode, value{i}]=MSXgetID(6,i,len);
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, len] = MSXgetIDlen(6,i);
                    [obj.errorCode, value{k}]=MSXgetID(6,i,len);
                    k=k+1;
                end
            end
        end
        
        %Parameters ID
        function value = getParametersIDMsx(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountParametersMsx
                    [obj.errorCode, len] = MSXgetIDlen(5,i);
                    [obj.errorCode, value{i}]=MSXgetID(5,i,len);
                end
                if ~obj.getCountParametersMsx
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, len] = MSXgetIDlen(5,i);
                    [obj.errorCode, value{k}]=MSXgetID(5,i,len);
                    k=k+1;
                end
            end
        end
        
        %Patterns ID
        function value = getPatternsIDMsx(obj,varargin)
            if isempty(varargin)
                for i=1:obj.getCountPatternsMsx
                    [obj.errorCode, len] = MSXgetIDlen(7,i);
                    [obj.errorCode, value{i}]=MSXgetID(7,i,len);
                end
                if ~obj.getCountPatternsMsx
                    value=0;
                end
            else
                k=1;
                for i=varargin{1}
                    [obj.errorCode, len] = MSXgetIDlen(7,i);
                    [obj.errorCode, value{k}]=MSXgetID(7,i,len);
                    k=k+1;
                end
            end
        end
        
        %Species Index
        function value = getSpeciesIndexMsx(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountSpeciesMsx;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = MSXgetindex(3,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                   [obj.errorCode, value] = MSXgetindex(3,varargin{1});
            end
        end
        
        %Constant Index
        function value = getConstantsIndexMsx(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountConstantsMsx;
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = MSXgetindex(6,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                   [obj.errorCode, value] = MSXgetindex(6,varargin{1});
            end
        end
        
        %Parameter Index
        function value = getParametersIndexMsx(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountParametersMsx;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = MSXgetindex(5,varargin{1}{j});
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                   [obj.errorCode, value] = MSXgetindex(5,varargin{1});
            end
        end
        
        %Pattern Index
        function value = getPatternIndexMsx(obj,varargin)
            if isempty(varargin)
                value=1:obj.getCountPatternsMsx;
                if ~length(value)
                    value=0;
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, len] = MSXgetIDlen(7,j);
                    [obj.errorCode, value{k}] = MSXgetID(7, obj.PatternIndexMsx,len);
                    if obj.errorCode
                        value{k}=0;
                    end
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, obj.PatternIndexMsx] = MSXgetindex(obj,7,varargin{1});
                [obj.errorCode, len] = MSXgetIDlen(7,obj.PatternIndexMsx);
                [obj.errorCode, value] = MSXgetID(7, obj.PatternIndexMsx,len);
                if obj.errorCode
                    value=0;
                end
            end
        end
        
        %MSXgetpatternlen
        function value = getPatternLengthMsx(obj,varargin)
            if isempty(varargin)
                tmpPatterns=1:obj.getCountPatternsMsx;
                if length(tmpPatterns)==0
                    value=0;
                end
                for i=tmpPatterns
                    [obj.errorCode, value(i)]=MSXgetpatternlen(i);
                end
            elseif isa(varargin{1},'cell')
                k=1;
                for j=1:length(varargin{1})
                    [obj.errorCode, value(k)] = MSXgetpatternlen(obj.getPatternIndexMsx(varargin{1}{j}));
                    k=k+1;
                end
            elseif isa(varargin{1},'char')
                [obj.errorCode, value] = MSXgetpatternlen(obj.getPatternIndexMsx(varargin{1}));
            elseif isa(varargin{1},'numeric')
                k=1;
                for i=varargin{1}
                    [obj.errorCode, value(k)]=MSXgetpatternlen(i);
                    k=k+1;
                end
            end
        end
        
        %MSXgetpatternvalue
        function value = getPatternMsx(obj) %Mass flow rate per minute of a chemical source
            tmpmaxlen=max(obj.getPatternLengthMsx);
            value=nan(obj.getCountPatternsMsx,tmpmaxlen);
            for i=1:obj.getCountPatternsMsx
                tmplength=obj.getPatternLengthMsx(i);
                for j=1:tmplength
                    [obj.errorCode, value(i,j)] = MSXgetpatternvalue(i, j);
                end
                if tmplength<tmpmaxlen
                    for j=(tmplength+1):tmpmaxlen
                        value(i,j)=value(i,j-tmplength);
                    end
                end
                    
            end
        end
        
        %MSXgetpatternvalue
        function value = getPatternValueMsx(obj,patternIndex, patternStep) %Mass flow rate per minute of a chemical source
            [obj.errorCode, value] = MSXgetpatternvalue(patternIndex, patternStep);
        end
        
        %MSXreport
        function getReportMsx(obj)
            [obj.errorCode]=MSXreport();
        end
               
        function value=getComputedQualityNodeMsx(obj)
            if obj.getCountSpeciesMsx==0
                value=0;
                return;
            end
            for i=1:obj.CountNodes
                % Obtain a hydraulic solution
                obj.solveCompleteHydraulicsMsx();
                % Run a step-wise water quality analysis
                % without saving results to file
                obj.initializeQualityAnalysisMsx(0);

                [t, tleft]=obj.stepQualityAnalysisTimeLeftMsx();

                % Retrieve species concentration at node
                k=1;
                while(tleft>0 && obj.errorCode==0)
                    [t, tleft]=obj.stepQualityAnalysisTimeLeftMsx();
                    value.Time(k,:)=t;
                    for j=1:obj.getCountSpeciesMsx
                        value.Quality{i}(k,:)=obj.getSpeciesConcentration(0, i, j);%node code0
                    end
                    k=k+1;
                end
            end
        end
        
        function value=getComputedQualityLinkMsx(obj)
            if obj.getCountSpeciesMsx==0
                value=0;
                return;
            end
            for i=1:obj.CountLinks
                % Obtain a hydraulic solution
                obj.solveCompleteHydraulicsMsx();
                % Run a step-wise water quality analysis
                % without saving results to file
                obj.initializeQualityAnalysisMsx(0);

                [t, tleft]=obj.stepQualityAnalysisTimeLeftMsx();

                % Retrieve species concentration at node
                k=1;
                while(tleft>0 && obj.errorCode==0)
                    [t, tleft]=obj.stepQualityAnalysisTimeLeftMsx();
                    value.Time(k,:)=t;
                    for j=1:obj.getCountSpeciesMsx
                        value.Quality{i}(k,:)=obj.getSpeciesConcentration(1, i, j);%node code0
                    end
                    k=k+1;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%% OPERATIONS %%%%%%%%%%%%%%%%%%%
        
        %MSXclose & MSXMatlabCleanup
        function unloadMsx(varargin)
            MSXclose;
            MSXMatlabCleanup;
        end
        
        %MSXsaveoutfile
        function saveQualityFileMsx(obj,outfname)
            [obj.errorCode]=MSXsaveoutfile(outfname);
        end
        
        %MSXusehydfile
        function useHydraulicFileMsx(obj,hydname)
            [obj.errorCode]=MSXusehydfile(hydname);
        end
        
        %MSXinit 
        function initializeQualityAnalysisMsx(obj,flag)
            [obj.errorCode] = MSXinit(flag);
        end
        
        %MSXstep
        function [t, tleft]=stepQualityAnalysisTimeLeftMsx(obj)
            [obj.errorCode, t, tleft] = MSXstep();
        end
        
        %MSXgetqual
        function value=getSpeciesConcentration(obj, type, index, species)
            [obj.errorCode, value] = MSXgetqual(type, index, species);
        end
        
        %MSXsavemsxfile
        function saveMsxFile(obj,msxname)
            [obj.errorCode] = MSXsavemsxfile(msxname);
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
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
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


% function [errcode,from,to] = ENgetalllinknodes()
%     global EN_SIZE;
%     fval1=int32(0);
%     fval2=int32(0);
%     p1=libpointer('int32Ptr',fval1);
%     p2=libpointer('int32Ptr',fval2);
%     from=int32(zeros(EN_SIZE.nlinks,1));
%     to=int32(zeros(EN_SIZE.nlinks,1));
%     for i=1:EN_SIZE.nlinks
%         [errcode]=calllib('epanet2','ENgetlinknodes',i,p1,p2);
%         if errcode 
%             ENerror(errcode); 
%         end
%         from(i)=get(p1,'Value');
%         to(i)=get(p2,'Value');
%     end
% end


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
    len=80;
    [e,errmsg] = calllib('epanet2','ENgeterror',errcode,errmsg,len);
    if e 
        ENerror(e); 
    end
end

function [errcode,flowUnits] = ENgetflowunits()
    [errcode, flowUnits]=calllib('epanet2','ENgetflowunits',0);
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


% function [nnodes,ntanks,nlinks,npats,ncurves,ncontrols,errcode] = ENgetnetsize()
%     [errcode,nnodes]=calllib('epanet2','ENgetcount',0,0);
%     if errcode 
%         ENerror(errcode); 
%     end
%     [errcode,ntanks]=calllib('epanet2','ENgetcount',1,0);
%     if errcode 
%         ENerror(errcode); 
%     end
%     [errcode,nlinks]=calllib('epanet2','ENgetcount',2,0);
%     if errcode 
%         ENerror(errcode); 
%     end
%     [errcode,npats]=calllib('epanet2','ENgetcount',3,0);
%     if errcode 
%         ENerror(errcode); 
%     end
%     [errcode,ncurves]=calllib('epanet2','ENgetcount',4,0);
%     if errcode 
%         ENerror(errcode); 
%     end
% 
%     [errcode,ncontrols]=calllib('epanet2','ENgetcount',5,0);
%     if errcode 
%         ENerror(errcode); 
%     end
% end

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
        ENDLLNAME='epanet2';
%     end;
    % Load library
    if libisloaded(ENDLLNAME)
        unloadlibrary(ENDLLNAME);
    else
        errstring =['Library ', ENDLLNAME, '.dll was not loaded'];
        disp(errstring);
    end;
end

function [errcode] = ENMatlabSetup(DLLname,Hname)
    %currentversion = 20012;
    % Load library
    ENDLLNAME=DLLname;
    ENHNAME=Hname;
    if ~libisloaded(ENDLLNAME)
%         unloadlibrary(ENDLLNAME);
        loadlibrary(ENDLLNAME,ENHNAME);
    end
    % Check version of EPANET DLL
    [errcode, version] = ENgetversion();
    %if version ~= currentversion
    %    errcode = 1;
    versionString = ['Current version of EPANET:',num2str(version)];
    disp(versionString);
    if errcode 
        ENerror(errcode); 
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
%     global EN_SIZE;
    
    repname='';
    binname='';
    
    errcode=calllib('epanet2','ENopen',inpname,repname,binname);
%     if errcode 
        while errcode~=0
            try
                errcode=calllib('epanet2','ENopen',inpname,repname,binname);
            catch err
            end
        end
%         ENerror(errcode); 
%     end
%     [nnodes,ntanks,nlinks,npats,ncurves,ncontrols,errcode] = ENgetnetsize();
%     if errcode 
%         ENerror(errcode); 
%     end
%     delete(repname,binname);
%     EN_SIZE = struct(...
%         'nnodes', nnodes,...
%         'ntanks', ntanks,...
%         'nlinks', nlinks,...
%         'npats',  npats,...
%         'ncurves',ncurves,...
%         'ncontrols',ncontrols);
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
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end


function [errcode] = ENsetlinkvalue(index, paramcode, value)
    [errcode]=calllib('epanet2','ENsetlinkvalue',index, paramcode, value);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetnodevalue(index, paramcode, value)
    [errcode]=calllib('epanet2','ENsetnodevalue',index, paramcode, value);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetoption(optioncode,value)
    [errcode]=calllib('epanet2','ENsetoption',optioncode,value);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetpattern(index, factors, nfactors)
    [errcode]=calllib('epanet2','ENsetpattern',index,factors,nfactors);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetpatternvalue(index, period, value)
    [errcode]=calllib('epanet2','ENsetpatternvalue',index, period, value);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end


function [errcode] = ENsetqualtype(qualcode,chemname,chemunits,tracenode)
    [errcode]=calllib('epanet2','ENsetqualtype',qualcode,chemname,chemunits,tracenode);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetreport(command)
    [errcode]=calllib('epanet2','ENsetreport',command);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsetstatusreport(statuslevel)
    [errcode]=calllib('epanet2','ENsetstatusreport',statuslevel);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [errcode] = ENsettimeparam(paramcode, timevalue)
    [errcode]=calllib('epanet2','ENsettimeparam',paramcode,timevalue);
    if errcode 
        ENerror(errcode); 
    end
    tmp=1;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
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


function [errcode, value] = ENgetpumptype(index)
    [errcode, value]=calllib('epanet2','ENgetpumptype',index, 0);
    if errcode 
        ENerror(errcode); 
    end
end


function [errcode, value] = ENgetheadcurve(index)
    [errcode, value]=calllib('epanet2','ENgetheadcurve',index, '');
    if errcode 
        ENerror(errcode); 
    end
end


% function [errcode, nValues, xValues, yValues] = ENgetcurve(curveIndex)
%     global 'epanet2';    
%     xValues=single(0);
%     yValues=single(0);
%     nValues=int32(0);
%     [errcode, nValues, xValues, yValues] =	calllib('epanet2','ENgetcurve', curveIndex, nValues, xValues, yValues);
% end


% ENplot 

function CoordinatesXY=ENplot(obj,varargin)

    % Initiality
    highlightnode=0;
    highlightlink=0;
    highlightnodeindex=[];
    highlightlinkindex=[];
    Node=char('no');
    Link=char('no');
    fontsize=10;
    
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
            otherwise
                warning('Invalid property found.');
                return
        end
    end

    cla
    % Get node names and x, y coordiantes
    CoordinatesXY = obj.getCoordinates;
    nodes = obj.getNodesInfo;
    if isa(highlightnode,'cell')       
        for i=1:length(highlightnode)
            n = strcmp(obj.getNodeID,highlightnode{i});
            if sum(n)==0
                warning('Undefined node with id "%s" in function call therefore the index is zero.', char(highlightnode{i})); 
            else
                highlightnodeindex(i) = strfind(n,1);
            end
        end
    end

    if isa(highlightlink,'cell') 
        for i=1:length(highlightlink)
            n = strcmp(obj.getLinkID,highlightlink{i});
            if sum(n)==0
                warning('Undefined link with id "%s" in function call therefore the index is zero.', char(highlightlink{i})); 
            else
                highlightlinkindex(i) = strfind(n,1);
            end
        end
    end

    % Coordinates for node FROM
    for i=1:nodes.CountNodes
        [x] = double(CoordinatesXY{1}(i));
        [y] = double(CoordinatesXY{2}(i));

        hh=strfind(highlightnodeindex,i);
        h(:,1)=plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','b',...
                      'MarkerFaceColor','b',...
                      'MarkerSize',5);
        legendString{1}= char('Junctions');

        % Plot Reservoirs
        if sum(strfind(nodes.NodeReservoirIndex,i))
            colornode = 'g';
            if length(hh)
                colornode = 'r';
            end
            h(:,2)=plot(x,y,'s','LineWidth',2,'MarkerEdgeColor','r',...
                      'MarkerFaceColor','g',...
                      'MarkerSize',13);
            plot(x,y,'s','LineWidth',2,'MarkerEdgeColor','r',...
                      'MarkerFaceColor',colornode,...
                      'MarkerSize',13);
                  
           legendString{2} = char('Reservoirs');
        end
        % Plot Tanks
        if sum(strfind(nodes.NodeTankIndex,i)) 
            colornode = 'k';
            if length(hh)
                colornode = 'r';
            end
            h(:,3)=plot(x,y,'p','LineWidth',2,'MarkerEdgeColor','r',...
              'MarkerFaceColor','k',...
              'MarkerSize',16);
          
            plot(x,y,'p','LineWidth',2,'MarkerEdgeColor','r',...
                      'MarkerFaceColor',colornode,...
                      'MarkerSize',16);

            legendString{3} = char('Tanks');
        end

        % Show Node id
        if (strcmp(lower(Node),'yes') && ~length(hh))
            text(x,y,nodes.NodesAll(i),'Fontsize',fontsize);%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
        end

        if length(hh) 
            plot(x, y,'o','LineWidth',2,'MarkerEdgeColor','r',...
                      'MarkerFaceColor','r',...
                      'MarkerSize',10)

            text(x,y,nodes.NodesAll(i),'Fontsize',fontsize)%'BackgroundColor',[.7 .9 .7],'Margin',margin/4);
        end
        hold on
    end
    links = obj.getLinksInfo;

    for i=1:links.CountLinks
        FromNode=strfind(strcmp(links.FromNode{i},nodes.NodesAll),1);
        ToNode=strfind(strcmp(links.ToNode{i},nodes.NodesAll),1);
        
        if FromNode
            x1 = double(CoordinatesXY{1}(FromNode));
            y1 = double(CoordinatesXY{2}(FromNode));
        end
        if ToNode
            x2 = double(CoordinatesXY{1}(ToNode));
            y2 = double(CoordinatesXY{2}(ToNode));
        end
        
        hh=strfind(highlightlinkindex,i);

%         h(:,4)=line([x1,x2],[y1,y2],'LineWidth',1);
        h(:,4)=line([x1 CoordinatesXY{3}{i} x2],[y1 CoordinatesXY{4}{i} y2],'LineWidth',1);
        
        legendString{4} = char('Pipes');
        % Plot Pumps
        if sum(strfind(links.LinkPumpIndex,i)) 
            colornode = 'b';
            if length(hh)
                colornode = 'r';
            end
            h(:,5)=plot((x1+x2)/2,(y1+y2)/2,'bv','LineWidth',2,'MarkerEdgeColor','b',...
                      'MarkerFaceColor','b',...
                      'MarkerSize',5);
            plot((x1+x2)/2,(y1+y2)/2,'bv','LineWidth',2,'MarkerEdgeColor',colornode,...
                      'MarkerFaceColor',colornode,...
                      'MarkerSize',5);
                  
           legendString{5} = char('Pumps');
        end

        % Plot Valves
        if sum(strfind(links.LinkValveIndex,i)) 
            h(:,6)=plot((x1+x2)/2,(y1+y2)/2,'b*','LineWidth',2,'MarkerEdgeColor','b',...
                      'MarkerFaceColor','b',...
                      'MarkerSize',7);
            legendString{6} = char('Valves');
        end

        % Show Link id
        if (strcmp(lower(Link),'yes') && ~length(hh))
            text((x1+x2)/2,(y1+y2)/2,links.LinksAll(i),'Fontsize',fontsize);
        end

        if length(hh) 
            line([x1,x2],[y1,y2],'LineWidth',2,'Color','g');
            text((x1+x2)/2,(y1+y2)/2,links.LinksAll(i),'Fontsize',fontsize);
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
    [xmax,~]=max(CoordinatesXY{1});
    [xmin,~]=min(CoordinatesXY{1});
    [ymax,~]=max(CoordinatesXY{2});
    [ymin,~]=min(CoordinatesXY{2});

%     xmax=yxmax(1); ymax=yxmax(2);
%     xmin=yxmin(1); ymin=yxmin(2);
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
    axis off
    whitebg('w');
end


function [vx,vy,vertx,verty] = getNodeCoord(obj)
    % Initialize 
    nodes=obj.getNodesInfo;
    links=obj.getLinksInfo;
    
    vx = NaN(nodes.CountNodes,1);
    vy = NaN(nodes.CountNodes,1);
    vertx = cell(links.CountLinks,1);
    verty = cell(links.CountLinks,1);
    nvert = zeros(links.CountLinks,1);

    % Open epanet input file
    [fid,message]=fopen(obj.PathFile,'rt');
    if fid < 0
        disp(message)
        return
    end

    sect = 0;i=1;t=1;
    % Read each line from input file.
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end

        % Get first token in the line
        tok = strtok(tline);

        % Skip blank lines and comments
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
%             % get the node index
%             [errcode,index] = ENgetnodeindex(char(A{1}));
%             if errcode ~=0 
%                 return; 
%             end
            vx(i) = A{2};
            vy(i) = A{3};
            i=i+1;
        % Vertices
        elseif sect == 2
            A = textscan(tline,'%s %f %f');
%             [errcode,index] = ENgetlinkindex(char(A{1}));
%             if errcode ~=0 
%                 return; 
%             end
            nvert(t) = nvert(t) + 1;
            vertx{t}(nvert(t)) = A{2};
            verty{t}(nvert(t)) = A{3};
            t=t+1;
        end
    end
end


%%%%%%%%%%%%%%% EPANET - MSX %%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj] = MSXMatlabSetup(obj,msxname)
    if ~libisloaded('epanetmsx')
       loadlibrary('epanetmsx','epanetmsx.h');
    end
    [obj.errorCode, obj.MsxPathFile, obj.MsxFile] = MSXopen(msxname);
    %%%%%%%%%%%%%%%%%%%% EPANET - MSX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [obj.TermsFormula,obj.PipesFormula,obj.TanksFormula] = GetFormulas(msxname);

    %MSXgetcount  
    [obj.errorCode, obj.CountSpeciesMsx] = MSXgetcount(3);
    [obj.errorCode, obj.CountConstantsMsx] = MSXgetcount(6);
    [obj.errorCode, obj.CountParametersMsx] = MSXgetcount(5);
    [obj.errorCode, obj.CountPatternsMsx] = MSXgetcount(7);

    obj.SpeciesIndexMsx=1:obj.CountSpeciesMsx;
    for i=1:obj.CountSpeciesMsx
        %MSXgetIDlen
        [obj.errorCode,len] = MSXgetIDlen(3,i);
        %MSXgetID 
        [obj.errorCode,obj.SpeciesNameIDMsx{i}] = MSXgetID(3,i,len);
        %MSXgetspecies
        [obj.errorCode,obj.SpeciesTypeMsx{i},obj.SpeciesUnitsMsx{i},obj.SpeciesAtolMsx(i),obj.SpeciesRtolMsx(i)] = MSXgetspecies(i);
    end

    %MSXgetconstant
    for i=1:obj.CountConstantsMsx
        [obj.errorCode, len] = MSXgetIDlen(6,i);
        [obj.errorCode, obj.ConstantNameIDMsx{i}] = MSXgetID(6,i,len);
        [obj.errorCode, obj.ConstantValueMsx(i)] = MSXgetconstant(i);
    end

    % MSXgetparameter
    for i=1:obj.CountParametersMsx
        [obj.errorCode, len] = MSXgetIDlen(5,i);
        [obj.errorCode, obj.ParameterNameIDMsx{i}] = MSXgetID(5,i,len); 
        %MSXgetindex
        [obj.errorCode, obj.ParametersIndexMsx(i)] = MSXgetindex(obj,5,obj.ParameterNameIDMsx{i}); 
    end
    for i=1:length(obj.CountTanks)
        for j=1:obj.CountParametersMsx
           [obj.errorCode, obj.TankParameterValueMsx{i}(j)] = MSXgetparameter(0,obj.CountTanks(i),j);   
        end
    end
    for i=1:obj.CountPipes
        for j=1:obj.CountParametersMsx
           [obj.errorCode, obj.PipeParameterValueMsx{i}(j)] = MSXgetparameter(1,i,j);   
        end
    end

    %MSXgetpatternlen
    for i=1:obj.CountPatternsMsx
        [obj.errorCode, len] = MSXgetIDlen(7,i);
        [obj.errorCode, obj.PatternIDMsx{i}] = MSXgetID(7,i,len);
        [obj.errorCode, obj.PatternIndexMsx(i)] = MSXgetindex(obj,7,obj.PatternIDMsx{i});
        [obj.errorCode, obj.PatternLengthsMsx(i)] = MSXgetpatternlen(i);
    end

    %MSXgetinitqual
    for i=1:obj.CountNodes
        for j=1:obj.CountSpeciesMsx
           [obj.errorCode, obj.NodeInitqualValueMsx{i}(j)] = MSXgetinitqual(0,i,j);   
        end
    end
    for i=1:obj.CountLinks
        for j=1:obj.CountSpeciesMsx
           [obj.errorCode, obj.LinkInitqualValueMsx{i}(j)] = MSXgetinitqual(1,i,j);   
        end
    end

    %MSXgetsource
    for i=1:obj.CountNodes
        for j=1:obj.CountSpeciesMsx 
           [obj.errorCode, obj.SourceTypeMsx{i}{j},obj.SourceLevelMsx{i}(j),obj.SourcePatternIndexMsx{i}(j)] = MSXgetsource(i,j);
           [obj.errorCode, len] = MSXgetIDlen(7,j);
           [obj.errorCode,obj.SourcePatternIDMsx{i}{j}] = MSXgetID(7,obj.SourcePatternIndexMsx{i}(j),len);
           [obj.errorCode, len] = MSXgetIDlen(3,j);
%                    [obj.errorCode,obj.SourceSpeciesNameIDMsx{i}{j}] = MSXgetID(3,j,len);
        end
           obj.SourceNodeIDMsx{i} = obj.NodeNameID(i);
    end
    obj.SourceAllMsx={obj.SourceTypeMsx,obj.SourceLevelMsx,obj.SourcePatternIndexMsx,obj.SourcePatternIDMsx,obj.SourceNodeIDMsx};        
end

function [errcode, pathfile, msxname] = MSXopen(msxname)
    pathfile = which(msxname);
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

function [errcode, index] = MSXgetindex(obj,varargin)
    index =0;
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

%%%%%%%
function [errcode] = MSXusehydfile(hydfname)
    [errcode] = ENsaveH();
    [errcode] = ENsavehydfile('test.hyd');
    [errcode]=calllib('epanetmsx','MSXusehydfile',hydfname);
    if errcode 
        MSXerror(errcode); 
    end
end

function [errcode, t, tleft] = MSXstep()
    t=0;
    tleft=0;
    [errcode,t,tleft]=calllib('epanetmsx','MSXstep',t,tleft);
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
    len=80;
    [e,errmsg] = calllib('epanetmsx','MSXgeterror',errcode,errmsg,len);
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


function  [Terms,Pipes,Tanks] = GetFormulas(msxname)
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

        % Skip blank lines and comments
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
        
    end
end
        
        
function [curvesID,X,Y,sectCurve]=CurveInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.PathFile,'rt');
    if fid < 0
        disp(message)
        return
    end

    curvesID={};
    X={};
    Y={};
    sect=0; i=1; sectCurve=0;
    % Read each line from msx file.
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank lines and comments
        if isempty(tok), continue, end
        if (tok(1) == ';'), continue, end
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
            a = textscan(tline,'%s %f %f');
            curvesID(i)=a{1};
            X(i)=a(2);
            Y(i)=a(3);
            i=i+1;
        end
    end          
end


function addCurve(obj,newCurveID,CurveX,CurveY)  
    load([pwd,'\Results\','tmpInp'],'tmp','-mat') 
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    % Check if new ID already exists
    [pCurveID,~,~,sectCurve]=CurveInfo(obj);
    sect=0;
    
    i=1;exists=0;
    while i<length(pCurveID)+1
        exists(i) = strcmp(newCurveID,char(pCurveID(i)));
        i=i+1;
    end

    if sum(exists)>0
        s = sprintf('Curve "%s" already exists.',newCurveID);
        warning(s);
        return
    end

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=1;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        info{t} = tline;
        t = t+1;
    end
    
    % write
    fid2 = fopen(obj.PathFile,'w');

    nn=0;yy=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

        u=1;
            while u < length(a)+1
                if strcmp(a{u},'[CURVES]') && sectCurve==0
                    fprintf(fid2,'[CURVES]');
                    sect=1; break;
                end

                spaces= 15 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                rr = regexp(a,'\w*[\w*]\w*','split');
                check_brackets = rr{:};
                ch1 = strcmp(check_brackets,'[');
                ch2 = strcmp(check_brackets,']');

                if (ch1(1)==1 && ch2(2)==1 && (sectCurve~=0 && nn==0))
                    if yy==0
                        if sect==0
                            fprintf(fid2,'[CURVES]\n;ID                X-Value            Y-Value\n');
                        end
                        fprintf(fid2,';PUMP:%sX-Value%sY-Value\n',sps{:},sps{:}); yy=1;
                    end
                    fprintf(fid2, '%s%s%d%s%d', newCurveID,sps{:},CurveX,sps{:},CurveY);
                    fprintf(fid2,'\r\n');
                    fprintf(fid2,'\r\n');
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2,'\r\n');
                    nn=1;

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
    
    [errcode] = Epanet('tempInpFile.inp');
    obj.saveInputFile('tempInpFile.inp'); % OK %
    movefile('tempInpFile.inp',[pwd,'\Results']);
end


function rmCurveID(obj,CurveID)

    obj.saveInputFile('tempInpFile.inp'); % OK %
    movefile('tempInpFile.inp',[pwd,'\Results']);
    
    % Check if id new already exists
    [pCurveID,~,~,sectCurve]=CurveInfo(obj);
    if length(pCurveID)==0 
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    i=1;
    while i<length(pCurveID)+1
        exists(i) = strcmp(CurveID,char(pCurveID(i)));
        i=i+1;
    end

    if (sum(exists)~=1)
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    value=obj.getLinksInfo;
    pp=obj.getPumpInfo;
    i=1;
    while (i<length(pp.PumpCurveID)+1)
        p = strcmp(CurveID,pp.PumpCurveID{i});

        if p==1
%             s = sprintf('Pump %s refers to undefined curve.',value.PumpsID{i});
            obj.removeLinkID(value.PumpsID{i})
%             warning(s);
        end

        i=i+1;
    end

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=1;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        info{t} = tline;
        t = t+1;
    end
    
    fid2 = fopen(obj.PathFile,'w');

    i=1;e=0;n=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

        u=1; 
        while u < length(a)+1

            spaces= 10 - length(a{u});
            k=1;sps={' '};
            while k < spaces+1
                sps = strcat(sps,{' '});
                k=k+1;
            end

            rr = regexp(a,'\w*[\w*]\w*','split');
            check_brackets = rr{:};
            ch1 = strcmp(check_brackets,'[');
            ch2 = strcmp(check_brackets,']');

            if strcmp(a{u},'[CURVES]')
                fprintf(fid2,'%s',a{u}); 
                n=1; 
            elseif ch1(1)==1 && ch2(2)==1 && n==1
                if (isempty(a{u})&& n==1) break; end
                e=1;
            end
            if strcmp(a{u},'[END]')  e=1; fprintf(fid2,'%s',a{u});break;   end

            if n==1 && e==0  
                if strcmp(a{u},'[CURVES]') break; end
                if isempty(a{u})                
                elseif strfind(a{u},';')
                    u = length(a)+1; break;
                else 
                        tt = strcmp(a{u},CurveID);  
                        if tt==1
                            u = length(a)+1; break;
                        else
                        fprintf(fid2,'%s%s',a{u},sps{:});
                        end
                end
            else
                if isempty(a{u})
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
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
%     [errcode] = Epanet('tempInpFile.inp');
%     obj.saveInputFile('tempInpFile.inp'); % OK %
%     movefile('tempInpFile.inp',[pwd,'\Results']);
end


function addLink(obj,typecode,newLink,fromNode,toNode,curveID,varargin)
    % Link type codes consist of the following constants:  
    % CVPIPE   0   Pipe with Check Valve 
    % PIPE     1   Pipe 
    % PUMP     2   Pump 
    % PRV      3   Pressure Reducing Valve 
    % PSV      4   Pressure Sustaining Valve 
    % PBV      5   Pressure Breaker Valve 
    % FCV      6   Flow Control Valve 
    % TCV      7   Throttle Control Valve 
    % GPV      8   General Purpose Valve    
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    
    % Initial PIPE
    % plength     value for length of new pipe
    % pdiameter   value for diameter of new pipe
    % proughness  value for roughness of new pipe 
    v=obj.getFlowUnitsHeadlossFormula;
    if v.UScustomary==1 
        plength=1000; %(ft)
        pdiameter=12; %(in)
        vdiameter=12; %valves
    else %SI metric  
        plength=304.8; %(m)
        pdiameter=304.8; %(mm)
        vdiameter=304.8; %valves
    end
    proughness=100;
    vsetting=0; %valves

    % Check if id new already exists
    Nodes = obj.getNodesInfo;
    if length(Nodes.NodesAll)==0 
        return
    end
    
    i=1;existsFrom=0;existsTo=0;
    while i<length(Nodes.NodesAll)+1
        existsFrom(i) = strcmp(fromNode,char(Nodes.NodesAll{i}));
        existsTo(i) = strcmp(toNode,char(Nodes.NodesAll{i}));
        i=i+1; 
    end
    if sum(existsFrom)~=1
        s = sprintf('There is no node "%s" in the network.',fromNode);
        warning(s);
        return
    elseif sum(existsTo)~=1
        s = sprintf('There is no node "%s" in the network.',toNode);
        warning(s);
        return
    end
    
    A = [0 1 2 3 4 5 6 7 8];
    code = strfind(A,typecode);
    if length(code)==0 
        warning('There is no such typecode(0-8)');
        return
    else 
        if typecode==0 type_valv = 'CVPIPE';  end
        if typecode==3 type_valv = 'PRV';     end
        if typecode==4 type_valv = 'PSV';     end
        if typecode==5 type_valv = 'PBV';     end
        if typecode==6 type_valv = 'FCV';     end
        if typecode==7 type_valv = 'TCV';     end
        if typecode==8 type_valv = 'GPV';     end
        if typecode~=1 && typecode~=2
            typecode=3;
        end
    end

    % Valve illegally connected to a tank or reservoir
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if typecode==3
        i=1;ifFromReservoir=0;ifToReservoir=0;
        while i<length(Nodes.ReservoirsID)+1
            ifFromReservoir(i) = strcmp(fromNode,char(Nodes.ReservoirsID{i}));
            ifToReservoir(i) = strcmp(toNode,char(Nodes.ReservoirsID{i}));
            i=i+1; 
        end
        i=1;ifFromTank=0;ifToTank=0;
        while i<length(Nodes.TanksID)+1
            ifFromTank(i) = strcmp(fromNode,char(Nodes.TanksID{i}));
            ifToTank(i) = strcmp(toNode,char(Nodes.TanksID{i}));
            i=i+1; 
        end
        if sum(ifFromReservoir)==1 || sum(ifFromTank)==1 || sum(ifToReservoir)==1 || sum(ifToTank)==1
            s = sprintf('Valve "%s" illegally connected to a tank.',newLink);
            warning(s);
            return
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    crvs = obj.getCurveInfo;
    % Check if newLink already exists
    Links = obj.getLinksInfo;

    i=1;
    while i<length(Links.LinksAll)+1
        exists_link = strcmp(newLink,char(Links.LinksAll(i)));
        i=i+1;

        if exists_link==1
            s = sprintf('Link %s already exists.',newLink);
            warning(s);
            return
        end

    end

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=1;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        info{t} = tline;
        t = t+1;
    end
    fid2 = fopen(obj.PathFile,'w');

    % Add pipe
    i=1;nn=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;
            while u < length(a)+1

                spaces= 18 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                t =  regexp(a{u}, '[(\w*)]','split');
                y=1;cnt=0;
                while y<length(t)+1
                     tt = isempty(t{y});
                     if tt==0
                         cnt=cnt+1;
                     end
                     y=y+1;
                end

                if (cnt==2 && strcmp(a{u},'[PIPES]') && nn==0 && typecode==1)
                    if Links.sectPipes==0  
                        fprintf(fid2,'[PIPES]');
                        fprintf(fid2,'\r\n');
                    end
 
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%d%s%d%s%d;',newLink,sps{:},fromNode,sps{:},...
                    toNode,sps{:},plength,sps{:},pdiameter,sps{:},proughness);

                elseif (cnt==2 && strcmp(a{u},'[PUMPS]') && nn==0 && typecode==2)
                    if Links.sectPumps==0  
                        fprintf(fid2,'[PIPES]');
                        fprintf(fid2,'\r\n');
                    end
 
                    if isempty(char(crvs.CurvesID))
                        s = sprintf('No head curve supplied for pump %s.',newLink);
                        warning(s);
                        warning('addNewCurve must be called after this function.');
                        fclose(fid2);
                        return
                    else
    %                     curve=input('Please enter the ID of curve:'); 
                          curve=curveID;
                    end
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%s%s%s;',newLink,sps{:},fromNode,sps{:},...
                    toNode,sps{:},'HEAD',sps{:},curve);
                
                elseif typecode==3 && strcmp(a{u},'[VALVES]')
                    if Links.sectValves==0  
                        fprintf(fid2,'[VALVES]');
                        fprintf(fid2,'\r\n');
                    end
                    fprintf(fid2,'%s',a{u});
                    fprintf(fid2, '\n%s%s%s%s%s%s%d%s%s%s%d;',newLink,sps{:},fromNode,sps{:},...
                    toNode,sps{:},vdiameter,sps{:},type_valv,sps{:},vsetting);
                    nn=1;
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
    
    [errcode] = Epanet('tempInpFile.inp');
    obj.saveInputFile('tempInpFile.inp'); % OK %
    movefile('tempInpFile.inp',[pwd,'\Results']);
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')

end
 
function value=LinksInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.PathFile,'rt');
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
    value.CountLinks=0;

    value.sectPipes=0;
    value.sectPumps=0;
    value.sectValves=0;
    sect=0; i=1;t=1;q=1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank lines and comments
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
         value.CountLinks= value.CountLinks+1;
    end
    value.LinksAll=[value.PipesID value.PumpsID value.ValvesID];
end

function value=NodesInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.PathFile,'rt');
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
    value.CountNodes=0;
    
    value.sectTanks=0;
    value.sectJunctions=0;
    value.sectReservoirs=0;
    sect=0; t=1;i=1;q=1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank lines and comments
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
        value.CountNodes=value.CountNodes+1;
    end

    value.NodesAll=[value.JunctionsID value.ReservoirsID value.TanksID];
end

function value=ControlsInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.PathFile,'rt');
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
        % Skip blank lines and comments
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
    [fid,message] = fopen(obj.PathFile,'rt');
    if fid < 0
        disp(message)
        return
    end
    value.FlowUnits={};
    value.Headloss={};
    value.sectOptions=0;
    value.SImetric=0;
    value.UScustomary=0;
    value.ConcentrationUnits={};
    
    sect=0;i=1;u=1;t=1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank lines and comments
        if isempty(tok), continue, end
        if (tok(1) == ';'), continue, end
        if (tok(1) == '[')
            % [CONTROLS] section
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
            if strcmp(upper(atline{1}),'UNITS')
                value.FlowUnits{i}=atline{2};
                u=u+1;
            end
            if strcmp(upper(atline{1}),'HEADLOSS')
                value.Headloss{t}=atline{2};
                t=t+1;
            end
            if strcmp(upper(atline{1}),'PRESSURE')
                value.PressureUnits=atline{2};
            end
            if strcmp(upper(atline{1}),'QUALITY')
                if ~strcmp(upper(atline{2}),'NONE')
                    value.ConcentrationUnits=atline{3};
                else
                    value.ConcentrationUnits='NONE';
                end
            end
            i=i+1;
        end

    end
    
%     US Customary - SI metric 
    switch char(value.FlowUnits)
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
        value.DemandsUnits=value.FlowUnits;
        value.DiameterPipesUnits='inches';
        value.DiameterTanksUnits='feet';
        value.EfficiencyUnits='percent'; 
        value.ElevationUnits='feet';        
        value.EmitterCoeffUnits='flow units @ 1 psi drop';
        value.EnergyUnits='kwatt-hours'; 
        value.FrictionFactorUnits='unitless'; 
        value.HeadUnits='feet';
        value.LengthUnits='feet';
        value.MinorLossCoeffUnits='unitless'; 
        value.PowerUnits='horsepower';
        value.ReactionCoeffBulkUnits='1/day (1st-order)'; 
        value.ReactionCoeffWallUnits='mass/sq-ft/day (0-order), ft/day (1st-order)'; 
        value.RoughnessCoeffUnits='millifeet(Darcy-Weisbach), unitless otherwise';
        value.SourceMassInjectionUnits='mass/minute';
        value.VelocityUnits='ft/sec';
        value.VolumeUnits='cubic feet';
        value.WaterAgeUnits='hours'; 
    else % SI Metric
        value.DemandsUnits=value.FlowUnits;
        value.DiameterPipesUnits='millimeters';
        value.DiameterTanksUnits='meters';
        value.EfficiencyUnits='percent'; 
        value.ElevationUnits='meters';        
        value.EmitterCoeffUnits='flow units @ 1 meter drop';
        value.EnergyUnits='kwatt-hours'; 
        value.FrictionFactorUnits='unitless'; 
        value.HeadUnits='meters';
        value.LengthUnits='meters';
        value.MinorLossCoeffUnits='unitless'; 
        value.PowerUnits='kwatts';
        value.ReactionCoeffBulkUnits='1/day (1st-order)'; 
        value.ReactionCoeffWallUnits='mass/sq-m/day(0-order), meters/day (1st-order)'; 
        value.RoughnessCoeffUnits='mm(Darcy-Weisbach), unitless otherwise';
        value.SourceMassInjectionUnits='mass/minute';
        value.VelocityUnits='meters/sec';
        value.VolumeUnits='cubic meters';
        value.WaterAgeUnits='hours';        
    end
end    

function value=PumpInfo(obj)
    % Open epanet input file
    [fid,message] = fopen(obj.PathFile,'rt');
    if fid < 0
        disp(message)
        return
    end

    value.PumpsID={};
    value.fromNodePumps={};
    value.toNodePumps={};
    sect=0; i=1; value.sectPump=0;
    % Read each line from msx file.
    while 1
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        % Get first token in the line
        tok = strtok(tline);
        % Skip blank lines and comments
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


function addNode(obj,typecode,newID,X,Y,varargin)
    % Addnode - Add node in the network.
    % Node type codes consist of the following constants:
    % EN_JUNCTION	0  Junction node
    % EN_RESERVOIR	1  Reservoir node
    % EN_TANK       2  Tank node

    % elevation value of elevation for new node
    % initqual  value for InitQual of new node(0.5,1..)
    
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    v=obj.getFlowUnitsHeadlossFormula;

    if typecode==1 || typecode==0 % junction & reservoir
        if v.UScustomary==1  %(ft)
            elevation=500;  
        else %SI metric   %(m)
            elevation=152.4;
        end        
        initqual=0;
    end

    if typecode==2
        % Initial TANK
        if v.UScustomary==1  %(ft)
            MaxLevel=20;
            Diameter=50;
            Initlevel=10;
            elevation=500;
        else %SI metric   %(m)
            MaxLevel=6.0960;
            Diameter=15.24;
            Initlevel=3.048;
            elevation=152.4;
        end
        initqual=0;
        MinLevel=0;
        MinVol=0;
    end
    
    % Check if id new already exists
    nodes = obj.getNodesInfo;
    if length(nodes.NodesAll)==0
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    i=1;
    while i<length(nodes.NodesAll)+1
        exists(i) = strcmp(newID,char(nodes.NodesAll(i)));
        i=i+1;
    end

    if (sum(exists)==1)
        s = sprintf('Node "%s" already exists.',newID);
        warning(s);
        return
    end

    % Get type of node
    if typecode==0 type_new = '[JUNCTIONS]'; end
    if typecode==1 type_new = '[RESERVOIRS]'; end
    if typecode==2 type_new = '[TANKS]'; end

    A = [0 1 2];
    code=strfind(A,typecode);

    if length(code)==0
        warning('There is no such typecode(0-2)!');
        return
    end

    % check section in inpname, [JUNCTIONS], [RESERVOIRS], [TANKS]
    stank_check=1;
    sreservoir_check=1;
    sjunction_check=1;
    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=1;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
%         a=regexp(tline,'\s*','split');
%         for i=1:length(a)
%             if isempty(a{i})
%             else
                info{t} = tline;
%                 t=t+1;
%             end
%         end
        t=t+1;
    end
    fid2 = fopen(obj.PathFile,'w');

    % Initiality
    i=1;qualch=0;
    qq=0;
    Coordch=0;
    onetime=1;
    gg=0;

    for t = 1:length(info)
        c = info{t};
        if ~isempty(c)
            a = regexp(c, '\s*','split');
        end
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;
            while u < length(a)+1

                % Spaces
                spaces1= 18 - length(a{u});
                k=1;sps1={' '};
                while k < spaces1+1
                    sps1 = strcat(sps1,{' '});
                    k=k+1;
                end

                % Find [brackets] cnt=2;
                t =  regexp(a{u}, '[(\w*)]','split');
                y=1;cnt=0;
                while y<length(t)+1
                    tt = isempty(t{y});
                    if tt==0
                        cnt=cnt+1;
                    end
                    y=y+1;
                end

                %%%%%%%% Quality Section %%%%%%%%
                if strcmp(a{u},'[QUALITY]')
                    fprintf(fid2,'[QUALITY]');
                    qualch=1;
                    break;
                end

                if (cnt==2 && qualch==1)
                    fprintf(fid2, '%s%s%d;', newID,sps1{:},initqual);
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
                        fprintf(fid2, '%s%s%d%s%d\n', newID,sps1{:},X,sps1{:},Y);
                    end                
                    gg=gg+1;
                end

%                 if isempty(obj.CoordinatesXY)
%                     if strcmp(a{u},'[END]')
%                         fprintf(fid2,'%s','[COORDINATES]');
%                         fprintf(fid2,'\r\n');
% 
%                         for qq=1:length(X)
%                             fprintf(fid2,'%s%s%d%s%d',char(NodesID(qq)),sps1{:},X(qq),sps1{:},Y(qq));
%                             fprintf(fid2,'\r\n');
%                         end
% 
%                         fprintf(fid2, '%s%s%d%s%d\n', newID,sps1{:},X,sps1{:},Y);
%                         fprintf(fid2,'%s',a{u});
%                         fprintf(fid2,'\r\n');
%                         fclose all;
%                         [errcode] = Epanet('tempInpFile.inp');
%                         obj.saveInputFile('tempInpFile.inp'); % OK %
%                         movefile('tempInpFile.inp',[pwd,'\Results']);
%                         return
%                     end
%                 end

                %%%%%%%% Nodes Section %%%%%%%%
                if (cnt==2 && (strcmp(a{u},'[TANKS]') || strcmp(a{u},'[JUNCTIONS]') || strcmp(a{u},'[RESERVOIRS]')))
                    if sjunction_check==0 && typecode==0  && strcmp(a{u},'[RESERVOIRS]')
                        fprintf(fid2,'[JUNCTIONS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s\n',newID,sps1{:},elevation,sps1{:},0,sps1{:});
                    end

                    if sreservoir_check==0 && typecode==1 && strcmp(a{u},'[TANKS]')
                        fprintf(fid2,'[RESERVOIRS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s\n', newID,sps1{:},elevation,sps1{:},'',sps1{:});
                    end

                    if stank_check==0 && typecode==2 && strcmp(a{u},'[PIPES]')
                        fprintf(fid2,'[TANKS]');
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d\n', newID,sps1{:},elevation,sps1{:},Initlevel,sps1{:},MinLevel,sps1{:},...
                        MaxLevel,sps1{:},Diameter,sps1{:},MinVol);
                    end

                    fprintf(fid2,'%s',a{u});

                    %%%%%%%% Jynctions Section %%%%%%%%
                    if typecode==0 && strcmp(a{u},'[JUNCTIONS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s', newID,sps1{:},elevation,sps1{:},0,sps1{:});
                    end
                    %%%%%%%% Reservoirs Section %%%%%%%%
                    if typecode==1 && strcmp(a{u},'[RESERVOIRS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s', newID,sps1{:},elevation,sps1{:},'',sps1{:});
                    end
                    %%%%%%%% Tanks Section %%%%%%%%
                    if typecode==2 && strcmp(a{u},'[TANKS]')
                        fprintf(fid2, '\n%s%s%d%s%d%s%d%s%d%s%d%s%d', newID,sps1{:},elevation,sps1{:},Initlevel,sps1{:},MinLevel,sps1{:},...
                        MaxLevel,sps1{:},Diameter,sps1{:},MinVol);
                    end

                elseif isempty(a{u})

                else

                    fprintf(fid2,'%s%s',a{u},sps1{:});
                end
                u=u+1;
            end

            %%%%%%%% Coordinates Section %%%%%%%%
            if gg~=0 && onetime==1
                % Correction Index
                if isempty(char(nodes.JunctionsID))
                    nodes.JunctionsID=[];
                end
                if isempty(nodes.ReservoirsID)
                    nodes.ReservoirsID=[];
                end
                if isempty(nodes.ReservoirsID)
                    nodes.NodesAll=[];
                end

                if (gg==length(nodes.JunctionsID)+length(nodes.NodeReservoirIndex)-1) && (typecode==1) || (gg==length(nodes.NodesAll)-1) && (typecode==2)
                    fprintf(fid2,'\r\n');
                    fprintf(fid2, '%s%s%d%s%d', newID,sps1{:},X,sps1{:},Y);
                    gg=0; onetime=0;
                end
            end

            if qualch==1 && qq==1
                qualch=0;
            end
            fprintf(fid2,'\n');
        end
    end

    fclose all;
%     [errcode] = Epanet('tempInpFile.inp');
%     obj.saveInputFile('tempInpFile.inp'); % OK %
%     movefile('tempInpFile.inp',[pwd,'\Results']);
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function addNewControl(obj,x,status,y_t_c,param,z,varargin)
    % Add control in the network. 
    % Syntax:  Addcontrol(inpname,x,status,y,param,z)
    %          LINK x status IF NODE y ABOVE/BELOW z
    %          
    %          Addcontrol(inpname,x,status,t)
    %          LINK x status AT TIME t
    % 
    %          Addcontrol(inpname,x,status,c,param)
    %          LINK x status AT CLOCKTIME c AM/PM
    % 
    % Inputs:  
    % inpname   name of an EPANET Input file.
    % x         a link ID label 
    % status    OPEN or CLOSED, a pump speed setting, or a control valve setting 
    % y         a node ID label 
    % param     ABOVE/BELOW or AM/PM
    % z         a pressure for a junction or a water level for a tank 
    % t         a time since the start of the simulation in decimal hours or in hours:minutes notation (string)
    % c         a 24-hour clock time (string)    % Examples: 
    % %1%
    % Addcontrol('Net1.inp','10','OPEN','10','ABOVE',100);
    % 
    % %2%
    % Addcontrol('Net1.inp','10','OPEN','10.00');
    % Removecontrol('Net1.inp',1,'10');
    % 
    % %3%
    % Addcontrol('Net1.inp','10','OPEN','12.00','AM');
    
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    
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
        nodes = obj.getNodesInfo;
        if length(char(nodes.NodesAll))==0 
            return
        end

        i=1;
        while i<length(char(nodes.NodesAll))+1
            exists(i) = strcmp(y_t_c,char(nodes.NodesAll(i)));
            i=i+1;
        end

        if (sum(exists)~=1)
            s = sprintf('There is no such object in the network.');
            warning(s);
            return
        end
    end

    % Check if id new already exists
    links = obj.getLinksInfo;
    if length(char(links.LinksAll))==0 
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    i=1;
    while i<length(char(links.LinksAll))+1
        exists(i) = strcmp(x,char(links.LinksAll(i)));
        i=i+1;
    end

    if (sum(exists)~=1)
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    type_n='[CONTROLS]';
    controls = obj.getControlsInfo;
    if controls.sectControls==0  type_n='[PATTERNS]'; end

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=0;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        t = t+1;
        info{t} = tline;
    end
    fid2 = fopen(obj.PathFile,'w');

    noo=0;s=0;yy=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;
            while u < length(a)+1
                if strcmp(a{u},type_n);
                    fprintf(fid2,'[CONTROLS]');
                    s=1; break;
                end

                spaces= 15 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                rr = regexp(a,'\w*[\w*]\w*','split');
                check_brackets = rr{:};
                ch1 = strcmp(check_brackets,'[');
                ch2 = strcmp(check_brackets,']');

                if (ch1(1)==1 && ch2(2)==1 && (s==1) && noo==0)
                    if yy==0
                        if controls.sectControls==0
                            fprintf(fid2,'[CONTROLS]\n');
                        end
                    end
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
                    fprintf(fid2,'%s%s',a{u},sps{:});
                    end
                end


                u=u+1;
            end
        end
        fprintf(fid2,'\n');
    end

    fclose all;
    [errcode] = Epanet('tempInpFile.inp');
    obj.saveInputFile('tempInpFile.inp'); % OK %
    movefile('tempInpFile.inp',[pwd,'\Results']);
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function [warn1] = rmLink(obj,LinkID)
    warn1=0;
    % Remove link from the network. 
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    % Check if id new already exists
    links = obj.getLinksInfo;
    if length(links.LinksAll)==0 
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    i=1;
    while i<length(links.LinksAll)+1
        exists(i) = strcmp(LinkID,char(links.LinksAll(i)));
        if exists(i)==1
            index_rmvlink=i;
        end
        i=i+1;
    end

    if (sum(exists)~=1)
        s = sprintf('There is no such object in the network.');
        warning(s);
        return
    end

    nodes = obj.getNodesInfo;
    from_node = links.FromNode(index_rmvlink);
    r = strcmp(nodes.NodesAll,from_node);
    if sum(r)==0, from_node=''; end
    to_node = links.ToNode(index_rmvlink);
    r = strcmp(nodes.NodesAll,to_node);
    if sum(r)==0, to_node=''; end

    % Remove control, code 1(LINK)
    obj.removeControlLinkID(LinkID);

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=0;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        t = t+1;
        info{t}=tline;
    end
    
    fid2 = fopen(obj.PathFile,'w');

    % section [JUNCTIONS]
    i=1;out=0;YY=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;x=0;xx=0;q=0;
            while u < length(a)+1

                spaces= 10 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                spaces1= 18 - length(a{u});
                k=1;sps1={' '};
                while k < spaces1+1
                    sps1 = strcat(sps1,{' '});
                    k=k+1;
                end

                if strcmp(a{u},'[PIPES]') YY=1;end
                if YY==1   
                    if isempty(a{u}) && (x==0)
                     u=u+1; x=1;xx=1;
                     if u==length(a)+1
                         break
                     end
                    end

                if strcmp(a{u},'[TAGS]') out=1; end
                if strcmp(a{u},'[STATUS]') out=1; end
                if strcmp(a{u},'[DEMANDS]') out=1; end
                if strcmp(a{u},'[PATTERNS]') out=1; end

                if strcmp(a{u},LinkID) && q~=1 && out==0
                    if xx==1 || strcmp(a{u},LinkID)
                    u=length(a)+1;
                    end
                else
                    q=1;
                    fprintf(fid2,'%s%s',a{u},sps{:});
                end
                else
                    if isempty(a{u})
                        u=u+1;
                        if u==length(a)+1
                           break
                        end
                    end
                    fprintf(fid2,'%s%s',a{u},sps{:});
                end
                u=u+1;
            end
        end

        fprintf(fid2,'\n');

    end

    fclose all;

    % Get nodes which delete with function Remove Node
    links = obj.getLinksInfo;

    i=1; 
    while i<length(links.LinksAll)+1
        t(i) = strcmp(from_node,char(links.FromNode(i)));
        tt(i) = strcmp(to_node,char(links.ToNode(i)));
        ttt(i) = strcmp(to_node,char(links.FromNode(i)));
        tttt(i) = strcmp(from_node,char(links.ToNode(i)));
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
            warning(s) 
        end
    end
    if sum(tt)+sum(ttt)==0
        if ~isempty(char(to_node)) 
            s = sprintf('Node %s disconnected.',char(to_node));
            warning(s) 
        end
    end
    %%%%%%%%
    if sum(t)+sum(tttt)==0 || sum(tt)+sum(ttt)==0                  
        if ~isempty(char(from_node)) || ~isempty(char(to_node))
            if ~sum(strcmp(from_node,nodes.ReservoirsID)) || ~sum(strcmp(from_node,nodes.TanksID))...
               || ~sum(strcmp(to_node,nodes.ReservoirsID)) || ~sum(strcmp(to_node,nodes.TanksID))
                warn1=0;
            else 
                warn1=1;
            end
        end
    end
    %%%%%%%%
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
%     [errcode] = Epanet('tempInpFile.inp');
%     obj.saveInputFile('tempInpFile.inp'); % OK %
%     movefile('tempInpFile.inp',[pwd,'\Results']);
end

function rmNode(obj,NodeID)
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    % Remove node from the network. 
    checklinks='';
    checknodes='';

    % Check if id new already exists
    nodes = obj.getNodesInfo;
    if length(nodes.NodesAll)==0 
        return;
    end

    i=1;
    while i<length(nodes.NodesAll)+1
        exists(i) = strcmp(NodeID,char(nodes.NodesAll(i)));
        i=i+1;
    end

    if (sum(exists)~=1)
        s = sprintf('There is no such object in the network.');
        warning(s);
        return;
    elseif length(char(nodes.ReservoirsID))+length(char(nodes.TanksID))<2
        warning('One or more errors in input file.')
        warning('This tank/reservoir has not removed.')
        return;
    end
    % Get links which delete with function Remove Link
    links = obj.getLinksInfo;
    checklinks_index=0;

    % Check NodeID, from from..to node 
    w=1;u=1;
    while (w<length(links.LinksAll)+1)
        t = strcmp(links.FromNode,NodeID);
        tt = strcmp(links.ToNode,NodeID);
        if t(w)==1
            id_from=links.ToNode(w);
            checknodes{u}=id_from;
            checklinks{u}=links.LinksAll(w);
            checklinks_index(u)=w;
            u=u+1;
        end
        if tt(w)==1
            id_to = links.FromNode(w);
            checknodes{u}=id_to;
            checklinks{u}=links.LinksAll(w);
            checklinks_index(u)=w;
            u=u+1;
        end
        w=w+1;
    end

    i=1;u=length(checklinks)+1;
    while i<length(links.ToNode)+1
        t = strcmp(nodes.NodesAll,links.FromNode(i));
        tt = strcmp(nodes.NodesAll,links.ToNode(i));
        if t==0  
            checklinks{u}=links.LinksAll(i);
            u=u+1;tt=1;
        end
        if tt==0
            checklinks{u}=links.LinksAll(i);
            u=u+1;
        end
        i=i+1;
    end

    % Remove control, code 0(NODE)
    obj.removeControlNodeID(NodeID);

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=0;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        t = t+1;
        info{t} = tline;
    end

    fid2 = fopen(obj.PathFile,'w');

    i=1;out=0;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c,'\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;x=0;xx=0;q=0;
            while u < length(a)+1

                spaces= 10 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                spaces1= 18 - length(a{u});
                k=1;sps1={' '};
                while k < spaces1+1
                    sps1 = strcat(sps1,{' '});
                    k=k+1;
                end

                if isempty(a{u}) && (x==0)
                     u=u+1; x=1;xx=1;
                     if u==length(a)+1
                         break
                     end
                end

                if strcmp(a{u},'[TANKS]') out=0; end
                if strcmp(a{u},'[PIPES]') out=1; end
                if strcmp(a{u},'[DEMANDS]') out=0; end %out=0; delete line
                if strcmp(a{u},'[QUALITY]') out=0; end
                if strcmp(a{u},'[SOURCES]') out=1; end 
                if strcmp(a{u},'[MIXING]') out=0; end
                if strcmp(a{u},'[COORDINATES]') out=0; end

                if strcmp(a{u},NodeID) && q~=1 && out==0
                    if xx==1 || strcmp(a{u},NodeID)
                        u=length(a)+1;
                    end
                else
                    q=1;
                    fprintf(fid2,'%s%s',a{u},sps{:});
                end

                u=u+1;
            end
        end
        fprintf(fid2,'\n');
    end

    fclose all;

    % Remove links
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if length(checklinks)
        for i=1:length(checklinks)
            warn1(i)=obj.removeLinkID(checklinks{i});
        end
    else
        warn1=1;
    end
    % Find who other id must be delete
    remove_link={''};
    remove_link_index = zeros(1,length(links.FromNode));

    if length(checklinks) 
        i=1;
        while i < (length(checklinks)+1)
            remove_link(i)=checklinks{i};
            remove_link_index(i)=i;
            s=sprintf('Removed link:%s',char(remove_link(i)));
            warning(s);
            i=i+1;
        end
    else
        warn1=sum(warn1)+1;
    end
    
    if sum(warn1)~=0
        [errcode] = Epanet('tempInpFile.inp');
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
end

function rmControl(obj,type,id)
    % Remove control from the network.
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    
    exists=0;
    exists1=0;

    controls = obj.getControlsInfo;
    
    if type==1
        if length(char(controls.linksID))==0 
            s = sprintf('There is no such object in the network.');
            warning(s);
            return
        end

        i=1;
        while i<length(char(controls.linksID))+1
            if type==1
                exists(i) = strcmp(controls.linksID{i},char(id));
            else
                warning('Type is NODE(0) or LINK(1)');
                return
            end
            i=i+1;
        end
    end

    if type==0
        if length(char(controls.nodesID))==0 
            s = sprintf('There is no such object in the network.');
            warning(s);
            return
        end

        i=1;
        while i<length(char(controls.nodesID))+1
            if type==0
                exists1(i) = strcmp(controls.nodesID{i},char(id));
            else
                warning('Type is NODE(0) or LINK(1)');
                return
            end
            i=i+1;
        end
    end

    if (sum(exists)==0) && (sum(exists1)==0)
        warning('There are no Controls in the network.');
        return
    end

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=0;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        t = t+1;
        info{t} = tline;
    end
    
    fid2 = fopen(obj.PathFile,'w');

    i=1;e=0;n=0;kk=1; 
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1; 
            while u < length(a)+1

                spaces= 10 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
                end

                rr = regexp(a,'\w*[\w*]\w*','split');
                check_brackets = rr{:};
                ch1 = strcmp(check_brackets,'[');
                ch2 = strcmp(check_brackets,']');

                if strcmp(a{u},'[CONTROLS]')
                    fprintf(fid2,'%s',a{u}); 
                    n=1; 
                elseif ch1(1)==1 && ch2(2)==1 && n==1
                    if (isempty(a{u})&& n==1) break; end
                    e=1;
                end
                if strcmp(a{u},'[END]')  e=1; fprintf(fid2,'%s',a{u});break;   end

                if n==1 && e==0 && kk==1
                    if strcmp(a{u},'[CONTROLS]') break; end
                    if isempty(a{u})
                    elseif strfind(a{u},';')
                        u = length(a)+1;break;
                    else 
                        %%%%%%%%
                        if type==1  
                            tt = strcmp(a{u+1},id); kk=0;
                            if tt==1
                                u = length(a)+1;break;
                            else
                            fprintf(fid2,'%s%s',a{u},sps{:});
                            end

                        elseif type==0
                            tt = strcmp(a{u+5},id); kk=0;
                            if tt==1
                                u = length(a)+1;break;
                            else
                            fprintf(fid2,'%s%s',a{u},sps{:});
                            end
                        end
                    end
                else
                    if isempty(a{u})
                    else 
                        fprintf(fid2,'%s%s',a{u},sps{:});
                    end
                end

                u=u+1;

            end
        end
            fprintf(fid2,'\n');kk=1;
    end

    if n==1
        while ~feof(fid)
            tline=fgetl(fid);
            if strcmp(controls.controlsInfo{index},tline)==0 
                fprintf(fid2,'%s',tline);
                fprintf(fid2,'\n');
            else
                fprintf(fid2,'\n');
            end
        end
    else
        warning('There was at least one error in input file.');
    end
    while ~feof(fid)
        tline=fgetl(fid);
        fprintf(fid2,'%s',tline);
        fprintf(fid2,'\n');
    end
    fclose all;
    if tmp==1
        [errcode] = Epanet('tempInpFile.inp');
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function Options(obj,newFlowUnits,headloss,varargin)
    % Notes:
    % Flow units codes are as follows:
    % CFS	cubic feet per second
    % GPM	gallons per minute
    % MGD	million gallons per day
    % IMGD	Imperial mgd
    % AFD	acre-feet per day
    % LPS	liters per second
    % LPM	liters per minute
    % MLD	million liters per day
    % CMH	cubic meters per hour
    % CMD	cubic meters per day
    
    load([pwd,'\Results\','tmpInp'],'tmp','-mat')
    if tmp==1
        obj.saveInputFile('tempInpFile.inp'); % OK %
        movefile('tempInpFile.inp',[pwd,'\Results']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    value=obj.getFlowUnitsHeadlossFormula;
    previousFlowUnits=value.FlowUnits;
    
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
        headloss=value.Headloss;
        nheadl=1;
    end
    nodes = obj.getNodesInfo;
    links = obj.getLinksInfo;
    controls=obj.getControlsInfo;
    curves=obj.getCurveInfo;

    % Open and read inpname
    fid = fopen(obj.PathFile,'r+');

    t=0;a=0;
    % Read all file and save in variable info
    while ~feof(fid)
        tline=fgetl(fid);
        t = t+1;
        info{t} = tline;
        c = info{t};
        a = regexp(c, '\s*','split');
    end
    fid2 = fopen(obj.PathFile,'w');
                    
    sections=[0 0 0 0 0 0 0 0 0];
    
    nn=0;i=1;pp=1;
    for t = 1:length(info)
        c = info{t};
        a = regexp(c, '\s*','split');
        y=1;
        while y < length(a)+1
            j(y) = isempty(a{y});
            y=y+1;
        end
        j = sum(j);
        if j == length(a)
            % skip
            i=i+1;
        elseif isempty(c)
            % skip
            i=i+1;
        else
            i=i+1;

            u=1;
            while u < length(a)+1
                if strcmp(a{u},'[JUNCTIONS]') && Units
                    fprintf(fid2,'[JUNCTIONS]');
                    sections=[1 0 0 0 0 0 0 0 0];
                    break;
                elseif strcmp(a{u},'[RESERVOIRS]') && Units
                    fprintf(fid2,'[RESERVOIRS]');
                    nn=0;pp=1;    
                    sections=[0 1 0 0 0 0 0 0 0];
                    break;
                elseif strcmp(a{u},'[TANKS]') && Units
                    fprintf(fid2,'[TANKS]');
                    nn=0;pp=1;    
                    sections=[0 0 1 0 0 0 0 0 0];
                    break;
                elseif strcmp(a{u},'[PIPES]') && Units
                    fprintf(fid2,'[PIPES]');
                    nn=0;pp=1;    
                    sections=[0 0 0 1 0 0 0 0 0];
                    break;   
                elseif strcmp(a{u},'[VALVES]') && Units
                    fprintf(fid2,'[VALVES]');
                    nn=0;pp=1;    
                    sections=[0 0 0 0 1 0 0 0 0];
                    break; 
                elseif strcmp(a{u},'[DEMANDS]') && ((Units || ~changes) && nheadl)
                    fprintf(fid2,'[DEMANDS]');
                    nn=0;pp=1;    
                    sections=[0 0 0 0 0 1 0 0 0];
                    break;  
%                 elseif strcmp(a{u},'[CURVES]') && ((Units || ~changes) && nheadl)
%                     fprintf(fid2,'[CURVES]');
%                     nn=0;pp=1;    
%                     sections=[0 0 0 0 0 0 1 0 0];
%                     break;                     
                elseif strcmp(a{u},'[CONTROLS]') && Units
                    fprintf(fid2,'[CONTROLS]');
                    nn=0;pp=1;    
                    sections=[0 0 0 0 0 0 0 1 0];
                    break;                      
                elseif strcmp(a{u},'[OPTIONS]') 
                    fprintf(fid2,'[OPTIONS]');
                    sections=[0 0 0 0 0 0 0 0 1];nn=0;
                    break;
                end
                
                spaces= 15 - length(a{u});
                k=1;sps={' '};
                while k < spaces+1
                    sps = strcat(sps,{' '});
                    k=k+1;
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
                                for mm=mm:4
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
                    
                % section [VALVES]
                elseif (sect==5) && (nn==0)  
                    mm=1;
                    if mm < length(a)+1
                        if isempty(a{mm})
                        % skip
                            mm=mm+1;
                        end
                        if pp<length(char(links.ValvesID))+1 
                            if strcmp(a{mm},links.ValvesID{pp})  
                                pp=pp+1;
                                for mm=mm:4
                                    fprintf(fid2,'%s%s',char(a{mm}),sps{:});
                                end
                                if changes==1
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*25.4),sps{:});
                                elseif changes==2
                                    fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.03937007874),sps{:});
                                end
                                for mm=(mm+1):(mm+4)
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
                elseif (sect==6) && (nn==0)  
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
%                 % section [CURVES]
%                 elseif (sect==7) && (nn==0)  
%                     mm=1;
%                     if mm < length(a)+1
%                         if isempty(a{mm})
%                         % skip
%                             mm=mm+1;
%                         end
%                         if pp<length(curves.CurvesID)+1 
%                             pp=pp+1;
%                             if strcmp(upper(a{mm+6}),';PUMP')  
%                                 fprintf(fid2,'%s%s',char(a{mm}),sps{:});
%                                 setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
%                                 if changes==1
%                                     fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.3048),sps{:});
%                                 elseif changes==2
%                                     fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.281),sps{:}); 
%                                 end
%                             end
%                         else
%                            nn=1;
%                            fprintf(fid2,'%s%s',char(a{1}),sps{:});
%                         end
%                     end
%                     break;
                    
               % section [CONTROLS]
                elseif (sect==8) && (nn==0)  
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
                                end
                            end
                        else
                           nn=1;
                           fprintf(fid2,'%s%s',char(a{1}),sps{:});
                        end
                    end
                    break;
                    
                % section [OPTIONS]
                elseif (sect==9) && (nn==0)
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
                           fprintf(fid2,'%s%s',char(a{1}),sps{:});
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
    [errcode] = Epanet('tempInpFile.inp');
    obj.saveInputFile('tempInpFile.inp');  
    movefile('tempInpFile.inp',[pwd,'\Results']);
    tmp=0;
    save([pwd,'\Results\','tmpInp'],'tmp','-mat')
end

function setflow(previousFlowUnits,newFlowUnits,fid2,a,sps,mm)
    if strcmp(previousFlowUnits,'GPM')
        switch newFlowUnits %(GPM)
            case 'CFS'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.228009e-3),sps{:});
            case 'MGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.00144),sps{:});
            case 'IMGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.19905e-3),sps{:});
            case 'AFD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.419191e-3),sps{:});                                    
            case 'LPS'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.0630902),sps{:});                                    
            case 'LPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.785412),sps{:});                                    
            case 'MLD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*5.450993e-3),sps{:});                                    
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
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.531466e-2),sps{:});
            case 'GPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*15.85032),sps{:});
            case 'MGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.282446e-2),sps{:});
            case 'IMGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.900533e-2),sps{:});                                    
            case 'AFD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*7.004562e-2),sps{:});                                    
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
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*5.885777e-4),sps{:});
            case 'GPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.264172),sps{:});
            case 'MGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.804078e-4),sps{:});
            case 'IMGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*3.167556e-4),sps{:});                                    
            case 'AFD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.1674272e-3),sps{:});                                    
            case 'LPS'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.666667e-2),sps{:});                                    
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
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*9.809635e-3),sps{:});
            case 'GPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.402868),sps{:});
            case 'MGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*6.340129e-3),sps{:});
            case 'IMGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*5.27926e-3),sps{:});                                    
            case 'AFD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.945712e-2),sps{:});                                    
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
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.087345e-4),sps{:});
            case 'GPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.1834528),sps{:});
            case 'MGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.64172e-4),sps{:});
            case 'IMGD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*2.199692e-4),sps{:});                                    
            case 'AFD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*8.107132e-4),sps{:});                                    
            case 'LPS'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*1.157407e-2),sps{:});                                    
            case 'LPM'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.6944444),sps{:});                                    
            case 'MLD'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*0.001),sps{:});                                    
            case 'CMH'
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})*4.166667e-2),sps{:});                                    
            otherwise
                fprintf(fid2,'%s%s',num2str(str2num(a{mm+1})),sps{:});   
        end                                    
    end
end
