%% EPANET-Matlab Class Test Part 2
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.

% Execute "addpath(genpath(pwd))" in main folder before running this, to load all EPANET
% functions

clc;
clear;
close all;clear class;

% Create EPANET object using the INP file
inpname='networks/example.inp'; %Net2_Rossman2000 example
inpname='networks/Net2_Rossman2000.inp'

%% MSX Functions
d=epanet(inpname);
d.loadMSXFile([inpname(1:end-4),'.msx'])
d

% New functions - Read MSX File
d.getMSXSolver
d.getMSXAreaUnits
d.getMSXRateUnits
d.getMSXRtol
d.getMSXAtol
d.getMSXTimeStep
d.getMSXCoupling
d.getMSXCompiler

d.setMSXTimeStep(3600) 

d.setMSXAreaUnitsFT2
d.setMSXAreaUnitsM2
d.setMSXAreaUnitsCM2

d.setMSXRateUnitsSEC
d.setMSXRateUnitsMIN
d.setMSXRateUnitsHR
d.setMSXRateUnitsDAY

d.setMSXSolverEUL
d.setMSXSolverRK5
d.setMSXSolverROS2
 
d.setMSXCouplingFULL
d.setMSXCouplingNONE

d.setMSXCompilerVC % depends on the C compiler you are using
d.setMSXCompilerGC % depends on the C compiler you are using
d.setMSXCompilerNONE % depends on the C compiler you are using

d.setMSXAtol(0.1)
d.setMSXRtol(0.2)

%% GET PARAMETERS
d.getMSXEquationsTerms
d.getMSXEquationsPipes
d.getMSXEquationsTanks
d.getMSXTimeStep
d.getMSXSpeciesCount
d.getMSXConstantsCount
d.getMSXParametersCount
d.getMSXPatternsCount
d.getMSXSpeciesNameID
d.getMSXSpeciesType
d.getMSXSpeciesUnits
d.getMSXSpeciesATOL  
d.getMSXSpeciesRTOL
d.getMSXSpeciesIndex
d.getMSXSpeciesIndex('AS5') 
d.getMSXConstantsNameID
d.getMSXConstantsValue
d.getMSXConstantsIndex
d.getMSXConstantsIndex('K2')
d.getMSXParametersNameID
d.getMSXParametersIndex
d.getMSXParametersTanksValue
d.getMSXParametersPipesValue
d.getMSXPatternsNameID
d.getMSXPatternsIndex
d.getMSXPatternsLengths
d.getMSXNodeInitqualValue
d.getMSXLinkInitqualValue
d.getMSXSources
d.getMSXSourceType
d.getMSXSourceLevel
d.getMSXSourcePatternIndex
d.getMSXPattern %Mass flow rate per minute of a chemical source
d.getMSXPatternValue(1,2) %Mass flow rate per minute of a chemical source
% % d.getMSXSpeciesConcentration
disp('Press any key to continue...')
pause


%% GET, SET TIME STEPS
d.setTimeHydraulicStep(3600)
d.setTimeQualityStep(3600)
d.getMSXComputedQualityNode(1)%index node
d.getMSXComputedQualityNode(1,1)%index node, index species
disp('Press any key to continue...')
pause


%% MSX PLOTS
lll=d.getMSXComputedQualityLink
figure;cmap=hsv(5);for i=1:d.getMSXSpeciesCount;plot(lll.Time,lll.Quality{1}{i},'Color',cmap(i,:));hold on; end; legend(d.MSXSpeciesNameID)

nnn=d.getMSXComputedQualityNode
figure;cmap=hsv(5);for i=1:d.getMSXSpeciesCount;plot(nnn.Time,nnn.Quality{1}{i},'Color',cmap(i,:));hold on; end; legend(d.MSXSpeciesNameID)

d.plotMSXSpeciesNodeConcentration(1,1:d.MSXSpeciesCount)
d.plotMSXSpeciesNodeConcentration(2,1:d.MSXSpeciesCount)
d.plotMSXSpeciesNodeConcentration(3,1:d.MSXSpeciesCount)
d.plotMSXSpeciesNodeConcentration(4,1:d.MSXSpeciesCount)
d.plotMSXSpeciesNodeConcentration(5,1:d.MSXSpeciesCount)
d.getMSXComputedQualityLink(1,1:d.MSXSpeciesCount)%index link, index species
d.plotMSXSpeciesLinkConcentration(1,1:d.MSXSpeciesCount)
disp('Press any key to continue...')
pause


%% Print Errors MSX
for e=[0,200,501:524]
    disp(d.getMSXError(e))
    % 0 bug, no error
    % 200 bug, cannot read EPANET-MSX file
    % 510 bug, could not open algebraic "equation" solver.
end
disp('Press any key to continue...')
pause

%% WRITE REPORT
% Solve for hydraulics & water quality (stored in memory)
d.solveMSXCompleteHydraulics
d.solveMSXCompleteQuality

% Write results to "./networks" folder
d.writeMSXReport %a specific water quality report file is named in the [REPORT] section of the MSX input file. 
% There is a possible bug: The report is not complete.
open([d.MSXTempFile(1:end-4),'.txt']);
disp('Press any key to continue...')
pause

% Different way to write the report (the "bug" above does not appear)
d.runMSXexe
open([d.MSXTempFile(1:end-4),'.txt']);

disp('Press any key to continue...')
pause


%% GET, ADD PATTERNS
d.addMSXPattern('testpat',[2 .3 .4 6 5 2 4]);
d.getMSXPatternsNameID
d.getMSXPatternsIndex
d.getMSXPatternsLengths  
disp('Press any key to continue...')
pause


%% GET, SET SOURCES
v=d.getMSXSources
node = 1;
spec=1;
type = 0;
level=0.2;
pat = 1;
d.setMSXSources(node, spec, type, level, pat)
v=d.getMSXSources
disp('Press any key to continue...')
pause


%% GET, SET CONSTANTS
d.getMSXConstantsValue     
value = [2 10 8];%index[1 2 3]
d.setMSXConstantsValue(value);
d.getMSXConstantsNameID
d.getMSXConstantsValue
d.getMSXConstantsIndex
disp('Press any key to continue...')
pause


%% GET, SET PARAMETERS TANKS/PIPES
d.getMSXParametersTanksValue
d.getMSXParametersPipesValue      
disp('Press any key to continue...')
pause

if d.getMSXParametersCount %runs for Net2
    d.setMSXParametersPipesValue(1,[1.5 2]) % pipeIndex,value 

    d.getMSXParametersPipesValue{1}        

    a=d.getNodeTankIndex
    d.getMSXParametersTanksValue{a(1)}
    d.setMSXParametersTanksValue(a(1),1,0.5) % tank_index, parameter_index, value
    d.getMSXParametersTanksValue{a(1)}  
    disp('Press any key to continue...')
    pause
end


%% QUALITY
values = d.getMSXLinkInitqualValue
nodeIndex=1; speciesIndex=1;
values{nodeIndex}(speciesIndex)=1000;%
d.setMSXLinkInitqualValue(values)     
d.getMSXLinkInitqualValue 
disp('Press any key to continue...')
pause

linkIndex=1; speciesIndex=1;
values = d.getMSXNodeInitqualValue
values{linkIndex}(speciesIndex)=1500;%
d.setMSXNodeInitqualValue(values)
d.getMSXNodeInitqualValue   

%% GET, SET PATTERN
d.setMSXPatternMatrix([.1 .2 .5 .2 1 .9]);
d.getMSXPattern

d.setMSXPatternValue(1,1,2);
d.getMSXPattern 

d.setMSXPattern(1,[1 0.5 0.8 2 1.5]);
d.getMSXPattern 
d.getMSXPatternValue(1,5) 

%% SAVE, USE FILES
d.saveMSXFile('testMSX.msx');                                                               
          
d.saveMSXQualityFile('testMSXQuality.bin')

d.saveHydraulicsOutputReportingFile
d.saveHydraulicFile('testMSXHydraulics.hyd')

d.useMSXHydraulicFile('testMSXHydraulics.hyd')

% % initializeMSXQualityAnalysis
% % stepMSXQualityAnalysisTimeLeft

d.unloadMSX

d.unload
%Delete s files (temporary files created by the C library
sfilesexist = dir('s*'); 
if (~isempty(sfilesexist)), delete s*, end;
delete('testMSX.msx','*.hyd','*.bin','*bat*','*_temp*')

fprintf('Test finished.\n')