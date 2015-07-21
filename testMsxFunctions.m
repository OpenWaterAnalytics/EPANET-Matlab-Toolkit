%% EPANET-Matlab Class Test Part 2
% This file is provided to ensure that all functions can be executed
% correctly.
% Press F10 for step-by-step execution. You may also use the breakpoints, 
% indicated with a short dash (-) on the left of each line number.
clc;
clear;
close all;clear class;

% Create EPANET object using the INP file
inpname='example.inp'; %Net2_Rossman2000 example

%% MSX Functions
if strcmp(computer('arch'),'win64') 
    version='epanet20013patchx64'; % epanet20012x64  epanet20013patchx64
elseif strcmp(computer('arch'),'win32')
    version='epanet20013patchx86'; % epanet20012x86  epanet20013patchx86
end
%d=epanet(inpname)
d=epanet(inpname,version);
d.msx([inpname(1:end-4),'.msx'])
d


%% GET PARAMETERS
d.getMsxEquationsTerms
d.getMsxEquationsPipes
d.getMsxEquationsTanks
d.getMsxTimeStep
d.getMsxSpeciesCount
d.getMsxConstantsCount
d.getMsxParametersCount
d.getMsxPatternsCount
d.getMsxSpeciesNameID
d.getMsxSpeciesType
d.getMsxSpeciesUnits
d.getMsxSpeciesATOL  
d.getMsxSpeciesRTOL
d.getMsxSpeciesIndex
d.getMsxSpeciesIndex('AS5') 
d.getMsxConstantsNameID
d.getMsxConstantsValue
d.getMsxConstantsIndex
d.getMsxConstantsIndex('K2')
d.getMsxParametersNameID
d.getMsxParametersIndex
d.getMsxParametersTanksValue
d.getMsxParametersPipesValue
d.getMsxPatternsNameID
d.getMsxPatternsIndex
d.getMsxPatternsLengths
d.getMsxNodeInitqualValue
d.getMsxLinkInitqualValue
d.getMsxSources
d.getMsxSourceType
d.getMsxSourceLevel
d.getMsxSourcePatternIndex
d.getMsxPattern %Mass flow rate per minute of a chemical source
d.getMsxPatternValue(1,5) %Mass flow rate per minute of a chemical source
% % d.getMsxSpeciesConcentration
disp('Press any key to continue...')
pause


%% GET, SET TIME STEPS
d.setTimeHydraulicStep(3600)
d.setTimeQualityStep(3600)
d.getMsxComputedQualityNode(1)%index node
d.getMsxComputedQualityNode(1,1)%index node, index species
disp('Press any key to continue...')
pause


%% MSX PLOTS
lll=d.getMsxComputedQualityLink
figure;cmap=hsv(5);for i=1:d.getMsxSpeciesCount;plot(lll.Time,lll.Quality{1}{i},'Color',cmap(i,:));hold on; end; legend(d.MsxSpeciesNameID)

nnn=d.getMsxComputedQualityNode
figure;cmap=hsv(5);for i=1:d.getMsxSpeciesCount;plot(nnn.Time,nnn.Quality{1}{i},'Color',cmap(i,:));hold on; end; legend(d.MsxSpeciesNameID)

d.MsxPlotConcentrationSpeciesOfNodes(1,1:d.MsxSpeciesCount)
d.MsxPlotConcentrationSpeciesOfNodes(2,1:d.MsxSpeciesCount)
d.MsxPlotConcentrationSpeciesOfNodes(3,1:d.MsxSpeciesCount)
d.MsxPlotConcentrationSpeciesOfNodes(4,1:d.MsxSpeciesCount)
d.MsxPlotConcentrationSpeciesOfNodes(5,1:d.MsxSpeciesCount)
d.getMsxComputedQualityLink(1,1:d.MsxSpeciesCount)%index link, index species
d.MsxPlotConcentrationSpeciesOfLinks(1,1:d.MsxSpeciesCount)
disp('Press any key to continue...')
pause


%% Errors Msx
for e=[0,200,501:524]
    d.getMsxError(e)
    % 0 bug, no error
    % 200 bug, cannot read EPANET-MSX file
    % 510 bug, could not open algebraic "equation" solver.
end
disp('Press any key to continue...')
pause

%% WRITE REPORT
% Solve for hydraulics & water quality
d.MsxSolveCompleteHydraulics
d.MsxSolveCompleteQuality
% Write results to the TestMsxReport file
d.MsxWriteReport %a specific water quality report file is named in the [REPORT] section of the MSX input file. %BUG
open([d.MsxTempFile(1:end-4),'.txt']);
disp('Press any key to continue...')
pause

%or
fid = fopen('ReportMsx.bat','w');
r = sprintf('epanetmsx %s %s %s',d.inputfile,d.MsxTempFile,[d.MsxTempFile(1:end-4),'.txt']); 
fprintf(fid,'%s \n',r);fclose all;
!ReportMsx.bat
open([d.MsxTempFile(1:end-4),'.txt'])
disp('Press any key to continue...')
pause


%% GET, ADD PATTERNS
d.MsxAddPattern('testpat',[2 .3 .4 6 5 2 4]);
d.getMsxPatternsNameID
d.getMsxPatternsIndex
d.getMsxPatternsLengths  
disp('Press any key to continue...')
pause


%% GET, SET SOURCES
v=d.getMsxSources
node = 1;
spec=1;
type = 0;
level=0.2;
pat = 1;
d.setMsxSources(node, spec, type, level, pat)
v=d.getMsxSources
disp('Press any key to continue...')
pause


%% GET, SET CONSTANTS
d.getMsxConstantsValue     
value = [2 10 8];%index[1 2 3]
d.setMsxConstantsValue(value);
d.getMsxConstantsNameID
d.getMsxConstantsValue
d.getMsxConstantsIndex
disp('Press any key to continue...')
pause


%% GET, SET PARAMETERS TANKS/PIPES
d.getMsxParametersTanksValue
d.getMsxParametersPipesValue      
disp('Press any key to continue...')
pause

if d.getMsxParametersCount 
    % d.setMsxParametersPipesValue(pipeIndex,value) 
    d.setMsxParametersPipesValue(1,[1.5 2]) 
    d.getMsxParametersPipesValue{1}        

    a=d.getNodeTankIndex
    d.setMsxParametersTanksValue(a(1),100)  
    d.getMsxParametersTanksValue{a(1)}  
    disp('Press any key to continue...')
    pause
end


%% QUALITY
values = d.getMsxLinkInitqualValue
nodeIndex=1; speciesIndex=1;
values{nodeIndex}(speciesIndex)=1000;%
d.setMsxLinkInitqualValue(values)     
d.getMsxLinkInitqualValue 
disp('Press any key to continue...')
pause

linkIndex=1; speciesIndex=1;
values = d.getMsxNodeInitqualValue
values{linkIndex}(speciesIndex)=1500;%
d.setMsxNodeInitqualValue(values)
d.getMsxNodeInitqualValue   

%% GET, SET PATTERN
d.setMsxPatternMatrix([.1 .2 .5 .2 1 .9]);
d.getMsxPattern

d.setMsxPatternValue(1,1,2);
d.getMsxPattern 

d.setMsxPattern(1,[1 0.5 0.8 2 1.5]);
d.getMsxPattern 
d.getMsxPatternValue(1,5) 

%% SAVE, USE FILES
d.MsxSaveFile('testMsx.msx');                                                               
          
d.MsxSaveQualityFile('testMsxQuality.bin')

d.saveHydraulicsOutputReportingFile
d.saveHydraulicFile('testMsxHydraulics.hyd')

d.MsxUseHydraulicFile('testMsxHydraulics.msx')

% % MsxInitializeQualityAnalysis
% % MsxStepQualityAnalysisTimeLeft

d.MsxUnload

d.unload
%Delete s files 
sfilesexist = dir('s*'); 
if (~isempty(sfilesexist)), delete s*, end;
delete('testMsx.msx','*.hyd','*.bin','*bat*','*_temp*',[d.inputfile(1:end-4),'.txt'])