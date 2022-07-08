% Using PARFOR
% More fast function with parfor is coming..

tic
start_toolkit;

% try
%     unloadlibrary('epanet2')
% catch
% end

d = epanet('Net1.inp'); 

clear H;clc;

number_scenarios = 100;

parfor i = 1:number_scenarios
    
% Uncomment section for MATLAB R2020 and previous versions.
%                                        
%     if isunix
%         loadlibrary(obj.LibEPANET, [obj.LibEPANETpath, obj.LibEPANET, '.h']);
%     else
%         loadlibrary([obj.LibEPANETpath, obj.LibEPANET], [obj.LibEPANETpath, obj.LibEPANET, '.h']);
%     end
%     d.loadEPANETFile(d.TempInpFile); 
    
    % set parameters
    elevations = d.getNodeElevations-d.getNodeElevations*rand(1)*.5;
    d.setNodeElevations(elevations*10);
    
    % Computed Hydraulics
    H{i} = d.getComputedHydraulicTimeSeries;
end
d.unload;
toc


