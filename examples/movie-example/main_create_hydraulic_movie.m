start_toolkit;

settings=[];
settings.inpname = 'Battle of the Calibration Networks System.inp'; %Net1 net2-cl2 Net3 example BWSN_Network_1
settings.duration = 15; % hours

settings.movname = [settings.inpname(1:end-4),'.avi'];
settings.fps = 8;
settings.quality = 100;
settings.fig = [];
settings.vsize = 5;  % Size of vertices in points (0 == omits verts)
settings.lwidth = 3; % Width of links in points
settings.tsize = 10;  % Size of tank/reservoir nodes
settings.hyd = 1; % code 1 for hydraulics, code 0 for quality

d=epanet(settings.inpname);
d.setTimeSimulationDuration(3600*settings.duration);

%% Get the simulation data using Epanet
NodeType = 'Pressure'; % Demand, Pressure
LinkType = 'Flow';   % Flow, Velocity
[V,L,T] = getHydraulicData(NodeType,LinkType,d);
        
[PData, SData] = movie_parameters(settings,V,L);

%% Write the Movie File
% NetworkMovie will display the first frame in a figure window and allow
% you to make adjustments to it (zoom, pan, etc.) before rendering the
% frames into an AVI movie file.
NetworkMovie(V,L,settings.fig,settings.movname,...
    settings.quality,settings.fps,PData,SData,d,NodeType,LinkType,settings.hyd);

%% Unload EPANET-MATLAB Toolkit  
d.unload;

%% Show the Movie
% You could display the movie in Matlab as follows, or just use 
% an external viewer...
implay(settings.movname);

%% Delete tmp msx file
if isfield(settings,'msxname')
    if exist(settings.msxname)==2, delete(settings.msxname); end
end