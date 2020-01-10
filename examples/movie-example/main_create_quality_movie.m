start_toolkit;

settings=[];
settings.inpname = 'Net1.inp'; %Net1 net2-cl2 Net3 example BWSN_Network_1
% settings.msxname = 'example.msx'; % add msx file or check for AGE
settings.bulk_specie_id = 'AGE'; % AGE CL2 NH2CL Animate the bulk specie
settings.wall_specie_id = ''; % AS5s Animate the wall specie
settings.duration = 48; % hours

settings.movname = [settings.inpname(1:end-4),'.avi'];
settings.fps = 8;
settings.quality = 100;
settings.fig = [];
settings.vsize = 5;  % Size of vertices in points (0 == omits verts)
settings.lwidth = 3; % Width of links in points
settings.tsize = 10;  % Size of tank/reservoir nodes

settings.hyd = 0; % code 1 for hydraulics, code 0 for quality

d=epanet(settings.inpname);
d.setTimeSimulationDuration(3600*settings.duration);
NodeType = settings.bulk_specie_id;

if strcmp(settings.bulk_specie_id,'AGE')
    out_msx = create_msx_file_water_age;%(input);
    d.writeMSXFile(out_msx)
    settings.msxname = out_msx.FILENAME;
    d.loadMSXFile(settings.msxname,d.LibEPANETpath)
elseif ~strcmp(settings.bulk_specie_id,'CL2') 
    d.loadMSXFile(settings.msxname);
else
    settings.msxname=[];
end

%   Get the simulation data using Epanet or Epanet-MSX
[V,L,T] = getQualityData(settings.bulk_specie_id,...
    settings.wall_specie_id,settings.msxname,d);
        
[PData, SData] = movie_parameters(settings,V,L);

%% Write the Movie File
% NetworkMovie will display the first frame in a figure window and allow
% you to make adjustments to it (zoom, pan, etc.) before rendering the
% frames into an AVI movie file.
NetworkMovie(V,L,settings.fig,settings.movname,...
    settings.quality,settings.fps,PData,SData,d,NodeType,[],settings.hyd);

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
