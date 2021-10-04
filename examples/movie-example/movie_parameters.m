function [PData, SData] = movie_parameters(settings,V,L)
%% Specify Movie Parameters
% These parameters are described in NetworkFrame()...
fig = settings.fig;         % Use a new figure window
movFname = settings.movname;% Movie file name
quality = settings.quality; % 0-100 movie quality (related to data compression)
fps = settings.fps;         % Frame rate - # to display per second

PData.c = 'jet';            % colormap - see 'help colormap'
PData.logtransv = 'n';      % Do not log transform the data
PData.vmin = 0;             % min vertex value for plot color mapping
PData.vmax = max(max(V));   % max vertex value
% PData.logtransl = 'n';    % We're not plotting link data so these are ignored
PData.lmin = 0;
PData.lmax = max(max(L));
PData.lwidth = settings.lwidth;	% Width of links in points
PData.vsize = settings.vsize;	% Size of vertices in points (0 == omits verts)
PData.tsize = settings.tsize; 	% Size of tank/reservoir nodes
PData.legend = 'v';             % Show a colorbar legend for vertex data

SData = [];                 % No special node symbols
% SData(1...).ivs = [];             % Node IDs to annotate with special symbols
% SData(1...).vsmarker = [];        % Marker symbols
% SData(1...).vsmarkersize = [];    % Marker symbol sizes
% SData(1...).vsmarkercolor = [];   % Marker symbol colors