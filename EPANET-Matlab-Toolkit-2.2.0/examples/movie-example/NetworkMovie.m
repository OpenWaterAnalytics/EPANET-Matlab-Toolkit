function [NData] = NetworkMovie(V,L,fig,movfname,quality,fps,...
    PData,SData,d,NodeType,LinkType,varargin)
%% Synopsis
% NetworkMovie uses NetworkFrame() to plot the network graph on figure
% handle fig with nodes and links colored by V,L, and store the animation
% as an *.avi movie.  For detailed descriptions of the PData and SData
% arguments, see NetworkFrame().
%
% After the first frame is plotted (using NetworkFrame()), control is
% passed to the standard input, and during this time you can use the
% standard Matlab figure UI controls to modify the figure view as you
% desire (e.g., zoom, pan).  After the network view is how you want it, the
% rest of the frames are generated consistent with that view.
%
%% Helpful Notes: 
% * Do not move or obscure the figure window after you position it
% and the figure within, when the movie frames are
% being rendered.  Otherwise you will have whatever junk covers the figure
% window in your movie.  If you actually resize the figure window once the
% movie frames are being rendered, you'll get an error message.
% * To increase your resolution try increasing the physical size of the
% figure window.  Remember that what you see will be what you get.
% * If a problem occurs when rendering the frames then Matlab seems to
% stubbornly hang onto the *.avi movie file handle despite, for example,
% issuing fclose('all').  The fix is 'clear all' and then you'll be able to
% either delete or rewrite the same *.avi filename without shutting down
% matlab first.
% * The exampleMovie.m (or exampleMovie.html) code shows a complete example
% of how to create a movie of Epanet or MSX simulation data, and for
% standard use you should be able to modify some of those parameters for
% your application.
%
%% Arguments
%  V           (nnodes x nframes) matrix of node values to plot
%  L           (nlinks x nframes) matrix of link values to plot
%              Either of V or L can be blank, in which case behavior is as
%              described in NetworkFrame().
%  fig         Existing figure handle to use.  If [] then new figure is 
%              created.
%  movfname    AVI filename to store movie.  If [] then no movie is made,
%              but the animation is displayed.
%  quality     0-100 quality of AVI movie
%  fps         frames per second of AVI movie
%  InpFname    Epanet input file name describing network geometry.
%  PData       PData structure of plot options as described in
%              NetworkFrame().
%  SData       SData structure of plot symbols as described in
%              NetworkFrame().
%
% Jim Uber
% 8/15/2005
% 4/10/2009 modified
%
% modified by Marios Kyriakou 25/09/2016

if ~isempty(varargin)
    hyd=varargin{1};
end

nframes = 0;
if isempty(V)
    colorV=false;
else
    colorV=true;
    nframes = length(V(1,:));
end
if isempty(L)
    colorL=false;
else
    colorL=true;
    nframes = length(L(1,:));
end

if colorV && colorL
    if length(V(1,:)) ~= length(L(1,:))
        disp 'Number of frames in V and L must agree.'
        return
    end
end

makeavi = ~isempty(movfname);

% Figure props
if isempty(fig)
    fig=figure;
    axis equal;
    axis off;
%    set(fig,'Color',[0 0 0]);
else
    figure(fig);
end
if makeavi
    mov=VideoWriter(movfname);%,'quality',quality,'fps',fps,'compression','Cinepak');
end

% Draw the first frame and stop to allow user adjustment of figure window
v = [];
l = [];
if colorV, v=V(:,1); end
if colorL, l=L(:,1); end
[NData,fig] = NetworkFrame(fig,v,l,PData,SData,[],d,NodeType,LinkType,hyd);

if makeavi
    mov.FrameRate=fps;
    mov.Quality=quality;
    open(mov)
    M=getframe(fig);
    writeVideo(mov,M);
else
    drawnow;
end

% Draw the rest
for i=2:nframes
    if colorV, v=V(:,i); end
    if colorL, l=abs(L(:,i)); end
    [NData,fig] = NetworkFrame(fig,v,l,PData,SData,NData,d,NodeType,LinkType,hyd);
    if makeavi
        M=getframe(fig);
        writeVideo(mov,M);
    else
        drawnow;
    end
end

if makeavi
    close(mov);
end
