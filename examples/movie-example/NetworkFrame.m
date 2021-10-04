function [NData,fig] = NetworkFrame(fig,vdata,ldata,...
    PData,SData,NData,d,NodeType,LinkType,varargin)
%
% NetworkFrame plots a static frame of a network graph described in Epanet
% input format, with nodes and/or links colored using the values in vdata
% and ldata.  On first call for a given network, the line and node objects
% are rendered, and the object handles stored so that subsequent calls are
% more efficient (e.g., for generating movie frames recording simulation
% results).
%
% Arguments
%
% fig
%           handle for existing figure, or [] for new figure.
% vdata         
%           vector of values used to color nodes, in Epanet node index
%           order.  Length of vdata must be the number of nodes for the
%           network.  If vdata = [] then nodes are not colored (black).
%           (note: the PData 'vsize' parameter can be set to zero, in which
%           case nodes will not be shown -- see below).
% ldata         
%           vector of values used to color links, in Epanet link index
%           order.  Length of ldata must be the number of links for the
%           network.  If ldata = [] then links are colored as the average
%           of their adjacent node values.
% InpFname
%           String containing name of Epanet format input file describing
%           the network geometry.
% PData
%           A structure containing information used for plotting, processed
%           during the first call for a given network (and ignored for
%           subsequent calls for the same network).  May be set to [] to
%           accept all defaults (described below), or values may be
%           specified for one or more of the following members.
%
%           c - colormap used to color nodes and links, as described for
%           function colormap(). Default is 'jet'.  See 'doc colormap' for
%           more information about colormaps supplied with Matlab.
%
%           logtransv - 'y' or 'n' to log10 transform the data in vdata
%           prior to plotting. Default is to leave data untransformed
%           ('n').
%
%           vmin,vmax - Minimum and maximum values of vdata that will be
%           colored - others are black. Default is vmin=min(vdata) and
%           vmax=max(vdata), for the node values vdata specified on first
%           call for a given network.
%
%           logtransl - 'y' or 'n' to log10 transform the data in ldata
%           prior to plotting. Default is to leave data untransformed
%           ('n').
%
%           lmin,lmax - Minimum and maximum values of ldata that will be
%           colored - others are black. Default is lmin=min(ldata) and
%           lmax=max(ldata), for the link values ldata specified on first
%           call for a given network.
%
%           lwidth - link line width, in points.  Default is lwidth=3 pts.
%
%           vsize - node size, in points.  Set to zero to hide all junction
%           nodes and plot only links.  Default is vsize=4 pts.
%
%           tsize - tank node size, in points.  Set to zero to hide tank
%           and reservoir nodes.  Default is tsize=4 pts.
%
%           legend - 'l', 'v', or 'n' to add a link colorbar legend on the
%           right side of the plot ('l'), a vertex node colorbar legend
%           ('v'), or no legend ('n'). Default is to exclude the legend
%           ('n').
%
%           DLLname - String containing name of Epanet DLL.  Default is
%           'epanet2'.
%           
%           Hname - String containing name of Epanet C header file.
%           Default is 'epanet2.h'.
% SData
%           A structure array containing information used to annotate the
%           network nodes with symbols.  Set to [] if the network is not
%           annotated. Valid members are as follows; all must be the same
%           length and corresponding entries describe the various
%           attributes describing the node symbols to be plotted:
%
%           SData(i).ivs -
%           ith Node ID to annotate with symbols.
%
%           SData(i).vsmarker -
%           Marker symbols corresponding to the node IDs defined in
%           SData(:).ivs.
%
%           SData(i).vsmarkersize -
%           Marker symbol sizes corresponding to the node IDs defined in 
%           SData(:).ivs.
%
%           SData(i).vsmarkercolor -
%           Marker symbol colors corresponding to the node IDs defined in 
%           SData(:).ivs.
% Ndata 
%           A structure containing various network data, set to [] on the
%           first use of NetworkFrame for a particular network, and
%           preserved across successive calls to NetworkFrame for the same
%           network (e.g., when using NetworkFrame to draw successive
%           frames in a movie using the same network but different node
%           values).  NData should not be modified between successive
%           calls.
%
%
% Jim Uber
% 8/15/2005
% modified 4/10/2009
%
% modified by Marios Kyriakou 25/09/2016

if ~isempty(varargin)
    hyd=varargin{1};
end

if isempty(fig)
    fig=figure;
    axis equal;
    axis off;
%    set(fig,'Color',[0 0 0]);
else
    figure(fig);
end

% Initial processing (first use only for a given network)
if isempty(NData)
    % Default Plot Parameters
    NData.c = 'default';
    NData.logtransv = 'n';
    NData.logtransl = 'n';
    NData.vmin = min(vdata);
    NData.vmax = max(vdata);
    NData.lmin = min(ldata);
    NData.lmax = max(ldata);
    NData.lwidth = 3;
    NData.vsize = 4;
    NData.tsize = 4;
    NData.legend = 'n';

    % User-specified Plot Parameters
    if isfield(PData,'c')         NData.c = PData.c; end
    if isfield(PData,'logtransv') NData.logtransv = PData.logtransv; end
    if isfield(PData,'logtransl') NData.logtransl = PData.logtransl; end
    if isfield(PData,'vmin')      NData.vmin = PData.vmin; end
    if isfield(PData,'vmax')      NData.vmax = PData.vmax; end
    if isfield(PData,'lmin')      NData.lmin = PData.lmin; end
    if isfield(PData,'lmax')      NData.lmax = PData.lmax; end
    if isfield(PData,'lwidth')    NData.lwidth = PData.lwidth; end
    if isfield(PData,'vsize')     NData.vsize = PData.vsize; end
    if isfield(PData,'tsize')     NData.tsize = PData.tsize; end
    if isfield(PData,'legend')    NData.legend = PData.legend; end

    % Colormaps
    colormap(NData.c);
    NData.cmap = colormap;

    % Open Epanet
    % Network Size
    NData.nnodes = d.NodeCount;
    NData.nlinks = d.LinkCount;

    % Network topology
    tmp = d.getLinkNodesIndex;
    NData.from = tmp(:,1);
    NData.to = tmp(:,2);
    
    Itank = false(NData.nnodes,1);
    tmp_nodetypes = d.getNodeTypeIndex;
    Itank(find(tmp_nodetypes~=0))=true;
    
    % Node indices for extra node symbol IDs
    if ~isempty(SData)
        Sindex=zeros(length(SData),1);
        for i=1:length(SData)
            Sindex(i) = d.getNodeIndex(SData(i).ivs);
        end
    end

    % Network Geometry (from processing text InpFname)
    vx = d.NodeCoordinates{1};
    vy = d.NodeCoordinates{2};
    vertx = d.NodeCoordinates{3};
    verty = d.NodeCoordinates{4};

    % Link handles
    x(1,:) = vx(NData.from);
    x(2,:) = vx(NData.to);
    y(1,:) = vy(NData.from);
    y(2,:) = vy(NData.to);
    for i=1:NData.nlinks
        NData.linkh(i) = ...
            line([x(1,i) vertx{i} x(2,i)],[y(1,i) verty{i} y(2,i)],'LineWidth',NData.lwidth);
    end

    % Node handles
    if NData.vsize > 0
        NData.nodeh=zeros(NData.nnodes,1);
        for i=1:NData.nnodes
            NData.nodeh(i)=line(vx(i),vy(i),'Marker','o','MarkerSize',NData.vsize);
        end
        NData.nodeh = NData.nodeh';
        pval = cell(size(NData.nodeh(Itank)))';
        pname(1)={'Marker'};
        pname(2)={'MarkerSize'};
        pval(:,1) = {'s'};
        pval(:,2) = {NData.tsize};
        set(NData.nodeh(Itank),pname,pval);
    end

    % Add colorbar legend
    logtrans='n';
    minvalL=NData.lmin;
    maxvalL=NData.lmax;
    if NData.logtransl=='y' 
        minvalL=log10(minval);
        maxvalL=log10(maxval);
        logtrans='y'; 
    end
    minval=NData.vmin;
    maxval=NData.vmax;
    if NData.logtransv=='y' 
        minval=log10(minval);
        maxval=log10(maxval);
        logtrans='y'; 
    end
        
    [mapsize,n] = size(NData.cmap);
    ytick = [1: floor(mapsize/4) : mapsize];   % ytick holds the colormap positions to label
    [~,ticksize] = size(ytick);
    if ytick(ticksize) ~= mapsize
        ticksize = ticksize + 1;
        ytick(ticksize) = mapsize;
    end

    unts_node=[]; unts_link=[];
    if hyd==1
        % Right hand side legend link
        NData.hvcL = colorbar('WestOutside');
        dv = (maxvalL - minvalL)/mapsize;
        labelvalue = minvalL + ((ytick - 1)*dv);
        if length(labelvalue)==length(ytick)
            labelvalue = [labelvalue labelvalue(end)+labelvalue(2)-labelvalue(1)];
        end
        labelstring = num2str( labelvalue', 5 );
        try
            set(NData.hvcL,'ticks',1:length(labelvalue));
            set(NData.hvcL,'ticklabels',cellstr(labelstring));
            set(gca, 'clim', [0.5 length(labelvalue)+0.5]);
        catch e
            labelvalue = [labelvalue labelvalue(end)+labelvalue(2)-labelvalue(1)];
            set(NData.hvcL,'yticklabel', labelvalue(2:end)');
        end
        
        try unts_node = eval(['d.Node',NodeType,'Units{1}']);
        catch e, unts_node = eval(['d.Node',NodeType,'Units']); end
        try unts_link = eval(['d.Link',LinkType,'Units{1}']);
        catch e, unts_link = eval(['d.Link',LinkType,'Units']); end
    else
        % Units for NodeType quality
        if ~isempty(d.MSXFile)
            sind = d.getMSXSpeciesIndex(NodeType);
            if sind, unts_node = d.MSXSpeciesUnits{sind}; end
        end
        if strcmp('CL2', NodeType), NodeType='Chlorine'; end
        if strcmp(d.QualityChemName, NodeType), unts_node = d.QualityChemUnits; end
    end
    % Right hand side legend node
    NData.hvc = colorbar('EastOutside');
    dv = (maxval - minval)/mapsize;
    labelvalue = minval + ((ytick - 1)*dv);
    if length(labelvalue)==length(ytick)
        labelvalue = [labelvalue labelvalue(end)+labelvalue(2)-labelvalue(1)];
    end
    labelstring = num2str( labelvalue', 5 );
    try
        set(NData.hvc,'ticks',1:length(labelvalue));
        set(NData.hvc,'ticklabels',cellstr(labelstring));
        set(gca, 'clim', [0.5 length(labelvalue)+0.5]);
    catch e
        labelvalue = [labelvalue labelvalue(end)+labelvalue(2)-labelvalue(1)];
        set(NData.hvc,'yticklabel', labelvalue(2:end)');
    end

    if isfield(NData,'hvcL'), ylabel(NData.hvcL, [LinkType,' (',unts_link,')'],'fontsize',12); end
    if isfield(NData,'hvc'), ylabel(NData.hvc, [NodeType,' (',unts_node,')'],'fontsize',12); end

    % Extra node symbols
    NData.snodeh=[];
    if ~isempty(SData)
        for i=1:length(Sindex)
            NData.snodeh(i)=line(vx(Sindex(i)),vy(Sindex(i)),'Marker',SData(i).vsmarker,...
                'MarkerSize',SData(i).vsmarkersize,...
                'MarkerFaceColor',SData(i).vsmarkercolor,...
                'MarkerEdgeColor',SData(i).vsmarkercolor);
        end
    end

end

% Check node and link colors
colornodes = false;
if ~isempty(vdata)
    colornodes = true;
    [m,n]=size(vdata);
    if n==1
        vdata=vdata';
    end
    if m*n ~= NData.nnodes
        disp 'vdata must be a vector of length equal to the number of nodes'
        return
    end
end
colorlinks = false;
if ~isempty(ldata)
    colorlinks = true;
    [m,n]=size(ldata);
    if n==1
        ldata=ldata';
    end
    if m*n ~= NData.nlinks
        disp 'ldata must be a vector of length equal to the number of links'
        return
    end
end

% Compute indices of nodes and links to color
vmin = NData.vmin;
vmax = NData.vmax;
Ivmin = vdata<vmin;
Ivmax = vdata>vmax;
Ivcolor = vdata>=vmin & vdata<=vmax;
lmin = NData.lmin;
lmax = NData.lmax;
Ilmin = ldata<lmin;
Ilmax = ldata>lmax;
Ilcolor = ldata>=lmin & ldata<=lmax;

% Transformation
if colornodes
    if NData.logtransv == 'y'
        vmin = log10(vmin);
        vmax = log10(vmax);
        vdata(Ivmin) = vmin;
        vdata(Ivmax) = vmax;
        vdata(Ivcolor) = log10(vdata(Ivcolor));
    else
        vdata(Ivmin) = vmin;
        vdata(Ivmax) = vmax;
    end
end
if colorlinks
    if NData.logtransl == 'y'
        lmin = log10(lmin);
        lmax = log10(lmax);
        ldata(Ilmin) = lmin;
        ldata(Ilmax) = lmax;
        ldata(Ilcolor) = log10(ldata(Ilcolor));
    else
        ldata(Ilmin) = lmin;
        ldata(Ilmax) = lmax;
    end
end

% COLOR LINKS
if colorlinks      % color by link data
    rgb = color(NData.cmap,ldata,lmin,lmax);       % Return the color from the map
    rgb(~Ilcolor,:) = rgb(~Ilcolor,:)*0;            % These links are black
    set(NData.linkh,{'Color'},num2cell(rgb,2));
elseif colornodes  % color by adjacent node data
    vavg = (vdata(NData.from) + vdata(NData.to))/2; % Average value of adjacent nodes
    rgb = color(NData.cmap,vavg,vmin,vmax);        % Return the color from the map
    Ilcolor = Ivcolor(NData.from)==1 & Ivcolor(NData.to)==1;
    rgb(~Ilcolor,:) = rgb(~Ilcolor,:)*0;            % These links are black
    set(NData.linkh,{'Color'},num2cell(rgb,2));
end

% COLOR JUNCTIONS
if NData.vsize > 0
    if colornodes
        rgb = color(NData.cmap,vdata,vmin,vmax);  % Return the color from the map
        rgb(~Ivcolor,:) = rgb(~Ivcolor,:)*0;  % nodes are black if no color
        set(NData.nodeh,{'MarkerFaceColor'},num2cell(rgb,2),{'MarkerEdgeColor'},num2cell(rgb,2));
    else
        rgb = zeros(NData.nnodes,3);          % nodes are black if no color
        set(NData.nodeh,{'MarkerFaceColor'},num2cell(rgb,2),{'MarkerEdgeColor'},num2cell(rgb,2));
    end
end

for i=1:length(NData.snodeh)
    set(NData.snodeh(i),'Marker',SData(i).vsmarker,...
        'MarkerSize',SData(i).vsmarkersize,...
        'MarkerFaceColor',SData(i).vsmarkercolor,...
        'MarkerEdgeColor',SData(i).vsmarkercolor);
end

function [rgb] = color(cmap,v,vmin,vmax)

[m,n] = size(cmap);                  % m is the number of colors in map
dv = (vmax-vmin)/m;                  % Divide interval into m bins
if dv > 0                            % Select the bin index where v falls
   i = max(ceil((v-vmin)/dv), 1);
   rgb = cmap(i,:);                  % Return the color from the map
else
   rgb = cmap(ones(size(v)),:);
end
