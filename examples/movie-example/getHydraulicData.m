%% Function getHydraulicData
function [V,L,T] = getHydraulicData(NodeType,LinkType,d)
%% Synopsis
%   [V,L,T,errcode] = getHydraulicData(NodeType,LinkType,inpFname) Returns
%   vertex (V) and link (L) data from running an Epanet hydraulic
%   simulation defined by the input file inpFname. The return value of V
%   contains the node data of type NodeType, while L contains the link data
%   of type LinkType. NodeType and LinkType are character strings as
%   defined in the Epanet toolkit documentation (i.e. 'EN_...').  The
%   ordering in V and/or L is node and/or link index order. The output
%   vector T contains the simulation times (in seconds) associated with
%   each frame.  If either of NodeType or LinkType are undefined then V or
%   L will be empty.
%
%   If defined, the output matrices V, L, and T are (nnodes x nframes),
%   (nlinks x nframes), and (nframes x 1) respectively.  The value of
%   nframes is defined by the values of REPORT START and REPORT STEP in the
%   Epanet input file.  Specifically, frames are stored only when Epanet
%   simulation time t satisfies: (t >= REPORT_START && mod(t,REPORT_STEP)
%   == 0).
% Jim Uber
%
% modified by Marios Kyriakou 29/09/2016

%% Check nodetype and linktype
if strcmp(lower(NodeType(1:4)),'pres')
    NodeType = 'Pressure';
end
if strcmp(lower(LinkType(1:4)),'flow')
    LinkType = 'Flow';
    linktype = 'Flows';
end
if strcmp(lower(NodeType(1:4)),'dema')
    NodeType = 'ActualDemand';
end
if strcmp(lower(LinkType(1:4)),'velo')
    LinkType = 'Velocity';
    linktype = LinkType;
end

%% Hydraulic Analysis
T = [];         % The time data
V = [];         % The vertex data
L = [];         % The link data
% Get the simulation time parameters
reportStart = d.getTimeReportingStart;
reportStep = d.getTimeReportingStep;

tleft = 1;      % Time left in simulation
iframe = 1;     % frame index
d.openHydraulicAnalysis;
d.initializeHydraulicAnalysis;
while (tleft > 0)
    t = d.runHydraulicAnalysis;
    if t >= reportStart && ~mod(t,reportStep)
        % Retrieve results at time t
        V = [V zeros(d.NodeCount,1)];
        for in=1:d.NodeCount
            V(in,iframe) = eval(['d.getNode',NodeType,'(in)']);
        end
        L = [L zeros(d.LinkCount,1)];
        for in=1:d.LinkCount
            L(in,iframe) = eval(['d.getLink',linktype,'(in)']);
        end
        T = [T; t];
        iframe = iframe + 1;
    end
    tleft = d.nextHydraulicAnalysisStep;
end
L=abs(L);
d.closeHydraulicAnalysis;

