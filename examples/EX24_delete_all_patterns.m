%% Delete all patterns from network file
% This example contains:
%   Load network.
%   Delete all patterns.
%   Save new file without patterns.
%   Unload library.

%%
% Clear
clear; close('all'); clc;
start_toolkit;

% Load network
d = epanet('Net1.inp');

% Delete all patterns
while ~isempty(d.getPatternIndex)
    indexPat = d.getPatternIndex;
    d.deletePattern(indexPat(1));
end

% Save new file without patterns
d.saveInputFile('No_patterns.inp');

%   Unload library
d.unload;