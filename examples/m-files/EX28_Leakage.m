%% FAVAD leakage (EPANET 2.3)
% Leakage model demo
 % Note: we can't just sum the leak rates at both end nodes together since 
 % in general the nodes can have leakage contributed by other connecting pipes.

clear; close all; clc;
start_toolkit;

%% Input file
INP = 'Net1.inp';
d = epanet(INP);

q_unit = d.getFlowUnits();
p_unit = d.getOptionsPressureUnits();
fprintf('[Units] Flow=%s | Pressure=%s\n', q_unit, p_unit);

%% Simulation settings
SIM_HOURS    = 24;
HYD_STEP_MIN = 60;

d.setTimeSimulationDuration(SIM_HOURS * 3600);
d.setTimeHydraulicStep(HYD_STEP_MIN * 60);

nLinks = d.getLinkCount();

%% Apply leakage to multiple links
% {LinkID, LC1, LC2}
leakLinks = { ...
    '21',  1.0,  0.10;
    '12',  0.8,  0.05;
    '31',  0.9,  0.08;
};

fprintf('\n===== APPLY LEAKAGE =====\n');

leakLinkIds     = leakLinks(:,1);
nLeakLinks      = size(leakLinks,1);
leakLinkIndices = zeros(nLeakLinks,1);

for k = 1:nLeakLinks
    linkID = leakLinks{k,1};
    LC1k   = leakLinks{k,2};
    LC2k   = leakLinks{k,3};

    linkIdx = d.getLinkIndex(linkID);
    if isempty(linkIdx) || linkIdx <= 0
        error('Link ID "%s" not found.', linkID);
    end

    leakLinkIndices(k) = linkIdx;

    d.setLinkLeakArea(linkIdx, LC1k);
    d.setLinkExpansionProperties(linkIdx, LC2k);

    fprintf('Applied leakage to link %s (index %d): LC1=%.6g, LC2=%.6g\n', ...
        linkID, linkIdx, LC1k, LC2k);
end

%% Run hydraulics and retrieve time series
res = d.getComputedHydraulicTimeSeries( ...
    'linkleakagerate', 'nodeleakageflow', 'emitterflow', ...
    'demanddelivered', 'demandrequested', 'time');

time_s = res.Time(:);
time_h = time_s / 3600;

leak_link_ts = res.LinkLeakageRate;   % [time x links] (L/s)
leak_node_ts = res.NodeLeakageFlow;   % [time x nodes] (L/s)

%% Core results: totals, averages, peak, volumes
totalLeak_ts_LPS = sum(leak_link_ts, 2, 'omitnan');   % L/s
avgTotalLeak_LPS = mean(totalLeak_ts_LPS, 'omitnan'); % L/s

if numel(time_s) >= 2
    totalLeakedVolume_m3 = trapz(time_s, totalLeak_ts_LPS) / 1000; % (L/s)*s -> L -> m^3
else
    totalLeakedVolume_m3 = NaN;
end

[peakTotalLeak_LPS, idxPeak] = max(totalLeak_ts_LPS);
tPeak_h = time_h(idxPeak);

fprintf('\n===== LEAKAGE SUMMARY =====\n');
fprintf('Average total leakage   : %.3f L/s\n', avgTotalLeak_LPS);
fprintf('Peak total leakage      : %.3f L/s at t = %.2f h\n', peakTotalLeak_LPS, tPeak_h);
fprintf('Total leaked volume     : %.2f m3/day (over %d h)\n', totalLeakedVolume_m3, SIM_HOURS);

%% Per-link mean leakage and peak snapshot
avgLeak_perLink_LPS = mean(leak_link_ts, 1, 'omitnan')';   % [links x 1]
leak_at_peak_LPS    = leak_link_ts(idxPeak, :)';           % [links x 1]

%% Links by mean leakage
topN = 3;

[sortedVals, sortedIx] = sort(avgLeak_perLink_LPS, 'descend');
topN = min(topN, nLinks);

topLinkIndex = sortedIx(1:topN);
topLinkLeak  = sortedVals(1:topN);

topLinkID = strings(topN,1);
for i = 1:topN
    idx = topLinkIndex(i);
    topLinkID(i) = string(d.getLinkNameID(idx));
end

fprintf('\nTop-%d links by MEAN leakage (L/s):\n', topN);
Ttop = table(topLinkIndex, topLinkID, topLinkLeak, ...
    'VariableNames', {'LinkIndex','LinkID','MeanLeak_LPS'});
disp(Ttop);

%% Consistency checks for ALL leaky links (nodes for each link)
fprintf('\n===== CONSISTENCY CHECKS (per leaky link) =====\n');

tol_L = 0.01;

checkTable = table('Size',[nLeakLinks 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'LinkID','LinkIndex','Node1','Node2','LeakLink_Total','LeakNodes_Total','AbsDiff'});

for k = 1:nLeakLinks
    linkID  = string(leakLinkIds{k});
    linkIdx = leakLinkIndices(k);

    nodes = d.getLinkNodesIndex(linkIdx);
    node1 = nodes(1);
    node2 = nodes(2);

    leak_nodes_total = sum(leak_node_ts(:, nodes), 'all', 'omitnan');
    leak_link_total  = sum(leak_link_ts(:, linkIdx), 'all', 'omitnan');
    absDiff          = abs(leak_link_total - leak_nodes_total);

    checkTable.LinkID(k)          = linkID;
    checkTable.LinkIndex(k)       = linkIdx;
    checkTable.Node1(k)           = node1;
    checkTable.Node2(k)           = node2;
    checkTable.LeakLink_Total(k)  = leak_link_total;
    checkTable.LeakNodes_Total(k) = leak_nodes_total;
    checkTable.AbsDiff(k)         = absDiff;

    fprintf('Link %s (idx %d) nodes (%d,%d): LinkSum=%.4f | NodeSum=%.4f | Diff=%.4f\n', ...
        linkID, linkIdx, node1, node2, leak_link_total, leak_nodes_total, absDiff);

    assert(absDiff < tol_L, 'Leak mismatch for link %s: link != nodes', linkID);
end

disp(checkTable);

%% Plots
figure('Name','Leakage Results','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) Total leakage over time
nexttile;
plot(time_h, totalLeak_ts_LPS, 'LineWidth', 1.5);
grid on;
xlabel('Time (h)');
ylabel('Total leakage (L/s)');
title('Total leakage time series');

hold on;
plot(time_h(idxPeak), peakTotalLeak_LPS, 'o', 'LineWidth', 1.5);
text(time_h(idxPeak), peakTotalLeak_LPS, ...
    sprintf('  Peak %.2f L/s @ %.2f h', peakTotalLeak_LPS, tPeak_h), ...
    'VerticalAlignment','bottom');
hold off;

% 2) Leakage over time for each leaky link
nexttile;
hold on;
for k = 1:nLeakLinks
    linkID  = string(leakLinkIds{k});
    linkIdx = leakLinkIndices(k);
    plot(time_h, leak_link_ts(:, linkIdx), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Link %s', linkID));
end
grid on;
xlabel('Time (h)');
ylabel('Leakage (L/s)');
title('Leakage on selected links');
legend('Location','best');
hold off;

% 3) Links bar (mean leakage)
nexttile;
bar(topLinkLeak);
grid on;
xlabel('Rank');
ylabel('Mean leakage (L/s)');
title('Links by mean leakage');
xticks(1:topN);
xticklabels(string(topLinkID));
xtickangle(0);

% 4) Sum of leakage at end nodes (per leaky link)
nexttile;
hold on;
for k = 1:nLeakLinks
    linkID  = string(leakLinkIds{k});
    linkIdx = leakLinkIndices(k);
    nodes   = d.getLinkNodesIndex(linkIdx);
    plot(time_h, sum(leak_node_ts(:, nodes), 2, 'omitnan'), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Nodes of %s', linkID));
end
grid on;
xlabel('Time (h)');
ylabel('End-node leakage (L/s)');
title('Leakage at end nodes of selected links');
legend('Location','best');
hold off;

%% Close EPANET object
d.unload;
