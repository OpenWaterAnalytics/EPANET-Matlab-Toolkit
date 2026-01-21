d = epanet('ky10.inp');

% Time settings
TIME_H     = 24 * 3600;    % 24 hours
HYD_STEP   = 15 * 60;      % 15 min
PAT_STEP   = 15 * 60;
RPT_STEP   = 15 * 60;

LINK_LEAK_AREA = 3e-7;
EMITTER_COEFF  = 0.002;

d.setTimeSimulationDuration(TIME_H);
d.setTimeHydraulicStep(HYD_STEP);
d.setTimePatternStep(PAT_STEP);
d.setTimeReportingStep(RPT_STEP);

linkID  = d.getLinkNameID(10);
pipe_ix = d.getLinkIndex(linkID);

nodeID  = d.getNodeNameID(10);
junc_ix = d.getNodeIndex(nodeID);

% Apply leakage settings
d.setLinkLeakArea(pipe_ix, LINK_LEAK_AREA);
d.setLinkExpansionProperties(pipe_ix, 0.5);   % example elasticity 
d.setNodeEmitterCoeff(junc_ix, EMITTER_COEFF);

p_unit = d.getOptionsPressureUnits();
q_unit = d.getFlowUnits();

fprintf('[Units] Flow=%s | Pressure=%s\n', q_unit, p_unit);
fprintf('[Params] LinkLeakArea=%.6g m^2 | EmitterCoeff=%.6g\n', LINK_LEAK_AREA, EMITTER_COEFF);

% Run hydraulics
d.openHydraulicAnalysis();
d.initializeHydraulicAnalysis();

t = 0.0;
while true
    d.runHydraulicAnalysis();

    % Leakages
    leak_pipe = d.getLinkLeakageRate(pipe_ix);
    leak_node = d.getNodeLeakageFlow(junc_ix);
    emitter   = d.getNodeEmitterFlow(junc_ix);

    % Demand
    dem_req   = d.getConsumerDemandRequested(junc_ix);
    dem_del   = d.getConsumerDemandDelivered(junc_ix);

    % Diagnostics
    Pj = d.getNodePressure(junc_ix);

    total_out  = max(dem_del + leak_pipe + leak_node + emitter, 1e-12);
    leak_share = 100.0 * (leak_pipe + leak_node + emitter) / total_out;

    fprintf(['t=%5.0fs | P=%7.3f %s | PipeLeak=%9.4f | NodeLeak=%9.4f | ' ...
             'Emitter=%8.4f | DemReq=%8.3f | DemDel=%8.3f | LeakShare=%5.1f%%\n'], ...
             t, Pj, p_unit, leak_pipe, leak_node, emitter, dem_req, dem_del, leak_share);

    dt = d.nextHydraulicAnalysisStep();
    if dt <= 0
        break
    end
    t = t + dt;
end

d.closeHydraulicAnalysis();
d.unload();