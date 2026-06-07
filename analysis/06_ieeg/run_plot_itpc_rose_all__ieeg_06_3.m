clear; clc;

% =========================================================================
% Force white background globally
% =========================================================================
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');

% =========================================================================
% Paths
% =========================================================================
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));

dirData = fullfile(repo, 'data', '06_ieeg', 'plot_data_for_itpc');

% =========================================================================
% Common configuration
% =========================================================================
cfgCommon = struct();

cfgCommon.bandOrder = ["Delta","Theta","Alpha","Beta"];

cfgCommon.nBins = 10;

% Use one shared histogram scale and one shared ITPC arrow scale.
% Here only ieeg_MTGR is included.
cfgCommon.histScaleMode  = "global_common";
cfgCommon.arrowScaleMode = "global_itpc";
cfgCommon.histRadiusFrac = 0.90;

% Slightly offset the two conditions within each angular bin.
cfgCommon.barAngleOffsetFrac = 0.080;
cfgCommon.barWidthFrac       = 0.70;

% Rose bars
cfgCommon.barAlpha      = 0.68;
cfgCommon.roseEdgeLW    = 0.95;
cfgCommon.roseEdgeAlpha = 1.00;

% Polar frame
cfgCommon.outerCircleColor = [0.10 0.10 0.10];
cfgCommon.outerCircleLW    = 0.90;
cfgCommon.outerCircleStyle = "-";

cfgCommon.crossColor       = [0.86 0.86 0.86];
cfgCommon.crossLW          = 0.55;
cfgCommon.crossStyle       = "--";

cfgCommon.innerCircleColor = [0.90 0.90 0.90];
cfgCommon.innerCircleLW    = 0.45;
cfgCommon.innerCircleStyle = ":";
cfgCommon.innerCircleFracs = [0.45, 0.70];

cfgCommon.showAngleLabels = false;

% ITPC arrows
cfgCommon.arrowMaxFrac = 0.95;
cfgCommon.arrowLW = 1.65;
cfgCommon.arrowHeadLengthFrac = 0.082;
cfgCommon.arrowHeadWidthFrac  = 0.052;

% Subject-level dots
cfgCommon.showSubDots = true;
cfgCommon.dotRadius = 1.018;
cfgCommon.subDotSize = 10.5;
cfgCommon.subDotAlpha = 0.88;
cfgCommon.subDotEdge = [1 1 1];
cfgCommon.subDotEdgeLW = 0.18;

% Fonts
cfgCommon.fontName = "Arial";

cfgCommon.angleFontSize = 8.5;
cfgCommon.angleFontWeight = "normal";

cfgCommon.panelTitleFontSize = 10.0;
cfgCommon.panelTitleFontWeight = "bold";

cfgCommon.itpcTextFontSize = 7.8;

cfgCommon.mainTitleFontSize = 12.2;
cfgCommon.mainTitleFontWeight = "bold";

cfgCommon.legendFontSize = 8.6;

% Figure size
cfgCommon.w_cm = 16.5;
cfgCommon.h_cm = 17.0;
cfgCommon.rlimMax = 0.125;

% =========================================================================
% Colors
% Outcome: green vs pink
% =========================================================================
line_outcome_v = hex2rgb_local("#83D768");
line_outcome_i = hex2rgb_local("#F98EAE");

edge_outcome_v = hex2rgb_local("#AFE28F");
edge_outcome_i = hex2rgb_local("#F6B7C9");

ribbon_outcome_v = hex2rgb_local("#EEF9E8");
ribbon_outcome_i = hex2rgb_local("#FDEBF1");

point_outcome_v = hex2rgb_local("#BDEB9E");
point_outcome_i = hex2rgb_local("#F8C3D3");

% =========================================================================
% Compute one global scale
% PostcentralL is excluded.
% =========================================================================
allGlobalStems = ["ieeg_MTGR"];

cfgScale = attach_one_global_scale(cfgCommon, allGlobalStems, dirData);

fprintf('\nGLOBAL scale for iEEG MTGR:\n');
fprintf('  max bin frequency = %.4f\n', cfgScale.globalMaxProb);
fprintf('  max ITPC          = %.4f\n\n', cfgScale.globalMaxITPC);

% =========================================================================
% iEEG outcome binding: MTGR only
% =========================================================================
cfg = cfgCommon;
cfg.globalMaxProb = cfgScale.globalMaxProb;
cfg.globalMaxITPC = cfgScale.globalMaxITPC;

cfg.condOrder = ["Voluntary_OutcomeOnsetReportSound", "Involuntary_OutcomeOnsetReportSound"];

cfg.condLabel = containers.Map( ...
    {'Voluntary_OutcomeOnsetReportSound', 'Involuntary_OutcomeOnsetReportSound'}, ...
    {'Voluntary', 'Involuntary'} ...
);

cfg.lineColor = containers.Map( ...
    {'Voluntary_OutcomeOnsetReportSound', 'Involuntary_OutcomeOnsetReportSound'}, ...
    {line_outcome_v, line_outcome_i} ...
);

cfg.ribbonColor = containers.Map( ...
    {'Voluntary_OutcomeOnsetReportSound', 'Involuntary_OutcomeOnsetReportSound'}, ...
    {ribbon_outcome_v, ribbon_outcome_i} ...
);

cfg.roseEdgeColor = containers.Map( ...
    {'Voluntary_OutcomeOnsetReportSound', 'Involuntary_OutcomeOnsetReportSound'}, ...
    {edge_outcome_v, edge_outcome_i} ...
);

cfg.pointColor = containers.Map( ...
    {'Voluntary_OutcomeOnsetReportSound', 'Involuntary_OutcomeOnsetReportSound'}, ...
    {point_outcome_v, point_outcome_i} ...
);

saveDir = fullfile(repo, 'figures', '06_ieeg', 'ITPC');

plot_itpc_rose_cartesian_ieeg("ieeg_MTGR", saveDir, dirData, "MTGR", cfg);

disp('Done.');

% =========================================================================
% Attach global histogram and ITPC scales
% =========================================================================
function cfg = attach_one_global_scale(cfg, regionStems, baseDir)

edges = linspace(0, 2*pi, cfg.nBins + 1);

cfg.globalMaxProb = get_global_max_prob_mixed(regionStems, baseDir, cfg.bandOrder, edges);
cfg.globalMaxITPC = get_global_max_itpc_mixed(regionStems, baseDir, cfg.bandOrder);

end

% =========================================================================
% Get the global maximum relative bin frequency
% =========================================================================
function globalMaxProb = get_global_max_prob_mixed(regionStems, baseDir, bandOrder, edges)

globalMaxProb = eps;

for rr = 1:numel(regionStems)
    anglesFile = fullfile(baseDir, string(regionStems(rr)) + "_angles.csv");

    if ~isfile(anglesFile)
        warning("Missing angles file: %s", anglesFile);
        continue;
    end

    A = readtable(anglesFile);
    A.condition = string(A.condition);
    A.band      = string(A.band);

    conds = unique(A.condition);

    for jj = 1:numel(bandOrder)
        b = bandOrder(jj);

        for cc = 1:numel(conds)
            c = conds(cc);

            idx = (A.condition == c) & (A.band == b);
            th = A.angle_rad(idx);

            counts = histcounts(th, edges);
            prob = counts / max(sum(counts), 1);

            if ~isempty(prob)
                globalMaxProb = max(globalMaxProb, max(prob));
            end
        end
    end
end

globalMaxProb = max(globalMaxProb, eps);

end

% =========================================================================
% Get the global maximum ITPC value
% =========================================================================
function globalMaxITPC = get_global_max_itpc_mixed(regionStems, baseDir, bandOrder)

globalMaxITPC = eps;

for rr = 1:numel(regionStems)
    itpcFile = fullfile(baseDir, string(regionStems(rr)) + "_itpc_mean.csv");

    if ~isfile(itpcFile)
        warning("Missing ITPC file: %s", itpcFile);
        continue;
    end

    I = readtable(itpcFile);
    I.condition = string(I.condition);
    I.band      = string(I.band);

    conds = unique(I.condition);

    for jj = 1:numel(bandOrder)
        b = bandOrder(jj);

        for cc = 1:numel(conds)
            c = conds(cc);

            row = I((I.condition == c) & (I.band == b), :);

            if ~isempty(row) && ismember("itpc_mean", string(row.Properties.VariableNames))
                vals = row.itpc_mean;
                vals = vals(isfinite(vals));

                if isempty(vals)
                    continue;
                end

                val = mean(vals);
                if isfinite(val)
                    globalMaxITPC = max(globalMaxITPC, val);
                end
            end
        end
    end
end

globalMaxITPC = max(globalMaxITPC, eps);

end

% =========================================================================
% Convert HEX color to RGB
% =========================================================================
function rgb = hex2rgb_local(hexColor)

hexColor = char(hexColor);

if startsWith(hexColor, '#')
    hexColor = hexColor(2:end);
end

rgb = sscanf(hexColor, '%2x%2x%2x', [1 3]) / 255;

end
