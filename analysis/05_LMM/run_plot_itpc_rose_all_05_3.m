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

dirData = fullfile(repo, 'data', '05_LMM', 'plot_data_for_itpc');

% =========================================================================
% Common configuration
% =========================================================================
cfgCommon = struct();

cfgCommon.bandOrder = ["Delta","Theta","Alpha","Beta"];

cfgCommon.nBins = 10;

% Use one shared histogram scale and one shared ITPC arrow scale
% across the final action and outcome figures.
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
% =========================================================================

% Action: blue vs orange
line_action_v = hex2rgb_local("#2A95E8");
line_action_i = hex2rgb_local("#FFAA3D");

edge_action_v = hex2rgb_local("#6FC1F6");
edge_action_i = hex2rgb_local("#FFC067");

ribbon_action_v = hex2rgb_local("#E8F5FF");
ribbon_action_i = hex2rgb_local("#FFF0D8");

point_action_v = hex2rgb_local("#7FC8F8");
point_action_i = hex2rgb_local("#FFC874");

% Outcome: green vs pink
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
% This excludes PreSMAL and excludes all probability-related figures.
% =========================================================================
allGlobalStems = [ ...
    "PreSMAR_action_ba", ...
    "SMAR_action_ba", ...
    "PostcentralL_outcome_bs", ...
    "MTGR_outcome_bs" ...
    ];

cfgScale = attach_one_global_scale(cfgCommon, allGlobalStems, dirData);

fprintf('\nGLOBAL scale for selected action + outcome figures:\n');
fprintf('  max bin frequency = %.4f\n', cfgScale.globalMaxProb);
fprintf('  max ITPC          = %.4f\n\n', cfgScale.globalMaxITPC);

% =========================================================================
% Action binding: PreSMAR / SMAR
% =========================================================================
cfg = cfgCommon;
cfg.globalMaxProb = cfgScale.globalMaxProb;
cfg.globalMaxITPC = cfgScale.globalMaxITPC;

cfg.condOrder = ["va_a","ia_a"];
cfg.condLabel = containers.Map({'va_a','ia_a'}, {'Voluntary','Involuntary'});

cfg.lineColor     = containers.Map({'va_a','ia_a'}, {line_action_v, line_action_i});
cfg.ribbonColor   = containers.Map({'va_a','ia_a'}, {ribbon_action_v, ribbon_action_i});
cfg.roseEdgeColor = containers.Map({'va_a','ia_a'}, {edge_action_v, edge_action_i});
cfg.pointColor    = containers.Map({'va_a','ia_a'}, {point_action_v, point_action_i});

saveDir = fullfile(repo, 'figures', '05_LMM', 'ITPC', 'action_binding');

plot_itpc_rose_cartesian("PreSMAR_action_ba", saveDir, dirData, "PreSMAR", cfg);
plot_itpc_rose_cartesian("SMAR_action_ba",    saveDir, dirData, "SMAR",   cfg);

% =========================================================================
% Outcome binding: PostcentralL / MTGR
% =========================================================================
cfg = cfgCommon;
cfg.globalMaxProb = cfgScale.globalMaxProb;
cfg.globalMaxITPC = cfgScale.globalMaxITPC;

cfg.condOrder = ["vs_s","is_s"];
cfg.condLabel = containers.Map({'vs_s','is_s'}, {'Voluntary','Involuntary'});

cfg.lineColor     = containers.Map({'vs_s','is_s'}, {line_outcome_v, line_outcome_i});
cfg.ribbonColor   = containers.Map({'vs_s','is_s'}, {ribbon_outcome_v, ribbon_outcome_i});
cfg.roseEdgeColor = containers.Map({'vs_s','is_s'}, {edge_outcome_v, edge_outcome_i});
cfg.pointColor    = containers.Map({'vs_s','is_s'}, {point_outcome_v, point_outcome_i});

saveDir = fullfile(repo, 'figures', '05_LMM', 'ITPC', 'outcome_binding');

plot_itpc_rose_cartesian("PostcentralL_outcome_bs", saveDir, dirData, "PostcentralL", cfg);
plot_itpc_rose_cartesian("MTGR_outcome_bs",         saveDir, dirData, "MTGR",        cfg);

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
% Get the global maximum ITPC value from exported TF-style ITPC tables
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

    if ~ismember("itpc_mean", string(I.Properties.VariableNames))
        error("itpc_mean.csv must contain an itpc_mean column.");
    end

    conds = unique(I.condition);

    for jj = 1:numel(bandOrder)
        b = bandOrder(jj);

        for cc = 1:numel(conds)
            c = conds(cc);

            row = I((I.condition == c) & (I.band == b), :);

            if ~isempty(row)
                val = row.itpc_mean(1);
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
