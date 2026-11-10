%% Read data
clear; clc;

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(fullfile(here, 'lib'));
W = readtable(fullfile(repo, 'analysis', '.cache', 'derived', '06_ieeg', 'windows_used.csv'));
figureTitle = "MTG-R";

dirData = fullfile(repo, 'analysis', '.cache', 'derived', '06_ieeg', 'ITPC');

%% Plot settings
cfgCommon = struct();

cfgCommon.bandOrder = ["Delta","Theta","Alpha","Beta"];

cfgCommon.nBins = 10;

cfgCommon.histScaleMode  = "global_common";
cfgCommon.arrowScaleMode = "global_itpc";
cfgCommon.histRadiusFrac = 0.90;

cfgCommon.barAngleOffsetFrac = 0.080;
cfgCommon.barWidthFrac       = 0.70;

cfgCommon.barAlpha      = 0.68;
cfgCommon.roseEdgeLW    = 0.65;
cfgCommon.roseEdgeAlpha = 1.00;

cfgCommon.outerCircleColor = [0.10 0.10 0.10];
cfgCommon.outerCircleLW    = 0.75;
cfgCommon.outerCircleStyle = "-";

cfgCommon.crossColor       = [0.86 0.86 0.86];
cfgCommon.crossLW          = 0.55;
cfgCommon.crossStyle       = "--";

cfgCommon.innerCircleColor = [0.90 0.90 0.90];
cfgCommon.innerCircleLW    = 0.45;
cfgCommon.innerCircleStyle = ":";
cfgCommon.innerCircleFracs = [0.45, 0.70];

cfgCommon.showAngleLabels = false;

cfgCommon.arrowMaxFrac = 0.95;
cfgCommon.arrowLW = 1.4;
cfgCommon.arrowHeadLengthFrac = .16;
cfgCommon.arrowHeadWidthFrac  = .11;

cfgCommon.showSubDots = true;
cfgCommon.dotRadius = 1.018;
cfgCommon.subDotSize = 4.0;
cfgCommon.subDotAlpha = 0.88;
cfgCommon.subDotEdge = [1 1 1];
cfgCommon.subDotEdgeLW = 0.18;

cfgCommon.fontName = "Arial";

cfgCommon.angleFontSize = 8.5;
cfgCommon.angleFontWeight = "normal";

cfgCommon.panelTitleFontSize = 6.0;
cfgCommon.panelTitleFontWeight = "bold";

cfgCommon.itpcTextFontSize = 6.0;

cfgCommon.mainTitleFontSize = 7.0;
cfgCommon.mainTitleFontWeight = "bold";

cfgCommon.legendFontSize = 7.0;

cfgCommon.panelWidthCm = 4.0;
cfgCommon.minWidthCm = 5.0;
cfgCommon.circleSizeCm = 2.72 * 1.12; % outer circle diameter from Fig5.ai
cfgCommon.rowHeightCm = 4.1;
cfgCommon.headerHeightCm = 1.0;
cfgCommon.sideMarginCm = .3;
cfgCommon.bottomMarginCm = .2;
cfgCommon.rlimMax = 0.125;

line_outcome_v = hex2rgb_local("#13AE67");
line_outcome_i = hex2rgb_local("#E65FC8");

edge_outcome_v = hex2rgb_local("#0BDA51");
edge_outcome_i = hex2rgb_local("#FF69B4");

ribbon_outcome_v = hex2rgb_local("#EEF9E8");
ribbon_outcome_i = hex2rgb_local("#FDEBF1");

point_outcome_v = hex2rgb_local("#0BDA51");
point_outcome_i = hex2rgb_local("#FF69B4");

allGlobalStems = ["ieeg_MTGR"];

%% Plot and save
cfgScale = attach_one_global_scale(cfgCommon, allGlobalStems, dirData);

fprintf('\nGLOBAL scale for iEEG MTGR:\n');
fprintf('  max bin frequency = %.4f\n', cfgScale.globalMaxProb);
fprintf('  max ITPC          = %.4f\n\n', cfgScale.globalMaxITPC);

cfg = cfgCommon;
cfg.outputStem = "ieeg_itpc_four_bands";
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

saveDir = fullfile(repo, 'figures', '06_ieeg');

plot_itpc_rose_cartesian_ieeg("ieeg_MTGR", saveDir, dirData, figureTitle, cfg);

disp('Done.');

function cfg = attach_one_global_scale(cfg, regionStems, baseDir)

edges = linspace(0, 2*pi, cfg.nBins + 1);

cfg.globalMaxProb = get_global_max_prob_mixed(regionStems, baseDir, cfg.bandOrder, edges);
cfg.globalMaxITPC = get_global_max_itpc_mixed(regionStems, baseDir, cfg.bandOrder);

end

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

function rgb = hex2rgb_local(hexColor)

hexColor = char(hexColor);

if startsWith(hexColor, '#')
    hexColor = hexColor(2:end);
end

rgb = sscanf(hexColor, '%2x%2x%2x', [1 3]) / 255;

end
