%% Read data
clear; clc;

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(fullfile(repo, 'analysis', '05_LMM'));

plotScope = "all";

dirData = fullfile(repo, 'analysis', '.cache', 'derived', '05_LMM', 'ITPC');
W = table();
for binding = ["action", "outcome"]
    C = readtable(fullfile(repo, 'results', '03_functional_regression', binding + "_binding.csv"), 'TextType', 'string');
    C = unique(C(:, {'region','primary_lower_ms','primary_upper_ms'}), 'rows', 'stable');
    C.binding = repmat(binding, height(C), 1);
    W = [W; C];
end
E = readtable(fullfile(dirData, 'windows_used.csv'), 'TextType', 'string');
W = sortrows(W, {'binding','region'});
E = sortrows(E, {'binding','region'});
assert(isequal(W.binding, E.binding) && isequal(W.region, E.region) && ...
       isequal(W.primary_lower_ms, E.primary_lower_ms) && ...
       isequal(W.primary_upper_ms, E.primary_upper_ms), 'ITPC windows differ; rerun analysis/05_LMM/02_itpc.py.');
for rr = 1:height(E)
    stem = erase(E.file(rr), '.csv');
    assert(isfile(fullfile(dirData, stem + '_angles.csv')) && ...
           isfile(fullfile(dirData, stem + '_itpc_mean.csv')), 'Missing ITPC export: %s', stem);
end

E.roi = erase(E.region, '-');
% Display all four bands for the manuscript regions.
E = E(ismember(E.roi,["PreSMAL","PrecentralR","MTGR"]),:);
selected = repelem(E(:,{'binding','roi','region','file'}),4,1);
selected.band = repmat(["Delta";"Theta";"Alpha";"Beta"],height(E),1);
plotRoot = fullfile(repo,'figures','05_LMM');
if ~exist(plotRoot,'dir'); mkdir(plotRoot); end
disp(selected);

%% Plot settings
cfgCommon = struct();
cfgCommon.returnFigureOnly = true;
sourceFigures = gobjects(0);
sourceLabels = strings(0);

cfgCommon.bandOrder = ["Delta","Theta","Alpha","Beta"];

cfgCommon.nBins = 10;

cfgCommon.histScaleMode  = "global_common";
cfgCommon.arrowScaleMode = "global_itpc";
cfgCommon.histRadiusFrac = 0.90;

cfgCommon.barAngleOffsetFrac = 0.080;
cfgCommon.barWidthFrac       = 0.70;

cfgCommon.barAlpha      = 0.678436;
cfgCommon.roseEdgeLW    = 1.00;
cfgCommon.roseEdgeAlpha = 1.00;

cfgCommon.outerCircleColor = [0.349 0.341 0.341];
cfgCommon.outerCircleLW    = 0.75;
cfgCommon.outerCircleStyle = "-";

cfgCommon.crossColor       = [0.243 0.227 0.224];
cfgCommon.crossLW          = 0.25;
cfgCommon.crossStyle       = "--";

cfgCommon.innerCircleColor = [0.243 0.227 0.224];
cfgCommon.innerCircleLW    = 0.50;
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
cfgCommon.subDotAlpha = 0.90;
cfgCommon.subDotEdge = [1 1 1];
cfgCommon.subDotEdgeLW = 0.18;

cfgCommon.fontName = "Arial";

cfgCommon.panelTitleFontSize = 6.0;
cfgCommon.panelTitleFontWeight = "bold";

cfgCommon.itpcTextFontSize = 5.5;

cfgCommon.mainTitleFontSize = 7.0;
cfgCommon.mainTitleFontWeight = "bold";

cfgCommon.legendFontSize = 7.0;

cfgCommon.panelWidthCm = 3.3;
cfgCommon.minWidthCm = 5.0;
cfgCommon.circleSizeCm = 2.10 * 1.12; % actual outer circle diameter = 21 mm
cfgCommon.rowHeightCm = 3.5;
cfgCommon.headerHeightCm = 1.0;
cfgCommon.sideMarginCm = .3;
cfgCommon.bottomMarginCm = .2;
cfgCommon.rlimMax = 0.125;

line_action_v = hex2rgb_local("#006AB6");
line_action_i = hex2rgb_local("#FF542D");
head_action_v = hex2rgb_local("#2222FF");

edge_action_v = hex2rgb_local("#1E90FF");
edge_action_i = hex2rgb_local("#FF7518");

ribbon_action_v = hex2rgb_local("#FFFFFF");
ribbon_action_i = hex2rgb_local("#FFF0D8");

point_action_v = hex2rgb_local("#1E90FF");
point_action_i = hex2rgb_local("#FF7518");

line_outcome_v = hex2rgb_local("#13AE67");
line_outcome_i = hex2rgb_local("#E65FC8");

edge_outcome_v = hex2rgb_local("#0BDA51");
edge_outcome_i = hex2rgb_local("#FF69B4");

ribbon_outcome_v = hex2rgb_local("#EEF9E8");
ribbon_outcome_i = hex2rgb_local("#FDEBF1");

point_outcome_v = hex2rgb_local("#0BDA51");
point_outcome_i = hex2rgb_local("#FF69B4");

%% Plot and save
cfgScale = attach_one_global_scale(cfgCommon, selected, dirData);

fprintf('\nGLOBAL scale for selected action + outcome figures:\n');
fprintf('  max bin frequency = %.4f\n', cfgScale.globalMaxProb);
fprintf('  max ITPC          = %.4f\n\n', cfgScale.globalMaxITPC);

cfg = cfgCommon;
cfg.globalMaxProb = cfgScale.globalMaxProb;
cfg.globalMaxITPC = cfgScale.globalMaxITPC;

cfg.condOrder = ["va_a","ia_a"];
cfg.condLabel = containers.Map({'va_a','ia_a'}, {'Voluntary','Involuntary'});

cfg.lineColor     = containers.Map({'va_a','ia_a'}, {line_action_v, line_action_i});
cfg.arrowHeadColor = containers.Map({'va_a','ia_a'}, {head_action_v, line_action_i});
cfg.roseEdgeAlphaByCondition = containers.Map({'va_a','ia_a'}, {0.90,0.80});
cfg.ribbonColor   = containers.Map({'va_a','ia_a'}, {ribbon_action_v, ribbon_action_i});
cfg.roseEdgeColor = containers.Map({'va_a','ia_a'}, {edge_action_v, edge_action_i});
cfg.pointColor    = containers.Map({'va_a','ia_a'}, {point_action_v, point_action_i});

saveDir = plotRoot;

jobs = E(E.binding == "action" & ismember(E.roi, selected.roi(selected.binding == "action")), :);
for rr = 1:height(jobs)
    bands = selected.band(selected.binding == "action" & selected.roi == jobs.roi(rr));
    cfg.bandOrder = cfgCommon.bandOrder(ismember(cfgCommon.bandOrder, bands));
    if plotScope == "all"
        cfg.outputStem = jobs.roi(rr) + "_ITPC_polar";
    end
    stem = erase(jobs.file(rr), '.csv');
    [~,sourceFigures(end+1)] = plot_itpc_rose_cartesian(stem,saveDir,dirData,jobs.region(rr),cfg);
    sourceLabels(end+1) = jobs.region(rr);
end

cfg = cfgCommon;
cfg.globalMaxProb = cfgScale.globalMaxProb;
cfg.globalMaxITPC = cfgScale.globalMaxITPC;

cfg.condOrder = ["vs_s","is_s"];
cfg.condLabel = containers.Map({'vs_s','is_s'}, {'Voluntary','Involuntary'});

cfg.lineColor     = containers.Map({'vs_s','is_s'}, {line_outcome_v, line_outcome_i});
cfg.arrowHeadColor = containers.Map({'vs_s','is_s'}, {line_outcome_v, line_outcome_i});
cfg.roseEdgeAlphaByCondition = containers.Map({'vs_s','is_s'}, {1.00,0.90});
cfg.ribbonColor   = containers.Map({'vs_s','is_s'}, {ribbon_outcome_v, ribbon_outcome_i});
cfg.roseEdgeColor = containers.Map({'vs_s','is_s'}, {edge_outcome_v, edge_outcome_i});
cfg.pointColor    = containers.Map({'vs_s','is_s'}, {point_outcome_v, point_outcome_i});

saveDir = plotRoot;

jobs = E(E.binding == "outcome" & ismember(E.roi, selected.roi(selected.binding == "outcome")), :);
for rr = 1:height(jobs)
    bands = selected.band(selected.binding == "outcome" & selected.roi == jobs.roi(rr));
    cfg.bandOrder = cfgCommon.bandOrder(ismember(cfgCommon.bandOrder, bands));
    if plotScope == "all"
        cfg.outputStem = jobs.roi(rr) + "_ITPC_polar";
    end
    stem = erase(jobs.file(rr), '.csv');
    [~,sourceFigures(end+1)] = plot_itpc_rose_cartesian(stem,saveDir,dirData,jobs.region(rr),cfg);
    sourceLabels(end+1) = jobs.region(rr);
end

[present,order] = ismember(["PreSMA-L","Precentral-R","MTG-R"],sourceLabels);
assert(all(present));
combine_itpc_figures(sourceFigures(order),["PreSMA-L","PreCG-R","MTG-R"],fullfile(plotRoot,'itpc_four_bands'));
disp('Done.');

function cfg = attach_one_global_scale(cfg, selected, baseDir)
edges = linspace(0, 2*pi, cfg.nBins + 1);
cfg.globalMaxProb = eps;
cfg.globalMaxITPC = eps;
for rr = 1:height(selected)
    stem = erase(selected.file(rr), '.csv');
    A = readtable(fullfile(baseDir,stem + "_angles.csv"),'TextType','string');
    I = readtable(fullfile(baseDir,stem + "_itpc_mean.csv"),'TextType','string');
    A = A(A.band == selected.band(rr), :);
    I = I(I.band == selected.band(rr), :);
    assert(~isempty(A) && ~isempty(I), 'Missing selected band: %s', selected.band(rr));
    conds = unique(A.condition);
    for cc = 1:numel(conds)
        th = A.angle_rad(A.condition == conds(cc));
        counts = histcounts(th, edges);
        prob = counts / max(sum(counts),1);
        cfg.globalMaxProb = max(cfg.globalMaxProb,max(prob));
        vals = I.itpc_mean(I.condition == conds(cc));
        assert(~isempty(vals) && all(isfinite(vals)), 'Missing selected condition ITPC.');
        cfg.globalMaxITPC = max(cfg.globalMaxITPC,mean(vals));
    end
end
end

function rgb = hex2rgb_local(hexColor)

hexColor = char(hexColor);

if startsWith(hexColor, '#')
    hexColor = hexColor(2:end);
end

rgb = sscanf(hexColor, '%2x%2x%2x', [1 3]) / 255;

end
