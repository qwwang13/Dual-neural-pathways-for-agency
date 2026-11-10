function outPdf = plot_itpc_rose_cartesian_ieeg(regionStem, saveDir, baseDir, figureTitle, cfg)

if ~exist(saveDir, "dir")
    mkdir(saveDir);
end

cfg = fill_default_cfg(cfg);

anglesFile = fullfile(baseDir, string(regionStem) + "_angles.csv");
itpcFile   = fullfile(baseDir, string(regionStem) + "_itpc_mean.csv");

A = readtable(anglesFile);
I = readtable(itpcFile);

A.condition = string(A.condition);
A.band      = string(A.band);

I.condition = string(I.condition);
I.band      = string(I.band);

if ~ismember("subject", string(A.Properties.VariableNames))
    error("angles.csv must contain a subject column.");
end

if iscell(A.subject) || iscategorical(A.subject)
    A.subject = string(A.subject);
end

condOrder = string(cfg.condOrder);
bandOrder = string(cfg.bandOrder);

if numel(condOrder) ~= 2
    error("cfg.condOrder must contain exactly two conditions.");
end

c1 = condOrder(1);
c2 = condOrder(2);

nBins = cfg.nBins;
edges = linspace(0, 2*pi, nBins + 1);

rlimMax  = cfg.rlimMax;
rHistMax = cfg.histRadiusFrac * rlimMax;

wrap2pi = @(x) mod(x, 2*pi);

if ~isfield(cfg, 'globalMaxProb') || isempty(cfg.globalMaxProb) || ~isfinite(cfg.globalMaxProb)
    error("Missing cfg.globalMaxProb. Please compute the global histogram scale in the run script first.");
end

if ~isfield(cfg, 'globalMaxITPC') || isempty(cfg.globalMaxITPC) || ~isfinite(cfg.globalMaxITPC)
    error("Missing cfg.globalMaxITPC. Please compute the global ITPC scale in the run script first.");
end

cfg.globalMaxProb = max(cfg.globalMaxProb, eps);
cfg.globalMaxITPC = max(cfg.globalMaxITPC, eps);

fig = figure( ...
    "Color", "w", "Visible", "off", "WindowStyle", "normal", ...
    "InvertHardcopy", "off", ...
    "Renderer", "painters");

nCols = min(2, numel(bandOrder));
nRows = ceil(numel(bandOrder)/nCols);
figW = max(cfg.minWidthCm, 2*cfg.sideMarginCm + nCols*cfg.panelWidthCm);
gridLeft = (figW-nCols*cfg.panelWidthCm)/2;
figH = cfg.headerHeightCm + nRows*cfg.rowHeightCm + cfg.bottomMarginCm;
set(fig, "Units", "centimeters", "Position", [2 2 figW figH], "Color", "w");
annotation(fig, "textbox", [.03, 1-.42/figH, .94, .30/figH], ...
    "String", string(figureTitle), "EdgeColor", "none", "Margin", 0, ...
    "HorizontalAlignment", "center", "FontName", cfg.fontName, ...
    "FontSize", cfg.mainTitleFontSize, "FontWeight", cfg.mainTitleFontWeight);

axH = gobjects(1, numel(bandOrder));
hLeg1 = gobjects(1);
hLeg2 = gobjects(1);

for j = 1:numel(bandOrder)

    b = bandOrder(j);

    row = ceil(j/nCols);
    col = mod(j-1,nCols)+1;
    x0 = gridLeft + (col-1)*cfg.panelWidthCm + (cfg.panelWidthCm-cfg.circleSizeCm)/2;
    y0 = cfg.bottomMarginCm + (nRows-row)*cfg.rowHeightCm + .12;
    ax = axes(fig, "Units", "centimeters", ...
        "Position", [x0,y0,cfg.circleSizeCm,cfg.circleSizeCm], ...
        "PositionConstraint", "innerposition");
    setappdata(ax,"fixedPositionCm",[x0,y0,cfg.circleSizeCm,cfg.circleSizeCm]);
    axH(j) = ax;

    hold(ax, "on");
    axis(ax, "equal");
    ax.Visible = "off";
    set(ax, "Color", [1 1 1]);

    axis(ax, [-1.12*rlimMax 1.12*rlimMax -1.12*rlimMax 1.12*rlimMax]);

    if j == 1
        hLeg1 = patch(ax, nan, nan, get_map_value(cfg.ribbonColor, c1), ...
            "FaceAlpha", cfg.barAlpha, ...
            "EdgeColor", get_map_value(cfg.roseEdgeColor, c1), ...
            "LineWidth", cfg.roseEdgeLW, ...
            "LineStyle", "-");

        hLeg2 = patch(ax, nan, nan, get_map_value(cfg.ribbonColor, c2), ...
            "FaceAlpha", cfg.barAlpha, ...
            "EdgeColor", get_map_value(cfg.roseEdgeColor, c2), ...
            "LineWidth", cfg.roseEdgeLW, ...
            "LineStyle", "-");
    end

    draw_polar_grid(ax, rlimMax, cfg);

    idx1 = (A.condition == c1) & (A.band == b);
    idx2 = (A.condition == c2) & (A.band == b);

    th1 = A.angle_rad(idx1);
    th2 = A.angle_rad(idx2);

    counts1 = histcounts(th1, edges);
    counts2 = histcounts(th2, edges);

    prob1 = counts1 / max(sum(counts1), 1);
    prob2 = counts2 / max(sum(counts2), 1);

    if strcmpi(string(cfg.histScaleMode), "global_common")
        h1 = prob1 / cfg.globalMaxProb * rHistMax;
        h2 = prob2 / cfg.globalMaxProb * rHistMax;
    else
        error("Unknown histScaleMode: %s", string(cfg.histScaleMode));
    end

    draw_rose_pair_short_on_top_offset(ax, h1, h2, edges, ...
        get_map_value(cfg.ribbonColor, c1), ...
        get_map_value(cfg.roseEdgeColor, c1), ...
        get_map_value(cfg.ribbonColor, c2), ...
        get_map_value(cfg.roseEdgeColor, c2), ...
        cfg.barAlpha, ...
        cfg.roseEdgeLW, ...
        cfg.roseEdgeAlpha, ...
        cfg.barAngleOffsetFrac, ...
        cfg.barWidthFrac);

    zUnit1 = subject_label_phase_vectors(A, idx1);
    zUnit2 = subject_label_phase_vectors(A, idx2);

    mu1 = group_mean_direction(zUnit1, wrap2pi);
    mu2 = group_mean_direction(zUnit2, wrap2pi);

    itpc1 = get_itpc_value(I, c1, b);
    itpc2 = get_itpc_value(I, c2, b);

    muSub1 = phase_vector_angles(zUnit1, wrap2pi);
    muSub2 = phase_vector_angles(zUnit2, wrap2pi);

    if isfinite(itpc1) && itpc1 > 0
        rArrow1 = cfg.arrowMaxFrac * rlimMax * itpc1 / cfg.globalMaxITPC;
        draw_clean_arrow(ax, mu1, rArrow1, get_map_value(cfg.lineColor, c1), cfg, rlimMax);
    end

    if isfinite(itpc2) && itpc2 > 0
        rArrow2 = cfg.arrowMaxFrac * rlimMax * itpc2 / cfg.globalMaxITPC;
        draw_clean_arrow(ax, mu2, rArrow2, get_map_value(cfg.lineColor, c2), cfg, rlimMax);
    end

    draw_outer_ring(ax, rlimMax, cfg);

    if cfg.showSubDots
        draw_subject_dots(ax, muSub1, rlimMax, get_map_value(cfg.pointColor, c1), cfg);
        draw_subject_dots(ax, muSub2, rlimMax, get_map_value(cfg.pointColor, c2), cfg);
    end

    draw_panel_title(ax, b, itpc1, itpc2, c1, c2, cfg);

    hold(ax, "off");
end

lgd = legend(axH(1), [hLeg1, hLeg2], ...
    {char(string(get_map_value(cfg.condLabel, c1))), ...
     char(string(get_map_value(cfg.condLabel, c2)))}, ...
    "Orientation", "horizontal", "Location", "none", ...
    "Box", "off", ...
    "FontName", cfg.fontName, ...
    "FontSize", cfg.legendFontSize);

set(lgd, "Units", "centimeters", "ItemTokenSize", [12 6], ...
    "Position", [cfg.sideMarginCm, figH-.98, figW-2*cfg.sideMarginCm, .35]);

for k = 1:numel(axH)
    set(axH(k),"Units","centimeters","Position",getappdata(axH(k),"fixedPositionCm"));
end
drawnow;

outPdf = fullfile(saveDir, "ieeg_itpc_four_bands.pdf");

set(fig,"PaperUnits","centimeters","PaperSize",[figW figH], ...
    "PaperPosition",[0 0 figW figH],"PaperPositionMode","manual");
print(fig,char(outPdf),'-dpdf','-painters');
close(fig);

end

function cfg = fill_default_cfg(cfg)

if ~isfield(cfg,'panelWidthCm'); cfg.panelWidthCm = 4.0; end
if ~isfield(cfg,'minWidthCm'); cfg.minWidthCm = 5.0; end
if ~isfield(cfg,'circleSizeCm'); cfg.circleSizeCm = 2.72*1.12; end
if ~isfield(cfg,'rowHeightCm'); cfg.rowHeightCm = 4.1; end
if ~isfield(cfg,'headerHeightCm'); cfg.headerHeightCm = 1.0; end
if ~isfield(cfg,'sideMarginCm'); cfg.sideMarginCm = .3; end
if ~isfield(cfg,'bottomMarginCm'); cfg.bottomMarginCm = .2; end

if ~isfield(cfg, 'histScaleMode'); cfg.histScaleMode = "global_common"; end
if ~isfield(cfg, 'histRadiusFrac'); cfg.histRadiusFrac = 0.90; end

if ~isfield(cfg, 'barAngleOffsetFrac'); cfg.barAngleOffsetFrac = 0.080; end
if ~isfield(cfg, 'barWidthFrac');       cfg.barWidthFrac       = 0.70; end

if ~isfield(cfg, 'barAlpha'); cfg.barAlpha = 0.68; end
if ~isfield(cfg, 'roseEdgeLW'); cfg.roseEdgeLW = 0.95; end
if ~isfield(cfg, 'roseEdgeAlpha'); cfg.roseEdgeAlpha = 1.00; end

if ~isfield(cfg, 'outerCircleColor'); cfg.outerCircleColor = [0.10 0.10 0.10]; end
if ~isfield(cfg, 'outerCircleLW'); cfg.outerCircleLW = 0.90; end
if ~isfield(cfg, 'outerCircleStyle'); cfg.outerCircleStyle = "-"; end

if ~isfield(cfg, 'crossColor'); cfg.crossColor = [0.86 0.86 0.86]; end
if ~isfield(cfg, 'crossLW'); cfg.crossLW = 0.55; end
if ~isfield(cfg, 'crossStyle'); cfg.crossStyle = "--"; end

if ~isfield(cfg, 'innerCircleColor'); cfg.innerCircleColor = [0.90 0.90 0.90]; end
if ~isfield(cfg, 'innerCircleLW'); cfg.innerCircleLW = 0.45; end
if ~isfield(cfg, 'innerCircleStyle'); cfg.innerCircleStyle = ":"; end
if ~isfield(cfg, 'innerCircleFracs'); cfg.innerCircleFracs = [0.45, 0.70]; end

if ~isfield(cfg, 'showAngleLabels'); cfg.showAngleLabels = false; end

if ~isfield(cfg, 'arrowScaleMode'); cfg.arrowScaleMode = "global_itpc"; end
if ~isfield(cfg, 'arrowMaxFrac'); cfg.arrowMaxFrac = 0.95; end
if ~isfield(cfg, 'arrowLW'); cfg.arrowLW = 1.4; end
if ~isfield(cfg, 'arrowHeadLengthFrac'); cfg.arrowHeadLengthFrac = .16; end
if ~isfield(cfg, 'arrowHeadWidthFrac'); cfg.arrowHeadWidthFrac = .11; end

if ~isfield(cfg, 'showSubDots'); cfg.showSubDots = true; end
if ~isfield(cfg, 'dotRadius'); cfg.dotRadius = 1.018; end
if ~isfield(cfg, 'subDotSize'); cfg.subDotSize = 10.5; end
if ~isfield(cfg, 'subDotAlpha'); cfg.subDotAlpha = 0.88; end
if ~isfield(cfg, 'subDotEdge'); cfg.subDotEdge = [1 1 1]; end
if ~isfield(cfg, 'subDotEdgeLW'); cfg.subDotEdgeLW = 0.18; end

if ~isfield(cfg, 'fontName'); cfg.fontName = "Arial"; end
if ~isfield(cfg, 'angleFontSize'); cfg.angleFontSize = 8.5; end
if ~isfield(cfg, 'angleFontWeight'); cfg.angleFontWeight = "normal"; end

if ~isfield(cfg, 'panelTitleFontSize'); cfg.panelTitleFontSize = 10.0; end
if ~isfield(cfg, 'panelTitleFontWeight'); cfg.panelTitleFontWeight = "bold"; end
if ~isfield(cfg, 'itpcTextFontSize'); cfg.itpcTextFontSize = 7.8; end

if ~isfield(cfg, 'mainTitleFontSize'); cfg.mainTitleFontSize = 12.2; end
if ~isfield(cfg, 'mainTitleFontWeight'); cfg.mainTitleFontWeight = "bold"; end
if ~isfield(cfg, 'legendFontSize'); cfg.legendFontSize = 8.6; end

if ~isfield(cfg, 'rlimMax'); cfg.rlimMax = 0.125; end

end

function draw_polar_grid(ax, rlimMax, cfg)

tt = linspace(0, 2*pi, 720);

for aa = [0, pi/2]
    plot(ax, [-rlimMax rlimMax] * cos(aa), ...
             [-rlimMax rlimMax] * sin(aa), ...
        "Color", cfg.crossColor, ...
        "LineWidth", cfg.crossLW, ...
        "LineStyle", cfg.crossStyle);
end

for k = 1:numel(cfg.innerCircleFracs)
    rr = cfg.innerCircleFracs(k) * rlimMax;

    plot(ax, rr*cos(tt), rr*sin(tt), ...
        "Color", cfg.innerCircleColor, ...
        "LineWidth", cfg.innerCircleLW, ...
        "LineStyle", cfg.innerCircleStyle);
end

end

function draw_outer_ring(ax, rlimMax, cfg)

tt = linspace(0, 2*pi, 720);

plot(ax, rlimMax*cos(tt), rlimMax*sin(tt), ...
    "Color", cfg.outerCircleColor, ...
    "LineWidth", cfg.outerCircleLW, ...
    "LineStyle", cfg.outerCircleStyle);

end

function draw_rose_pair_short_on_top_offset(ax, h1, h2, edges, ...
    face1, edge1, face2, edge2, barAlpha, edgeLW, edgeAlpha, offsetFrac, widthFrac)

nBins = numel(h1);

for k = 1:nBins

    r1 = h1(k);
    r2 = h2(k);

    if r1 <= 0 && r2 <= 0
        continue;
    end

    off1 = +offsetFrac;
    off2 = -offsetFrac;

    if r1 >= r2
        draw_one_rose_bar_offset(ax, r1, edges, k, face1, edge1, ...
            barAlpha, edgeLW, edgeAlpha, off1, widthFrac);

        draw_one_rose_bar_offset(ax, r2, edges, k, face2, edge2, ...
            barAlpha, edgeLW, edgeAlpha, off2, widthFrac);
    else
        draw_one_rose_bar_offset(ax, r2, edges, k, face2, edge2, ...
            barAlpha, edgeLW, edgeAlpha, off2, widthFrac);

        draw_one_rose_bar_offset(ax, r1, edges, k, face1, edge1, ...
            barAlpha, edgeLW, edgeAlpha, off1, widthFrac);
    end
end

end

function draw_one_rose_bar_offset(ax, r, edges, k, faceColor, edgeColor, ...
    barAlpha, edgeLW, edgeAlpha, offsetFrac, widthFrac)

if r <= 0
    return;
end

t1 = edges(k);
t2 = edges(k + 1);
binWidth = t2 - t1;

tCenter = (t1 + t2) / 2 + offsetFrac * binWidth;

halfWidth = 0.5 * widthFrac * binWidth;
halfWidth = min(halfWidth, 0.48 * binWidth);

ang = linspace(tCenter - halfWidth, tCenter + halfWidth, 60);

x = [0, r*cos(ang), 0];
y = [0, r*sin(ang), 0];

p = patch(ax, x, y, faceColor, ...
    "FaceAlpha", barAlpha, ...
    "EdgeColor", edgeColor, ...
    "LineWidth", edgeLW, ...
    "LineStyle", "-", ...
    "LineJoin", "round");

try
    p.EdgeAlpha = edgeAlpha;
catch
end

end

function draw_subject_dots(ax, muSub, rlimMax, faceColor, cfg)

if isempty(muSub)
    return;
end

rDot = cfg.dotRadius * rlimMax;
x = rDot * cos(muSub);
y = rDot * sin(muSub);

h = scatter(ax, x, y, cfg.subDotSize, ...
    "Marker", "o", ...
    "MarkerFaceColor", faceColor, ...
    "MarkerEdgeColor", cfg.subDotEdge, ...
    "LineWidth", max(cfg.subDotEdgeLW, 0.01));

try
    h.MarkerFaceAlpha = cfg.subDotAlpha;
    h.MarkerEdgeAlpha = min(0.95, cfg.subDotAlpha + 0.05);
catch
end

end

function draw_clean_arrow(ax, theta, rArrow, color, cfg, rlimMax)

if ~(isfinite(rArrow) && rArrow > 0)
    return;
end

headLen = cfg.arrowHeadLengthFrac * rlimMax;
headWid = cfg.arrowHeadWidthFrac  * rlimMax;

headLen = min(headLen, 0.45 * rArrow);

tip = [rArrow*cos(theta), rArrow*sin(theta)];
baseCenter = tip - headLen * [cos(theta), sin(theta)];

perp = [-sin(theta), cos(theta)];

p1 = baseCenter + 0.5 * headWid * perp;
p2 = baseCenter - 0.5 * headWid * perp;

shaftEnd = tip - 0.72 * headLen * [cos(theta), sin(theta)];

plot(ax, [0, shaftEnd(1)], [0, shaftEnd(2)], ...
    "Color", color, ...
    "LineWidth", cfg.arrowLW);

patch(ax, [tip(1), p1(1), p2(1)], ...
          [tip(2), p1(2), p2(2)], ...
      color, ...
      "EdgeColor", "none", ...
      "FaceAlpha", 1.0);

end

function zUnit = subject_label_phase_vectors(A, idx)

zUnit = [];

if ~any(idx)
    return;
end

if ismember("label", string(A.Properties.VariableNames))
    groupKeys = unique(A(idx, ["subject", "label"]), "rows");
else
    groupKeys = unique(A(idx, "subject"), "rows");
end

zUnit = complex(nan(height(groupKeys), 1), nan(height(groupKeys), 1));

for s = 1:height(groupKeys)

    sid = groupKeys.subject(s);

    if isnumeric(A.subject)
        subIdx = A.subject == sid;
    else
        subIdx = string(A.subject) == string(sid);
    end

    if ismember("label", string(A.Properties.VariableNames))
        lab = groupKeys.label(s);

        if isnumeric(A.label)
            labelIdx = A.label == lab;
        else
            labelIdx = string(A.label) == string(lab);
        end
    else
        labelIdx = true(height(A), 1);
    end

    th_s = A.angle_rad(idx & subIdx & labelIdx);

    if isempty(th_s)
        continue;
    end

    zUnit(s) = mean(exp(1i * th_s));
end

keep = isfinite(real(zUnit)) & isfinite(imag(zUnit));
zUnit = zUnit(keep);

end

function muSub = phase_vector_angles(zUnit, wrap2pi)

muSub = [];

if isempty(zUnit)
    return;
end

keep = isfinite(real(zUnit)) & isfinite(imag(zUnit)) & abs(zUnit) > eps;
muSub = wrap2pi(angle(zUnit(keep)));

end

function mu = group_mean_direction(zUnit, wrap2pi)

mu = 0;

if isempty(zUnit)
    return;
end

zGroup = mean(zUnit);

if isfinite(abs(zGroup)) && abs(zGroup) > eps
    mu = wrap2pi(angle(zGroup));
end

end

function itpc = get_itpc_value(I, cond, band)

row = I((I.condition == cond) & (I.band == band), :);

if ~isempty(row) && ismember("itpc_mean", string(row.Properties.VariableNames))
    vals = row.itpc_mean;
    vals = vals(isfinite(vals));

    if isempty(vals)
        itpc = NaN;
    else
        itpc = mean(vals);
    end
else
    itpc = NaN;
end

end

function draw_panel_title(ax, bandName, itpc1, itpc2, c1, c2, cfg)

text(ax, .5, 1.16, string(bandName), "Units", "normalized", ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom", ...
    "FontName", cfg.fontName, "FontSize", cfg.panelTitleFontSize, ...
    "FontWeight", cfg.panelTitleFontWeight, "Clipping", "off", "Color", [0 0 0]);
text(ax, .5, 1.05, sprintf("Vol. %.3f | Invol. %.3f", itpc1, itpc2), ...
    "Units", "normalized", "HorizontalAlignment", "center", ...
    "VerticalAlignment", "bottom", "FontName", "Times New Roman", ...
    "FontSize", cfg.itpcTextFontSize, "Interpreter", "none", ...
    "Clipping", "off", "Color", [0 0 0]);

end

function val = get_map_value(mp, key)

keyStr  = string(key);
keyChar = char(keyStr);

try
    val = mp(keyChar);
    return;
catch
end

try
    val = mp(keyStr);
    return;
catch
end

try
    val = mp(key);
    return;
catch
end

error("Key not found in containers.Map: %s", keyChar);

end

function out = short_condition_label(x)

x = string(x);
xl = lower(x);

if strcmpi(x, "Voluntary")
    out = "V";
elseif strcmpi(x, "Involuntary")
    out = "I";
elseif contains(xl, "involuntary")
    out = "I";
elseif contains(xl, "voluntary")
    out = "V";
else
    out = x;
end

end
