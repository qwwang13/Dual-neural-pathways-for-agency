function outPdf = plot_itpc_rose_cartesian(regionStem, saveDir, baseDir, figureTitle, cfg)

anglesFile = fullfile(baseDir, string(regionStem) + "_angles.csv");
itpcFile   = fullfile(baseDir, string(regionStem) + "_itpc_mean.csv");

A = readtable(anglesFile);
I = readtable(itpcFile);

A.condition = string(A.condition);
A.band      = string(A.band);
I.condition = string(I.condition);
I.band      = string(I.band);

condOrder = string(cfg.condOrder);
bandOrder = string(cfg.bandOrder);

condLabel = cfg.condLabel;
ribbonColor = cfg.ribbonColor;
lineColor   = cfg.lineColor;
pointColor  = cfg.pointColor;

nBins  = cfg.nBins;
edges  = linspace(0, 2*pi, nBins+1);

barAlpha = cfg.barAlpha;

outerCircleColor = cfg.outerCircleColor;
outerCircleLW    = cfg.outerCircleLW;
outerCircleStyle = cfg.outerCircleStyle;

crossColor       = cfg.crossColor;
crossLW          = cfg.crossLW;
crossStyle       = cfg.crossStyle;

innerCircleColor = cfg.innerCircleColor;
innerCircleLW    = cfg.innerCircleLW;
innerCircleStyle = cfg.innerCircleStyle;
innerCircleFrac  = cfg.innerCircleFrac;

arrowGain = cfg.arrowGain;
arrowLW   = cfg.arrowLW;
headSize  = cfg.headSize;

subDotSize = cfg.subDotSize;
subDotEdge = cfg.subDotEdge;

angleFontSize   = cfg.angleFontSize;
angleFontWeight = cfg.angleFontWeight;

w_cm = cfg.w_cm;
h_cm = cfg.h_cm;

rlimMax = cfg.rlimMax;
rInner  = innerCircleFrac * rlimMax;

fig = figure("Color","w");
set(fig, "Units","centimeters");
set(fig, "Position", [2 2 w_cm h_cm]);
set(fig, "Color","w");

t = tiledlayout(2,4,"TileSpacing","compact","Padding","compact");

wrap2pi = @(x) mod(x, 2*pi);
axH = gobjects(numel(condOrder), numel(bandOrder));

for i = 1:numel(condOrder)
    c = condOrder(i);

    for j = 1:numel(bandOrder)
        b = bandOrder(j);

        ax = nexttile(t);
        axH(i,j) = ax;

        hold(ax,"on");
        axis(ax,"equal");

        pad = 1;
        axis(ax,[-pad*rlimMax pad*rlimMax+cfg.xpad_right -rlimMax rlimMax]);
        ax.Visible = "off";

        for aa = [0, pi/2]
            plot(ax, [-rlimMax rlimMax]*cos(aa), [-rlimMax rlimMax]*sin(aa), ...
                "Color", crossColor, "LineWidth", crossLW, "LineStyle", crossStyle);
        end

        tt = linspace(0,2*pi,500);
        plot(ax, rlimMax*cos(tt), rlimMax*sin(tt), ...
            "Color", outerCircleColor, "LineWidth", outerCircleLW, "LineStyle", outerCircleStyle);

        plot(ax, rInner*cos(tt), rInner*sin(tt), ...
            "Color", innerCircleColor, "LineWidth", innerCircleLW, "LineStyle", innerCircleStyle);

        labR = 1.08 * rlimMax;
        txtCol = [0.1 0.1 0.1];
        text(ax, labR, 0,    "0°",   "HorizontalAlignment","left",  "VerticalAlignment","middle", ...
            "FontSize",angleFontSize,"FontWeight",angleFontWeight,"Color",txtCol);
        text(ax, 0,    labR, "90°",  "HorizontalAlignment","center","VerticalAlignment","bottom", ...
            "FontSize",angleFontSize,"FontWeight",angleFontWeight,"Color",txtCol);
        text(ax,-labR, 0,    "180°", "HorizontalAlignment","right", "VerticalAlignment","middle", ...
            "FontSize",angleFontSize,"FontWeight",angleFontWeight,"Color",txtCol);
        text(ax, 0,   -labR, "270°", "HorizontalAlignment","center", "VerticalAlignment","top", ...
            "FontSize",angleFontSize,"FontWeight",angleFontWeight,"Color",txtCol);

        idx = (A.condition == c) & (A.band == b);
        thAll = A.angle_rad(idx);

        if ~isempty(thAll)
            counts  = histcounts(thAll, edges);
            heights = counts / max(sum(counts),1);

            for k = 1:nBins
                r = heights(k);
                if r <= 0; continue; end
                t1 = edges(k); t2 = edges(k+1);
                ang = linspace(t1, t2, 40);

                x = [0, r*cos(ang), 0];
                y = [0, r*sin(ang), 0];

                patch(ax, x, y, ribbonColor(c), ...
                    "FaceAlpha", barAlpha, ...
                    "EdgeColor", "none");
            end
        end

        row = I((I.condition == c) & (I.band == b), :);
        if ~isempty(row); itpc = row.itpc_mean(1); else; itpc = NaN; end

        mu_sub = [];
        if ~isempty(thAll)
            subIDs = unique(A.subject(idx));
            mu_sub = nan(numel(subIDs),1);
            for s = 1:numel(subIDs)
                sid = subIDs(s);
                th_s = A.angle_rad(idx & (A.subject == sid));
                if isempty(th_s); continue; end
                R_s = mean(exp(1i*th_s));
                mu_sub(s) = wrap2pi(angle(R_s));
            end
            mu_sub = mu_sub(isfinite(mu_sub));
        end

        if ~isempty(mu_sub)
            xdot = rlimMax * cos(mu_sub);
            ydot = rlimMax * sin(mu_sub);
            scatter(ax, xdot, ydot, subDotSize, ...
                "MarkerFaceColor", pointColor(c), ...
                "MarkerEdgeColor", subDotEdge);
        end

        mu = 0;
        if ~isempty(mu_sub)
            Rg = mean(exp(1i*mu_sub));
            mu = wrap2pi(angle(Rg));
        elseif ~isempty(thAll)
            Rp = mean(exp(1i*thAll));
            mu = wrap2pi(angle(Rp));
        end

        if isfinite(itpc) && itpc > 0
            rArrow = arrowGain * itpc * rlimMax;
            xEnd = rArrow * cos(mu);
            yEnd = rArrow * sin(mu);

            quiver(ax, 0, 0, xEnd, yEnd, 0, ...
                "Color", lineColor(c), ...
                "LineWidth", arrowLW, ...
                "MaxHeadSize", headSize);
        end

        if isfinite(itpc)
            ttl = sprintf("%s-band ITPC=%.3f", b, itpc);
        else
            ttl = sprintf("%s-band ITPC=NA", b);
        end
        text(ax, 0, cfg.title_y, ttl, ...
            "HorizontalAlignment","center", "VerticalAlignment","bottom", ...
            "FontSize", cfg.titleFontSize, "FontWeight","normal", "Color",[0.1 0.1 0.1]);

        hold(ax,"off");
    end
end

drawnow;

for i = 1:numel(condOrder)
    leftAx  = axH(i,1);
    rightAx = axH(i,end);

    pL = leftAx.Position;
    pR = rightAx.Position;

    xLeft  = pL(1);
    xRight = pR(1) + pR(3);
    yTop   = pL(2) + pL(4);

    y = yTop + cfg.rowLabel_y_pad;
    annotation("textbox", [xLeft, y, (xRight-xLeft), 0.04], ...
        "String", condLabel(condOrder(i)), ...
        "EdgeColor","none", ...
        "FontSize", cfg.rowLabelFontSize, ...
        "FontWeight","bold", ...
        "HorizontalAlignment","center", ...
        "VerticalAlignment","middle", ...
        "Color",[0.1 0.1 0.1]);
end

outPdf = fullfile(saveDir, string(figureTitle) + "_ITPC_polar.pdf");
exportgraphics(gcf, outPdf, "ContentType", "vector");

end