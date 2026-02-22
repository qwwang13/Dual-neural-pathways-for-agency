function fig = plot3d_tc_from_csv(csvPath, roiOrder, outPdf, opts)
% plot3d_tc_from_csv
% Read long-format CSV (Time, region, Value), reorder ROIs, plot 3D time-course,
% export vector PDF.
%
% Required:
%   csvPath  - path to CSV file
%   roiOrder - cell array of ROI names in desired order
%   outPdf   - output PDF path
% Optional opts fields (all optional):
%   .titleStr, .xLabel, .yLabel, .zLabel
%   .xTicks, .zTicks, .xLim, .zLim
%   .viewAzEl = [az el]
%   .baseZ = 0
%   .reverseZ = true
%   .lineWidth = 1.8
%   .faceAlpha = 0.32
%   .colors (Ny x 3) or empty -> auto
%
% Notes:
% - Missing values are treated as NaN (NOT zero), so true zeros are preserved.
% - Only exports PDF (vector).

arguments
    csvPath (1,:) char
    roiOrder cell
    outPdf (1,:) char
    opts.titleStr (1,:) char = ''
    opts.xLabel (1,:) char = 'Time (ms)'
    opts.yLabel (1,:) char = 'ROI'
    opts.zLabel (1,:) char = 'Amplitude (z-score)'
    opts.xTicks double = []
    opts.zTicks double = []
    opts.xLim double = []
    opts.zLim double = []
    opts.viewAzEl double = [50 25]
    opts.baseZ double = 0
    opts.reverseZ logical = true
    opts.lineWidth double = 1.8
    opts.faceAlpha double = 0.32
    opts.colors double = []
end

% ---------- read ----------
T = readtable(csvPath);
time = T.Time;
roi  = string(T.region);
val  = T.Value;

% ---------- pick & order ROIs that exist ----------
roiExist = unique(roi);
roiSorted = strings(0,1);
for i = 1:numel(roiOrder)
    if any(roiExist == string(roiOrder{i}))
        roiSorted(end+1,1) = string(roiOrder{i}); %#ok<AGROW>
    end
end
Ny = numel(roiSorted);
if Ny == 0
    error("No ROI in roiOrder found in file: %s", csvPath);
end

% ---------- build matrix z (Ny x Nx) ----------
timeUnique = unique(time);
timeUnique = sort(timeUnique);
Nx = numel(timeUnique);

[tfR, rIdx] = ismember(roi, roiSorted);
[tfT, tIdx] = ismember(time, timeUnique);

keep = tfR & tfT & ~isnan(val);
rIdx = rIdx(keep);
tIdx = tIdx(keep);
v    = val(keep);

% mean-aggregate if duplicates exist
z = accumarray([rIdx, tIdx], v, [Ny, Nx], @mean, NaN);

% ---------- colors ----------
if isempty(opts.colors)
    opts.colors = default_palette(max(Ny,12));
end
colors = opts.colors(1:Ny, :);

% ---------- plot ----------
fig = figure('Position',[100 100 950 850], 'Color','white');
hold on;

baseZ = opts.baseZ;

for i = 1:Ny
    zi = z(i, :);
    valid = ~isnan(zi);

    xi = timeUnique(valid);
    zi = zi(valid);

    if numel(xi) > 1
        plot3(xi, i*ones(size(xi)), zi, 'Color', colors(i,:), 'LineWidth', opts.lineWidth);

        x_fill = [xi(:); flipud(xi(:))]';
        y_fill = [i*ones(size(xi(:))); i*ones(size(xi(:)))]';
        z_fill = [zi(:); baseZ*ones(size(zi(:)))]';
        fill3(x_fill, y_fill, z_fill, colors(i,:), 'FaceAlpha', opts.faceAlpha, 'EdgeColor','none');
    end
end

hold off;

xlabel(opts.xLabel, 'FontSize', 13, 'FontWeight','bold');
ylabel(opts.yLabel, 'FontSize', 13, 'FontWeight','bold');
zlabel(opts.zLabel, 'FontSize', 13, 'FontWeight','bold');

yticks(1:Ny);
yticklabels(cellstr(roiSorted));

set(gca, 'FontSize', 12, 'FontName','Arial');
grid off;

view(opts.viewAzEl(1), opts.viewAzEl(2));

% ticks/limits
if ~isempty(opts.xTicks), xticks(opts.xTicks); end
if ~isempty(opts.zTicks), zticks(opts.zTicks); end
if ~isempty(opts.xLim), xlim(opts.xLim); else, xlim([min(timeUnique) max(timeUnique)]); end
ylim([0.5, Ny+0.5]);
if ~isempty(opts.zLim), zlim(opts.zLim); end

if opts.reverseZ
    set(gca, 'ZDir','reverse');
end

if ~isempty(opts.titleStr)
    ht = title(opts.titleStr, 'FontSize', 13, 'FontWeight','bold');
    set(ht, 'VerticalAlignment','bottom');
    set(ht, 'Position', [0.5, 2, (opts.zLim(1) + 0.2*(opts.zLim(2)-opts.zLim(1)))]);
end

% ---------- export PDF ----------
outDir = fileparts(outPdf);
if ~exist(outDir, 'dir'), mkdir(outDir); end
exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','white');

end

% ================= helper: palette =================
function C = default_palette(n)
base = [
    0.806, 0.106, 0.489;
    0.151, 0.071, 0.558;
    0.141, 0.573, 0.898;
    0.175, 0.876, 0.918;
    0.247, 0.796, 0.416;
    1.00,  0.70,  0.20;
    1.00,  0.50,  0.30;
    1.00,  0.40,  0.60;
    0.90,  0.30,  0.80;
    0.70,  0.30,  0.90;
    0.50,  0.30,  1.00;
    0.890, 0.161, 0.222;
];
if n <= size(base,1)
    C = base(1:n,:);
else
    extra = lines(n - size(base,1));
    C = [base; extra];
end
end
