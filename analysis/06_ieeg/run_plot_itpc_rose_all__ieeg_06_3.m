clear; clc;

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));

dirData = fullfile(repo,'data', '06_ieeg', 'plot_data_for_itpc');

cfgCommon = struct();
cfgCommon.bandOrder = ["Delta","Theta","Alpha","Beta"];

cfgCommon.nBins = 20;
cfgCommon.barAlpha = 1;

cfgCommon.outerCircleColor = [0.15 0.15 0.15];
cfgCommon.outerCircleLW    = 0.9;
cfgCommon.outerCircleStyle = ":";

cfgCommon.crossColor       = [0.45 0.45 0.45];
cfgCommon.crossLW          = 1;
cfgCommon.crossStyle       = "--";

cfgCommon.innerCircleColor = [0.65 0.65 0.65];
cfgCommon.innerCircleLW    = 0.9;
cfgCommon.innerCircleStyle = ":";
cfgCommon.innerCircleFrac  = 0.50;

cfgCommon.arrowLW   = 2;
cfgCommon.headSize  = 0.8;

cfgCommon.subDotSize = 12;
cfgCommon.subDotEdge = "none";

cfgCommon.angleFontSize = 10;
cfgCommon.angleFontWeight = "bold";

cfgCommon.w_cm = 18;
cfgCommon.h_cm = 12;

cfgCommon.rlimMax = 0.125;

cfgCommon.xpad_right = 0.05;

cfgCommon.titleFontSize = 10;
cfgCommon.title_y = 0.18;

cfgCommon.rowLabel_y_pad = 0.1;
cfgCommon.rowLabelFontSize = 10;

ribbon_va = [255 210 157] / 255;
ribbon_ia = [189 215 231] / 255;
line_va   = [251 133   0] / 255;
line_ia   = [  8  81 156] / 255;

point_va = [188  57   8] / 255;
point_ia = [ 31 119 180] / 255;

point_strength = 0.75;
mixWhite = @(rgb, s) s*rgb + (1-s)*[1 1 1];
point_va = mixWhite(point_va, point_strength);
point_ia = mixWhite(point_ia, point_strength);


% ---------------- outcome binding ----------------
cfg = cfgCommon;
cfg.condOrder = ["Voluntary_OutcomeOnsetReportSound","Involuntary_OutcomeOnsetReportSound"];
cfg.condLabel = containers.Map(["Voluntary_OutcomeOnsetReportSound","Involuntary_OutcomeOnsetReportSound"], ["Voluntary","Involuntary"]);
cfg.ribbonColor = containers.Map(["Voluntary_OutcomeOnsetReportSound","Involuntary_OutcomeOnsetReportSound"], {ribbon_va, ribbon_ia});
cfg.lineColor   = containers.Map(["Voluntary_OutcomeOnsetReportSound","Involuntary_OutcomeOnsetReportSound"], {line_va,   line_ia});
cfg.pointColor  = containers.Map(["Voluntary_OutcomeOnsetReportSound","Involuntary_OutcomeOnsetReportSound"], {point_va,  point_ia});
cfg.arrowGain = 2.2;

saveDir = fullfile(repo,'figures', '06_ieeg', 'ITPC');
plot_itpc_rose_cartesian_ieeg("ieeg_PostcentralL", saveDir, dirData, "PostcentralL", cfg);
plot_itpc_rose_cartesian_ieeg("ieeg_MTGR",         saveDir, dirData, "MTGR",        cfg);
