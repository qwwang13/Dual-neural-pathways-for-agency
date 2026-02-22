clear; clc;

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));

dirData = fullfile(repo,'data', '05_LMM', 'plot_data_for_itpc');

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

% ---------------- action binding ----------------
cfg = cfgCommon;
cfg.condOrder = ["va_a","ia_a"];
cfg.condLabel = containers.Map(["va_a","ia_a"], ["Voluntary","Involuntary"]);
cfg.ribbonColor = containers.Map(["va_a","ia_a"], {ribbon_va, ribbon_ia});
cfg.lineColor   = containers.Map(["va_a","ia_a"], {line_va,   line_ia});
cfg.pointColor  = containers.Map(["va_a","ia_a"], {point_va,  point_ia});
cfg.arrowGain = 2.2;

saveDir = fullfile(repo,'figures', '05_LMM', 'ITPC', "action_binding");
plot_itpc_rose_cartesian("PreSMAR_action_ba", saveDir, dirData, "PreSMAR", cfg);
plot_itpc_rose_cartesian("PreSMAL_action_ba", saveDir, dirData, "PreSMAL", cfg);
plot_itpc_rose_cartesian("SMAR_action_ba",    saveDir, dirData, "SMAR",   cfg);

% ---------------- outcome binding ----------------
cfg = cfgCommon;
cfg.condOrder = ["vs_s","is_s"];
cfg.condLabel = containers.Map(["vs_s","is_s"], ["Voluntary","Involuntary"]);
cfg.ribbonColor = containers.Map(["vs_s","is_s"], {ribbon_va, ribbon_ia});
cfg.lineColor   = containers.Map(["vs_s","is_s"], {line_va,   line_ia});
cfg.pointColor  = containers.Map(["vs_s","is_s"], {point_va,  point_ia});
cfg.arrowGain = 2.2;

saveDir = fullfile(repo,'figures', '05_LMM', 'ITPC', "outcome_binding");
plot_itpc_rose_cartesian("PostcentralL_outcome_bs", saveDir, dirData, "PostcentralL", cfg);
plot_itpc_rose_cartesian("MTGR_outcome_bs",         saveDir, dirData, "MTGR",        cfg);

% ---------------- 75/50 without outcome ----------------
cfg = cfgCommon;
cfg.condOrder = ["a75_an","a50_an"];
cfg.condLabel = containers.Map(["a75_an","a50_an"], ["75% without outcome","50% without outcome"]);
cfg.ribbonColor = containers.Map(["a75_an","a50_an"], {ribbon_va, ribbon_ia});
cfg.lineColor   = containers.Map(["a75_an","a50_an"], {line_va,   line_ia});
cfg.pointColor  = containers.Map(["a75_an","a50_an"], {point_va,  point_ia});
cfg.arrowGain = 2.2;

saveDir = fullfile(repo,'figures', '05_LMM', 'ITPC',"75_50_without");
plot_itpc_rose_cartesian("AGR_75_50_ba",        saveDir, dirData, "AGR",        cfg);
plot_itpc_rose_cartesian("DLPFCR_75_50_ba",     saveDir, dirData, "DLPFCR",     cfg);
plot_itpc_rose_cartesian("PrecuneusR_75_50_ba", saveDir, dirData, "PrecuneusR", cfg);

% ---------------- 50 with/without outcome ----------------
cfg = cfgCommon;
cfg.condOrder = ["a50_ay","a50_an"];
cfg.condLabel = containers.Map(["a50_ay","a50_an"], ["50% with outcome","50% without outcome"]);
cfg.ribbonColor = containers.Map(["a50_ay","a50_an"], {ribbon_va, ribbon_ia});
cfg.lineColor   = containers.Map(["a50_ay","a50_an"], {line_va,   line_ia});
cfg.pointColor  = containers.Map(["a50_ay","a50_an"], {point_va,  point_ia});
cfg.arrowGain = 1;

saveDir = fullfile(repo,'figures', '05_LMM', 'ITPC',"50_with_without");
plot_itpc_rose_cartesian("STGR_50_ba", saveDir, dirData, "STGR", cfg);