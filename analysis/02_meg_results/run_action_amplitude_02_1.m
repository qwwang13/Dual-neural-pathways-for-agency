% run_action_amplitude.m
clear; clc;

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));

dirData = fullfile(repo,'data', '02_meg_results');
dirOut  = fullfile(repo,'figures', '02_meg_results');

roi6  = {'AGR','PreSMAL','PreSMAR','SMAR','DLPFCR','InsulaR'};
roi12 = {'AGR','PreSMAL','PreSMAR','SMAR','DLPFCR','InsulaR','PrecuneusR','MTGR','PostcentralL','SPGL','PrecentralR','PostcentralR'};

common = struct();
common.xTicks   = -200:200:1000;
common.viewAzEl = [50 25];
common.baseZ    = 0;
common.reverseZ = true;

% ---- 6 ROI ----
opts = common; opts.zLim = [-3 5]; opts.zTicks = -3:1:5;
plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_action_va_6.csv'), roi6, ...
    fullfile(dirOut,'action_binding_va_6.pdf'), ...
    zLim=[-3 5],titleStr="Voluntary" );

plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_action_ia_6.csv'), roi6, ...
    fullfile(dirOut,'action_binding_ia_6.pdf'), ...
    zLim=[-3 5],titleStr="Involuntary" );

% ---- 12 ROI ----
% opts = common; opts.zLim = [-3 5]; opts.zTicks = -3:1:5;
plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_action_va_12.csv'), roi12, ...
    fullfile(dirOut,'action_binding_va_12.pdf'), ...
    zLim=[-3 5],titleStr="Voluntary"); 

plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_action_ia_12.csv'), roi12, ...
    fullfile(dirOut,'action_binding_ia_12.pdf'), ...
    zLim=[-3 5],titleStr="Involuntary"  ); 
