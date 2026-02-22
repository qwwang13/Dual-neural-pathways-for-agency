% run_probability_amplitude.m
clear; clc;

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));

dirData = fullfile(repo,'data', '02_meg_results');
dirOut  = fullfile(repo,'figures', '02_meg_results');

roi12 = {'AGR','PreSMAL','PreSMAR','SMAR','DLPFCR','InsulaR','PrecuneusR','MTGR','PostcentralL','SPGL','PrecentralR','PostcentralR'};

common = struct();
common.xTicks   = -200:200:1500;
common.viewAzEl = [50 25];
common.baseZ    = 0;
common.reverseZ = true;


plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_a50n_12.csv'), roi12, ...
    fullfile(dirOut,'prob_50_without_12.pdf'), ...
    zLim=[-3 3],titleStr=sprintf('50%% probability\n  (without outcome)')); 

plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_a50y_12.csv'), roi12, ...
    fullfile(dirOut,'prob_50_with_12.pdf'), ...
    zLim=[-3 3],titleStr=sprintf('50%% probability\n  (with outcome)')) 

plot3d_tc_from_csv( ...
    fullfile(dirData,'ta_a75n_12.csv'), roi12, ...
    fullfile(dirOut,'prob_75_without_12.pdf'), ...
    zLim=[-3 3],titleStr=sprintf('75%% probability\n  (without outcome)')) ;

