clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
%filter_name = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl';
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_timebin'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);

%%
subplot(1,4,1);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time bin index'; 
session_list_plot = bpGlobal.rolo.session_list.switching;
fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions)

%%
subplot(1,4,2);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'oblique';
plotOptions.xlabelStr = 'Time bin index'; 
session_list_plot = bpGlobal.rolo.session_list.switching;
fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions)
%%
subplot(1,4,3);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time bin index'; 
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions)

%%
subplot(1,4,4);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'oblique';
plotOptions.xlabelStr = 'Time bin index'; 
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions)