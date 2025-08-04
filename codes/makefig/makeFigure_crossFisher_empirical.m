clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
load(fullfile(saveFolder, 'results_SubsampleOrganized_perCohr_fisherInfo_all_sessions'));
results_all_perCohr = results_cross_sizeControl_perCohr;
results_all_perCohr = get_sample_CI_cross(results_all_perCohr);
%% 1. bar plot of I_real and I_cross

save_folder = '../../figures/figures_final/fisher_info';
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;

figure
ftsize = 16;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
session_list_plot = bpGlobal.rolo.session_list.switching;
fig.plot_bar_Info(results_all, session_list_plot, plotOptions); title('Monkey R')
subplot(1,2,2)
%fig.plot_bar_Info(results_all, session_list_gremlin_good,plotOptions); title('Monkey G Good sessions');
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_bar_Info(results_all, session_list_plot,plotOptions); title('Monkey G');
save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar.svg'));
print(save_name, '-dsvg');


figure
ftsize = 16;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
session_list_plot = bpGlobal.rolo.session_list.switching;
fig.plot_bar_Info(results_all_perCohr, session_list_plot, plotOptions); title('Monkey R')
subplot(1,2,2)
%fig.plot_bar_Info(results_all, session_list_gremlin_good,plotOptions); title('Monkey G Good sessions');
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_bar_Info(results_all_perCohr, session_list_plot,plotOptions); title('Monkey G');
save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar_percohr.svg'));
print(save_name, '-dsvg');
%% 2. bar plot of I_redundacy, I_redundacy_cross

plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;

figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
session_list_plot = bpGlobal.rolo.session_list.switching;
fig.plot_bar_deltaInfo(results_all, session_list_plot, plotOptions); title('Monkey R');
set(gca,'fontsize',20)

subplot(1,2,2)
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_bar_deltaInfo(results_all, session_list_plot, plotOptions); title('Monkey G');
set(gca,'fontsize',20)

figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
session_list_plot = bpGlobal.rolo.session_list.switching;
fig.plot_bar_deltaInfo(results_all_perCohr, session_list_plot, plotOptions); title('Monkey R');
set(gca,'fontsize',20)

subplot(1,2,2)
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_bar_deltaInfo(results_all_perCohr, session_list_plot, plotOptions); title('Monkey G');
set(gca,'fontsize',20)
%% 3. scatter plot of I_cross and I_real
