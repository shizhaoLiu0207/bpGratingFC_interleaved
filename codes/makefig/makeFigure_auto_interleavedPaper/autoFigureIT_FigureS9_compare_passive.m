clear all
clc
close all
%%
global   bpGlobal ftsize
ftsize                                     = 14;
bpGratingFCGlobal();
versionName = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl';
[session_list, results_neural] = data_prep(versionName);
%%
figure
set(gcf,'unit','inches','position',[0,0,12,4])
save_folder = '../../../figures/figures_auto_interleavedpaper';

save_name = fullfile(save_folder,'Figure_S9_passive_task_Iredundancy_sizeControl.svg');
tex_bar_name = fullfile(save_folder,'Figure_S9_passive_task_Iredundancy_sizeControl.tex');
%% monkey R cardinal
ax_1 = subplot(1,4,1);
set(ax_1,'position',get(ax_1,'position')+[-0.05 0.02 0.01 -0.15]);
plotOptions.nTimebin = 4;
plotOptions.dataPlot = 'combine_delta';
plotOptions.taskType = 'cardinal';
plotOptions.plotType = 'bar';
sessionlist_plot_late = session_list.rolo_switching;
sessionlist_plot_early = [];
stats_info = fig.plot_task_passive_compare(results_neural, sessionlist_plot_late,sessionlist_plot_early, plotOptions);
ylabel('$I_\textrm{redundancy}$','FontSize',ftsize,'Interpreter','latex')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');

monkeyR_cardinal_bar_text = sprintf('Monkey R cardinal: task performing: paired t-test, $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n',...
    stats_info.df_late, stats_info.tval_late, stats_info.p_late);


%% monkey R oblique
ax_2 = subplot(1,4,2);
set(ax_2,'position',get(ax_2,'position')+[-0.03 0.02 0.01 -0.15]);
plotOptions.nTimebin = 4;
plotOptions.dataPlot = 'combine_delta';
plotOptions.taskType = 'oblique';
plotOptions.plotType = 'bar';
sessionlist_plot_late = session_list.rolo_switching;
sessionlist_plot_early = [];
stats_info = fig.plot_task_passive_compare(results_neural, sessionlist_plot_late,sessionlist_plot_early, plotOptions);
ylabel('','FontSize',ftsize,'Interpreter','latex')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');

monkeyR_oblique_bar_text = sprintf('Monkey R oblique: task performing: paired t-test, $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n',...
    stats_info.df_late, stats_info.tval_late, stats_info.p_late);
ylim([0,0.020])
%% monkey G cardinal
ax_3 = subplot(1,4,3);
set(ax_3,'position',get(ax_3,'position')+[0.02 0.02 0.01 -0.15]);
plotOptions.nTimebin = 4;
plotOptions.dataPlot = 'combine_delta';
plotOptions.taskType = 'cardinal';
plotOptions.plotType = 'bar';
sessionlist_plot_late = session_list.gremlin_switching;
sessionlist_plot_early = [];
stats_info = fig.plot_task_passive_compare(results_neural, sessionlist_plot_late,sessionlist_plot_early, plotOptions);
ylabel('$I_\textrm{redundancy}$','FontSize',ftsize,'Interpreter','latex')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');

monkeyG_cardinal_bar_text = sprintf('Monkey G cardinal: task performing: paired t-test, $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n',...
    stats_info.df_late, stats_info.tval_late, stats_info.p_late);

%% monkey G oblique
ax_4 = subplot(1,4,4);
set(ax_4,'position',get(ax_4,'position')+[0.04 0.02 0.01 -0.15]);
plotOptions.nTimebin = 4;
plotOptions.dataPlot = 'combine_delta';
plotOptions.taskType = 'oblique';
plotOptions.plotType = 'bar';
sessionlist_plot_late = session_list.gremlin_switching;
sessionlist_plot_early = [];
stats_info = fig.plot_task_passive_compare(results_neural, sessionlist_plot_late,sessionlist_plot_early, plotOptions);
ylabel('','FontSize',ftsize,'Interpreter','latex')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');

monkeyG_oblique_bar_text = sprintf('Monkey G oblique: task performing: paired t-test, $\\t(%d) = %.2f$, $p = \\num{%.2e}$; \n',...
    stats_info.df_late, stats_info.tval_late, stats_info.p_late);

%% add anotations
annotation('textbox',[0.22,0.95,0.15,0.04],'string','Monkey R','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.70,0.95,0.15,0.04],'string','Monkey G','fontsize',16,'FontWeight','bold','EdgeColor','none')
%% save figure as svg
print(save_name,'-dsvg','-vector')

%% write statistic report
fid = fopen(tex_bar_name,'wt');

fwrite(fid, ['t test \n',...
    monkeyR_cardinal_bar_text, monkeyR_oblique_bar_text, monkeyG_oblique_bar_text, monkeyG_cardinal_bar_text]);

fclose(fid)