clear all
clc
close all
%%

global   bpGlobal ftsize
ftsize                                     = 14;
bpGratingFCGlobal();
versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
[session_list, results_neural]  = data_prep(versionName);
%%
figure
set(gcf,'unit','inches','position',[0,0,12,4])
save_folder = '../../../figures/figures_auto_interleavedpaper';

save_name = fullfile(save_folder,'Figure_S8_withinTrial_Iredundancy_sizeControl.svg');
tex_name = fullfile(save_folder,'Figure_S8_withinTrial_Iredundancy_sizeControl.tex');
%% monkey R cardinal
ax_1 = subplot(1,4,1);
set(ax_1,'position',get(ax_1,'position')+[-0.05 0.15 0.01 -0.20]);
plotOptions.includeOri = 0;
plotOptions.nTimebin   = 8;
plotOptions.dataPlot   = 'combine_delta'; 

sessionlist_plot_late = session_list.rolo_switching;
sessionlist_plot_early = [];
plotOptions.taskType = 'cardinal';

stats_info = fig.plot_timebins_outside(results_neural, sessionlist_plot_late, sessionlist_plot_early, plotOptions);

ylabel('$I_\textrm{redundancy}$','FontSize',ftsize+2,'Interpreter','latex')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
add_stats_annotations(stats_info, ax_1, bpGlobal.color_list.color_cardinal);

monkeyR_cardinal_text = sprintf(['Monket R cardinal, linear regression. \n' ...
    '$\\beta = \\num{%.4e}, p = \\num{%.4e}$'], ...
    stats_info.beta_late, stats_info.p_late);
%% monkey R oblique
ax_2 = subplot(1,4,2);
set(ax_2,'position',get(ax_2,'position')+[-0.03 0.15 0.01 -0.20]);
plotOptions.includeOri = 0;
plotOptions.nTimebin   = 8;
plotOptions.dataPlot   = 'combine_delta'; 

sessionlist_plot_late = session_list.rolo_switching;
sessionlist_plot_early = [];
plotOptions.taskType = 'oblique';

stats_info = fig.plot_timebins_outside(results_neural, sessionlist_plot_late, sessionlist_plot_early, plotOptions);

ylabel('')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
add_stats_annotations(stats_info, ax_2, bpGlobal.color_list.color_oblique);

monkeyR_oblique_text = sprintf(['Monket R oblique, linear regression. \n' ...
    'Late-learning: $\\beta = \\num{%.4e}, p = \\num{%.4e}$'], ...
     stats_info.beta_late, stats_info.p_late);
%% monkey G cardinal
ax_3 = subplot(1,4,3);
set(ax_3,'position',get(ax_3,'position')+[0.02 0.15 0.01 -0.20]);
plotOptions.includeOri = 0;
plotOptions.nTimebin   = 8;
plotOptions.dataPlot   = 'combine_delta'; 

sessionlist_plot_late = session_list.gremlin_switching;
sessionlist_plot_early = [];
plotOptions.taskType = 'cardinal';

stats_info = fig.plot_timebins_outside(results_neural, sessionlist_plot_late, sessionlist_plot_early, plotOptions);

ylabel('$I_\textrm{redundancy}$','FontSize',ftsize+2,'Interpreter','latex')
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
add_stats_annotations(stats_info, ax_3, bpGlobal.color_list.color_cardinal);

monkeyG_cardinal_text = sprintf(['Monket G cardinal, linear regression. \n' ...
    '$\\beta = \\num{%.4e}, p = \\num{%.4e}$'], ...
    stats_info.beta_late, stats_info.p_late);

%% monkey G oblique
ax_4 = subplot(1,4,4);
set(ax_4,'position',get(ax_4,'position')+[0.04 0.15 0.01 -0.20]);
plotOptions.includeOri = 0;
plotOptions.nTimebin   = 8;
plotOptions.dataPlot   = 'combine_delta'; 

sessionlist_plot_late = session_list.gremlin_switching;
sessionlist_plot_early = [];
plotOptions.taskType = 'oblique';

stats_info = fig.plot_timebins_outside(results_neural, sessionlist_plot_late, sessionlist_plot_early, plotOptions);

ylabel('')
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
add_stats_annotations(stats_info, ax_4, bpGlobal.color_list.color_oblique);

monkeyG_oblique_text = sprintf(['Monket G oblique, linear regression. \n' ...
    'Late-learning: $\\beta = \\num{%.4e}, p = \\num{%.4e}$'], ...
     stats_info.beta_late, stats_info.p_late);

%% add anotations
annotation('textbox',[0.22,0.97,0.15,0.04],'string','Monkey R','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.70,0.97,0.15,0.04],'string','Monkey G','fontsize',16,'FontWeight','bold','EdgeColor','none')
%% save figure as svg
print(save_name,'-dsvg','-vector')

%% write statistic report
fid = fopen(tex_name,'wt');

fwrite(fid, ['Linear regression statistics: \n',...
    monkeyR_cardinal_text, monkeyR_oblique_text,  monkeyG_cardinal_text, monkeyG_oblique_text]);

fclose(fid)
%%
function add_stats_annotations(stats_info, ax_plot, plotColor)
    global ftsize
    late_stats_string = sprintf('$\\beta = %.2e^{%s}$', stats_info.beta_late,  fig.p2star(stats_info.p_late));
   
    ax_pos = get(ax_plot,'position');
    annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.23, ax_pos(2)-0.23],'LineWidth',2,'Color',plotColor)
    
    annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.23,0.15,0.04],'string',late_stats_string,'fontsize',ftsize,'EdgeColor','none','Interpreter','latex')

end