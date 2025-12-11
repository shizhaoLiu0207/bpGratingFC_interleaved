clear all
clc
close all

global   bpGlobal ftsize
ftsize                                     = 14;
bpGratingFCGlobal();


versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';

[session_list, results_neural] = data_prep(versionName);
%%
figure
set(gcf,'unit','inches','position',[0,0,12,8])
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name = fullfile(save_folder,'Figure_S7_FisherInfo_sizeControl_nobehav.svg');
tex_name = fullfile(save_folder,'Figure_S7_FisherInfo_sizeControl_report.tex');
%% 1. Fisher information timecourse - Monkey R
ax_1 = subplot(2,4,[1,2]);
set(ax_1,'position',get(ax_1,'position')+[-0.04 0 0.03 0]);
plotOptions.dataName = 'combine_delta';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.uncertaintyName = 'CI';
%plotOptions.plotBehav = 1;
plotOptions.plotBehav = 0;
plotOptions.plotZeroLine = 1;
stats_info = fig_it.plot_timecourse_multiple_interleaved(results_neural,'Ro',session_list.rolo_switching,plotOptions);
%ylim([-0.003,0.09])
ylim([-0.003,0.08])
set(gca,'TickDir','out','fontsize',ftsize)
ylabel('Linear Fisher information','fontsize',ftsize+2)
title('Monkey R', 'fontsize',ftsize+2)

if plotOptions.plotBehav
    text(1,0.07, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df_fisher_behav_cardinal, stats_info.r_fisher_behav_cardinal, fig.p2star(stats_info.p_fisher_behav_cardinal)),...
        'Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'Interpreter','latex')
    text(27,0.07, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df_fisher_behav_oblique, stats_info.r_fisher_behav_oblique, fig.p2star(stats_info.p_fisher_behav_oblique)),...
        'Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'Interpreter','latex')
end
monkeyR_behav_Fisher_text = sprintf(['Correlation between neural and behavioral fisher information: \n' ...
        'Monkey R, Cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$; ' ...
        'Monkey R, Oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n'], ...
stats_info.df_fisher_behav_cardinal, stats_info.r_fisher_behav_cardinal, stats_info.p_fisher_behav_cardinal, ...
stats_info.df_fisher_behav_oblique, stats_info.r_fisher_behav_oblique, stats_info.p_fisher_behav_oblique);

%% 2. Fisher information timecourse - Monkey G
ax_2 = subplot(2,4,[3,4]);
set(ax_2,'position',get(ax_2,'position')+[0.02 0 0.03 0]);
plotOptions.dataName = 'combine_delta';
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.uncertaintyName = 'CI';
%plotOptions.plotBehav = 1;
plotOptions.plotBehav = 0;
plotOptions.plotZeroLine = 1;
stats_info = fig_it.plot_timecourse_multiple_interleaved(results_neural,'Gr',session_list.gremlin_switching,plotOptions);
ylim([-0.001,0.015])
set(gca,'TickDir','out','fontsize',ftsize)
ylabel('Linear Fisher information','fontsize',ftsize+2)
title('Monkey G', 'fontsize',ftsize+2)

if  plotOptions.plotBehav
    text(52,0.016, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df_fisher_behav_cardinal, stats_info.r_fisher_behav_cardinal, fig.p2star(stats_info.p_fisher_behav_cardinal)),...
        'Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'Interpreter','latex')
    text(5,0.016, sprintf('$r(%d) = %.2f^{%s}$', stats_info.df_fisher_behav_oblique, stats_info.r_fisher_behav_oblique, fig.p2star(stats_info.p_fisher_behav_oblique)),...
        'Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'Interpreter','latex')
end
monkeyG_behav_Fisher_text = sprintf(['Correlation between neural and behavioral fisher information: \n' ...
    'Monkey G, Cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$; ' ...
    'Monkey G, Oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n'], ...
    stats_info.df_fisher_behav_cardinal, stats_info.r_fisher_behav_cardinal, stats_info.p_fisher_behav_cardinal, ...
    stats_info.df_fisher_behav_oblique, stats_info.r_fisher_behav_oblique, stats_info.p_fisher_behav_oblique);

%% 3. Scatter (Real/shuffle vs learning index) - Monkey R cardinal
ax_3 = subplot(2,4,5);
set(ax_3,'position',get(ax_3,'position')+[-0.04 0.04 0.01 -0.05]);
plotOptions =  struct();
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.taskPlot = 'cardinal';
[stats_info,h] = fig_it.plot_scatterplot_real_shuffle_interleaved(results_neural,'Ro',plotOptions);
set(gca,'TickDir','out','fontsize',ftsize)
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
ylabel('Linear Fisher information','fontsize',ftsize+2)
xlabel('Learning Index')



real_stats_string = sprintf('$I_\\textrm{real}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_real, fig.p2star(stats_info.p_real));
shuffle_stats_string = sprintf('$I_\\textrm{shuffle}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_shuffle, fig.p2star(stats_info.p_shuffle));

ax_pos = get(ax_3,'position');
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.09, ax_pos(2)-0.09],'LineWidth',2)
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.12, ax_pos(2)-0.12],'LineWidth',2,'LineStyle','--')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.12,0.15,0.04],'string',real_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_cardinal,'Interpreter','latex')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.15,0.15,0.04],'string',shuffle_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_cardinal,'Interpreter','latex')


monkeyR_cardinal_real_text = sprintf('Monkey R, with $I_\\textrm{real}$, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_real, stats_info.p_real);
monkeyR_cardinal_shuffle_text = sprintf('Monkey R, with $I_\\textrm{shuffle}$, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_shuffle, stats_info.p_shuffle);
%% 4. Scatter (Real/shuffle vs learning index) - Monkey R oblique
ax_4 = subplot(2,4,6);
set(ax_4,'position',get(ax_4,'position')+[-0.03 0.04 0.01 -0.05]);
plotOptions =  struct();
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.taskPlot = 'oblique';
[stats_info,h]  =  fig_it.plot_scatterplot_real_shuffle_interleaved(results_neural,'Ro',plotOptions);
set(gca,'TickDir','out','fontsize',ftsize)
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
ylabel('','fontsize',ftsize+2)
xlabel('Learning Index')

real_stats_string = sprintf('$I_\\textrm{real}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_real, fig.p2star(stats_info.p_real));
shuffle_stats_string = sprintf('$I_\\textrm{shuffle}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_shuffle, fig.p2star(stats_info.p_shuffle));
ax_pos = get(ax_4,'position');
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.09, ax_pos(2)-0.09],'LineWidth',2)
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.12, ax_pos(2)-0.12],'LineWidth',2,'LineStyle','--')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.12,0.15,0.04],'string',real_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_oblique,'Interpreter','latex')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.15,0.15,0.04],'string',shuffle_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_oblique,'Interpreter','latex')

monkeyR_oblique_real_text = sprintf('Monkey R, with $I_\\textrm{real}$, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_real, stats_info.p_real);
monkeyR_oblique_shuffle_text = sprintf('Monkey R, with $I_\\textrm{shuffle}$, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_shuffle, stats_info.p_shuffle);



%% 5. Scatter  (Real/shuffle vs learning index) - Monkey G cardinal
ax_5 = subplot(2,4,7);
set(ax_5,'position',get(ax_5,'position')+[0.01 0.04 0.01 -0.05]);
plotOptions =  struct();
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.taskPlot = 'cardinal';
[stats_info,h] =  fig_it.plot_scatterplot_real_shuffle_interleaved(results_neural,'Gr',plotOptions);
set(gca,'TickDir','out','fontsize',ftsize)
title('Cardinal Task','Color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize,'FontWeight','bold');
ylabel('Linear Fisher information','fontsize',ftsize+2)
xlabel('Learning Index')

real_stats_string = sprintf('$I_\\textrm{real}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_real, fig.p2star(stats_info.p_real));
shuffle_stats_string = sprintf('$I_\\textrm{shuffle}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_shuffle, fig.p2star(stats_info.p_shuffle));
ax_pos = get(ax_5,'position');
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.09, ax_pos(2)-0.09],'LineWidth',2)
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.12, ax_pos(2)-0.12],'LineWidth',2,'LineStyle','--')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.12,0.15,0.04],'string',real_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_cardinal,'Interpreter','latex')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.15,0.15,0.04],'string',shuffle_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_cardinal,'Interpreter','latex')


monkeyG_cardinal_real_text = sprintf('Monkey G, with $I_\\textrm{real}$, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e$} \n' , ...
    stats_info.df, stats_info.r_real, stats_info.p_real);
monkeyG_cardinal_shuffle_text = sprintf('Monkey G, with $I_\\textrm{shuffle}$, cardinal task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_shuffle, stats_info.p_shuffle);

%% 6. Scatter (Real/shuffle vs learning index) - Monkey G oblique
ax_6 = subplot(2,4,8);
set(ax_6,'position',get(ax_6,'position')+[0.03 0.04 0.01 -0.05]);
plotOptions =  struct();
plotOptions.tBin_toPlot = 0;
plotOptions.sessionType = 'mainTask';
plotOptions.taskPlot = 'oblique';
[stats_info,h] =  fig_it.plot_scatterplot_real_shuffle_interleaved(results_neural,'Gr',plotOptions);
set(gca,'TickDir','out','fontsize',ftsize)
title('Oblique Task','Color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize,'FontWeight','bold');
ylabel('','fontsize',ftsize+2)


real_stats_string = sprintf('$I_\\textrm{real}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_real, fig.p2star(stats_info.p_real));
shuffle_stats_string = sprintf('$I_\\textrm{shuffle}: r(%d) = %.2f^{%s}$', stats_info.df, stats_info.r_shuffle, fig.p2star(stats_info.p_shuffle));
ax_pos = get(ax_6,'position');
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.09, ax_pos(2)-0.09],'LineWidth',2)
annotation('line',[ax_pos(1)-0.016, ax_pos(1) + 0.008],[ax_pos(2)-0.12, ax_pos(2)-0.12],'LineWidth',2,'LineStyle','--')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.12,0.15,0.04],'string',real_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_oblique,'Interpreter','latex')
annotation('textbox',[ax_pos(1) + 0.01,ax_pos(2)-0.15,0.15,0.04],'string',shuffle_stats_string,'fontsize',ftsize,'EdgeColor','none','Color',bpGlobal.color_list.color_oblique,'Interpreter','latex')

monkeyG_oblique_real_text = sprintf('Monkey G, with $I_\\textrm{real}$, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_real, stats_info.p_real);
monkeyG_oblique_shuffle_text = sprintf('Monkey G, with $I_\\textrm{shuffle}$, oblique task: $\\rho(%d) = %.2f$, $p = \\num{%.2e}$ \n' , ...
    stats_info.df, stats_info.r_shuffle, stats_info.p_shuffle);
%%  add some annotations

annotation('textbox',[0.23,0.47,0.15,0.04],'string','Monkey R','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.70,0.47,0.15,0.04],'string','Monkey G','fontsize',16,'FontWeight','bold','EdgeColor','none')


%% save figure as svg
print(save_name,'-dsvg','-vector')
%% write statistic report
fid = fopen(tex_name,'wt');

fwrite(fid, [monkeyR_behav_Fisher_text, monkeyG_behav_Fisher_text,...
    monkeyR_cardinal_real_text, monkeyR_oblique_real_text, monkeyG_oblique_real_text, monkeyG_cardinal_real_text,...
     monkeyR_cardinal_shuffle_text, monkeyR_oblique_shuffle_text, monkeyG_oblique_shuffle_text, monkeyG_cardinal_shuffle_text]);

fclose(fid);

