clear all
clc
close all
%%% Automaticaly make the behavior figure with minimum subsequenct edit
%%% with Inkscape
%%
global   bpGlobal ftsize
ftsize = 14;
bpGratingFCGlobal();
%%
figure
set(gcf, 'unit','inches','position',[0,0,12,9]);
save_folder = '/Users/liushizhao/Documents/projects/bpGratingEx/figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'figure1_behav.svg');
tex_name    = fullfile(save_folder ,'figure1_behav_stats.tex');
%% 1. pychometric curves - monkey R
ax_1 = subplot(3,4,1); 
set(ax_1,'position',get(ax_1,'position')+[-0.06 0.03 0.01 -0.03]);
subjectCode = 'Ro';
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;

[h_c,~] = fig.plot_psychometric_curve_sessions(sessionlist_cardinal_plot,'C',subjectCode); hold on
[~,h_o] = fig.plot_psychometric_curve_sessions(sessionlist_oblique_plot,'O',subjectCode);

N_cardinal = numel(sessionlist_cardinal_plot);
N_oblique   = numel(sessionlist_oblique_plot);
text(-15,0.7,sprintf('N = %d',N_cardinal),'color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize);
text(5,0.3,sprintf('N = %d',N_oblique),'color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize);
xlabel('Signal level (%)','fontsize', 14)
%ylabel('Percent of 0^{\circ}/45^{\circ} choice')
title('Monkey R','fontsize',ftsize);

text(-0.2, 1.2, 'd', ...
         'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');
%% 2. pychometric curves - monkey G
ax_2 = subplot(3,4,2); 
subjectCode = 'Gr';
set(ax_2,'position',get(ax_2,'position')+[-0.05 0.03 0.01 -0.03]);
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;

[h_c,~] = fig.plot_psychometric_curve_sessions(sessionlist_cardinal_plot,'C',subjectCode); hold on
[~,h_o] = fig.plot_psychometric_curve_sessions(sessionlist_oblique_plot,'O',subjectCode);

N_cardinal = numel(sessionlist_cardinal_plot);
N_oblique   = numel(sessionlist_oblique_plot);
text(-15,0.7,sprintf('N = %d',N_cardinal),'color',bpGlobal.color_list.color_cardinal, 'fontsize', ftsize);
text(5,0.3,sprintf('N = %d',N_oblique),'color',bpGlobal.color_list.color_oblique, 'fontsize', ftsize);
xlabel('Signal level (%)','fontsize', 14)
title('Monkey G','fontsize',ftsize);
%% 3. spatial psychometric kernels - monkey R
ax_3 = subplot(3,4,3); 
set(ax_3,'position',get(ax_3,'position')+[-0.01 0.03 0.01 -0.03]);
subjectCode = 'Ro';
plotBias = 1;
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;
load('/Users/liushizhao/Documents/projects/bpGratingEx/results/behavior/Rolo_psyKernel_table_final');
yyaxis left
[h_c] = fig.plot_spatial_psychometricKernel(psyKernel_table,sessionlist_cardinal_plot,'C',plotBias); hold on
yyaxis left
[~,h_o] = fig.plot_spatial_psychometricKernel(psyKernel_table,sessionlist_oblique_plot,'O',plotBias);
yyaxis right 
ylabel('')
yyaxis left

set(gca,'fontsize', 11,'tickdir','out')
ylabel('Psychometric weight','fontsize',14)
xlabel('Orientation','fontsize', 14)
title('Monkey R','fontsize',ftsize);

text(-0.2, 1.2, 'e', ...
         'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');
%% 4. spatial psychometric kernels - monkey G
ax_4 = subplot(3,4,4); 
set(ax_4,'position',get(ax_4,'position')+[-0.01 0.03 0.01 -0.03]);
subjectCode = 'Gr';
plotBias = 1;
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;
load('/Users/liushizhao/Documents/projects/bpGratingEx/results/behavior/Gremlin_psyKernel_table');
yyaxis left
[h_c] = fig.plot_spatial_psychometricKernel(psyKernel_table,sessionlist_cardinal_plot,'C',plotBias); hold on
yyaxis left
[~,h_o] = fig.plot_spatial_psychometricKernel(psyKernel_table,sessionlist_oblique_plot,'O',plotBias);

set(gca,'fontsize', 11,'tickdir','out')
xlabel('Orientation','fontsize', 14)
yyaxis right
ylabel('Bias')
yyaxis left 
ylabel('')

title('Monkey G','fontsize',ftsize);
%% 5. Time course of learning index - monkey R
ax_4  = subplot(3,4,[5,6]);
set(ax_4,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;


load('/Users/liushizhao/Documents/projects/bpGratingEx/results/behavior/Rolo_psyKernel_table_final');

session_all = [sessionlist_cardinal_plot;sessionlist_oblique_plot];
psyKernel_table(~ismember({psyKernel_table(:).sessionName},session_all)) = [];


plotOptions.CI_level = 68;
plotOptions.use_real_date = 0;
plotOptions.dataPlot = 'match_amplitude';
plotOptions.x_axis = 'part_idx';
[~,~,y_c] = fig.plot_behav_psyKernel_timecourse(psyKernel_table,'Cardinal',sessionlist_cardinal_plot,plotOptions); hold on
[~,~,y_o] = fig.plot_behav_psyKernel_timecourse(psyKernel_table,'Oblique',sessionlist_oblique_plot,plotOptions);

[~,p,~,stats] = ttest(y_c,y_o);
monkeyR_compare_task_string = sprintf('Monkey R, $t(%d) = %.2f$, $p = \\num{%.2e}$ \n', stats.df, stats.tstat, p);

[rho, pval] = corr([1:numel(y_c)]', y_c', 'Type','Spearman');
monkeyR_cardinal_corr_string = sprintf('Monkey R, cardinal, $r = %.2f$, $p = \\num{%.2e}$ \n', rho, pval);
[rho, pval] = corr([1:numel(y_o)]', y_o', 'Type','Spearman');
monkeyR_oblique_corr_string = sprintf('Monkey R, oblique, $r = %.2f$, $p = \\num{%.2e}$ \n', rho, pval);


ylim([-0.01,0.17])
idx_plot = find(ismember({psyKernel_table(:).sessionName},session_all));
xlimit = [min(idx_plot),max(idx_plot)];
xlim(xlimit);
line([0,xlimit(2)+1],[0,0],'color','black','linestyle','--','linewidth',2)
set(gca,'TickDir','out')
xlabel('Session number');
ylabel('Learning index');

title('Monkey R')

text(-0.1, 1.1, 'f', ...
         'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');
%% 6. Time course of learning index - monkey G
ax_5  = subplot(3,4,[7,8]);
set(ax_5,'position',get(ax_5,'position')+[0 0 0.03 0]);
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;

load('/Users/liushizhao/Documents/projects/bpGratingEx/results/behavior/Gremlin_psyKernel_table');

session_all = [sessionlist_cardinal_plot;sessionlist_oblique_plot];
psyKernel_table(~ismember({psyKernel_table(:).sessionName},session_all)) = [];


plotOptions.CI_level = 68;
plotOptions.use_real_date = 0;
plotOptions.dataPlot = 'match_amplitude';
plotOptions.x_axis = 'part_idx';
[h_c, ~, y_c] = fig.plot_behav_psyKernel_timecourse(psyKernel_table,'Cardinal',sessionlist_cardinal_plot,plotOptions); hold on
[h_o, ~, y_o] = fig.plot_behav_psyKernel_timecourse(psyKernel_table,'Oblique',sessionlist_oblique_plot,plotOptions);


[~,p,~,stats] = ttest(y_c,y_o);
monkeyG_compare_task_string = sprintf('Monkey G, $t(%d) = %.2f$, $p = \\num{%.2e}$ \n', stats.df, stats.tstat, p);

[rho, pval] = corr([1:numel(y_c)]', y_c', 'Type','Spearman');
monkeyG_cardinal_corr_string = sprintf('Monkey G, cardinal, $r = %.2f$, $p = \\num{%.2e}$ \n', rho, pval);
[rho, pval] = corr([1:numel(y_o)]', y_o', 'Type','Spearman');
monkeyG_oblique_corr_string = sprintf('Monkey G, oblique, $r = %.2f$, $p = \\num{%.2e}$ \n', rho, pval);

ylim([-0.01,0.17])
idx_plot = find(ismember({psyKernel_table(:).sessionName},session_all));
xlimit = [min(idx_plot),max(idx_plot)];
xlim(xlimit);
line([0,xlimit(2)+1],[0,0],'color','black','linestyle','--','linewidth',2)

set(gca,'TickDir','out')
xlabel('Session number');
ylabel('Learning index');

title('Monkey G')


%% 7. Within and cross prediction - monkey R
ax_6  = subplot(3,4,9);
%set(ax_,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
set(ax_6,'position',get(ax_6,'position')+[-0.06 -0.03 0.01 -0.03]);
ax_7 = subplot(3,4,10);
set(ax_7,'position',get(ax_7,'position')+[-0.05 -0.03 0.01 -0.03]);
plotOptions = struct();
plotOptions.figurestyle = 'scatterplot';
sessionStr_list = bpGlobal.rolo.session_list.switching;
dataPath =  '../../../results/behavior/psyKernel_interleaved_crossvalidate/Ro';
stats_info = fig.plot_cross_predict(sessionStr_list,dataPath,plotOptions, ax_6, ax_7);

axes(ax_6)
text(-1.15,-0.55,sprintf('$t(%d) = %.2f^{%s}$',stats_info.t_cardinal.df,stats_info.t_cardinal.tstat,fig.p2star(stats_info.t_cardinal.p)),...
    'Interpreter','latex','Color',bpGlobal.color_list.color_cardinal,'fontsize',ftsize);

text(-0.3,1.3, 'g', ...
         'Units', 'normalized', 'FontSize', 28, 'FontWeight', 'bold');

text(0.9,1.25, 'Monkey R', ...
         'Units', 'normalized', 'FontSize', ftsize+2, 'FontWeight', 'bold');

axes(ax_7)
ylabel('')
text(-1,-0.5,sprintf('$t(%d) = %.2f^{%s}$',stats_info.t_oblique.df,stats_info.t_oblique.tstat,fig.p2star(stats_info.t_oblique.p)),...
    'Interpreter','latex','Color',bpGlobal.color_list.color_oblique,'fontsize',ftsize);

monkeyR_cross_pred_string = sprintf(['Monkey R, cardinal, $t(%d) = %.2f$, $p = \\num{%.2e}$; \n' ...
    'Monkey R, oblique, $t(%d) = %.2f$, $p = \\num{%.2e}$ \n'],...
    stats_info.t_cardinal.df, stats_info.t_cardinal.tstat, stats_info.t_cardinal.p,...
    stats_info.t_oblique.df, stats_info.t_oblique.tstat, stats_info.t_oblique.p);


%% 8. Within and cross prediction - monkey G

ax_8  = subplot(3,4,11);
%set(ax_,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
set(ax_8,'position',get(ax_8,'position')+[-0.01 -0.03 0.01 -0.03]);
ax_9 = subplot(3,4,12);
set(ax_9,'position',get(ax_9,'position')+[0.01 -0.03 0.01 -0.03]);
plotOptions = struct();
plotOptions.figurestyle = 'scatterplot';
sessionStr_list = bpGlobal.gremlin.session_list.interleaved_real;
dataPath =  '../../../results/behavior/psyKernel_interleaved_crossvalidate/Gr';
stats_info = fig.plot_cross_predict(sessionStr_list,dataPath,plotOptions, ax_8, ax_9);


axes(ax_8)
text(-1.9,-0.55,sprintf('$t(%d) = %.2f^{%s}$',stats_info.t_cardinal.df,stats_info.t_cardinal.tstat,fig.p2star(stats_info.t_cardinal.p)),...
    'Interpreter','latex','Color',bpGlobal.color_list.color_cardinal,'fontsize',ftsize);
text(0.9,1.25, 'Monkey G', ...
         'Units', 'normalized', 'FontSize', ftsize+2, 'FontWeight', 'bold');
axes(ax_9)
ylabel('')
text(-1,-0.5,sprintf('$t(%d) = %.2f^{%s}$',stats_info.t_oblique.df,stats_info.t_oblique.tstat,fig.p2star(stats_info.t_oblique.p)),...
    'Interpreter','latex','Color',bpGlobal.color_list.color_oblique,'fontsize',ftsize);

monkeyG_cross_pred_string = sprintf(['Monkey G, cardinal, $t(%d) = %.2f$, $p = \\num{%.2e}$; \n' ...
    'Monkey G, oblique, $t(%d) = %.2f$, $p = \\num{%.2e}$ \n'],...
    stats_info.t_cardinal.df, stats_info.t_cardinal.tstat, stats_info.t_cardinal.p,...
    stats_info.t_oblique.df, stats_info.t_oblique.tstat, stats_info.t_oblique.p);
%%
print(save_name,'-dsvg','-vector')
%%
fid = fopen(tex_name,'wt');

fwrite(fid, ['Correlation between learning index and time: ', ...
    monkeyR_cardinal_corr_string, monkeyR_oblique_corr_string, ...
    monkeyG_cardinal_corr_string, monkeyG_oblique_corr_string,...
    'Compare learning index of two tasks: ',...
    monkeyR_compare_task_string, monkeyG_compare_task_string, ...
    'Compare within and cross prediction:',...
    monkeyR_cross_pred_string, monkeyG_cross_pred_string]);

fclose(fid);
