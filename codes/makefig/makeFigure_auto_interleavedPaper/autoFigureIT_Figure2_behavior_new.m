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
set(gcf, 'unit','inches','position',[0,0,12,6]);
save_folder = '/Users/liushizhao/Documents/projects/bpGratingEx/figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'v2_Figure_2_exp_behav_interleaved.svg');
tex_name    = fullfile(save_folder ,'v2_Figure_2_exp_behav_interleaved.tex');
%% 1. pychometric curves - monkey R
ax_1 = subplot(2,4,1); 
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
ax_2 = subplot(2,4,2); 
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
ax_3 = subplot(2,4,3); 
set(ax_3,'position',get(ax_3,'position')+[-0.01 0.03 0.01 -0.03]);
subjectCode = 'Ro';
plotBias = 1;
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;
load('../../../results/behavior/Rolo_psyKernel_table_final');
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
ax_4 = subplot(2,4,4); 
set(ax_4,'position',get(ax_4,'position')+[-0.01 0.03 0.01 -0.03]);
subjectCode = 'Gr';
plotBias = 1;
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;
load('../../../results/behavior/Gremlin_psyKernel_table');
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


%% 7. Within and cross prediction - monkey R
ax_6  = subplot(2,4,5);
%set(ax_,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
set(ax_6,'position',get(ax_6,'position')+[-0.06 -0.03 0.01 -0.03]);
ax_7 = subplot(2,4,6);
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

ax_8  = subplot(2,4,7);
%set(ax_,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
set(ax_8,'position',get(ax_8,'position')+[-0.01 -0.03 0.01 -0.03]);
ax_9 = subplot(2,4,8);
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
