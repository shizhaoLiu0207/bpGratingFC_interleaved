clear all
clc
close all
%%% Automaticaly make the supplementary behavior figure: learning index and
%%% temporal kernel
%%
global   bpGlobal ftsize
ftsize = 14;
bpGratingFCGlobal();
%%
figure
set(gcf, 'unit','inches','position',[0,0,12,7]);
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'v2_Figure_S1_suppBehav.svg');
%tex_name    = fullfile(save_folder ,'figure1_behav_stats_new.tex');

%% 1. Temporal psychometric kernel - monkey R, cardinal
ax_1 = subplot(2,4,1); 
set(ax_1,'position',get(ax_1,'position')+[-0.06 0.03 0.01 -0.06]);
load('../../../results/behavior/Rolo_psyKernel_table_final');

timebin = [100:100:1600];
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
fig.plot_temporal_psychometricKernel(psyKernel_table,sessionlist_cardinal_plot,'C',timebin); 
title('Cardinal','Color',bpGlobal.color_list.color_cardinal)
%% 2. Temporal psychometric kernel - monkey R, oblique
ax_2 = subplot(2,4,2); 
subjectCode = 'Gr';
set(ax_2,'position',get(ax_2,'position')+[-0.05 0.03 0.01 -0.06]);

sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;
fig.plot_temporal_psychometricKernel(psyKernel_table,sessionlist_oblique_plot,'O',timebin); 
ylabel('')
title('Oblique','Color',bpGlobal.color_list.color_oblique)
%% 3. Temporal psychometric kernel - monkey G, cardinal
ax_3 = subplot(2,4,3); 
set(ax_3,'position',get(ax_3,'position')+[-0.01 0.03 0.01 -0.06]);

load('../../../results/behavior/Gremlin_psyKernel_table');
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_temporal_psychometricKernel(psyKernel_table,sessionlist_cardinal_plot,'C',timebin); 
title('Cardinal','Color',bpGlobal.color_list.color_cardinal)
%% 4. Temporal psychometric kernel - monkey G, oblique
ax_4 = subplot(2,4,4); 
set(ax_4,'position',get(ax_4,'position')+[-0.01 0.03 0.01 -0.06]);
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;
fig.plot_temporal_psychometricKernel(psyKernel_table,sessionlist_oblique_plot,'O',timebin); 
ylabel('')
title('Oblique','Color',bpGlobal.color_list.color_oblique)
%% 5. Time course of learning index - monkey R
ax_4  = subplot(2,4,[5,6]);
set(ax_4,'position',get(ax_4,'position')+[-0.06 0 0.03 0]);
sessionlist_cardinal_plot = bpGlobal.rolo.session_list.switching;
sessionlist_oblique_plot  = bpGlobal.rolo.session_list.switching;


load('../../../results/behavior/Rolo_psyKernel_table_final');
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

%% 6. Time course of learning index - monkey G
ax_5  = subplot(2,4,[7,8]);
set(ax_5,'position',get(ax_5,'position')+[0 0 0.03 0]);
sessionlist_cardinal_plot = bpGlobal.gremlin.session_list.interleaved_real;
sessionlist_oblique_plot  = bpGlobal.gremlin.session_list.interleaved_real;

load('../../../results/behavior/Gremlin_psyKernel_table');

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

%% add annotations
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0.2,0.95,0.15,0.04],'string','Monkey R','fontsize',16,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.68,0.95,0.15,0.04],'string','Monkey G','fontsize',16,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0,0.985,0.1,0.04],'string','a','fontsize',30,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0,0.48,0.1,0.04],'string','b','fontsize',30,'FontWeight','bold','EdgeColor','none')
%%
print(save_name,'-dsvg','-vector')