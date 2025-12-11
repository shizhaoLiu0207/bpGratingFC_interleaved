clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
%filter_name = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl';
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
load(fullfile(saveFolder, 'results_SubsampleOrganized_perCohr_fisherInfo_all_sessions'));
results_all_perCohr = results_cross_sizeControl_perCohr;
results_all_perCohr = get_sample_CI_cross(results_all_perCohr);
%% 1. bar plot of I_real and I_cross
plot_per_cohr = false;
save_folder = '../../figures/figures_final/fisher_info_bar';

if plot_per_cohr
    results_plot = results_all_perCohr;
    save_name = fullfile(save_folder, sprintf('cross_fisher_info_shuffle_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('cross_fisher_info_shuffle_bar_perCohr.tex'));
    ylim_rolo = 0.038;
    ylim_gremlin = 0.008;
    ytext_rolo = 0.036;
    ytext_gremlin = 0.0075;
else
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('cross_fisher_shuffle_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('cross_fisher_shuffle_info_bar.tex'));
    ylim_rolo = 0.050;
    ylim_gremlin = 0.009;
    ytext_rolo = 0.045;
    ytext_gremlin = 0.0085;
end
figure
set(gcf,'unit','inches','position',[0,0,12,5])
ax_1 = subplot(3,2,[1,3]); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);

plotOptions = struct();
plotOptions.plot_data = 'shuffle';
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_Info(results_plot, session_list_plot, plotOptions); 

stats_real_cross_string_monkeyR_cardinal = sprintf('Monkey R cardinal, $I_\\textrm{real}$ v.s. $I_\\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_real_cross_string_monkeyR_oblique = sprintf('Monkey R oblique, $I_\\textrm{real}$ v.s. $I_\\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);

ylim([0, ylim_rolo])
text(0.5,ytext_rolo,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_rolo,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey R')

ax_2 = subplot(3,2,5); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0 0.03 -0.02]);
plotOptions = struct();
plotOptions.plot_data = 'shuffle';
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 16;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
stats_info_diff_info_rolo = fig_it.plot_diff_errorbar(results_plot, session_list_plot, plotOptions);
stats_diff_info_string_monkeyR_cardinal = sprintf('Monkey R cardinal, Diff between $I_\\textrm{real}$ and $I_\\textrm{cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_info_rolo.mean_cardinal,stats_info_diff_info_rolo.std_cardinal);
stats_diff_info_string_monkeyR_oblique =  sprintf('Monkey R oblique, Diff between $I_\\textrm{real}$ and $I_\\textrm{cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_info_rolo.mean_oblique,stats_info_diff_info_rolo.std_oblique);
ylim([-10,30])

ax_3 = subplot(3,2,[2,4]); hold on
set(ax_3,'position',get(ax_3,'position')+[0 0.04 0.02 -0.03]);
plotOptions = struct();
plotOptions.plot_data = 'shuffle';
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
%fig.plot_bar_Info(results_all, session_list_gremlin_good,plotOptions); title('Monkey G Good sessions');
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_Info(results_plot, session_list_plot,plotOptions);

stats_real_cross_string_monkeyG_cardinal = sprintf('Monkey G cardinal, $I_\\textrm{real}$ v.s. $I_\\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_real_cross_string_monkeyG_oblique = sprintf('Monkey G oblique, $I_\\textrm{real}$ v.s. $I_\\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);

ylim([0,ylim_gremlin])
text(0.5,ytext_gremlin,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_gremlin,'Oblique','color','blue','FontSize',14,'FontWeight','bold')

title('Monkey G');

ax_4 = subplot(3,2,6); hold on
set(ax_4,'position',get(ax_4,'position')+[0 0 0.02 -0.02]);
plotOptions = struct();
plotOptions.plot_data = 'shuffle';
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 16;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
stats_info_diff_info_gremlin = fig_it.plot_diff_errorbar(results_plot, session_list_plot, plotOptions);
stats_diff_info_string_monkeyG_cardinal = sprintf('Monkey G cardinal, Diff between $I_\\textrm{real}$ and $I_\\textrm{cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_info_gremlin.mean_cardinal, stats_info_diff_info_gremlin.std_cardinal);
stats_diff_info_string_monkeyG_oblique =  sprintf('Monkey G oblique, Diff between $I_\\textrm{real}$ and $I_\\textrm{cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_info_gremlin.mean_oblique, stats_info_diff_info_gremlin.std_oblique);
ylim([-10,30])

print(save_name, '-dsvg');

fid = fopen(tex_name,'wt');
fwrite(fid,[stats_real_cross_string_monkeyR_cardinal, stats_real_cross_string_monkeyR_oblique,...
            stats_real_cross_string_monkeyG_cardinal, stats_real_cross_string_monkeyG_oblique,...
            stats_diff_info_string_monkeyR_cardinal, stats_diff_info_string_monkeyR_oblique,...
            stats_diff_info_string_monkeyG_cardinal, stats_diff_info_string_monkeyG_oblique]);
fclose(fid);