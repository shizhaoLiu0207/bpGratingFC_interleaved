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
    save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('cross_fisher_info_bar_perCohr.tex'));
    ylim_rolo = 0.038;
    ylim_gremlin = 0.008;
    ytext_rolo = 0.036;
    ytext_gremlin = 0.0075;
else
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('cross_fisher_info_bar.tex'));
    ylim_rolo = 0.038;
    ylim_gremlin = 0.007;
    ytext_rolo = 0.035;
    ytext_gremlin = 0.0065;
end
figure
set(gcf,'unit','inches','position',[0,0,10,5])
ax_1 = subplot(3,2,[1,3]); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);

plotOptions = struct();
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



%% 2. bar plot of I_redundacy, I_redundacy_cross in percent
plot_per_cohr = false;
save_folder = '../../figures/figures_final/fisher_info_bar';


if plot_per_cohr % & plotOptions.plotPercent
    results_plot = results_all_perCohr;
    save_name = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar_perCohr.tex'));
    ylim_rolo = 60;
    ylim_gremlin_min = -10;
    ylim_gremlin_max = 80;
    ytext_rolo = 50;
    ytext_gremlin = 70;

    y_diff_lim_rolo = 20;
    y_diff_lim_gremlin = 80;
elseif ~plot_per_cohr % & plotOptions.plotPercent
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar.tex'));
    ylim_rolo = 60;
    ylim_gremlin_min = -10;
    ylim_gremlin_max = 70;
    ytext_rolo = 50;
    ytext_gremlin = 60;

    y_diff_lim_rolo = 40
    y_diff_lim_gremlin = 40;
end

figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,10,5])
ax_1 = subplot(3,2,[1,3]); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_plot, session_list_plot, plotOptions);

stats_redundacy_within_cross_string_monkeyR_cardinal = sprintf(['Monkey R cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);


stats_redundacy_within_cross_string_monkeyR_oblique = sprintf(['Monkey R oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);


ylim([0, ylim_rolo])
text(0.5,ytext_rolo,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_rolo,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey R');


ax_2 = subplot(3,2,5); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0 0.03 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 16;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
stats_info_diff_delta_rolo = fig_it.plot_diff_errorbar(results_plot, session_list_plot, plotOptions);
ylim([-10,y_diff_lim_rolo])

stats_diff_delta_string_monkeyR_cardinal = sprintf('Monkey R cardinal, Diff between $I_\\textrm{redundancy,within}$ and $I_\\textrm{redundancy,cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_delta_rolo.mean_cardinal,stats_info_diff_delta_rolo.std_cardinal);
stats_diff_delta_string_monkeyR_oblique =  sprintf('Monkey R oblique, Diff between $I_\\textrm{redundancy,within}$ and $I_\\textrm{redundancy,cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_delta_rolo.mean_oblique,stats_info_diff_delta_rolo.std_oblique);


ax_3 = subplot(3,2,[2,4]); hold on
set(ax_3,'position',get(ax_3,'position')+[0 0.04 0.02 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_plot, session_list_plot, plotOptions); 
stats_redundacy_within_cross_string_monkeyG_cardinal = sprintf(['Monkey G cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);


stats_redundacy_within_cross_string_monkeyG_oblique = sprintf(['Monkey G oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);
ylim([ylim_gremlin_min,ylim_gremlin_max])
text(0.5,ytext_gremlin,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_gremlin,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey G');

ax_4 = subplot(3,2,6); hold on
set(ax_4,'position',get(ax_4,'position')+[0 0 0.02 -0.02]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 16;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
stats_info_diff_delta_gremlin = fig_it.plot_diff_errorbar(results_plot, session_list_plot, plotOptions);
ylim([-10,y_diff_lim_gremlin])
stats_diff_delta_string_monkeyG_cardinal = sprintf('Monkey R cardinal, Diff between $I_\\textrm{redundancy,within}$ and $I_\\textrm{redundancy,cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_delta_gremlin.mean_cardinal,stats_info_diff_delta_gremlin.std_cardinal);
stats_diff_delta_string_monkeyG_oblique =  sprintf('Monkey R oblique, Diff between $I_\\textrm{redundancy,within}$ and $I_\\textrm{redundancy,cross}$: $Mean = %.2f$, $s.t.d = %.2f$\n',...
                                        stats_info_diff_delta_gremlin.mean_oblique,stats_info_diff_delta_gremlin.std_oblique);


print(save_name, '-dsvg');
fid = fopen(tex_name,'wt');
fwrite(fid,[stats_redundacy_within_cross_string_monkeyR_cardinal, stats_redundacy_within_cross_string_monkeyR_oblique,...
            stats_redundacy_within_cross_string_monkeyG_cardinal, stats_redundacy_within_cross_string_monkeyG_oblique,...
            stats_diff_delta_string_monkeyR_cardinal, stats_diff_delta_string_monkeyR_oblique,...
            stats_diff_delta_string_monkeyG_cardinal, stats_diff_delta_string_monkeyG_oblique]);
fclose(fid);

%% 3. bar plot of I_redundacy, I_redundacy_cross 
plot_per_cohr = true;
save_folder = '../../figures/figures_final/fisher_info_bar';
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotPercent = false;
plotOptions.ftsize = 16;
if plot_per_cohr 
    results_plot = results_all_perCohr;
    save_name = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar_perCohr.tex'));
    ylim_rolo = 0.015;
     ylim_gremlin_min = 0;
    ylim_gremlin_max = 0.0025;
    ytext_rolo = 0.013;
    ytext_gremlin = 0.0022;
elseif ~plot_per_cohr 
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar.tex'));
    ylim_rolo = 0.015;
    ylim_gremlin_min = 0;
    ylim_gremlin_max = 0.0035;
    ytext_rolo = 0.013;
    ytext_gremlin = 0.003;
end

figure
set(gcf,'unit','inches','position',[0,0,10,4])
ax_1 = subplot(1,2,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.01 0.03 0.04 -0.03]);

plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 14;
plotOptions.plotPercent = false;
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_plot, session_list_plot, plotOptions); 

stats_redundacy_within_cross_string_monkeyR_cardinal = sprintf('Monkey R cardinal, $I_\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_redundacy_within_cross_string_monkeyR_oblique = sprintf('Monkey R oblique, $I_\textrm{redundacy}$, within v.s. cross:  $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);
ylim([0, ylim_rolo])
text(0.5,ytext_rolo,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_rolo,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey R');



ax_2 = subplot(1,2,2); hold on
set(ax_2,'position',get(ax_2,'position')+[0.01 0.03 0.04 -0.03]);
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_plot, session_list_plot, plotOptions); 
stats_redundacy_within_cross_string_monkeyG_cardinal = sprintf('Monkey G cardinal, $I_\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_redundacy_within_cross_string_monkeyG_oblique = sprintf('Monkey G oblique, $I_\textrm{redundacy}$, within v.s. cross:  $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);
ylim([ylim_gremlin_min,ylim_gremlin_max])
text(0.5,ytext_gremlin,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_gremlin,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey G');



print(save_name, '-dsvg');
fid = fopen(tex_name,'wt');
fwrite(fid,[stats_redundacy_within_cross_string_monkeyR_cardinal, stats_redundacy_within_cross_string_monkeyR_oblique,...
            stats_redundacy_within_cross_string_monkeyG_cardinal, stats_redundacy_within_cross_string_monkeyG_oblique]);
fclose(fid);