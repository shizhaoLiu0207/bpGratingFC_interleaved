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
plot_per_cohr = false;
save_folder = '../../figures/figures_final/fisher_info_bar';
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;

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
ftsize = 16;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_Info(results_plot, session_list_plot, plotOptions); 

stats_real_cross_string_monkeyR_cardinal = sprintf('Monkey R cardinal, $I_\textrm{real}$ v.s. $I_\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_real_cross_string_monkeyR_oblique = sprintf('Monkey R oblique, $I_\textrm{real}$ v.s. $I_\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);

ylim([0, ylim_rolo])
text(0.5,ytext_rolo,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_rolo,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey R')

subplot(1,2,2)
%fig.plot_bar_Info(results_all, session_list_gremlin_good,plotOptions); title('Monkey G Good sessions');
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_Info(results_plot, session_list_plot,plotOptions);

stats_real_cross_string_monkeyG_cardinal = sprintf('Monkey G cardinal, $I_\textrm{real}$ v.s. $I_\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_real_cross_string_monkeyG_oblique = sprintf('Monkey G oblique, $I_\textrm{real}$ v.s. $I_\textrm{cross}$: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);

ylim([0,ylim_gremlin])
text(0.5,ytext_gremlin,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_gremlin,'Oblique','color','blue','FontSize',14,'FontWeight','bold')

title('Monkey G');



print(save_name, '-dsvg');

fid = fopen(tex_name,'wt');
fwrite(fid,[stats_real_cross_string_monkeyR_cardinal, stats_real_cross_string_monkeyR_oblique,...
            stats_real_cross_string_monkeyG_cardinal, stats_real_cross_string_monkeyG_oblique]);
fclose(fid);



%% 2. bar plot of I_redundacy, I_redundacy_cross
plot_per_cohr = true;
save_folder = '../../figures/figures_final/fisher_info_bar';

plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotPercent = true;

if plot_per_cohr & ~plotOptions.plotPercent
    results_plot = results_all_perCohr;
    save_name = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar_perCohr.tex'));
    ylim_rolo = 0.015;
    ylim_gremlin = 0.0025;
    ytext_rolo = 0.013;
    ytext_gremlin = 0.0022;
elseif ~plot_per_cohr & ~plotOptions.plotPercent
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_within_cross_info_bar.tex'));
    ylim_rolo = 0.015;
    ylim_gremlin = 0.0035;
    ytext_rolo = 0.013;
    ytext_gremlin = 0.003;
elseif plot_per_cohr & plotOptions.plotPercent
    results_plot = results_all_perCohr;
    save_name = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar_perCohr.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar_perCohr.tex'));
      ylim_rolo = 60;
    ylim_gremlin = 70;
    ytext_rolo = 50;
    ytext_gremlin = 60;
elseif ~plot_per_cohr & plotOptions.plotPercent
    results_plot = results_all;
    save_name = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar.svg'));
    tex_name  = fullfile(save_folder, sprintf('redundacy_percent_within_cross_info_bar.tex'));
     ylim_rolo = 60;
    ylim_gremlin = 70;
    ytext_rolo = 50;
    ytext_gremlin = 60;
end

figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
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


subplot(1,2,2)
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_plot, session_list_plot, plotOptions); 
stats_redundacy_within_cross_string_monkeyG_cardinal = sprintf('Monkey G cardinal, $I_\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value);
stats_redundacy_within_cross_string_monkeyG_oblique = sprintf('Monkey G oblique, $I_\textrm{redundacy}$, within v.s. cross:  $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value);
ylim([-10,ylim_gremlin])
text(0.5,ytext_gremlin,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.5,ytext_gremlin,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Monkey G');


print(save_name, '-dsvg');
fid = fopen(tex_name,'wt');
fwrite(fid,[stats_redundacy_within_cross_string_monkeyR_cardinal, stats_redundacy_within_cross_string_monkeyR_oblique,...
            stats_redundacy_within_cross_string_monkeyG_cardinal, stats_redundacy_within_cross_string_monkeyG_oblique]);
fclose(fid);

