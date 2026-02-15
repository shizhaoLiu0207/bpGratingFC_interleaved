clear all
clc
close all
%%%% This scripts makes the new version of Figure 3
%%%% First row: 2 bars of information redundancy (same v.s. diff) Across
%%%% two tasks. Three columns are: monkey 1, monkey 2 and model
%%%% Second row: 4 bars, separate each task


global   bpGlobal  ftsize
bpGratingFCGlobal();
%% load empirical fisher information results

filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);



%% load simulated fisher information results
b_PF            = 0;
cardinal_delta  = 0.08;
oblique_delta   = 0.08;
cardinal_prior  = 1;
oblique_prior   = 1;
session_name_str_afterlearning = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);  
data_name_afterlearning =['synthData_use_interleaved_',session_name_str_afterlearning];
session_name_afterlearning = ['Model_',session_name_str_afterlearning];

fisher_folder = '../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';
fisher_afterlearning = load(fullfile(fisher_folder, session_name_afterlearning));
results_fisher_afterlearning = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_afterlearning.dat_fisher_cross);
results_fisher_afterlearning = get_sample_CI_cross(results_fisher_afterlearning);

%%
figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,12,6])
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'v2_Figure_3_informationRedundancy_percent.svg');
tex_name    = fullfile(save_folder ,'v2_Figure_3_informationRedundancy_percent.tex');
%% 1. Monkey R, two bars
ax_1 = subplot(2,3,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = bpGlobal.rolo.session_list.switching;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyR_avg = sprintf(['Monkey R avg, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);


ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex', 'Position', [-3,-3,1]);

title('Monkey R')
ylim([0, 40])
%% 2. Monkey G, two bars

ax_2 = subplot(2,3,2); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyG_avg = sprintf(['Monkey G avg, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);
ylabel('')
title('Monkey G')
ylim([0, 45])
%% 3. Model, two bars
ax_3 = subplot(2,3,3); hold on
set(ax_3,'position',get(ax_3,'position')+[0.04 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
fig_it.plot_bar_cross_deltaInfo(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions); 
ylabel('')
ylim([40, 100])
title('Bayesian model')
%% 4. Monkey R, four bars
ax_4 = subplot(2,3,4); hold on
set(ax_4,'position',get(ax_4,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 4;
session_list_plot = bpGlobal.rolo.session_list.switching;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
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
text(0.3,45,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,45,'Oblique','color','blue','FontSize',14,'FontWeight','bold')

ylabel('')
ylim([0,55])
title('Monkey R')
%% 5. Monkey G, four bars
ax_5 = subplot(2,3,5); hold on
set(ax_5,'position',get(ax_5,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 4;
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);

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

text(0.3,60,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,60,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
ylabel('')
ylim([0,65])
title('Monkey G')
%% 6. Model, four bars

ax_6 = subplot(2,3,6); hold on
set(ax_6,'position',get(ax_6,'position')+[0.04 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 4;
fig_it.plot_bar_cross_deltaInfo(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions); 
ylabel('')
ylim([40, 100])
title('Bayesian model')
text(0.3,90,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,90,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
%% add annotations
delete(findall(gcf,'type','annotation'))

annotation('textbox',[0.04,0.985,0.1,0.04],'string','a','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65,0.985,0.1,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.04,0.5,0.1,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65,0.5,0.1,0.04],'string','d','fontsize',40,'FontWeight','bold','EdgeColor','none')
%% save figure
print(save_name, '-dsvg');
%% save stats

fid = fopen(tex_name,'wt');
fwrite(fid,[ 'Redundancy in percentage: \n'...
            stats_redundacy_within_cross_string_monkeyR_avg, stats_redundacy_within_cross_string_monkeyG_avg,...
            stats_redundacy_within_cross_string_monkeyR_cardinal, stats_redundacy_within_cross_string_monkeyR_oblique,...
            stats_redundacy_within_cross_string_monkeyG_cardinal, stats_redundacy_within_cross_string_monkeyG_oblique]);
fclose(fid);
