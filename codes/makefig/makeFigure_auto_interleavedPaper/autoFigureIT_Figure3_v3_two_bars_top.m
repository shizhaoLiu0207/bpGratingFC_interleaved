clear all
clc
%close all
%%

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
set(gcf,'Units','inches','position',[0,0,12,5])
save_folder = '../../../figures/interleaved_paperFigures_V3';
save_name   = fullfile(save_folder,'v3_Figure_3_informationRedundancy_percent_twobars_top.svg');
tex_name    = fullfile(save_folder ,'v3_Figure_3_informationRedundancy_percent_twobars_top.tex');
%% 1. Monkey R, two bars
ax_1 = subplot(1,3,1); hold on
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


ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex');

title('Monkey R')
ylim([0, 40])
%% 2. Monkey G, two bars

ax_2 = subplot(1,3,2); hold on
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
ax_3 = subplot(1,3,3); hold on
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
%% add annotations
delete(findall(gcf,'type','annotation'))

annotation('textbox',[0.04,0.985,0.1,0.04],'string','a','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.65,0.985,0.1,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')

%annotation('textbox',[0.04,0.5,0.1,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')
%% save figure
print(save_name, '-dsvg');
%% save stats

fid = fopen(tex_name,'wt');
fwrite(fid,[ 'Redundancy in percentage: \n'...
            stats_redundacy_within_cross_string_monkeyR_avg, stats_redundacy_within_cross_string_monkeyG_avg]);
fclose(fid);
