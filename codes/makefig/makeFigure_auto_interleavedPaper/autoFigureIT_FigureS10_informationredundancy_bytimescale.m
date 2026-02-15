clear all
clc
close all
global   bpGlobal  ftsize
bpGratingFCGlobal();

rolo_block_size = bpGlobal.rolo.session_list.mini_block;
gremlin_block_size = bpGlobal.gremlin.session_list.interleaved_blockSize;

rolo_mini_block_list = rolo_block_size(cell2mat(rolo_block_size(:,2)) < 200, 1);

rolo_random_list = setdiff(bpGlobal.rolo.session_list.switching, rolo_block_size(:,1));

gremlin_block_list = gremlin_block_size(cell2mat(gremlin_block_size(:,2)) > 0,1);
gremlin_random_list = setdiff(bpGlobal.gremlin.session_list.interleaved_real, gremlin_block_list);

%% load empirical results
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
%%
figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,12,6])
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'v2_Figure_S10_informationRedundancy_percent_bytimescale.svg');
tex_name    = fullfile(save_folder ,'v2_Figure_S10_informationRedundancy_percent_bytimescale.tex');

%% 1. Monkey R, mini block
ax_1 = subplot(2,2,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = rolo_mini_block_list;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyR_block_avg = sprintf(['Monkey R mini-block, avg, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);


ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex', 'Position', [-2,0,1]);

title('Monkey R, mini-block')
ylim([0, 40])
%% 2. Monkey R, trial-by-trial
ax_2 = subplot(2,2,2); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = rolo_random_list;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyR_trial_avg = sprintf(['Monkey R avg trial-by-trial, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);
ylabel('')
title('Monkey R, trial-by-trial')
ylim([0, 45])
%% 3. Monkey G,  block
ax_3 = subplot(2,2,3); hold on
set(ax_3,'position',get(ax_3,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = gremlin_block_list;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyG_block_avg = sprintf(['Monkey G mini-block, avg, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);


ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex', 'Position', [-3,-3,1]);

title('Monkey G, blockwise')
ylim([0, 50])
%% 3. Monkey G, trial-by-trial
ax_4 = subplot(2,2,4); hold on
set(ax_4,'position',get(ax_4,'position')+[-0.02 0.04 0.03 -0.03]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
session_list_plot = gremlin_random_list;

[~, ~, stats_info_avg]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
stats_redundacy_within_cross_string_monkeyG_trial_avg = sprintf(['Monkey G avg trial-by-trial, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_avg.df, stats_info_avg.tstat, stats_info_avg.p_value,...
                                        stats_info_avg.mu_within, stats_info_avg.std_within,...
                                        stats_info_avg.mu_cross, stats_info_avg.std_cross);
ylabel('')
title('Monkey G, trial-by-trial')
ylim([0, 50])
%% direct compare the neural signature, block vs trial-by-trial

session_list_plot = rolo_mini_block_list;
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);
diff_delta_monkeyR_block = [results_all(idx).diff_delta_percent_median];


session_list_plot = rolo_random_list;
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);
diff_delta_monkeyR_trial = [results_all(idx).diff_delta_percent_median];


session_list_plot = gremlin_block_list;
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);
diff_delta_monkeyG_block = [results_all(idx).diff_delta_percent_median];


session_list_plot = gremlin_random_list;
idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
    ismember({results_all(:).sessionStr}, session_list_plot);
diff_delta_monkeyG_trial = [results_all(idx).diff_delta_percent_median];

[~, p_val, ~, stats] = ttest2(diff_delta_monkeyR_block, diff_delta_monkeyR_trial);
stats_compare_timescale_monkeyR = sprintf('Monkey R, compare neural signature across time scales: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats.df, stats.tstat, p_val);

[~, p_val, ~, stats] = ttest2(diff_delta_monkeyG_block, diff_delta_monkeyG_trial);
stats_compare_timescale_monkeyG = sprintf('Monkey G, compare neural signature across time scales: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        stats.df, stats.tstat, p_val);
%% add annotations
delete(findall(gcf,'type','annotation'))

annotation('textbox',[0.05,0.985,0.1,0.04],'string','a','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.5,0.985,0.1,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.05,0.5,0.1,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.5,0.5,0.1,0.04],'string','d','fontsize',40,'FontWeight','bold','EdgeColor','none')

%% save figure
print(save_name, '-dsvg')

%% save stats
fid = fopen(tex_name,'wt');

fwrite(fid, ['Redundancy in percentage: \n', ...
            stats_redundacy_within_cross_string_monkeyR_block_avg,...
            stats_redundacy_within_cross_string_monkeyR_trial_avg,...
            stats_redundacy_within_cross_string_monkeyG_block_avg, ...
            stats_redundacy_within_cross_string_monkeyG_trial_avg, ...
            stats_compare_timescale_monkeyR, stats_compare_timescale_monkeyG]);

fclose(fid);
