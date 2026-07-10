clear all
clc
close all
%%%% This scripts makes the new version of Figure 3
%%%% First row: 2 bars of information redundancy (same v.s. diff) Across
%%%% two tasks. Three columns are: monkey 1, monkey 2 and model
%%%% Second row: 4 bars, separate each task


global   bpGlobal  ftsize
bpGratingFCGlobal();
save_folder         = '../../figures/figures_final/fisher_info_bar' ;
plotOptions.group   = 'high_low'; %'high_low', 'consistent_opposite','two_opposite'
plotOptions.source  = 'passiveViewing'; % 'mainTask', 'passiveViewing';
plotOptions.dosize_control = true;
switch plotOptions.group
    case 'high_low'
        filter_name_top = sprintf('all_trials_highC_highO_%s', plotOptions.source);
        filter_name_bottom = sprintf('all_trials_lowC_lowO_%s', plotOptions.source);

        if plotOptions.dosize_control
            filter_name_top = [filter_name_top,'_sizeControl'];
            filter_name_bottom = [filter_name_bottom, '_sizeControl'];
        end
        if  plotOptions.dosize_control
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_high_low_%s_sizeControl_averaged.svg', plotOptions.source));
        else
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_high_low_%s_averaged.svg', plotOptions.source));
        end
        
        titleStr_top = 'highC-highO';
        titleStr_bottom = 'lowC-lowO';

    case 'consistent_opposite'
        filter_name_top = sprintf('all_trials_highC_highO_lowC_lowO_%s', plotOptions.source);
        filter_name_bottom = sprintf('all_trials_highC_lowO_lowC_highO_%s', plotOptions.source);
        if plotOptions.dosize_control
            filter_name_top = [filter_name_top,'_sizeControl'];
            filter_name_bottom = [filter_name_bottom, '_sizeControl'];
        end
        if plotOptions.dosize_control
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_consistent_opposite_%s_sizeControl_averaged.svg', plotOptions.source));
        else
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_consistent_opposite_%s_averaged.svg', plotOptions.source));
        end

        titleStr_top = 'consistent';
        titleStr_bottom = 'opposite';
    case 'two_opposite'
        filter_name_top = sprintf('all_trials_highC_lowO_%s', plotOptions.source);
        filter_name_bottom = sprintf('all_trials_lowC_highO_%s', plotOptions.source);

        if plotOptions.dosize_control
            filter_name_top = [filter_name_top,'_sizeControl'];
            filter_name_bottom = [filter_name_bottom, '_sizeControl'];
        end
        if  plotOptions.dosize_control
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_two_opposite_%s_sizeControl_averaged.svg', plotOptions.source));
        else
            save_name = fullfile(save_folder, sprintf('informationRedundancy_percent_two_opposite_%s_averaged.svg', plotOptions.source));
        end
        
        titleStr_top = 'highC-lowO';
        titleStr_bottom = 'lowC-highO';
end

plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 2;
%% 
figure
set(gcf,'unit','inches','position',[0,0,10,8]);

%% load empirical fisher information results

filter_name = filter_name_top;
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);

%%% 1. Monkey R, four bars
ax_1 = subplot(2,2,1); hold on
set(ax_1,'position',get(ax_4,'position')+[-0.02 0.04 0.03 -0.03]);

session_list_plot = bpGlobal.rolo.session_list.switching;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
% stats_redundacy_within_cross_string_monkeyR_cardinal = sprintf(['Monkey R cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
%                                         stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
%                                         stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);
% 
% stats_redundacy_within_cross_string_monkeyR_oblique = sprintf(['Monkey R oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
%                                         stats_info_oblique.mu_within, stats_info_oblique.std_within,...
%                                         stats_info_oblique.mu_cross, stats_info_oblique.std_cross);


ylabel('')
ylim([0,40])
title(sprintf('Monkey R, %s',titleStr_top))
%% 5. Monkey G, four bars
ax_2 = subplot(2,2,2); hold on
%set(ax_5,'position',get(ax_5,'position')+[-0.02 0.04 0.03 -0.03]);
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);

% stats_redundacy_within_cross_string_monkeyG_cardinal = sprintf(['Monkey G cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
%                                         stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
%                                         stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);
% 
% 
% stats_redundacy_within_cross_string_monkeyG_oblique = sprintf(['Monkey G oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
%                                         stats_info_oblique.mu_within, stats_info_oblique.std_within,...
%                                         stats_info_oblique.mu_cross, stats_info_oblique.std_cross);


ylabel('')
ylim([-10,40])
title(sprintf('Monkey G, %s',titleStr_top))


%%
filter_name = filter_name_bottom;
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);

%%% 1. Monkey R, four bars
ax_1 = subplot(2,2,3); hold on
%set(ax_4,'position',get(ax_4,'position')+[-0.02 0.04 0.03 -0.03]);

session_list_plot = bpGlobal.rolo.session_list.switching;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);
% stats_redundacy_within_cross_string_monkeyR_cardinal = sprintf(['Monkey R cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
%                                         stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
%                                         stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);
% 
% stats_redundacy_within_cross_string_monkeyR_oblique = sprintf(['Monkey R oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
%                                         stats_info_oblique.mu_within, stats_info_oblique.std_within,...
%                                         stats_info_oblique.mu_cross, stats_info_oblique.std_cross);


ylabel('')
ylim([0,300])
title(sprintf('Monkey R, %s',titleStr_bottom))
%% 5. Monkey G, four bars
ax_2 = subplot(2,2,4); hold on
%set(ax_5,'position',get(ax_5,'position')+[-0.02 0.04 0.03 -0.03]);

session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;

[stats_info_cardinal, stats_info_oblique]  = fig_it.plot_bar_cross_deltaInfo(results_all, session_list_plot, plotOptions);

% stats_redundacy_within_cross_string_monkeyG_cardinal = sprintf(['Monkey G cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
%                                         stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
%                                         stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);
% 
% 
% stats_redundacy_within_cross_string_monkeyG_oblique = sprintf(['Monkey G oblique, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
%                                         'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
%                                         'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
%                                         stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
%                                         stats_info_oblique.mu_within, stats_info_oblique.std_within,...
%                                         stats_info_oblique.mu_cross, stats_info_oblique.std_cross);

ylabel('')
ylim([0,65])
title(sprintf('Monkey G, %s',titleStr_bottom))


%%
print(save_name,'-dsvg');
