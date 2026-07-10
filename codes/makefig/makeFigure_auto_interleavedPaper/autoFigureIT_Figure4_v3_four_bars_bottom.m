clear all
clc
close all
%%%% This scripts makes the new version of Figure 3
%%%% First row: 2 bars of information redundancy (same v.s. diff) Across
%%%% two tasks. Three columns are: monkey 1, monkey 2 and model
%%%% Second row: 4 bars, separate each task


global   bpGlobal  ftsize
bpGratingFCGlobal();



filter_name_high = 'all_trials_highC_highO_passiveViewing_sizeControl';
filter_name_low = 'all_trials_lowC_lowO_passiveViewing_sizeControl';




plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
plotOptions.ftsize = 16;
plotOptions.plotPercent = true;
plotOptions.numBar = 4;

%% load data
%%%%% high sensitivity data
saveFolder = sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name_high);
load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all_high = results_cross_sizeControl;
results_all_high = get_sample_CI_cross(results_all_high);

%%%%% low sensitivity data
saveFolder = sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name_low);
load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all_low = results_cross_sizeControl;
results_all_low = get_sample_CI_cross(results_all_low);
%% 
figure
set(gcf,'unit','inches','position',[0,0,12,5]);
save_folder = '../../../figures/interleaved_paperFigures_V3';
save_name   = fullfile(save_folder,'v3_Figure_4_informationRedundancy_percent_fourbars_bottom.svg');
tex_name    = fullfile(save_folder ,'v3_Figure_4_informationRedundancy_percent_fourbars_bottom.tex');
%% 1. Monkey R, high sensitivity
ax_1 = subplot(1,4,1); hold on
set(ax_1,'position',get(ax_1,'position')+[-0.05 0.04 0.02 -0.15]);
session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_all_high, session_list_plot, plotOptions);

stats_redundacy_within_cross_string_monkeyR_cardinal_high = sprintf(['Monkey R, high sensitivity units, cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);

stats_redundacy_within_cross_string_monkeyR_oblique_high = sprintf(['Monkey R, high sensitivity units, oblique,  $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);

ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex');
ylim([0,35])

text(0.3,30,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,30,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('High senstivity')
%% 2. Monkey R, low sensitivity
ax_2 = subplot(1,4,2); hold on
set(ax_2,'position',get(ax_2,'position')+[-0.03 0.04 0.02 -0.15]);

session_list_plot = bpGlobal.rolo.session_list.switching;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_all_low, session_list_plot, plotOptions);


stats_redundacy_within_cross_string_monkeyR_cardinal_low = sprintf(['Monkey R, low sensitivity units, cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);

stats_redundacy_within_cross_string_monkeyR_oblique_low = sprintf(['Monkey R, low sensitivity units, oblique,  $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);

ylabel('');
ylim([0,300])
text(0.3,275,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,275,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Low senstivity')
%% 3. Monkey G, high sensitivity
ax_3 = subplot(1,4,3); hold on
set(ax_3,'position',get(ax_3,'position')+[0.01 0.04 0.02 -0.15]);

session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_all_high, session_list_plot, plotOptions);


stats_redundacy_within_cross_string_monkeyG_cardinal_high = sprintf(['Monkey G, high sensitivity units, cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);

stats_redundacy_within_cross_string_monkeyG_oblique_high = sprintf(['Monkey G, high sensitivity units, oblique,  $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);

ylabel('$I_\textrm{redundancy}$ (Percent)','FontSize',18,'Interpreter','Latex');
ylim([-10,40])
text(0.3,35,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,35,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('High senstivity')
%% 4. Monkey G, low sensitivity

ax_4 = subplot(1,4,4); hold on
set(ax_4,'position',get(ax_4,'position')+[0.05 0.04 0.02 -0.15]);

session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
[stats_info_cardinal, stats_info_oblique] = fig_it.plot_bar_cross_deltaInfo(results_all_low, session_list_plot, plotOptions);


stats_redundacy_within_cross_string_monkeyG_cardinal_low = sprintf(['Monkey G, low sensitivity units, cardinal, $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_cardinal.df, stats_info_cardinal.tstat, stats_info_cardinal.p_value,...
                                        stats_info_cardinal.mu_within, stats_info_cardinal.std_within,...
                                        stats_info_cardinal.mu_cross, stats_info_cardinal.std_cross);

stats_redundacy_within_cross_string_monkeyG_oblique_low = sprintf(['Monkey G, low sensitivity units, oblique,  $I_\\textrm{redundacy}$, within v.s. cross: $t(%d) = %.2f$, $p = \\num{%.2e}$ \n',...
                                        'Within-redundancy: $Mean = %.2f$, $std = %.2f$' ...
                                        'Cross-redundancy: $Mean = %.2f$, $std = %.2f$ \n'], ...
                                        stats_info_oblique.df, stats_info_oblique.tstat, stats_info_oblique.p_value,...
                                        stats_info_oblique.mu_within, stats_info_oblique.std_within,...
                                        stats_info_oblique.mu_cross, stats_info_oblique.std_cross);

ylabel('');
ylim([0,70])
text(0.3,65,'Cardinal','color','red','FontSize',14,'FontWeight','bold')
text(5.3,65,'Oblique','color','blue','FontSize',14,'FontWeight','bold')
title('Low senstivity')
%% add annotations
delete(findall(gcf,'type','annotation'))

annotation('textbox',[0.23,0.95,0.20,0.04],'string','Monkey R','fontsize',18,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.71,0.95,0.20,0.04],'string','Monkey G','fontsize',18,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.01,0.88,0.1,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.49,0.88,0.1,0.04],'string','d','fontsize',40,'FontWeight','bold','EdgeColor','none')


%% save figure
print(save_name, '-dsvg');
%% save stats

fid = fopen(tex_name,'wt');
fwrite(fid,[ 'Redundancy in percentage: \n'...
            stats_redundacy_within_cross_string_monkeyR_cardinal_high, stats_redundacy_within_cross_string_monkeyR_oblique_high,...
            stats_redundacy_within_cross_string_monkeyR_cardinal_low, stats_redundacy_within_cross_string_monkeyR_oblique_low,...
            stats_redundacy_within_cross_string_monkeyG_cardinal_high, stats_redundacy_within_cross_string_monkeyG_cardinal_high,...
            stats_redundacy_within_cross_string_monkeyG_cardinal_low, stats_redundacy_within_cross_string_monkeyG_oblique_low]);
fclose(fid);