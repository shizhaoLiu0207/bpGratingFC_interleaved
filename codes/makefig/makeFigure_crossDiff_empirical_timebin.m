clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
%filter_name = 'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl';
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_timebin'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
%%
save_folder = '../../figures/figures_final/';
save_name = fullfile(save_folder, 'within_trial_empirical_monkeyG_cardinal.svg');
tex_name  = fullfile(save_folder, 'within_trial_empirical_monkeyG_cardinal.tex');
figure
set(gcf,'Units','inches','position',[0,0,12,8])
%%

ax_1 = subplot(2,2,1);
set(ax_1,'position',get(ax_1,'position')+[-0.02 0.04 0.03 -0.03 ]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time within trial (ms)';
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
stats_info = fig_it.plot_info_withintrial(results_all, session_list_plot, plotOptions);
title('')
add_two_stats_annotations(stats_info,ax_1, bpGlobal.color_list.color_cardinal);

monkeyR_cardinal_fisher_text = sprintf(['Monket G cardinal, $\\Ir$: $\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n' ...
    ' $\\Icross$: $\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n'], ...
    stats_info.beta_within, stats_info.beta_p_within, stats_info.beta_cross, stats_info.beta_p_cross);
%%

ax_2 = subplot(2,2,2);
set(ax_2,'position',get(ax_2,'position')+[0.02 0.04 0.03 -0.03 ]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time within trial (ms)';
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
stats_info = fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions);
title('')
add_stats_annotations(stats_info, ax_2, bpGlobal.color_list.color_cardinal);

monkeyR_cardinal_info_diff_text = sprintf(['Monket G cardinal, difference between $\\Ir$ and $\\Icross$:' ...
    '$\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n'], ...
    stats_info.beta, stats_info.beta_p);

%%

ax_3 = subplot(2,2,3);
set(ax_3,'position',get(ax_3,'position')+[-0.02 0.04 0.03 -0.03 ]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta_percent';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time within trial (ms)';
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
stats_info = fig_it.plot_info_withintrial(results_all, session_list_plot, plotOptions);
title('')
add_two_stats_annotations(stats_info,ax_3, bpGlobal.color_list.color_cardinal);
monkeyR_cardinal_redundancy_text = sprintf(['Monket G cardinal, $\\Ideltawithin$: $\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n' ...
    ' $\\Ideltacross$: $\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n'], ...
    stats_info.beta_within, stats_info.beta_p_within, stats_info.beta_cross, stats_info.beta_p_cross);
%%
ax_4 = subplot(2,2,4);
set(ax_4,'position',get(ax_4,'position')+[0.02 0.04 0.03 -0.03 ]);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time within trial (ms)';
session_list_plot = bpGlobal.gremlin.session_list.interleaved_real;
stats_info = fig_it.plot_diff_withintrial(results_all, session_list_plot, plotOptions);
title('')
add_stats_annotations(stats_info, ax_4, bpGlobal.color_list.color_cardinal);

monkeyR_cardinal_delta_diff_text = sprintf(['Monket G cardinal, difference between $\\Ideltawithin$ and $\\Ideltacross$:' ...
    '$\\beta = \\num{%.4e}, p = \\num{%.4e}$ \n'], ...
    stats_info.beta, stats_info.beta_p);

%% add annotations
annotation('textbox',[0.02,0.97,0.15,0.04],'string','a','fontsize',36,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.50,0.97,0.15,0.04],'string','b','fontsize',36,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.02,0.47,0.15,0.04],'string','c','fontsize',36,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.50,0.47,0.15,0.04],'string','d','fontsize',36,'FontWeight','bold','EdgeColor','none')
%%
print(save_name,'-dsvg','-vector')



fid = fopen(tex_name,'wt');

fwrite(fid, [monkeyR_cardinal_fisher_text, monkeyR_cardinal_info_diff_text,  ...
    monkeyR_cardinal_redundancy_text, monkeyR_cardinal_delta_diff_text]);

fclose(fid)

%% helper functions

function add_stats_annotations(stats_info, ax_plot, plotColor)
    global ftsize
    late_stats_string = sprintf('Diff: $\\beta = %.2e^{%s}$', stats_info.beta,  fig.p2star(stats_info.beta_p));
   
    ax_pos = get(ax_plot,'position');
    annotation('line',[ax_pos(1)+0.02, ax_pos(1) + 0.05],[ax_pos(2)-0.11, ax_pos(2)-0.11],'LineWidth',2,'Color',plotColor)
    
    annotation('textbox',[ax_pos(1) + 0.06,ax_pos(2)-0.13,0.15,0.04],'string',late_stats_string,'fontsize',ftsize,'EdgeColor','none','Interpreter','latex')

end


function add_two_stats_annotations(stats_info, ax_plot, plotColor)
    global ftsize
    within_stats_string = sprintf('Within: $\\beta = %.2e^{%s}$', stats_info.beta_within,  fig.p2star(stats_info.beta_p_within));
    cross_stats_string = sprintf('Cross: $\\beta = %.2e^{%s}$', stats_info.beta_cross,  fig.p2star(stats_info.beta_p_cross));
    
    ax_pos = get(ax_plot,'position');

    annotation('line',[ax_pos(1)+0.01, ax_pos(1) + 0.04],[ax_pos(2)-0.11, ax_pos(2)-0.11],'LineWidth',2,'Color',plotColor)
    annotation('textbox',[ax_pos(1) + 0.04,ax_pos(2)-0.13,0.15,0.04],'string',within_stats_string,'fontsize',ftsize,'EdgeColor','none','Interpreter','latex')

    annotation('line',[ax_pos(1)+0.19, ax_pos(1) + 0.22],[ax_pos(2)-0.11, ax_pos(2)-0.11],'LineWidth',2,'Color',plotColor,'LineStyle','--')
    annotation('textbox',[ax_pos(1) + 0.23,ax_pos(2)-0.13,0.15,0.04],'string',cross_stats_string,'fontsize',ftsize,'EdgeColor','none','Interpreter','latex')

end