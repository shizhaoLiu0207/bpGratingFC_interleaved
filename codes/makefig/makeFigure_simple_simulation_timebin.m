clear all
clc
close all
%%%%% Make figure of results of a simplied model simulation
global  bpGlobal
bpGratingFCGlobal();
color_C = bpGlobal.color_list.color_cardinal; color_O = bpGlobal.color_list.color_oblique;

%%
save_name = '../../results/neural/organized_simple_simulation_timebin';
doOrganize = 0;
if doOrganize
    b_PF            = 0;
    cardinal_delta  = 0.08;
    oblique_delta   = 0.08;
    cardinal_prior  = 1;
    oblique_prior   = 1;
    session_name_str_afterlearning = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);  
    
    cardinal_delta  = 0;
    oblique_delta   = 0;
    
    session_name_str_beforelearning = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);  
    
    data_name_afterlearning =['synthData_use_interleaved_',session_name_str_afterlearning];
    session_name_afterlearning = ['Model_',session_name_str_afterlearning(1:end-4),'_timebin'];
    
    data_name_beforelearning =['synthData_use_interleaved_',session_name_str_beforelearning];
    session_name_beforelearning = ['Model_',session_name_str_beforelearning(1:end-4),'_timebin'];
    
    fisher_folder = '../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';
    
    
    fisher_afterlearning = load(fullfile(fisher_folder, session_name_afterlearning));
    fisher_beforelearning = load(fullfile(fisher_folder, session_name_beforelearning));
    
    % results_fisher_afterlearning = organize_sample_fisher(fisher_afterlearning.dat_fisher_cross);
    % results_fisher_beforelearning = organize_sample_fisher(fisher_beforelearning.dat_fisher_cross);
    
    results_fisher_afterlearning = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_afterlearning.dat_fisher_cross);
    results_fisher_beforelearning = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_beforelearning.dat_fisher_cross);
    
    
    results_fisher_afterlearning = get_sample_CI_cross(results_fisher_afterlearning);
    results_fisher_beforelearning = get_sample_CI_cross(results_fisher_beforelearning);
    
    save(save_name,'results_fisher_afterlearning','results_fisher_beforelearning');
else
    load(save_name)
end
%% within trial Diff(I_cross, I_real)

figure
set(gcf,'Units','inches','Position',[0,0,12,6])
ymin = -30;
ymax = 30;
subplot(2,2,1)
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'info';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time bin index'; 
fig_it.plot_diff_withintrial(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title({'Before learning';'\color{red}Cardinal'},'FontSize',16)

subplot(2,2,3)
plotOptions.task    = 'oblique';
fig_it.plot_diff_withintrial(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title('\color{blue}Oblique','FontSize',16)

subplot(2,2,2)
plotOptions.task    = 'cardinal';
fig_it.plot_diff_withintrial(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title({'After learning';'\color{red}Cardinal'},'FontSize',16)

subplot(2,2,4)
plotOptions.task    = 'oblique';
fig_it.plot_diff_withintrial(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title('\color{blue}Oblique','FontSize',16)

save_folder = '../../figures/figures_final/model_fisher_ideal';
save_name = fullfile(save_folder,'within_trial_diff_fisher.svg');
print(save_name,'-dsvg','-vector')
%%


figure
set(gcf,'Units','inches','Position',[0,0,12,6])
ymin = -6;
ymax = 6;
subplot(2,2,1)
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.ftsize = 14;
plotOptions.plotdata = 'delta';
plotOptions.markersize = 6;
plotOptions.task    = 'cardinal';
plotOptions.xlabelStr = 'Time bin index'; 
fig_it.plot_diff_withintrial(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title({'Before learning';'\color{red}Cardinal'},'FontSize',16)

subplot(2,2,3)
plotOptions.task    = 'oblique';
fig_it.plot_diff_withintrial(results_fisher_beforelearning, results_fisher_beforelearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title('\color{blue}Oblique','FontSize',16)

subplot(2,2,2)
plotOptions.task    = 'cardinal';
fig_it.plot_diff_withintrial(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title({'After learning';'\color{red}Cardinal'},'FontSize',16)

subplot(2,2,4)
plotOptions.task    = 'oblique';
fig_it.plot_diff_withintrial(results_fisher_afterlearning, results_fisher_afterlearning(1).sessionStr, plotOptions)
ylim([ymin, ymax])
title('\color{blue}Oblique','FontSize',16)

save_folder = '../../figures/figures_final/model_fisher_ideal';
save_name = fullfile(save_folder,'within_trial_diff_redundancy.svg');
print(save_name,'-dsvg','-vector')