clear all
clc
close all
%%% This script automatically makes new figure 4. Simulation results with
%%% modified parameters
global   bpGlobal  ftsize
bpGratingFCGlobal();

figure
set(gcf,'Units','inches','Position',[0,0,12,12])
save_folder = '../../../figures/figures_auto_interleavedpaper';
save_name   = fullfile(save_folder,'v2_Figure_4_informationRedundancy_simulation_nonIdeal.svg');
%% 1. effect of delta

b_PF                = 0.00;
cardinal_delta      = 0.08;
oblique_delta       = 0.05;
cardinal_prior      = 1;
oblique_prior       = 1;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

ax_top = subplot(6,3,[2,5]);
ax_bottom = subplot(6,3,8);
set(ax_top,'position',get(ax_top,'position')+[0.02 0.03 0.01 -0.05]);
set(ax_bottom,'position',get(ax_bottom,'position')+[0.02 0.03 0.01 0]);
plot_simulation_session(results_fisher, ax_top, ax_bottom);
axes(ax_top)
title(sprintf('b_{PF} = %d \n \\color{red}Cardinal: \\delta = %.2f, prior = %.1f \n  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold');
%% 2. effect of prior

b_PF                = 0.00;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 0.5;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

ax_top = subplot(6,3,[3,6]);
ax_bottom = subplot(6,3,9);
set(ax_top,'position',get(ax_top,'position')+[0.05 0.03 0.01 -0.05]);
set(ax_bottom,'position',get(ax_bottom,'position')+[0.05 0.03 0.01 0]);
plot_simulation_session(results_fisher, ax_top, ax_bottom);
axes(ax_top)
title(sprintf('b_{PF} = %d \n \\color{red}Cardinal: \\delta = %.2f, prior = %.1f \n  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold');
%% 3. effect of b_PF
b_PF                = 0.8;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 1;
[data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior);

ax_top = subplot(6,3,[10,13]);
ax_bottom = subplot(6,3,16);
set(ax_top,'position',get(ax_top,'position')+[-0.01 -0.03 0.01 -0.05]);
set(ax_bottom,'position',get(ax_bottom,'position')+[-0.01 -0.03 0.01 0]);
plot_simulation_session(results_fisher, ax_top, ax_bottom);
axes(ax_top)
title(sprintf('b_{PF} = %.1f \n \\color{red}Cardinal: \\delta = %.2f, prior = %.1f \n  \\color{blue}Oblique: \\delta = %.2f, prior = %.1f',...
            b_PF, cardinal_delta, cardinal_prior, oblique_delta, oblique_prior),'interpreter','tex','fontweight','bold');
%% 4. effect of b_sample
b_PF                = 0;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 1;

session_name_str    = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior); 
session_name = ['Model_',session_name_str];

filter_name = 'nNeuron_256_bPF_original_0_bPF_sample_2_tao_cardinal_0_tao_oblique_0_nSample_128_random_1000.mat';
fisher_struct = load(sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/%s/individual_sessions_cross/%s',...
        filter_name, session_name ));

results_fisher = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_struct.dat_fisher_cross);
results_fisher = get_sample_CI_cross(results_fisher);

ax_top = subplot(6,3,[11,14]);
ax_bottom = subplot(6,3,17);
set(ax_top,'position',get(ax_top,'position')+[0.02 -0.03 0.01 -0.05]);
set(ax_bottom,'position',get(ax_bottom,'position')+[0.02 -0.03 0.01 0]);
plot_simulation_session(results_fisher, ax_top, ax_bottom);
axes(ax_top)
title({'b_{PF} = 0'; 'b_{sample} = 2'},'interpreter','tex','fontweight','bold');
%% 5. effect of tao_cardinal
b_PF                = 0;
cardinal_delta      = 0.08;
oblique_delta       = 0.08;
cardinal_prior      = 1;
oblique_prior       = 1;

session_name_str    = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior); 
session_name = ['Model_',session_name_str];

filter_name = 'nNeuron_256_bPF_original_0_bPF_sample_0_tao_cardinal_2_tao_oblique_0_nSample_128_random_1000.mat';
fisher_struct = load(sprintf('../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/%s/individual_sessions_cross/%s',...
        filter_name, session_name ));

results_fisher = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_struct.dat_fisher_cross);
results_fisher = get_sample_CI_cross(results_fisher);

ax_top = subplot(6,3,[12,15]);
ax_bottom = subplot(6,3,18);
set(ax_top,'position',get(ax_top,'position')+[0.05 -0.03 0.01 -0.05]);
set(ax_bottom,'position',get(ax_bottom,'position')+[0.05 -0.03 0.01 0]);

plot_simulation_session(results_fisher, ax_top, ax_bottom);
axes(ax_top)
title({'b_{PF} = 0'; '\tau_{cardinal} = 2'},'interpreter','tex','fontweight','bold');


%%
delete(findall(gcf,'type','annotation'))
annotation('textbox',[0.38,0.95,0.3,0.04],'string','b','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.68,0.95,0.3,0.04],'string','c','fontsize',40,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[0.05,0.45,0.3,0.04],'string','d','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.38,0.45,0.3,0.04],'string','e','fontsize',40,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.68,0.45,0.3,0.04],'string','f','fontsize',40,'FontWeight','bold','EdgeColor','none')
%%
print(save_name, '-dsvg');
%%

function plot_simulation_session(results_fisher, ax_top, ax_bottom)
plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.dottest = false;
    plotOptions.plotShuffle = false;
    plotOptions.ftsize = 14;
    plotOptions.plotPercent = true;
    axes(ax_top)
    fig_it.plot_bar_cross_deltaInfo(results_fisher, results_fisher(1).sessionStr, plotOptions); 
    ylim([40, 100])
    
     plotOptions = struct();
    plotOptions.errorbar = 'CI_sample';
    plotOptions.ftsize = 14;
    plotOptions.plotdata = 'delta';
    plotOptions.markersize = 6;
    axes(ax_bottom)
    fig_it.plot_diff_errorbar(results_fisher, results_fisher(1).sessionStr, plotOptions)
    ylim([-5,25])
end

function [data_struct,results_fisher] = get_data(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior)

    session_name_str    = util_it.para_to_namestr(b_PF, cardinal_delta, oblique_delta, cardinal_prior, oblique_prior); 
    
    data_name = ['synthData_use_interleaved_',session_name_str];
    session_name= ['Model_',session_name_str];
    
    data_folder = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/synthData_use_interleaved/real_interleaved/batch_1';
    fisher_folder = '../../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/subset_32_256_random_1000/individual_sessions_cross';
    
    data_struct = load(fullfile(data_folder, data_name));
    fisher_struct = load(fullfile(fisher_folder, session_name));
    
    
    results_fisher = util_it.run_organize_cross_fisherinfo_sizeControl(fisher_struct.dat_fisher_cross);
    results_fisher = get_sample_CI_cross(results_fisher);

end