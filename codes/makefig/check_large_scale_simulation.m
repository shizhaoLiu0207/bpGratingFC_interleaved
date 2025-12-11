clear all
clc
close all

global   bpGlobal  
bpGratingFCGlobal();
%% 1. load simulation results
filter_folder   = '../../results/filtered_neuron_synthetic';

%%%%% check simulation features using a uniform distribution
versionName = 'sampled_subset_uniform_nTotal_128_nSample_64_b_PF_0_00_random_250';

save_folder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/%s',versionName);

load(fullfile(save_folder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData'));

results_all_simulation = get_sample_CI_cross(results_cross_sizeControl);

nSession = numel(results_all_simulation);
[cardinal_delta, cardinal_prior, oblique_delta, oblique_prior] = deal(zeros(nSession,1));
%%%% key measure
%%% (I_cross - I_real) / I_cross; cardinal and oblique
%%%% I_{redundancy, within, pct} - I_{redundancy, cross, pct}; cardinal and
%%%% oblique
[diff_fisher_cardinal, diff_fisher_oblique, diff_redundancy_cardinal, diff_redundancy_oblique] = deal(zeros(nSession,1));
[redundancy_cardinal_within, redundancy_cardinal_cross, redundancy_oblique_within, redundancy_oblique_cross] =  deal(zeros(nSession,1));
for n = 1:nSession
    sessionStr = results_all_simulation(n).sessionStr;
    tokens = regexp(sessionStr, 'Model_bPF_([\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)','tokens');
    
    extracted_params        = tokens{1}; % Extract matched tokens
    cardinal_delta_str      = strrep(extracted_params{2}, '_', '.'); 
    cardinal_prior_str      = strrep(extracted_params{3}, '_', '.'); 
    oblique_delta_str       = strrep(extracted_params{4}, '_', '.'); 
    oblique_prior_str       = strrep(extracted_params{5}, '_', '.'); 

    cardinal_delta(n)       = str2double(cardinal_delta_str);
    cardinal_prior(n)       = str2double(cardinal_prior_str);
    oblique_delta(n)        = str2double(oblique_delta_str);
    oblique_prior(n)        = str2double(oblique_prior_str);
    

    

    redundancy_cardinal_within(n)   = results_all_simulation(n).delta_percent_cardinal_cardinal_median;
    redundancy_cardinal_cross(n)    = results_all_simulation(n).delta_percent_cardinal_oblique_median;
    redundancy_oblique_within(n)    = results_all_simulation(n).delta_percent_oblique_oblique_median;
    redundancy_oblique_cross(n)     = results_all_simulation(n).delta_percent_oblique_cardinal_median;

    diff_fisher_cardinal(n) = results_all_simulation(n).diff_info_percent_cardinal_median;
    diff_fisher_oblique(n)  = results_all_simulation(n).diff_info_percent_oblique_median;

    diff_redundancy_cardinal(n)     = results_all_simulation(n).diff_delta_percent_cardinal_median;
    diff_redundancy_oblique(n)      = results_all_simulation(n).diff_delta_percent_oblique_median;
end

diff_cardinal_oblique_fisher = diff_fisher_cardinal - diff_fisher_oblique;
diff_cardinal_oblique_redundancy = diff_redundancy_cardinal - diff_redundancy_oblique;

diff_cardinal_oblique_delta = cardinal_delta - oblique_delta;
diff_cardinal_oblique_prior = cardinal_prior - oblique_prior;





%% visualize the heat map
feature_target_str_list = {'diff_fisher_cardinal';'diff_redundancy_cardinal';...
                            'diff_fisher_oblique';'diff_redundancy_oblique'};
feature_title_list = { '\color{red}Diff. Information_{cardinal}';'\color{red}Diff. Redundancy_{cardinal}';...
                        '\color{blue}Diff. Information_{oblique}';'\color{blue}Diff. Redundancy_{oblique}'};
for k = 1:4
    plotOptions.ftsize = 14;
    plotOptions.feature_title = feature_title_list{k};
    eval(sprintf('feature_target = %s;',feature_target_str_list{k}));
    figure
    set(gcf,'Units','inches','position',[0,0,12,5])
    fig_it.make_heatmap_4D(feature_target, cardinal_delta,cardinal_prior, oblique_delta,oblique_prior, plotOptions)
    save_folder = '../../figures/figures_final/model_fisher_largescale';
    save_name = fullfile(save_folder,sprintf('heatmap_4D_%s.svg',feature_target_str_list{k}));
    print(save_name,'-dsvg')
end
%% visualize the heat map between difference in variables
f = diff_cardinal_oblique_redundancy;
p1 = diff_cardinal_oblique_delta;
p2 = diff_cardinal_oblique_prior;

p1_name = '\color{red}\delta_{cardinal} \color{black}- \color{blue}\delta_{oblique}';
p2_name =  '\color{red}prior_{cardinal} \color{black}- \color{blue}prior_{oblique}';
figure
set(gcf,'Units','inches','position',[0,0,5,4])
fig_it.make_heatmap(p1,p2,f,p1_name,p2_name);
set(gca,'fontsize',14)
title('\color{red}Diff. Redundancy_{cardinal} - \color{blue}Diff. Redundancy_{oblique}','Interpreter','tex',...
    'FontWeight','bold')

save_name = fullfile(save_folder,'heatmap_diff_diff_redundancy.svg');
print(save_name,'-dsvg','-vector')



f = diff_cardinal_oblique_fisher;
p1 = diff_cardinal_oblique_delta;
p2 = diff_cardinal_oblique_prior;

p1_name = '\color{red}\delta_{cardinal} \color{black}- \color{blue}\delta_{oblique}';
p2_name =  '\color{red}prior_{cardinal} \color{black}- \color{blue}prior_{oblique}';
figure
set(gcf,'Units','inches','position',[0,0,5,4])
fig_it.make_heatmap(p1,p2,f,p1_name,p2_name);
set(gca,'fontsize',14)
title('\color{red}Diff. Information_{cardinal} - \color{blue}Diff. Information_{oblique}','Interpreter','tex',...
    'FontWeight','bold')

save_name = fullfile(save_folder,'heatmap_diff_diff_information.svg');
print(save_name,'-dsvg')


%%

