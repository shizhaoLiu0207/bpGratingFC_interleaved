clear all
clc
close all

global   bpGlobal  
bpGratingFCGlobal();
%% 1. load simulation results
filter_folder   = '../../results/filtered_neuron_synthetic';

b_PF            = 0.8;
subject_code    = 'rolo';


b_PF_str        = strrep(sprintf('%.2f',b_PF),'.','_');
versionName = sprintf('sampled_subset_empirical_%s_nTotal_128_nSample_64_b_PF_%s_random_250',subject_code, b_PF_str);

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
%% 2. load empirical data
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
save_folder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);

load(fullfile(save_folder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all_empirical = get_sample_CI_cross(results_cross_sizeControl);

switch subject_code
    case 'rolo'
        idx_session = contains({results_all_empirical(:).sessionStr},'Ro');
    case 'gremlin'
        idx_session = contains({results_all_empirical(:).sessionStr},'Gr');
end

diff_redundancy_cardinal_empirical  = cell2mat({results_all_empirical(idx_session).diff_delta_percent_cardinal_median});
diff_redundancy_oblique_empirical   = cell2mat({results_all_empirical(idx_session).diff_delta_percent_oblique_median});

diff_fisher_cardinal_empirical      = cell2mat({results_all_empirical(idx_session).diff_info_percent_cardinal_median});
diff_fisher_oblique_empirical       = cell2mat({results_all_empirical(idx_session).diff_info_percent_oblique_median});

diff_cardinal_oblique_redundancy_empirical = diff_redundancy_cardinal_empirical - diff_redundancy_oblique_empirical;
%% visualize the heat map


feature_target_str_list = {'diff_redundancy_cardinal';'diff_redundancy_oblique'};

feature_title_list = {'\color{red}Diff. Redundancy_{cardinal}';'\color{blue}Diff. Redundancy_{oblique}'};

for k = 1:2
    plotOptions.ftsize = 12;
    plotOptions.feature_title = feature_title_list{k};
    plotOptions.doPeak = true;
    eval(sprintf('plotOptions.f_obs  =  mean(%s_empirical);', feature_target_str_list{k}));
    eval(sprintf('feature_target = %s;',feature_target_str_list{k}));
    figure
    set(gcf,'Units','inches','position',[0,0,12,5])
    fig_it.make_heatmap_4D(feature_target, cardinal_delta,cardinal_prior, oblique_delta,oblique_prior, plotOptions)
    save_folder = '../../figures/figures_final/model_fisher_largescale';
    save_name = fullfile(save_folder,sprintf('heatmap_4D_%s_%s_matched.svg',feature_target_str_list{k}, subject_code));
    print(save_name,'-dsvg')
end

%% visualize the heat map between difference in variables
f = diff_cardinal_oblique_redundancy;
p1 = diff_cardinal_oblique_delta;
p2 = diff_cardinal_oblique_prior;

p1_name = '\color{red}\delta_{cardinal} \color{black}- \color{blue}\delta_{oblique}';
p2_name =  '\color{red}prior_{cardinal} \color{black}- \color{blue}prior_{oblique}';
figure
plotOptions = struct();
plotOptions.ftsize = 12;

plotOptions.doPeak = false;
plotOptions.f_obs = mean(diff_cardinal_oblique_redundancy_empirical);
set(gcf,'Units','inches','position',[0,0,5,4])
fig_it.make_heatmap(p1,p2,f,p1_name,p2_name, plotOptions);

set(gca,'fontsize',14)
title('\color{red}Diff. Redundancy_{cardinal} - \color{blue}Diff. Redundancy_{oblique}','Interpreter','tex',...
    'FontWeight','bold')

save_name = fullfile(save_folder,sprintf('heatmap_diff_diff_redundancy_%s_matched.svg',subject_code));
print(save_name,'-dsvg','-vector')