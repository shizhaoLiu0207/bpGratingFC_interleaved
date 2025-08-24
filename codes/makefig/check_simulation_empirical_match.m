clear all
clc
close all

global   bpGlobal  
bpGratingFCGlobal();
%% 1. load simulation results
filter_folder   = '../../results/filtered_neuron_synthetic';

b_PF            = 0.8;
subject_code    = 'gremlin';


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

diff_redundancy_cardinal_empirical  = cell2mat({results_all_empirical(idx_session).diff_info_percent_cardinal_median});
diff_redundancy_oblique_empirical   = cell2mat({results_all_empirical(idx_session).diff_info_percent_oblique_median});

diff_fisher_cardinal_empirical      = cell2mat({results_all_empirical(idx_session).diff_info_percent_cardinal_median});
diff_fisher_oblique_empirical       = cell2mat({results_all_empirical(idx_session).diff_info_percent_oblique_median});

%% match between simulation and empirical data - average


dist_feature = sqrt((diff_redundancy_cardinal - mean(diff_redundancy_cardinal_empirical)).^ 2 + ...
            (diff_redundancy_oblique - mean(diff_redundancy_oblique_empirical)).^ 2 );

    
[~,idx_sort] = sort(dist_feature,'ascend');
[matched_feature, best_params] = deal(cell(1,5));
target_feature = [mean(diff_redundancy_cardinal_empirical), mean(diff_redundancy_oblique_empirical)];
for n = 1:5

    
    matched_feature{n} = [diff_redundancy_cardinal(idx_sort(n)),diff_redundancy_oblique(idx_sort(n))];
    
    best_params{n} = [cardinal_delta(idx_sort(n)), cardinal_prior(idx_sort(n)), oblique_delta(idx_sort(n)), oblique_prior(idx_sort(n))];
end

%% match between simulation and each session
nSession = numel(diff_redundancy_cardinal_empirical);
[cardinal_delta_best, oblique_delta_best] = deal(zeros(nSession,1));
for n = 1:nSession
    dist_feature = sqrt((diff_redundancy_cardinal - diff_redundancy_cardinal_empirical(n)).^ 2 + ...
            (diff_redundancy_oblique - diff_redundancy_oblique_empirical(n)).^ 2 );

    [~,idx_sort] = sort(dist_feature,'ascend');
    cardinal_delta_best(n)  = cardinal_delta(idx_sort(1));
    oblique_delta_best(n)   = oblique_delta(idx_sort(1));
end

figure;
plot(cardinal_delta_best,'-o','color','red'); hold on
plot(oblique_delta_best,'-o','color','blue')