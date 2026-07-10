clear all
clc
close all
global bpGlobal
bpGratingFCGlobal();

%%%% This script aims to check psychometric kernel of the model with the
%%%% parameters that best explain the neural data.

run_version = 'version_3';
switch run_version
    case 'version_1'
        %%%%%%%%%%% Version 1: only use the very best set of parameters for each
        %%%%%%%%%%% animal. Note: these parameters were based on only one repeat,
        %%%%%%%%%%% so might not be very reliable
        
        %%%% monkey R: delta_cardinal = 0.06, prior_cardinal = 0.60, 
        %%%% delta_oblique = 0.04, prior_oblique = 0.80;
        %%%% monkey G: delta_cardinal = 0.08, prior_cardinal = 1.00,
        %%%% delta_oblique = 0.04, prior_oblique = 0.60
        monkeyR_delta_cardinal(1) = 0.06;
        monkeyR_prior_cardinal(1) = 0.60;
        monkeyR_delta_oblique(1) = 0.04;
        monkeyR_prior_oblique(1) = 0.80;
        
        monkeyG_delta_cardinal(1) = 0.08;
        monkeyG_prior_cardinal(1) = 1.00;
        monkeyG_delta_oblique(1) = 0.04;
        monkeyG_prior_oblique(1) = 0.60;
    case 'version_2'
        %%% choose the best 10 sets of parameters for each animal
        load('../../results/model_match_neural_rolo_moreRepeats.mat');
        load('../../results/model_match_neural_gremlin_moreRepeats.mat');

        nChoose = 10;

        [~, index_rolo] = sort(neural_match_result_rolo.dist_feature, 'ascend');
        monkeyR_delta_cardinal  = neural_match_result_rolo.delta_cardinal(index_rolo(1:nChoose));
        monkeyR_prior_cardinal  = neural_match_result_rolo.prior_cardinal(index_rolo(1:nChoose));
        monkeyR_delta_oblique   = neural_match_result_rolo.delta_oblique(index_rolo(1:nChoose));
        monkeyR_prior_oblique   = neural_match_result_rolo.prior_oblique(index_rolo(1:nChoose));

        [~, index_gremlin] = sort(neural_match_result_gremlin.dist_feature, 'ascend');
        monkeyG_delta_cardinal  = neural_match_result_gremlin.delta_cardinal(index_gremlin(1:nChoose));
        monkeyG_prior_cardinal  = neural_match_result_gremlin.prior_cardinal(index_gremlin(1:nChoose));
        monkeyG_delta_oblique   = neural_match_result_gremlin.delta_oblique(index_gremlin(1:nChoose));
        monkeyG_prior_oblique   = neural_match_result_gremlin.prior_oblique(index_gremlin(1:nChoose));
    case 'version_3'
        %%% symmetric parameters as baseline
        delta_list = [0.04:0.01:0.08];
        prior_list = [0.5:0.1:1];
        [baseline_delta_cardinal, baseline_delta_oblique, baseline_prior_cardinal, baseline_prior_oblique] = ...
            deal(zeros(numel(delta_list) * numel(prior_list) ,1));
        k = 1;
        for i = 1:numel(delta_list)
            for j = 1:numel(prior_list)
                baseline_delta_cardinal(k) = delta_list(i);
                baseline_delta_oblique(k) = delta_list(i);
                baseline_prior_cardinal(k) = prior_list(j);
                baseline_prior_oblique(k) = prior_list(j);
                k = k+1;

            end
        end


end
% figure;
% subplot(2,1,1); histogram(neural_match_result_rolo.dist_feature,[0:1:40]);
% subplot(2,1,2); histogram(neural_match_result_gremlin.dist_feature,[0:1:40]);

%%%%%%


%% Generate synthetic data with pre-set parameters

doThis = 0;
if doThis
%%% on macbook
home_folder                     = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';
%%%% on linux
if ~exist(home_folder)
    home_folder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/data_for_psykernel';
end


    
runOptions.number_samples_per_evidence         = 6;
runOptions.stimulus_regime                     = 'dynamic-switching-signal-blocked';

%%%%% four "sessions" for each set of parameter
%%%%% monkey R
nSession = numel(monkeyR_delta_cardinal);
for n = 1:nSession

    runOptions.image_task_list                     = {'cardinal';'oblique'};

    runOptions.delta_list                          = [monkeyR_delta_cardinal(n), monkeyR_delta_oblique(n)];
    runOptions.prior_task_list                     = {[monkeyR_prior_cardinal(n), 1 - monkeyR_prior_cardinal(n)];...
                                                       [1 - monkeyR_prior_oblique(n), monkeyR_prior_oblique(n)]};
    
    runOptions.nRepeats                            = 5000;
    
    runOptions.stimulus_contrast_list              = {[0,6],[6,0]};
    runOptions.nRepeats_list                       = ones(size(runOptions.stimulus_contrast_list)) * runOptions.nRepeats;
    
    runOptions.run_ori_energy           = true;
    runOptions.n_ori_bin                = 12; 
    runOptions.task_mode                = 'single';
    %%%%%%%% 01/06/2025, change back to non-clamp prior
    runOptions.clamp_prior              = false;
    runOptions.save_folder              = fullfile(home_folder, sprintf('best_parameters_monkeyR_set_%d_contrast6',n));
    runOptions.nSession                 = 5; %%% run 5 sessions per condition
    
    generate_data_model_psyKernel(runOptions);
end

%%%%% monkey G
nSession = numel(monkeyG_delta_cardinal);
for n = 1:nSession

    runOptions.image_task_list                     = {'cardinal';'oblique'};

    runOptions.delta_list                          = [monkeyG_delta_cardinal(n), monkeyG_delta_oblique(n)];
    runOptions.prior_task_list                     = {[monkeyG_prior_cardinal(n), 1 - monkeyG_prior_cardinal(n)];...
                                                       [1 - monkeyG_prior_oblique(n), monkeyG_prior_oblique(n)]};
    
    runOptions.nRepeats                            = 5000;
    
    runOptions.stimulus_contrast_list              = {[0,6],[6,0]};
    runOptions.nRepeats_list                       = ones(size(runOptions.stimulus_contrast_list)) * runOptions.nRepeats;
    
    runOptions.run_ori_energy           = true;
    runOptions.n_ori_bin                = 12; 
    runOptions.task_mode                = 'single';
    %%%%%%%% 01/06/2025, change back to non-clamp prior
    runOptions.clamp_prior              = false;
    runOptions.save_folder              = fullfile(home_folder, sprintf('best_parameters_monkeyG_set_%d_contrast6',n));
    runOptions.nSession                 = 5; %%% run 5 sessions per condition
    
    generate_data_model_psyKernel(runOptions);
end
end
%% generate data with the baseline (symmetrical) parameters

doThis = 1;
if doThis
    %%% on macbook
    home_folder                     = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';
    %%%% on linux
    if ~exist(home_folder)
        home_folder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/data_for_psykernel';
    end
    
        
    runOptions.number_samples_per_evidence         = 6;
    runOptions.stimulus_regime                     = 'dynamic-switching-signal-blocked';
    
    nSession = numel(baseline_delta_cardinal);
    
    for n = 1:nSession
    
        runOptions.image_task_list                     = {'cardinal';'oblique'};
    
        runOptions.delta_list                          = [baseline_delta_cardinal(n), baseline_delta_oblique(n)];
        runOptions.prior_task_list                     = {[baseline_prior_cardinal(n), 1 - baseline_prior_cardinal(n)];...
                                                           [1 - baseline_prior_oblique(n), baseline_prior_oblique(n)]};
        
        runOptions.nRepeats                            = 5000;
        
        runOptions.stimulus_contrast_list              = {[0,6],[6,0]};
        runOptions.nRepeats_list                       = ones(size(runOptions.stimulus_contrast_list)) * runOptions.nRepeats;
        
        runOptions.run_ori_energy           = true;
        runOptions.n_ori_bin                = 12; 
        runOptions.task_mode                = 'single';
        %%%%%%%% 01/06/2025, change back to non-clamp prior
        runOptions.clamp_prior              = false;
        runOptions.save_folder              = fullfile(home_folder, sprintf('baseline_parameters_set_%d_contrast6',n));
        runOptions.nSession                 = 5; %%% run 5 sessions per condition
        
        generate_data_model_psyKernel(runOptions);
    end
end


