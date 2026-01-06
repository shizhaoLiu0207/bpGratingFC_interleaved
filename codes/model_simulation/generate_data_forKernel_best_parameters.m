clear all
clc
close all
%%
%%% on macbook
home_folder                     = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';
%%%% on linux
if ~exist(home_folder)
    home_folder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/data_for_psykernel';
end


%%%% monkey R: delta_cardinal = 0.06, prior_cardinal = 0.60, 
%%%% delta_oblique = 0.04, prior_oblique = 0.80;
%%%% monkey G: delta_cardinal = 0.08, prior_cardinal = 1.00,
%%%% delta_oblique = 0.04, prior_oblique = 0.60
monkeyR_delta_cardinal = 0.06;
monkeyR_prior_cardinal = 0.60;
monkeyR_delta_oblique = 0.04;
monkeyR_prior_oblique = 0.80;

monkeyG_delta_cardinal = 0.08;
monkeyG_prior_cardinal = 1.00;
monkeyG_delta_oblique = 0.04;
monkeyG_prior_oblique = 0.60;

runOptions.number_samples_per_evidence         = 6;
runOptions.stimulus_regime                     = 'dynamic-switching-signal-blocked';

runOptions.delta_list                          = [monkeyR_delta_cardinal, monkeyR_delta_oblique, ...
                                        monkeyG_delta_cardinal, monkeyG_delta_oblique];
runOptions.prior_task_list                     = {[monkeyR_prior_cardinal, 1 - monkeyR_prior_cardinal];...
                                       [1 - monkeyR_prior_oblique, monkeyR_prior_oblique];...
                                       [monkeyG_prior_cardinal, 1 - monkeyG_prior_cardinal];...
                                       [1 - monkeyG_prior_oblique, monkeyG_prior_oblique]};

runOptions.image_task_list                     = {'cardinal';'oblique';'cardinal';'oblique'};


runOptions.nRepeats                            = 5000;

runOptions.stimulus_contrast_list              = {[0,6],[6,0]};
runOptions.nRepeats_list                       = ones(size(runOptions.stimulus_contrast_list)) * runOptions.nRepeats;

runOptions.run_ori_energy           = true;
runOptions.n_ori_bin                = 12; 
runOptions.task_mode                = 'single';
%%%%%%%% 01/06/2025, change back to non-clamp prior
runOptions.clamp_prior              = false;
runOptions.save_folder              = fullfile(home_folder, 'best_parameters_contrast6');
runOptions.nSession                 = 5; %%% run 5 sessions per condition

generate_data_model_psyKernel(runOptions);
