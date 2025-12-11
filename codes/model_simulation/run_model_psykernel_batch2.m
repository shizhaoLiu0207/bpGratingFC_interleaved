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



number_samples_per_evidence         = 6;
stimulus_regime                     = 'dynamic-switching-signal-blocked';

prior_task_list                     = {[1, 0];[0.7, 0.3];[0.5, 0.5]; [0.3,0.7]};
image_task_list                     = {'cardinal';'cardinal';'cardinal';'cardinal'};


nRepeats                            = 1000;
run_ori_energy                      = true;
n_ori_bin                           = 12;

stimulus_contrast_list              = {[0,6],[6,0]};
nRepeats_list                       = ones(size(stimulus_contrast_list)) * nRepeats;

runOptions.run_ori_energy           = true;
runOptions.n_ori_bin                = 12; 
runOptions.task_mode                = 'single';
runOptions.clamp_prior              = true;
runOptions.save_folder              = fullfile(home_folder, 'moreRepeats');
runOptions.nSession                 = 20; %%% run 10 sessions per condition

generate_data_model_psyKernel(prior_task_list, image_task_list, stimulus_contrast_list, nRepeats_list, runOptions)