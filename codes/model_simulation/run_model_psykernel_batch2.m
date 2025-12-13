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

doGenerate = 0;
if doGenerate 




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

generate_data_model_psyKernel(prior_task_list, image_task_list, stimulus_contrast_list, nRepeats_list, runOptions);
end
%%
doCompute = true;
version_name = 'moreRepeats_contrast_6';
save_folder = fullfile(home_folder, version_name);
energy_use = 'proj';
save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_%s_energy_%s',version_name, energy_use);

if doCompute
    
    data_list = dir(fullfile(save_folder,'*.mat'));
    data_list(strcmp({data_list(:).name},'test_data.mat')) = [];
    
    psyKernel_model = struct();
    for n = 1:numel(data_list)
        fprintf('Computing session %d/%d \n',n,numel(data_list));
        load(fullfile(save_folder, data_list(n).name));
        
        [map_params,orientationsDEG] = util_it.compute_psyKernel_model(data_use, energy_use);
      
        tokens = regexp(data_list(n).name, '([a-zA-Z0-9]+)_image_task_([a-zA-Z0-9]+)_prior_cardinal_([-]?[\d_]+)_prior_oblique_([\d_]+)_session_([\d_]+)', 'tokens');
        % Convert extracted strings back to numbers
        extracted_params        = tokens{1}; % Extract matched tokens
        task_mode_str           = extracted_params{1};
        imagetask_str           = extracted_params{2};
        prior_cardinal_str      = strrep(extracted_params{3}, '_', '.');
        prior_oblique_str       = strrep(extracted_params{4}, '_', '.');
        session_num_str         = extracted_params{5};
        % bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
        % delta_str               = strrep(extracted_params{3}, '_', '.');
        % taskprior_str           = strrep(extracted_params{4}, '_', '.');
    
        psyKernel_model(n).mode             = task_mode_str;
        psyKernel_model(n).image_task       = imagetask_str;
        psyKernel_model(n).prior_cardinal   = str2double(prior_cardinal_str);
        psyKernel_model(n).prior_oblique    = str2double(prior_oblique_str);
        psyKernel_model(n).sessionNum       = str2double(session_num_str);
        psyKernel_model(n).w_ori            = map_params.w_ori;
        psyKernel_model(n).w_time           = map_params.w_time;
        psyKernel_model(n).orientationsDEG  = orientationsDEG;
    
        psyKernel_model(n).energy_use       = energy_use;
    end
    save(save_name,'psyKernel_model');
end
%% visualization results
version_name = 'moreRepeats_contrast_6';

energy_use = 'proj';
save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_%s_energy_%s',version_name, energy_use);
load(save_name);
%%% Alignment between the ideal kernel as a function of repeats
kernel_type = 'spatial';
figure
subplot(1,2,1)
image_task              = 'cardinal';
prior_task_list{1}      = [1,0];
mode_list{1}            = 'Single';
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
%%% Plot kernels of different priors together to check if makes sense
%%
figure
prior_task_list      = {[1,0]};
mode_list            = {'Single'};
fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)