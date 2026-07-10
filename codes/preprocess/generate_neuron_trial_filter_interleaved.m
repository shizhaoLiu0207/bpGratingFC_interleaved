clear all
clc
close all
%%%% This script generares neuron_trial filter by combining existing neuron
%%%% filter and trial filter
%% 
% neuron_filter_name_list =  {'highC_highO_passiveViewing_sizeControl';'highC_lowO_passiveViewing_sizeControl';'lowC_highO_passiveViewing_sizeControl';'lowC_lowO_passiveViewing_sizeControl';...
%                 'highC_highO_mainTask_sizeControl';'highC_lowO_mainTask_sizeControl';'lowC_highO_mainTask_sizeControl';'lowC_lowO_mainTask_sizeControl';...
%                 'highC_lowO_lowC_highO_passiveViewing_sizeControl';'highC_highO_lowC_lowO_passiveViewing_sizeControl';...
%                 'highC_lowO_lowC_highO_mainTask_sizeControl';'highC_highO_lowC_lowO_mainTask_sizeControl'};
%neuron_filter_name      = 'coef1_hVis2_FR1_interleaved_forEye_sizeControl'; 
neuron_filter_name_list = {'highC_highO_passiveViewing';'lowC_lowO_passiveViewing';...
     'highC_highO_mainTask';'lowC_lowO_mainTask'};
for i = 1:numel(neuron_filter_name_list)
    neuron_filter_name = neuron_filter_name_list{i};
    trial_filter_name       = 'all_trials';
    
    combined_filter_name    = [trial_filter_name,'_', neuron_filter_name];
    
    neuron_filter_folder    = fullfile('../../results/filter_neuron/', neuron_filter_name);
    trial_filter_folder     = fullfile('../../results/filter_trial/', trial_filter_name);
    combined_filter_folder  = fullfile('../../results/filter_trial_neuron',combined_filter_name);
    mkdir(combined_filter_folder);
    
    combined_saveName   = fullfile(combined_filter_folder,sprintf('filtered_trials_neurons_%s_%s',trial_filter_name, neuron_filter_name));
    combined_txtName    = fullfile(combined_filter_folder,sprintf('criterion_%s_%s.txt',trial_filter_name, neuron_filter_name));
    
    %%%%%%%%%%%%%%%%%%%%% combine the struct variable %%%%%%%%%%%%%%%%%%%%%%
    neuron_filter    = load(fullfile(neuron_filter_folder, sprintf('filtered_neurons_%s',neuron_filter_name)));
    trial_filter     = load(fullfile(trial_filter_folder,  sprintf('filtered_trials_%s',trial_filter_name)));
    filtered_neurons = neuron_filter.filtered_neurons;
    filtered_trials  = trial_filter.filtered_trials;
    
    % if numel(filtered_neurons) ~= numel(filtered_trials)
    %     error('Number of sessions in filtered neuron struct and filtered trial struct does not match');
    % end
    
    filtered_trials_neurons = filtered_neurons;
    f_list = fieldnames(filtered_trials);
    for n = 1 : numel(filtered_trials_neurons)
        sessionStr = filtered_neurons(n).sessionStr;
    
           
        idx = strcmp({filtered_trials(:).sessionStr}, sessionStr(1:10));
        filtered_trials_neurons(n).trialInd_kept = filtered_trials(idx).trialInd_kept;
        filtered_trials_neurons(n).nTrial_total = filtered_trials(idx).nTrial_toal;
    
    
    end
    
    keep_options.trialOptions =  trial_filter.keep_options.trialOptions;
    try
    keep_options.neuronOptions = neuron_filter.keep_options.neuronOptions;
    catch
        keep_options.neuronOptions = neuron_filter.keep_options.unitOptions;
    end
    
    save(combined_saveName,'filtered_trials_neurons','keep_options');
end