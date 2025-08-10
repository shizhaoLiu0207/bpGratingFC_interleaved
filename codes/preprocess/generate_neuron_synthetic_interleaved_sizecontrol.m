clear all
clc
close all
%%
save_folder         = '../../results/filtered_neuron_synthetic/subsample';
save_folder_sizecontrol  = '../../results/filtered_neuron_synthetic/subsample_random';
filename_list = dir(fullfile(save_folder,'*.mat'));

for n = 1:numel(filename_list)
    load(fullfile(save_folder, filename_list(n).name));

    tokens = regexp(filename_list(n).name, 'nNeuron_([\d_]+)_b*', 'tokens');
    nNeuron_total           = str2double(tokens{1}{1});
    nNeuron_sample          = numel(idx_sample);
    %nNeuron_sub_list       = [1,8,32];
    nNeuron_sub             = 32;
  
    nSample       = 1000;

    filtered_neurons.sessionStr              = sprintf('Model_interleaved_nNeuron%d',nNeuron_total);
    filtered_neurons.subjectCode             = 'Model';
    filtered_neurons.dateStr                 = sprintf('%d',1);
    filtered_neurons.nNeuron_kept            = nNeuron_sample;
    filtered_neurons.nNeuron_total           = nNeuron_total;
    combination_all = zeros(nSample, nNeuron_sub);
    for t = 1:nSample
        if nNeuron_sub == 1
            combination_all(t,:) = t;
        else
            combination_all(t,:) = randsample(idx_sample, nNeuron_sub); 
        end
    end
    filtered_neurons.neuronIdx_kept          = combination_all;
    filtered_neurons.unitOptions             = filename_list(n).name;

    save_name = fullfile(save_folder_sizecontrol,sprintf('%s_random_%d.mat',filename_list(n).name, nSample));
    save(save_name,'filtered_neurons')
end

