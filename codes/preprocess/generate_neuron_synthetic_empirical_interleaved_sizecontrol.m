clear all
clc
close all
%%
save_folder         = '../../results/filtered_neuron_synthetic/';
%save_folder_sizecontrol  = '../../results/filtered_neuron_synthetic/';
filename_list{1} = 'sampled_subset_empirical_gremlin_nTotal_128_nSample_64_b_PF_0_00';
filename_list{2} = 'sampled_subset_empirical_gremlin_nTotal_128_nSample_64_b_PF_0_80';

filename_list{3} = 'sampled_subset_empirical_rolo_nTotal_128_nSample_64_b_PF_0_00';
filename_list{4} = 'sampled_subset_empirical_rolo_nTotal_128_nSample_64_b_PF_0_80';

%%%% add a unsampled condition
filename_list{5} = 'sampled_subset_uniform_nTotal_128_nSample_64_b_PF_0_80';
for n = 1:numel(filename_list)
    if n < 5
        load(fullfile(save_folder, filename_list{n}));
    
        tokens = regexp(filename_list{n}, 'sampled_subset_empirical_([a-zA-Z0-9]+)_nTotal_([\d_]+)_nSample_([\d_]+)_b_PF_([\d_]+)', 'tokens');
        subject_code             = tokens{1}{1};
    
        switch subject_code
            case 'rolo'
                idx_sample = results_rolo.idx_sample;
            case 'gremlin'
                idx_sample = results_gremlin.idx_sample;
        end
    else
        idx_sample = [1:2:128]; % uniformly sample 64 from 128
    end
    nNeuron_total           = numel(idx_sample);
    %nNeuron_sub_list       = [1,8,32];
    nNeuron_sub             = 32;
  
    nSample       = 250;

    filtered_neurons.sessionStr              = sprintf('Model_interleaved');
    filtered_neurons.subjectCode             = 'Model';
    filtered_neurons.dateStr                 = sprintf('%d',1);
    filtered_neurons.nNeuron_kept            = nNeuron_sub;
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
    filtered_neurons.unitOptions             = filename_list{n};

    save_name = fullfile(save_folder,sprintf('%s_random_%d.mat',filename_list{n}, nSample));
    save(save_name,'filtered_neurons')
end

