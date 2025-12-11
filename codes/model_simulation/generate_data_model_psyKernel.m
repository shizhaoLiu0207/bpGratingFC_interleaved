function  generate_data_model_psyKernel(prior_task_list, image_task_list, stimulus_contrast_list, nRepeats_list, runOptions)
assert(numel(stimulus_contrast_list) == numel(nRepeats_list), 'Contrast and nRepeats number do not match');
assert(numel(prior_task_list) == numel(image_task_list), 'Task prior and image task numbers do not match');

run_ori_energy  = runOptions.run_ori_energy;
n_ori_bin       = 12; 
task_mode       = runOptions.task_mode;
clamp_prior     = runOptions.clamp_prior;
nSession        = runOptions.nSession;
save_folder     = runOptions.save_folder;

switch task_mode
    case 'single'
        run_model_mode = 'test-interleaved';
    case 'dual'
        run_model_mode = 'test-dualTask';
end


number_samples_per_evidence = 6;
stimulus_regime             = 'dynamic-switching-signal-blocked';

%%%%% generate data of single task mode
nCond_single    = numel(prior_task_list);
for i  = 1:nCond_single
    prior_task = prior_task_list{i};
    image_task = image_task_list{i};

    prior_task_cardinal = prior_task(1);
    prior_task_oblique  = prior_task(2);
    prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
    prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');
    for t = 1: nSession
        save_name  = sprintf('Single_image_task_%s_prior_cardinal_%s_prior_oblique_%s_session_%d.mat', image_task, prior_cardinal_str, prior_oblique_str, t);
        if isfile(fullfile(save_folder,save_name))
            continue
        end
        for n  = 1:numel(stimulus_contrast_list)
            stimulus_contrast   = stimulus_contrast_list{n};
            nRepeats            = nRepeats_list(n);

            P                   = S_Exp_Para(run_model_mode,'I.stimulus_regime',stimulus_regime,...
                'G.number_samples_per_evidence',number_samples_per_evidence,...
                'I.run_ori_energy', run_ori_energy,...
                'I.n_ori_bin', n_ori_bin,...
                'I.stimulus_contrast',stimulus_contrast,...
                'G.prior_task',prior_task,...
                'I.image_task',image_task,...
                'S.number_repetitions', nRepeats,...
                'G.clamp_prior', clamp_prior);

           
         
         

            dat = S_Experiment(P);

            %%%% extract choice and energy
            data_use(n).contrast = stimulus_contrast(2) - stimulus_contrast(1);
            data_use(n).image_task = image_task;
            data_use(n).prior_task = prior_task;
            data_use(n).decision = (dat.O(:,3,end) > 0.5) + 1;
            data_use(n).signal = dat.Signal;
            data_use(n).oriEnergy = dat.oriEnergy;
            data_use(n).orientation_bins_center = dat.orientation_bins_center / pi * 180;
            data_use(n).P = P;
            %data_use(n).n_zero_signal = P.I.n_zero_signal;

        end
        save(fullfile(save_folder, save_name),'data_use');
    end
end