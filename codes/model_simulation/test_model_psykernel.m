clear all
clc
close all
%%
save_folder = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';

doGenerate = 1;
if doGenerate
    number_samples_per_evidence = 6;
    stimulus_regime             = 'dynamic-switching-signal-blocked';
    %stimulus_regime             = 'dynamic-delayed';
    run_ori_energy              = true;
    n_ori_bin                   = 12;
    stimulus_contrast_list      = {[6,0],[3,0],[0,0],[0,3],[0,6]};
    prior_task                  = [1, 0];
    image_task                  = 'cardinal'; 
    nRepeats                    = 128;
    
    for n = 1:numel(stimulus_contrast_list)
        stimulus_contrast = stimulus_contrast_list{n};
        P = S_Exp_Para('test-interleaved','I.stimulus_regime',stimulus_regime,...
                    'G.number_samples_per_evidence',number_samples_per_evidence,...
                    'I.run_ori_energy', run_ori_energy,...
                    'I.n_ori_bin', n_ori_bin,...
                    'I.stimulus_contrast',stimulus_contrast,...
                    'G.prior_task',prior_task,...
                    'I.image_task',image_task,...
                    'S.number_repetitions', nRepeats);
    
        dat = S_Experiment(P);
    
        %%%% extract choice and energy
        data_use(n).contrast = stimulus_contrast(2) - stimulus_contrast(1);
        data_use(n).image_task = image_task;
        data_use(n).prior_task = prior_task;
        data_use(n).decision = (dat.O(:,3,end) > 0.5) + 1; 
        data_use(n).signal = dat.Signal;
        data_use(n).oriEnergy = dat.oriEnergy;
        data_use(n).orientation_bins_center = dat.orientation_bins_center / pi * 180;
        data_use(n).phi_x = P.G.phi_x;
        data_use(n).n_zero_signal = P.I.n_zero_signal;
    end
    save(fullfile(save_folder,'test_data'),'data_use')
end
%% use the same code for the empirical data to compute
load(fullfile(save_folder,'test_data'));

energy_use = 'proj';
% normalized ori energy of all trials  nTrial x nOri x nFrame
switch energy_use
    case 'fft'
        n_zero_signal = data_use(1).n_zero_signal;
        X_zero_trials = data_use(cell2mat({data_use(:).contrast}) == 0).oriEnergy;
        X_zero_trials = X_zero_trials(:,:,n_zero_signal+1:end);
        muZero = mean(X_zero_trials(:));
        stdZero = std(X_zero_trials(:));
        
        X_all_trials =  {data_use(:).oriEnergy};
        X_all_trials = cat(1,X_all_trials{:});
        X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;
        kernelOptions.orientationsDEG = data_use(1).orientation_bins_center;
    case 'proj'
        n_zero_signal = data_use(1).n_zero_signal;
        X_zero_trials = data_use(cell2mat({data_use(:).contrast}) == 0).signal;
        X_zero_trials = X_zero_trials(:,:,n_zero_signal+1:end);
        muZero = mean(X_zero_trials(:));
        stdZero = std(X_zero_trials(:));
        
        X_all_trials =  {data_use(:).signal};
        X_all_trials = cat(1,X_all_trials{:});
        X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;

        phi_x = data_use(1).phi_x / pi * 180;

        kernelOptions.orientationsDEG   = [0:15:165];
        
        for i = 1:numel(kernelOptions.orientationsDEG )
             [~, idx(i)] =  min(abs(kernelOptions.orientationsDEG(i) - phi_x));
        end
        X_all_trials_norm = X_all_trials_norm(:,idx,:);

        

end


choice_all_trials   = {data_use(:).decision};
choice_all_trials   = cat(1,choice_all_trials{:});
switch data_use(1).image_task
    case 'cardinal'
        choice_list         = [0,90];
    case 'oblique'
        choice_list         = [45,135];
end
index               = [1:numel(choice_all_trials)]; % use all trials 
choice_all_trials   = choice_list(choice_all_trials); %%%% tranlate index to orientation


%%%%% parameters same as what we used for the empirical data
fittingKernelParams.hypers  = [0.5 6 60 0.15 0.06];
fittingKernelParams.nLapse = 2;

map_params.cardinal =  psyKernel.compute_SPTP_kernel(X_all_trials_norm,choice_all_trials,...
                    kernelOptions,choice_list,index,fittingKernelParams);
%%
%%% use a simple way to get the kernels?

w_all = squeeze(mean(X_all_trials_norm(choice_all_trials == choice_list(2),:,:),1) ...
    - mean(X_all_trials_norm(choice_all_trials == choice_list(1),:,:),1));


%%
figure;
subplot(2,2,1)
plot(kernelOptions.orientationsDEG, map_params.cardinal.w_ori);
subplot(2,2,2)
plot(map_params.cardinal.w_time);

subplot(2,2,3)
plot(kernelOptions.orientationsDEG, mean(w_all,2));
subplot(2,2,4)
idx_pref = choice_all_trials == choice_list(2);
idx_anti = choice_all_trials == choice_list(1);
ixp = 1;
ixa = 6;
prefpref=mean(X_all_trials_norm(idx_pref,ixp,:));
prefanti=mean(X_all_trials_norm(idx_pref,ixa,:));
antipref=mean(X_all_trials_norm(idx_anti,ixp,:));
antianti=mean(X_all_trials_norm(idx_anti,ixa,:));
kernel=prefpref-prefanti-antipref+antianti;
plot(squeeze(kernel))