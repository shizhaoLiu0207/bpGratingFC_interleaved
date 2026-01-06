function [map_params, orientationsDEG, X_all_trials_unnorm] = compute_psyKernel_model(data_use, energy_use)
%%%%% parameters same as what we used for the empirical data
fittingKernelParams.hypers  = [0.5 6 60 0.15 0.06];
fittingKernelParams.nLapse  = 2;
% fittingKernelParams.hypers  = [0.5 6 60 2 0];
% fittingKernelParams.nLapse  = 0;

[X_all_trials_norm, choice_all_trials, ~, X_all_trials_unnorm, orientationsDEG] = util_it.extract_X_psyKernel(data_use, energy_use);




index               = [1:numel(choice_all_trials)]; % use all trials 
kernelOptions.orientationsDEG = orientationsDEG;

choice_list = unique(choice_all_trials);
map_params =  psyKernel.compute_SPTP_kernel(X_all_trials_norm,choice_all_trials,...
                    kernelOptions,choice_list,index,fittingKernelParams);


end