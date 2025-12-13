function [map_params, orientationsDEG, X_all_trials_return, contrast_list] = compute_psyKernel_model(data_use, energy_use)
%%%%% parameters same as what we used for the empirical data
% fittingKernelParams.hypers  = [0.5 6 60 0.15 0.06];
% fittingKernelParams.nLapse  = 2;
fittingKernelParams.hypers  = [0.5 6 60 2 0];
fittingKernelParams.nLapse  = 0;

% normalized ori energy of all trials  nTrial x nOri x nFrame
switch energy_use
    case 'fft'
        load('../../results/model_oriEnergy_zeroContrast.mat','zero_energy_fft');
        muZero = mean(zero_energy_fft(:));
        stdZero = std(zero_energy_fft(:));

        n_zero_signal = data_use(1).P.I.n_zero_signal;
      
        X_all_trials_return =  {data_use(:).oriEnergy};
        X_all_trials = cat(1,X_all_trials_return{:});
        X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;
        orientationsDEG = data_use(1).orientation_bins_center;

        contrast_list = cell2mat({data_use(:).contrast});
    case 'proj'
        load('../../results/model_oriEnergy_zeroContrast.mat','zero_energy_proj');
        muZero = mean(zero_energy_proj(:));
        stdZero = std(zero_energy_proj(:));
        %%%%% Need to use the square of the projections
        for n  = 1:numel(data_use)
            data_use(n).energy_proj = data_use(n).signal .^ 2;
        end
        n_zero_signal =  data_use(1).P.I.n_zero_signal;
      
        
        X_all_trials_return =  {data_use(:).energy_proj};
        %%%% downsampling orientation
        phi_x = data_use(1).P.G.phi_x / pi * 180;
        orientationsDEG   = [0:15:165];
        idx = zeros(numel(orientationsDEG, 1));
        for i = 1:numel(orientationsDEG)
             [~, idx(i)] =  min(abs(orientationsDEG(i) - phi_x));
        end
        
        for n = 1:numel(X_all_trials_return)
            X_all_trials_return{n} = X_all_trials_return{n}(:,idx,:);
        end

        X_all_trials = cat(1,X_all_trials_return{:});
        X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;

        contrast_list = cell2mat({data_use(:).contrast});
      
        

        
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




kernelOptions.orientationsDEG = orientationsDEG;
map_params =  psyKernel.compute_SPTP_kernel(X_all_trials_norm,choice_all_trials,...
                    kernelOptions,choice_list,index,fittingKernelParams);


end