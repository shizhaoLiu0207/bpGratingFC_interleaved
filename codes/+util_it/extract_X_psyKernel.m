function [X_all_trials_norm, choice_all_trials, contrast_all_trials, X_all_trials_unnorm, orientationsDEG] = extract_X_psyKernel(data_use, energy_use)
% normalized ori energy of all trials  nTrial x nOri x nFrame
for n = 1:numel(data_use)
    data_use(n).contrast_trials = ones(size(data_use(n).decision)) * data_use(n).contrast;
end
switch energy_use
    case 'fft'
        load('../../results/model_oriEnergy_zeroContrast.mat','zero_energy_fft');
        muZero = mean(zero_energy_fft(:));
        stdZero = std(zero_energy_fft(:));

        n_zero_signal = data_use(1).P.I.n_zero_signal;
      
        X_all_trials_unnorm =  {data_use(:).oriEnergy};
        X_all_trials = cat(1,X_all_trials_unnorm{:});
        X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;
        orientationsDEG = data_use(1).orientation_bins_center;

        contrast_all_trials = cell2mat({data_use(:).contrast_trials});
    case 'proj'
        load('../../results/model_oriEnergy_zeroContrast.mat','zero_energy_proj');
        muZero = mean(zero_energy_proj(:));
        stdZero = std(zero_energy_proj(:));
        %%%%% Need to use the square of the projections
        for n  = 1:numel(data_use)
            data_use(n).energy_proj = data_use(n).signal .^ 2;
        end
        n_zero_signal =  data_use(1).P.I.n_zero_signal;
      
        
        X_all_trials_unnorm =  {data_use(:).energy_proj};
        %%%% downsampling orientation
        try
            phi_x = data_use(1).P.G.phi_x / pi * 180;
        catch
            phi_x  = data_use(1).P.phi_x;
        end
        orientationsDEG   = [0:15:165];
        idx = zeros(numel(orientationsDEG, 1));
        for i = 1:numel(orientationsDEG)
             [~, idx(i)] =  min(abs(orientationsDEG(i) - phi_x));
        end
        
        for n = 1:numel(X_all_trials_unnorm)
            X_all_trials_unnorm{n} = X_all_trials_unnorm{n}(:,idx,:);
        end

        X_all_trials = cat(1,X_all_trials_unnorm{:});
        if ndims(X_all_trials) == 3
            X_all_trials = X_all_trials(:,:,n_zero_signal+1:end);
        end
        X_all_trials_norm = (X_all_trials - muZero) / stdZero;

        contrast_all_trials = cell2mat({data_use(:).contrast_trials});
      
        

        
end

choice_all_trials   = {data_use(:).decision};
choice_all_trials   = cat(1,choice_all_trials{:});
switch data_use(1).image_task
    case 'cardinal'
        choice_list         = [0,90];
    case 'oblique'
        choice_list         = [45,135];
end

choice_all_trials   = choice_list(choice_all_trials); %%%% tranlate index to orientation

end