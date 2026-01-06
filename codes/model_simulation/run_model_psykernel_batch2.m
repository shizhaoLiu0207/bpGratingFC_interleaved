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

prior_task_list                     = {[0, 1];[0.3,0.7];[0.5,0.5];[0.7,0.3]};
image_task_list                     = {'oblique';'oblique';'oblique';'oblique'};


nRepeats                            = 1000;
run_ori_energy                      = true;
n_ori_bin                           = 12;

stimulus_contrast_list              = {[0,6],[6,0]};
nRepeats_list                       = ones(size(stimulus_contrast_list)) * nRepeats;

runOptions.run_ori_energy           = true;
runOptions.n_ori_bin                = 12; 
runOptions.task_mode                = 'single';
runOptions.clamp_prior              = true;
runOptions.save_folder              = fullfile(home_folder, 'moreRepeats_contrast6');
runOptions.nSession                 = 20; %%% run 10 sessions per condition

generate_data_model_psyKernel(prior_task_list, image_task_list, stimulus_contrast_list, nRepeats_list, runOptions);
end
%%
doCompute = 0;
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

        %%%% number of repeats and list of coherence
        psyKernel_model(n).cohr_list        = unique(abs(cell2mat({data_use(:).contrast})));
        for k = 1:numel(psyKernel_model(n).cohr_list)
            idx = cell2mat({data_use(:).contrast}) == psyKernel_model(n).cohr_list;
            psyKernel_model(n).num_repeats(k) = numel(data_use(idx).decision);
        end
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
mode                    = 'Single';
figure
subplot(2,4,1)
image_task              = 'cardinal';
prior_task              = [1,0];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,2)
image_task              = 'cardinal';
prior_task              = [0.7,0.3];
[r_cardinal, r_oblique] =fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));



subplot(2,4,3)
image_task              = 'cardinal';
prior_task              = [0.5,0.5];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,4)
image_task              = 'cardinal';
prior_task              = [0.3,0.7];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,5)
image_task              = 'oblique';
prior_task              = [0, 1];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,6)
image_task              = 'oblique';
prior_task              = [0.3,0.7];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,7)
image_task              = 'oblique';
prior_task              = [0.5,0.5];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));


subplot(2,4,8)
image_task              = 'oblique';
prior_task              = [0.7,0.3];
[r_cardinal, r_oblique] = fig_it.plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task, mode);
title(sprintf('Prior_{cardinal} = %.1f,  Prior_{oblique} = %.1f \n r_{cardinal} = %.2f, r_{oblique} = %.2f', ...
    prior_task(1), prior_task(2), r_cardinal, r_oblique));

%%  correlation between the template as a function of sessions (each session has 1000 trials in this batch)
orientationsDEG = psyKernel_model(1).orientationsDEG;
ideal_kernel_cardinal            = sin(pi * (orientationsDEG - 45)/90);
ideal_kernel_oblique             = sin(pi * (orientationsDEG - 90)/90);


prior_cardinal  = 1;
prior_oblique   = 0;
image_task      = 'cardinal';
mode            = 'Single'; 

idx = find(cell2mat({psyKernel_model(:).prior_cardinal}) == prior_cardinal &...
      cell2mat({psyKernel_model(:).prior_oblique})  == prior_oblique & ...
      strcmp({psyKernel_model(:).image_task}, image_task) & ...
      strcmp({psyKernel_model(:).mode}, mode));


w_ori_all = {psyKernel_model(idx).w_ori};
w_ori_all = cat(2,w_ori_all{:});


nSession = size(w_ori_all,2);
[r_cardinal, r_oblique] = deal(cell(nSession,1));
for nSample = 1:nSession
    nBootstrap = 100;

     if nchoosek(nSession,nSample) <= nBootstrap
            combination_all = nchoosek([1:nSession],nSample); 
        elseif nchoosek(nSession,nSample) <= 5 * nBootstrap
            combination_all = nchoosek([1:nSession],nSample); 
            combination_all = combination_all(1:nBootstrap,:); % each row is one combination
        else
            combination_all = zeros(nBootstrap,nSample);
            for k = 1:nBootstrap
                 combination_all(k,:) = randsample(nSession,nSample,'false');
            end
     end
     nBootstrap_real = size(combination_all,1);
     [r_cardinal{nSample}, r_oblique{nSample}] = deal(zeros(nBootstrap_real, 1));
     for k = 1:nBootstrap_real
        w_ori_avg = mean(w_ori_all(:,combination_all(k,:)),2);
        r_cardinal{nSample}(k) = corr(w_ori_avg, ideal_kernel_cardinal');
        r_oblique{nSample}(k) = corr(w_ori_avg, ideal_kernel_oblique');
     end

end

figure
plot([1:nSession], cellfun(@mean,r_cardinal),'color','red'); hold on
errorbar([1:nSession], cellfun(@mean,r_cardinal),cellfun(@std,r_cardinal),'color','red')
plot(cellfun(@mean,r_oblique),'color','blue');
errorbar([1:nSession], cellfun(@mean,r_oblique),cellfun(@std,r_oblique),'color','blue')
grid on
