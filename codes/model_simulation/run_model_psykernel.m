clear all
clc
close all

data_mode = 'allLevels'; % zeroOnly or allLevels
energy_use = 'fft';

%% Generate data to use
switch data_mode
    case 'allLevels'
        save_folder                 = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel';
    case 'zeroOnly'
        save_folder                 = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel_zeroOnly';
    case 'clampprior'
        save_folder                 = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel_clampprior';
    case 'clampprior_zeroOnly'
        save_folder                 = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/data_for_psykernel_clampprior_zeroOnly';

end
if ~isfolder(save_folder)
    mkdir(save_folder);
end
doGenerate = 0;
if doGenerate
    number_samples_per_evidence = 6;
    stimulus_regime             = 'dynamic-switching-signal-blocked';
    
    nSession                    = 10; %%% run 10 sessions per condition
    
    run_ori_energy              = true;
    n_ori_bin                   = 12;
    switch data_mode
        case {'allLevels','clampprior'}
            stimulus_contrast_list      = {[12,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,12]};
            nRepeats_list               = [50, 100, 150, 150, 150, 100, 50];
        case {'zeroOnly','clampprior_zeroOnly'}
            stimulus_contrast_list = {[0,0]};
            nRepeats_list = [500];
    end
    
    assert(numel(stimulus_contrast_list) == numel(nRepeats_list), 'Contrast and nRepeats number do not match');
    
    prior_task_single_list      = {[1, 0]; [0, 1]; [0.7, 0.3]; [0.7,0.3]; [0.5, 0.5]};
    image_task_single_list      = {'cardinal';'oblique';'cardinal';'oblique';'cardinal'};
    assert(numel(prior_task_single_list) == numel(image_task_single_list), 'Task prior and image task numbers do not match (single)');
    
    
    prior_task_dual_list        = {[0.8, 0.2]; [0.8, 0.5]; [0.2, 0.5]}; 
    image_task_dual_list        = {'cardinal'; 'cardinal'; 'oblique'}; 
    assert(numel(prior_task_dual_list) == numel(image_task_dual_list), 'Task prior and image task numbers do not match (dual)');
    
    
   
    
    %%%%% generate data of dual task mode
    nCond_dual      = numel(prior_task_dual_list);
    
    for i  = 1:nCond_dual
        prior_task = prior_task_dual_list{i};
        image_task = image_task_dual_list{i};
    
        prior_task_cardinal = prior_task(1);
        prior_task_oblique  = prior_task(2);
        prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
        prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');
        for t = 1:nSession
            save_name  = sprintf('Dual_image_task_%s_prior_cardinal_%s_prior_oblique_%s_session_%d.mat', image_task, prior_cardinal_str, prior_oblique_str, t); 
            if isfile(fullfile(save_folder,save_name))
                continue
            end
            for n  = 1:numel(stimulus_contrast_list)
                stimulus_contrast   = stimulus_contrast_list{n};
                nRepeats            = nRepeats_list(n);
        
                P                   = S_Exp_Para('test-dualTask','I.stimulus_regime',stimulus_regime,...
                                        'G.number_samples_per_evidence',number_samples_per_evidence,...
                                        'I.run_ori_energy', run_ori_energy,...
                                        'I.n_ori_bin', n_ori_bin,...
                                        'I.stimulus_contrast',stimulus_contrast,...
                                        'G.prior_task',prior_task,...
                                        'I.image_task',image_task,...
                                        'S.number_repetitions', nRepeats);
                if strcmp(data_mode,'clampprior')
                    P.G.clamp_prior = true;
                end        
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
end
%% compute psychometric kernels of these data
doCompute = 0;


switch data_mode
    case 'allLevels'
        save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_%s.mat', energy_use);
    case 'zeroOnly'
        save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_zeroOnly_%s.mat', energy_use);
        do_plot_energy = 0;
    case 'clampprior'
        save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_clampprior_%s.mat', energy_use);
    case 'clampprior_zeroOnly'
        save_name = sprintf('../../results/model_behavior/psyKernel_model_single_dual_test_clampprior_zeroOnly_%s.mat', energy_use);

end
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
else
    load(save_name)
end
%% plot ori energy
do_plot_energy = 0;
if do_plot_energy
    load(fullfile(save_folder, data_list(1).name));
    figure
    subplot(1,2,1)
    energy_use = 'proj';
    [~,orientationsDEG, X_all_trials_return, contrast_list] = util_it.compute_psyKernel_model(data_use, energy_use);
    plot_ori_energy(X_all_trials_return, orientationsDEG, contrast_list, data_use(1).P.I.n_zero_signal);
    title(energy_use)
    subplot(1,2,2)
    energy_use = 'fft';
    [~,orientationsDEG, X_all_trials_return, contrast_list] = util_it.compute_psyKernel_model(data_use, energy_use);
    plot_ori_energy(X_all_trials_return, orientationsDEG, contrast_list, data_use(1).P.I.n_zero_signal);
    title(energy_use)
end
%% visualization results
%%% two ideal cases, single mode


kernel_type = 'spatial';


 

figure;
set(gcf, 'Units','normalized','Position',[0,0,1,1])
subplot(2,5,1)
image_task              = 'cardinal';
prior_task_list{1}      = [1,0];
mode_list{1}            = 'Single';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,2)
image_task              = 'oblique';
prior_task_list{1}      = [0,1];
mode_list{1}            = 'Single';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,3)
image_task              = 'cardinal';
prior_task_list{1}      = [0.7,0.3];
mode_list{1}            = 'Single';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,4)
image_task              = 'oblique';
prior_task_list{1}      = [0.7,0.3];
mode_list{1}            = 'Single';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,5)
image_task              = 'cardinal';
prior_task_list{1}      = [0.5,0.5];
mode_list{1}            = 'Single';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));


subplot(2,5,6)
image_task              = 'cardinal';
prior_task_list{1}      = [0.8,0.2];
mode_list{1}            = 'Dual';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,7)
image_task              = 'cardinal';
prior_task_list{1}      = [0.8,0.5];
mode_list{1}            = 'Dual';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));

subplot(2,5,8)
image_task              = 'oblique';
prior_task_list{1}      = [0.2,0.5];
mode_list{1}            = 'Dual';
plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
title(sprintf('Prior_{cardinal} = %.1f, \n Prior_{oblique} = %.1f', prior_task_list{1}(1), prior_task_list{1}(2)));


annotation('textbox',[0.05,0.85,0.1,0.1],'string','Single','fontsize',24,'fontweight','bold','LineStyle','none')
annotation('textbox',[0.05,0.45,0.1,0.1],'string','Dual','fontsize',24,'fontweight','bold','LineStyle','none')
%%
function plot_ori_energy(X_ori_energy, orientationsDEG, contrast_list,n_zero)
    for n = 1:numel(X_ori_energy)
        X = mean(X_ori_energy{n}(:,:,n_zero+1:end),3);
        avg = mean(X, 1);
        sdv = std(X, [], 1);
        h(n) = plot(orientationsDEG,avg, 'LineWidth', 2); hold on
    
        fill([orientationsDEG, fliplr(orientationsDEG)],...
            [avg - sdv, fliplr(avg + sdv)],h(n).Color,'FaceAlpha',0.5);
        lgd_str{n} = sprintf('contrast = %d',contrast_list(n));
    end
    box off
    xlabel('OrientationDEG')
    ylabel('Orientation Energy')
    set(gca,'fontsize', 16)
    legend(h,lgd_str)
end

function plot_psyKernel_model(psyKernel_model, kernel_type, image_task, prior_task_list, mode_list)
orientationsDEG = psyKernel_model(1).orientationsDEG;
 switch image_task
        case 'cardinal'
            plot_color = 'red';
            ideal_kernel            = 0.4*sin(pi * (orientationsDEG - 45)/90);
        case 'oblique'
            plot_color = 'blue';
            ideal_kernel            = 0.4*sin(pi * (orientationsDEG - 90)/90);
 end

switch kernel_type
    case 'spatial'
            hold on
            line([0,180],[0,0],'linestyle','--','color','black')
            plot(orientationsDEG,ideal_kernel,'linestyle','--','linewidth',2,'color',plot_color);
            line([0,0],[-0.6,0.6],'linestyle','--','color','black'); 
            line([90,90],[-0.6,0.6],'linestyle','--','color','black');
            line([45,45],[-0.6,0.6],'linestyle','--','color','black');
            line([135,135],[-0.6,0.6],'linestyle','--','color','black');
            xlim([-10,190])
            ylim([-0.5,0.5])
            xlabel('Orientation');
            ylabel('Spatial kernel')
    case 'temporal'
        hold on
end

assert(numel(prior_task_list) == numel(mode_list), 'Number of prior and number of mode do not match')
for n = 1:numel(prior_task_list)
    prior_cardinal = prior_task_list{n}(1);
    prior_oblique = prior_task_list{n}(2);

    idx = find(cell2mat({psyKernel_model(:).prior_cardinal}) == prior_cardinal &...
          cell2mat({psyKernel_model(:).prior_oblique})  == prior_oblique & ...
          strcmp({psyKernel_model(:).image_task}, image_task) & ...
          strcmp({psyKernel_model(:).mode}, mode_list{n}));
    
   
    
    w_ori_all = {psyKernel_model(idx).w_ori};
    w_ori_all = cat(2,w_ori_all{:});
    
    w_time_all = {psyKernel_model(idx).w_time};
    w_time_all = cat(2,w_time_all{:});
    
    
    w_ori_avg = mean(w_ori_all,2);
    w_ori_sem = std(w_ori_all, [], 2) / sqrt(size(w_ori_all, 2));
    
    w_time_avg = mean(w_time_all,2);
    w_time_sem = std(w_time_all, [], 2) / sqrt(size(w_time_all, 2));
    
   
    switch kernel_type
        case 'spatial'
            plot(orientationsDEG, w_ori_avg,'Color',plot_color,'LineWidth',2); hold on
            fill([orientationsDEG, fliplr(orientationsDEG)], [w_ori_avg' - w_ori_sem', fliplr(w_ori_avg' + w_ori_sem')],...
                plot_color, 'FaceAlpha',0.5);

        case 'temporal'
            x = [1:numel(w_time_avg)];
            
            plot(x, w_time_avg,'Color',plot_color,'LineWidth',2);  hold on
            fill([x, fliplr(x)], [w_time_avg' - w_time_sem', fliplr(w_time_avg' + w_time_sem')],...
                plot_color, 'FaceAlpha',0.5);
    end
    set(gca,'fontsize',18);
    box off
end
end