clear all
clc
close all


%% extract data needed for psykernel
doExtract = 0;

organizeData_saveFolder ='/Volumes/T7/dataTransfer/BayesianModel_interleaved_simulation/psyKernel_use_interleaved_moreRepeats_static';
if doExtract

    dataFolder = '/Volumes/T7/dataTransfer/BayesianModel_interleaved_simulation/run_interleaved_bPF_0_80';


    energy_use = 'proj';
    
    filename_list  = dir(fullfile(dataFolder,'*.mat'));
    nFile = numel(filename_list);
    [taskprior, delta, contrast_signed, repeats] = deal(zeros(nFile, 1));
    imagetask  = cell(nFile,1);
    for n = 1:nFile
    % Extract numbers using regular expressions
        tokens = regexp(filename_list(n).name, 'synthetic_data_([a-zA-Z0-9]+)_bPF_([-]?[\d_]+)_taskprior_([\d_]+)_delta_([\d_]+)_contrast_([\d_]+)_([\d_]+)_rep_([\d_]+)', 'tokens');
        % Convert extracted strings back to numbers
        extracted_params        = tokens{1}; % Extract matched tokens
        imagetask_str           = extracted_params{1};
        bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
        taskprior_str           = strrep(extracted_params{3}, '_', '.');
        delta_str               = strrep(extracted_params{4}, '_', '.');
        contrast_1_str          = extracted_params{5};
        contrast_2_str          = extracted_params{6};
        nRep_str                = extracted_params{7};
        
        imagetask{n}                = imagetask_str;
    %    b_PF(n)                    = str2double(bPF_str);
        taskprior(n)                = str2double(taskprior_str);
        delta(n)                    = str2double(delta_str);
        contrast_signed(n)          = str2double(contrast_2_str) - str2double(contrast_1_str);
        repeats(n)                     = str2double(nRep_str); 
    end
    
    imagetask_list          = unique(imagetask); 
    %bPF_list                = unique(b_PF);
    taskprior_list          = unique(taskprior);
    delta_list              = unique(delta); 
    contrast_signed_list    = unique(contrast_signed); 
    repeat_list               = unique(repeats); 
    
    nTask                   = numel(imagetask_list); 
    nPrior                  = numel(taskprior_list);
    %nPF                     = numel(bPF_list); 
    nDelta                  = numel(delta_list); 
    nRep                    = numel(repeat_list);


    if ~isfolder(organizeData_saveFolder)
        mkdir(organizeData_saveFolder)
    end
    timeBin_list{1} = [4:99];
    % T_total = numel(timeBin_list{1});
    % nTimebin = 8;
    % for n = 1:nTimebin
    %     binSize = T_total / nTimebin;
    %     timeBin_list{n+1} = [4 + (n-1) * binSize: n*binSize + 3];
    % end
    
    for i = 1:nTask
        for j = 1:nPrior
            for n = 1:nDelta
                for m = 1:5
                    %%%%% find all sessions of this parameter combination
                    idx = find(strcmp(imagetask, imagetask_list{i}) & ...
                          taskprior == taskprior_list(j) & ...
                          delta == delta_list(n) & ...
                          repeats == repeat_list(m) );
    
                      
                      delta_str             = strrep(sprintf('%.2f', delta_list(n)), '.', '_');
                      task_prior            = strrep(sprintf('%.2f', taskprior_list(j)), '.', '_');
                    
                      session_name          = sprintf('synthData_use_imagetask_%s_bPF_0_80_delta_%s_taskprior_%s_rep_%d',...
                                                imagetask_list{i}, delta_str, task_prior, repeats(m));
                      save_name = fullfile(organizeData_saveFolder, session_name);
                      %%%%% load and extract useful information for each
                      %%%%% session
                      data_psyKernel_use = struct();
                      d = 1;
                      for k = 1:numel(idx)
                          for t = 1:numel(timeBin_list)
                            load(fullfile(dataFolder, filename_list(idx(k)).name));
    
    
    
                   
                            data_psyKernel_use(d).sessionStr             = session_name;  
                            data_psyKernel_use(d).image_task             = dat.Projection.image_task;
                            data_psyKernel_use(d).b_PF                   = dat.Projection.b_PF;
                            data_psyKernel_use(d).learning_strength      = dat.Projection.delta;
                           
                            data_psyKernel_use(d).prior_task             = dat.Projection.prior_task;
                            data_psyKernel_use(d).cardinal_prior         = data_psyKernel_use(d).prior_task(1);
                            data_psyKernel_use(d).oblique_prior          = data_psyKernel_use(d).prior_task(2);
                            data_psyKernel_use(d).timebin                = timeBin_list{t};
                            data_psyKernel_use(d).timebinIndex           = t-1;
                            data_psyKernel_use(d).X_response             = sum(dat.X(:,:,timeBin_list{t}),3);
                            data_psyKernel_use(d).signal                 = dat.Signal;
                            data_psyKernel_use(d).decision               = (dat.O(:,3,end) > 0.5) + 1;
      
                            % stimulus contrast
                            data_psyKernel_use(d).stimulus_contrast      = dat.Projection.stimulus_contrast;
                            data_psyKernel_use(d).contrast               = data_psyKernel_use(d).stimulus_contrast(2)- data_psyKernel_use(d).stimulus_contrast(1);
                            % neurons's preferred orientation
                            data_psyKernel_use(d).phi_x                  = dat.Projection.phi_x;
                            data_psyKernel_use(d).P                      = dat.Projection;
                            d = d+1;
                          end
                      end
                      save(save_name, 'data_psyKernel_use')
                     
                    
                end
            end
        end
    end
end
%% compute kernel
doCompute = 1;
if doCompute
    data_list = dir(fullfile(organizeData_saveFolder, 'synthData*.mat'));
    energy_use = 'proj';
    save_name = '../../results/model_behavior/psyKernel_model_single_largescale_moreRepeats_staticImage';
    for n = 1:numel(data_list)
        fprintf('Running session %d/%d \n', n, numel(data_list))
        load(fullfile(organizeData_saveFolder, data_list(n).name));
        [map_params,orientationsDEG] = util_it.compute_psyKernel_model(data_psyKernel_use, energy_use);
        
    %synthData_use_imagetask_oblique_bPF_0_80_delta_0_05_taskprior_0_60_rep_4
            tokens = regexp(data_list(n).name, 'synthData_use_imagetask_([a-zA-Z0-9]+)_bPF_0_80_delta_([-]?[\d_]+)_taskprior_([-]?[\d_]+)_rep_([\d_]+)', 'tokens');
            % Convert extracted strings back to numbers
            extracted_params        = tokens{1}; % Extract matched tokens
           
            imagetask_str           = extracted_params{1};
            delta_str               = strrep(extracted_params{2}, '_', '.');
            prior_str               = strrep(extracted_params{3}, '_', '.');
            session_num_str         = extracted_params{4};
            % bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
            % delta_str               = strrep(extracted_params{3}, '_', '.');
            % taskprior_str           = strrep(extracted_params{4}, '_', '.');
        
            psyKernel_model(n).mode             = 'Single';
            psyKernel_model(n).image_task       = imagetask_str;
            psyKernel_model(n).delta            = str2double(delta_str); 
            psyKernel_model(n).prior            = str2double(prior_str);
            
            psyKernel_model(n).repeat           = str2double(session_num_str);
    
            psyKernel_model(n).w_ori            = map_params.w_ori;
            psyKernel_model(n).w_time           = map_params.w_time;
            psyKernel_model(n).orientationsDEG  = orientationsDEG;
        
            psyKernel_model(n).energy_use       = energy_use;
    
            % %%%% number of repeats and list of coherence
            % psyKernel_model(n).cohr_list        = unique(abs(cell2mat({data_psyKernel_use(:).contrast})));
            % for k = 1:numel(psyKernel_model(n).cohr_list)
            %     idx = cell2mat({data_psyKernel_use(:).contrast}) == psyKernel_model(n).cohr_list;
            %     psyKernel_model(n).num_repeats(k) = numel(data_psyKernel_use(idx).decision);
            % end
         
    
    end
    save(save_name, 'psyKernel_model');
end
%% visualization
load('../../results/model_behavior/psyKernel_model_single_largescale_moreRepeats_staticImage');
orientationsDEG = psyKernel_model(1).orientationsDEG;
image_task = 'oblique';

%idx_all = strcmp({psyKernel_model(:).image_task}, image_task);

delta = [psyKernel_model(:).delta];
prior = [psyKernel_model(:).prior];

delta_list = unique(delta);
prior_list = unique(prior);

nDelta = numel(delta_list);
nPrior = numel(prior_list);

for i = 1:nDelta
    for j = 1:nPrior
        idx  = strcmp({psyKernel_model(:).image_task}, image_task) & ...
                delta == delta_list(i) & prior == prior_list(j);
        w_ori = [psyKernel_model(idx).w_ori];

        subplot(nDelta, nPrior, (i-1)*nPrior + j)
        plot(orientationsDEG, mean(w_ori,2));

    
    end
end