clear all
clc
close all
%%
versionName_list = {'all_trials_coef1_hVis2_FR1_interleaved_sizeControl'};
for i_version = 1:numel(versionName_list)
    versionName = versionName_list{i_version};% 'all_trials_coef1_hVis2_FR1_sizeControl';
    saveFolder =  sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', versionName);
    saveName = fullfile(saveFolder, 'results_sizeControl_combined_allSessions');
    sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s/results_sizeControl_combined_allSessions', versionName);
    
    
    results_folder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s/individual_sessions/',versionName);
    
    flist = dir(fullfile(results_folder,'*.mat'));
    
    organized_folder = fullfile(saveFolder, 'results_cohr_combined_session');
    %results_sizeControl_combined = struct();
  
    for n = 1:numel(flist)
        k = 1;

        saveName_session = fullfile(organized_folder, sprintf('results_cohrCombined_%s', flist(n).name));
        fprintf('Organizing session %d/%d\n',n,numel(flist))
        load(fullfile(results_folder, flist(n).name));
        if isempty(dat_fisher)
            continue
        end
        nSample = max(cell2mat({dat_fisher(:).i_subSample}));
        %%% find all combinations of sessionType (mainTask/orientation) and
        %%% taskType (cardinal/oblique)
        [comb_field1, comb_field2] =  ndgrid( unique({dat_fisher(:).sessionType}),  unique({dat_fisher(:).task}));
        type_combinations = [comb_field1(:), comb_field2(:)];
    
        timeWin_list = unique(cell2mat({dat_fisher(:).timeWinIndex}));
        for i_t = 1:numel(timeWin_list)
            for i_type = 1:size(type_combinations,1) % loop over all possible types of combinations
                idx_base = strcmp({dat_fisher(:).sessionType}, type_combinations{i_type, 1}) & ...
                            strcmp({dat_fisher(:).task}, type_combinations{i_type, 2}) & ...
                            cell2mat({dat_fisher(:).timeWinIndex}) == timeWin_list(i_t);
                if sum(idx_base) == 0 % skip if this combination does not really appear in the struct
                    continue
                end
                % dat_fisher_combine = struct();
                parfor i = 1:nSample
                    idx = idx_base & cell2mat({dat_fisher(:).i_subSample}) == i;
                    options = struct();
                 
                
                    [dat_fisher_combine(i), options] = combine_coherence_fisherInfo(dat_fisher(idx),options);
                   
                end
        
                if isempty(fieldnames(dat_fisher_combine(1)))
                    clear dat_fisher_combine
                    continue
                end
        
                results_sizeControl_combined(k).sessionStr     = dat_fisher_combine(1).sessionStr;
                results_sizeControl_combined(k).taskType       = dat_fisher_combine(1).taskType;
                results_sizeControl_combined(k).sessionType    = dat_fisher_combine(1).sessionType;
                results_sizeControl_combined(k).timeWin        = dat_fisher_combine(1).timeWin;
                results_sizeControl_combined(k).timeWinIndex   = dat_fisher_combine(1).timeWinIndex;
                results_sizeControl_combined(k).nUnit          = dat_fisher_combine(1).nUnit;
        
        
                results_sizeControl_combined(k).combine_delta_sample = cell2mat({dat_fisher_combine(:).combine_delta});
                results_sizeControl_combined(k).combine_fisher_sample = cell2mat({dat_fisher_combine(:).combine_fisher});
                results_sizeControl_combined(k).combine_shuffle_sample = cell2mat({dat_fisher_combine(:).combine_shuffle});
        
                results_sizeControl_combined(k).combine_delta = mean(results_sizeControl_combined(k).combine_delta_sample);
                results_sizeControl_combined(k).combine_delta_std  = std(results_sizeControl_combined(k).combine_delta_sample);
                results_sizeControl_combined(k).combine_fisher = mean(results_sizeControl_combined(k).combine_fisher_sample);
                results_sizeControl_combined(k).combine_fisher_std  = std(results_sizeControl_combined(k).combine_fisher_sample);
                results_sizeControl_combined(k).combine_shuffle = mean(results_sizeControl_combined(k).combine_shuffle_sample);
                results_sizeControl_combined(k).combine_shuffle_std  = std(results_sizeControl_combined(k).combine_shuffle_sample);
                k  =  k + 1;
                clear dat_fisher_combine
            end
        end
       save(saveName_session, 'results_sizeControl_combined');
       clear results_sizeControl_combined
    end
    
    

    %results_sizeControl_combined = get_sample_CI(results_sizeControl_combined);
    
    %save(saveName,'results_sizeControl_combined');

end
%%
versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder =  sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', versionName);
saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');


individual_folder = fullfile(saveFolder,'results_cohr_combined_session');

file_list = dir(fullfile(individual_folder,'*.mat'));

%results_sizeControl_combined = struct();
for n = 1:numel(file_list)
    results_sizeControl_individual = load(fullfile(individual_folder, file_list(n).name));

    if n == 1
        results_sizeControl_combined = results_sizeControl_individual.results_sizeControl_combined;
    else
    results_sizeControl_combined = [results_sizeControl_combined,...
        results_sizeControl_individual.results_sizeControl_combined];
    end

end

save(saveName,'results_sizeControl_combined');