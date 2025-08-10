function results_cross_sizeControl = run_organize_cross_fisherinfo_sizeControl(dat_input)

if isfield(dat_input, 'sessionStr')
   dat_fisher_cross = dat_input;
   session_list = unique({dat_fisher_cross(:).sessionStr});
   nSession = numel(session_list);
   load_flag = false;
elseif isfield(dat_input, 'folder')
    nSession = numel(file_name_list);
    load_flag = true;
end
t = 1;
for n = 1:nSession
    fprintf('Organizing session %d/%d\n',n,nSession)
    if load_flag
        load(fullfile(file_name_list(n).folder, file_name_list(n).name));
    
        if isempty(dat_fisher_cross)
            continue
        end
         nSample = max(cell2mat({dat_fisher_cross(:).i_subSample}));
         sessionType_list = unique({dat_fisher_cross(:).sessionType});
         timebin_list = unique(cell2mat({dat_fisher_cross(:).timeWinIndex}));
    else
        idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n});
        nSample = max(cell2mat({dat_fisher_cross(idx_base).i_subSample}));
        sessionType_list = unique({dat_fisher_cross(idx_base).sessionType});
        timebin_list = unique(cell2mat({dat_fisher_cross(idx_base).timeWinIndex}));
    end
   
    
    for m = 1:numel(sessionType_list)
   
        for k = 1:numel(timebin_list)
            if load_flag
                idx_base = cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1 & ...
                    strcmp({dat_fisher_cross(:).sessionType}, sessionType_list{m});
            else
                idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n}) & ...
                        cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1 & ...
                        strcmp({dat_fisher_cross(:).sessionType}, sessionType_list{m});
            end

             for i = 1:nSample
                    idx = idx_base & cell2mat({dat_fisher_cross(:).i_subSample}) == i;
                    options = struct();
                  
                    [dat_fisher_cross_combine(i), options] = combine_coherence_crossfisherInfo(dat_fisher_cross(idx),options);
                end
        
                if isempty(fieldnames(dat_fisher_cross_combine(1)))
                    clear dat_fisher_cross_combine
                    continue
                end


                results_cross_sizeControl(t).sessionStr     = dat_fisher_cross_combine(1).sessionStr;
                results_cross_sizeControl(t).sessionType    = dat_fisher_cross_combine(1).sessionType;
                results_cross_sizeControl(t).timeWin        = dat_fisher_cross_combine(1).timeWin;
                results_cross_sizeControl(t).timeWinIndex   = dat_fisher_cross_combine(1).timeWinIndex;
                results_cross_sizeControl(t).nUnit          = dat_fisher_cross_combine(1).nUnit;
   


                results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal});
                results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique});
               
                results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal});
                results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique});

                results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal_shuffle});
                results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique_shuffle});
               
                results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal_shuffle});
                results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique_shuffle});

    
                results_cross_sizeControl(t).fisher_cardinal_cardinal_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample);
                results_cross_sizeControl(t).fisher_oblique_oblique_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample);
                results_cross_sizeControl(t).fisher_oblique_cardinal_median  = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample);
                results_cross_sizeControl(t).fisher_cardinal_oblique_median  = median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample);
    
                results_cross_sizeControl(t).fisher_cardinal_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample);
                results_cross_sizeControl(t).fisher_oblique_oblique_shuffle_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample);

                results_cross_sizeControl(t).fisher_oblique_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample);
                results_cross_sizeControl(t).fisher_cardinal_oblique_shuffle_median =  median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample);
               

                t = t+1;

        end
    end
   
end

end