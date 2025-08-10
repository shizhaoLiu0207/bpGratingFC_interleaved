function results_cross_sizeControl_perCohr = run_organize_cross_fisherinfo_perCohr_sizeControl(data_list)


nSession = numel(data_list);
t = 1;
for n = 1:nSession
    fprintf('Organizing session %d/%d\n',n,nSession)
    load(fullfile(data_list(n).folder, data_list(n).name));

    if isempty(dat_fisher_cross)
        continue
    end

    sessionType_list = unique({dat_fisher_cross(:).sessionType});

    nSample = max(cell2mat({dat_fisher_cross(:).i_subSample}));
    timebin_list = unique(cell2mat({dat_fisher_cross(:).timeWinIndex}));

    for m = 1:numel(sessionType_list)

        for k = 1:numel(timebin_list)
            idx_base = cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1 & strcmp({dat_fisher_cross(:).sessionType}, sessionType_list{m});
            cohr_list = unique(cell2mat({dat_fisher_cross(idx_base).coherence_level}));


            for i = 1:numel(cohr_list)
                idx = idx_base & cell2mat({dat_fisher_cross(:).coherence_level}) == cohr_list(i);
                dat_fisher_cross_sub = dat_fisher_cross(idx);

                results_cross_sizeControl_perCohr(t).sessionStr     = dat_fisher_cross_sub(1).sessionStr;
                results_cross_sizeControl_perCohr(t).sessionType    = dat_fisher_cross_sub(1).sessionType;
                results_cross_sizeControl_perCohr(t).timeWin        = dat_fisher_cross_sub(1).timeWin;
                results_cross_sizeControl_perCohr(t).timeWinIndex   = dat_fisher_cross_sub(1).timeWinIndex;
                results_cross_sizeControl_perCohr(t).nUnit          = dat_fisher_cross_sub(1).N;
                results_cross_sizeControl_perCohr(t).coherence_level = dat_fisher_cross_sub(1).coherence_level;





                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_bc});
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_bc});

                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_bc});
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_bc});

                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_shuffle_bc});
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_shuffle_bc});

                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_shuffle_bc});
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_shuffle_bc});


                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_median = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_median   = median(results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_median  = median(results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_sample);
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_median  = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_sample);

                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_median = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_median   = median(results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_sample);

                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_median = median(results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_sample);
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_median =  median(results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_sample);



                t = t+1;

            end
        end
    end
end

end