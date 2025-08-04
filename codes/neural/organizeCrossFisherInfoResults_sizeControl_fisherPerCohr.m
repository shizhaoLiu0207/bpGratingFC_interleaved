clear all
clc
close all

%%
versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', versionName);

dataName = fullfile(saveFolder,  sprintf('fisherInfo_cross_direct_all_sessions_%s',versionName));
load(dataName);
saveName = fullfile(saveFolder,'results_SubsampleCombined_perCohr_fisherInfo_all_sessions');

session_list = unique({dat_fisher_cross(:).sessionStr});
nSession = numel(session_list);
t = 1;
for n = 1:nSession
    fprintf('Organizing session %d/%d\n',n,nSession)
    idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n});
    nSample = max(cell2mat({dat_fisher_cross(idx_base).i_subSample}));
    timebin_list = unique(cell2mat({dat_fisher_cross(idx_base).timeWinIndex}));

    for k = 1:numel(timebin_list)
        idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n})  & cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1;
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
save(saveName,'results_cross_sizeControl_perCohr');
