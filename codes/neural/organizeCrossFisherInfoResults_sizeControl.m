clear all
clc
close all

%%
versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', versionName);

dataName = fullfile(saveFolder,  sprintf('fisherInfo_cross_direct_all_sessions_%s',versionName));
load(dataName);
saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');

results_cross_sizeControl = util_it.run_organize_cross_fisherinfo_sizeControl(dat_fisher_cross);

save(saveName,'results_cross_sizeControl');

% session_list = unique({dat_fisher_cross(:).sessionStr});
% nSession = numel(session_list);
% t = 1;
% for n = 1:nSession
%     fprintf('Organizing session %d/%d\n',n,nSession)
%     idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n});
%     nSample = max(cell2mat({dat_fisher_cross(idx_base).i_subSample}));
%     timebin_list = unique(cell2mat({dat_fisher_cross(idx_base).timeWinIndex}));
% 
% 
%         for k = 1:numel(timebin_list)
% 
%             idx_base = strcmp({dat_fisher_cross(:).sessionStr}, session_list{n})  & cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1;
% 
% 
%              for i = 1:nSample
%                     idx = idx_base & cell2mat({dat_fisher_cross(:).i_subSample}) == i;
%                     options = struct();
% 
%                     [dat_fisher_cross_combine(i), options] = combine_coherence_crossfisherInfo(dat_fisher_cross(idx),options);
%                 end
% 
%                 if isempty(fieldnames(dat_fisher_cross_combine(1)))
%                     clear dat_fisher_cross_combine
%                     continue
%                 end
% 
% 
%                 results_cross_sizeControl(t).sessionStr     = dat_fisher_cross_combine(1).sessionStr;
%                 results_cross_sizeControl(t).sessionType    = dat_fisher_cross_combine(1).sessionType;
%                 results_cross_sizeControl(t).timeWin        = dat_fisher_cross_combine(1).timeWin;
%                 results_cross_sizeControl(t).timeWinIndex   = dat_fisher_cross_combine(1).timeWinIndex;
%                 results_cross_sizeControl(t).nUnit          = dat_fisher_cross_combine(1).nUnit;
% 
% 
% 
%                 results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal});
%                 results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique});
% 
%                 results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal});
%                 results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique});
% 
%                 results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal_shuffle});
%                 results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique_shuffle});
% 
%                 results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal_shuffle});
%                 results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique_shuffle});
% 
% 
%                 results_cross_sizeControl(t).fisher_cardinal_cardinal_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample);
%                 results_cross_sizeControl(t).fisher_oblique_oblique_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample);
%                 results_cross_sizeControl(t).fisher_oblique_cardinal_median  = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample);
%                 results_cross_sizeControl(t).fisher_cardinal_oblique_median  = median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample);
% 
%                 results_cross_sizeControl(t).fisher_cardinal_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample);
%                 results_cross_sizeControl(t).fisher_oblique_oblique_shuffle_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample);
% 
%                 results_cross_sizeControl(t).fisher_oblique_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample);
%                 results_cross_sizeControl(t).fisher_cardinal_oblique_shuffle_median =  median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample);
% 
% 
%                 t = t+1;
% 
%         end
% 
% 
% end

