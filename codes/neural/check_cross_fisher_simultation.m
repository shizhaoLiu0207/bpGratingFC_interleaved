clear all
clc
close all
%%
global   bpGlobal ftsize
ftsize = 16;
bpGratingFCGlobal();
saveFolder = '../../results/neural/fisherInfo_direct/fisherInfo_direct_model_versionControl';

versionName = 'subset_32_random';

% load(fullfile(saveFolder, versionName, 'results_SubsampleCombined_perCohr_fisherInfo_all_sessions_syntheticData'));
% results_all = results_cross_sizeControl_perCohr;
load(fullfile(saveFolder, versionName, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData'));
results_all = results_cross_sizeControl;
nStage = 4;

for i = 1:nStage
  
        session_list{i} = sprintf('Model_session%d', i);

end
session_list = session_list([1,3,4]);
nStage = numel(session_list);
results_all = get_sample_CI_cross(results_all);
%% bar plot of each type of fisher information
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = 0;
plotOptions.plotShuffle = 0;
save_folder = '../../figures/figures_informal/cross_fisher';
figure
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
%fig.plot_bar_Info(results_cross_sizeControl_perCohr, session_list(i)); title(sprintf('Learning stage %d',i));
fig.plot_bar_Info(results_all, session_list(1), plotOptions); title(sprintf('Before learning'));
set(gca,'fontsize',20)
ylim([0,0.1])


subplot(1,2,2)
%fig.plot_bar_Info(results_cross_sizeControl_perCohr, session_list(i)); title(sprintf('Learning stage %d',i));
fig.plot_bar_Info(results_all, session_list(3), plotOptions); title(sprintf('After learning'));
set(gca,'fontsize',20)
save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar_synthetic.svg'));
print(save_name, '-dsvg');
ylim([0,0.1])
%% plot I_redundancy (within and cross task) 
plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.dottest = 0;
figure
for i = 1:nStage
    subplot(1,3,i)
    fig.plot_bar_deltaInfo(results_all, session_list(i), plotOptions); title(sprintf('Learning stage %d',i));
    ylim([0, 0.25])
    set(gca,'fontsize', 24)
end
%% plot I_redundancy as a function of time within a trial


plotOptions = struct();
plotOptions.errorbar = 'CI_sample';
plotOptions.nTimebin = 8;
fig.plot_bar_deltaInfo_withinTrial(results_all, session_list(end), plotOptions);

% %% plot fisher information as a function of learning stage
% nStage = 4; % four learning stages
% nWin = 9; % number of time windows. 0: the whole trial, 1-8: 8 timebins
% [fisher_real, fisher_shuffle, fisher_cross, fisher_cross_corr, fisher_shuffle_cross] = deal(cell(nStage, nWin));
% [cross_diff, cross_diff_corr, cross_diff_shuffle] = deal(cell(nStage, nWin));
% for n  = 1:nStage
%     for t = 1:nWin
%         idx = contains({dat_fisher_cross(:).sessionStr},sprintf('Model_session%d', n)) & ...
%                 strcmp({dat_fisher_cross(:).sessionType}, 'mainTask') & cell2mat({dat_fisher_cross(:).timeWinIndex}) == t-1;
% 
%         % only use cardinal, because oblique is totally symmetric
%         fisher_real{n,t}            = cell2mat({dat_fisher_cross(idx).fisher_cardinal_cardinal_bc});
%         fisher_shuffle{n,t}         = cell2mat({dat_fisher_cross(idx).fisher_cardinal_cardinal_shuffle_bc});
%         fisher_cross{n,t}           = cell2mat({dat_fisher_cross(idx).fisher_cardinal_oblique_bc});
%         fisher_shuffle_cross{n,t}   = cell2mat({dat_fisher_cross(idx).fisher_cardinal_oblique_shuffle_bc});
% 
%         cross_diff{n,t}  = fisher_cross{n,t} - fisher_real{n,t};
% 
%         cross_diff_shuffle{n,t}  =fisher_shuffle_cross{n,t} - fisher_cross{n,t};
% 
%     end
% end
% nItem = 48;
% figure;
% errorbar([1:nStage], cellfun(@mean, cross_diff(:,1)), cellfun(@std, cross_diff(:,1)) ./ sqrt(nItem) ,'LineWidth',2); hold on
% errorbar([1:nStage], cellfun(@mean, cross_diff_shuffle(:,1)), cellfun(@std, cross_diff_shuffle(:,1)) ./ sqrt(nItem), 'LineWidth',2);
% xlim([0.5, nStage + 0.5])
% set(gca,'fontsize', 18, 'xtick',[1:nStage]);
% xlabel('Learning stage'); ylabel('Linear Fisher information')
% 
% legend('I_{cross} - I_{within}', 'I_{cross,shuffle} - I_{shuffle}')
% %%
% 
% figure
% for i = 1:nStage
%     subplot(2,2,i)
%     fig.plot_bar_deltaInfo(dat_fisher_cross, session_list(i,:)); title(sprintf('Learning stage %d',i));
% end
% 
% 
% 
% %%
% figure
% for i = 1:nStage
%     subplot(2,2,i)
%     fig.plot_bar_crossInfo(dat_fisher_cross, session_list(i,:)); title(sprintf('Learning stage %d',i));
% end

