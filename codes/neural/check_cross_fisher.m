clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', filter_name);


% load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
% results_all = results_cross_sizeControl;
load(fullfile(saveFolder, 'results_SubsampleOrganized_perCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl_perCohr;
results_all = get_sample_CI_cross(results_all);

%% add behavioral performance
psyKernel_table_rolo         = load('../../results/behavior/Rolo_psyKernel_table_final.mat');
psyKernel_table_gremlin      = load('../../results/behavior/Gremlin_psyKernel_table');
%%%% all interleaved session
session_list_rolo = bpGlobal.rolo.session_list.switching;
session_list_gremlin = bpGlobal.gremlin.session_list.interleaved_real;
for n = 1:numel(results_all)
    sessionStr = results_all(n).sessionStr;
     %%% add learning index values
    if contains(sessionStr,'Ro')
        predllPath = '../../results/behavior/psyKernel_interleaved_crossvalidate/Ro';
        idx = strcmp({psyKernel_table_rolo.psyKernel_table(:).sessionName},sessionStr);
        results_all(n).cardinal_match_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_match_amplitude;
        results_all(n).oblique_match_amplitude   =  psyKernel_table_rolo.psyKernel_table(idx).oblique_match_amplitude;
        results_all(n).cardinal_template_match  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_template_match;
        results_all(n).oblique_template_match  = psyKernel_table_rolo.psyKernel_table(idx).oblique_template_match;
        results_all(n).cardinal_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_amplitude;
        results_all(n).oblique_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).oblique_amplitude;
        
        [fitCpredC_ll,fitOpredO_ll,fitOpredC_ll, fitCpredO_ll, fitOpredC_ll_final, fitCpredO_ll_final] = ...
            read_predll_result({sessionStr}, predllPath); 
        results_all(n).fitCpredC_ll = fitCpredC_ll;
        results_all(n).fitOpredO_ll = fitOpredO_ll;
        results_all(n).fitOpredC_ll = fitOpredC_ll;
        results_all(n).fitCpredO_ll = fitCpredO_ll;
        results_all(n).fitOpredC_ll_final = fitOpredC_ll_final;
        results_all(n).fitCpredO_ll_final = fitCpredO_ll_final;
    elseif contains(sessionStr, 'Gr')
        predllPath = '../../results/behavior/psyKernel_interleaved_crossvalidate/Gr';
        idx = strcmp({psyKernel_table_gremlin.psyKernel_table(:).sessionName},sessionStr);
        results_all(n).cardinal_match_amplitude  = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_match_amplitude;
        results_all(n).oblique_match_amplitude   =  psyKernel_table_gremlin.psyKernel_table(idx).oblique_match_amplitude;
        results_all(n).cardinal_template_match   = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_template_match;
        results_all(n).oblique_template_match    = psyKernel_table_gremlin.psyKernel_table(idx).oblique_template_match;
        results_all(n).cardinal_amplitude        = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_amplitude;
        results_all(n).oblique_amplitude         = psyKernel_table_gremlin.psyKernel_table(idx).oblique_amplitude;

        [fitCpredC_ll,fitOpredO_ll,fitOpredC_ll, fitCpredO_ll, fitOpredC_ll_final,fitCpredO_ll_final] = ...
            read_predll_result({sessionStr}, predllPath); 
        results_all(n).fitCpredC_ll = fitCpredC_ll;
        results_all(n).fitOpredO_ll = fitOpredO_ll;
        results_all(n).fitOpredC_ll = fitOpredC_ll;
        results_all(n).fitCpredO_ll = fitCpredO_ll;
        results_all(n).fitOpredC_ll_final = fitOpredC_ll_final;
        results_all(n).fitCpredO_ll_final = fitCpredO_ll_final;

    end
end
%%


%%% all good interleaved session
session_list_rolo_good = bpGlobal.rolo.session_list.switching_good;
session_list_gremlin_good = bpGlobal.gremlin.session_list.interleaved_good;

%%% separate by block size
session_list_rolo_block = bpGlobal.rolo.session_list.mini_block(:,1);
session_list_rolo_trial = setdiff(session_list_rolo, session_list_rolo_block);

gremlin_blocksize = bpGlobal.gremlin.session_list.interleaved_blockSize;
session_list_gremlin_trial = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) == 0, 1);
session_list_gremlin_block = gremlin_blocksize(cell2mat(gremlin_blocksize(:,2)) > 0, 1);



%%% good session, separate by block size
session_list_rolo_good_block = intersect(session_list_rolo_good, session_list_rolo_block);
session_list_rolo_good_trial = intersect(session_list_rolo_good, session_list_rolo_trial);

session_list_gremlin_good_block = intersect(session_list_gremlin_good, session_list_gremlin_block);
session_list_gremlin_good_trial = intersect(session_list_gremlin_good, session_list_gremlin_trial);

%% deltaI of each session, with same Cov and I_redundancy with cross Cov
save_folder = '../../figures/figures_informal/cross_fisher';
figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,8,6])
fig.plot_crossfisher_session(results_all, session_list_rolo)
sgtitle('Monkey R','fontsize',18,'fontweight','bold');
save_name = fullfile(save_folder, 'crossFisher_timecourse_monkeyR.svg');
print(save_name, '-dsvg');

figure
set(gcf,'Units','inches','position',[0,0,8,6])
fig.plot_crossfisher_session(results_all, session_list_gremlin)
sgtitle('Monkey G','fontsize',18,'fontweight','bold');
save_name = fullfile(save_folder, 'crossFisher_timecourse_monkeyG.svg');
print(save_name, '-dsvg');
%% plot_bar_Info: real, real_cross, shuffle, shuffle_cross  
save_folder = '../../figures/figures_informal/cross_fisher';
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
plotOptions.plotShuffle = false;
figure
ftsize = 20;
set(gcf,'Units','inches','position',[0,0,10,4])
subplot(1,2,1)
%fig.plot_bar_Info(results_all, session_list_rolo_good, plotOptions); title('MonkeyR Good sessions');
fig.plot_bar_Info(results_all, session_list_rolo, plotOptions); title('Monkey R')
subplot(1,2,2)
%fig.plot_bar_Info(results_all, session_list_gremlin_good,plotOptions); title('Monkey G Good sessions');
fig.plot_bar_Info(results_all, session_list_gremlin,plotOptions); title('Monkey G');
save_name = fullfile(save_folder, sprintf('cross_fisher_info_bar.svg'));
print(save_name, '-dsvg');
%% plot_bar_deltaInfo: I_reduandacy with same Cov and I_redundancy with cross Cov
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = true;
figure
ftsize = 14;
set(gcf,'Units','inches','position',[0,0,8,6])
subplot(2,1,1)
fig.plot_bar_deltaInfo(results_all, session_list_rolo, plotOptions); title('Monkey R');
set(gca,'fontsize',20)
% subplot(2,2,1)
% fig.plot_bar_deltaInfo(results_all, session_list_rolo, plotOptions); title('All sessions');
% subplot(2,2,2)
% fig.plot_bar_deltaInfo(results_all, session_list_rolo_good, plotOptions); title('Good sessions');
% subplot(2,2,3)
% fig.plot_bar_deltaInfo(results_all, session_list_rolo_block, plotOptions); title('Blockwise sessions')
% subplot(2,2,4)
% fig.plot_bar_deltaInfo(results_all, session_list_rolo_trial, plotOptions); title('Trial-by-trial sessions')
% % subplot(3,2,5)
% % fig.plot_bar_deltaInfo(dat_fisher_cross, session_list_rolo_good_block);title('Good blockwise')
% % subplot(3,2,6)
% % fig.plot_bar_deltaInfo(dat_fisher_cross, session_list_rolo_good_trial); title('Good trial-by-trial')
%  sgtitle('Monkey R','fontsize', 20)
subplot(2,1,2)
fig.plot_bar_deltaInfo(results_all, session_list_gremlin, plotOptions); title('Monkey G');
set(gca,'fontsize',20)
%save_name = fullfile(save_folder, sprintf('delta_cross_fisher_info_bar.svg'));
%print(save_name, '-dsvg');
% subplot(2,2,1)
% fig.plot_bar_deltaInfo(results_all, session_list_gremlin, plotOptions); title('All sessions');
% subplot(2,2,2)
% fig.plot_bar_deltaInfo(results_all, session_list_gremlin_good, plotOptions); title('Good sessions');
% subplot(2,2,3)
% fig.plot_bar_deltaInfo(results_all, session_list_gremlin_block, plotOptions); title('Blockwise sessions')
% subplot(2,2,4)
% fig.plot_bar_deltaInfo(results_all, session_list_gremlin_trial, plotOptions); title('Trial-by-trial sessions')
% % subplot(3,2,5)
% % fig.plot_bar_deltaInfo(dat_fisher_cross, session_list_gremlin_good_block);title('Good blockwise')
% % subplot(3,2,6)
% % fig.plot_bar_deltaInfo(dat_fisher_cross, session_list_gremlin_good_trial); title('Good trial-by-trial')
% sgtitle('Monkey G','fontsize', 20)

%% over time within a trial
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.nTimebin = 8;
fig.plot_bar_deltaInfo_withinTrial(results_all, session_list_rolo_good, plotOptions);
sgtitle('Monkey R','fontsize', 20)
fig.plot_bar_deltaInfo_withinTrial(results_all, session_list_gremlin_good, plotOptions);
sgtitle('Monkey G','fontsize', 20)