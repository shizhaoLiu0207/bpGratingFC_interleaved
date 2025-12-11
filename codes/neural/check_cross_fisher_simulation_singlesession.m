clear all
clc
%close all
global   bpGlobal  ftsize
bpGratingFCGlobal();
ftsize = 16;
figFolder = '../../figures/figures_informal/Bayesian_model_simulation';
%%
dataFolder = '../../results/neural/fisherInfo_direct/fisherInfo_direct_modelInterleaved_versionControl/subset_32_random_1000';

%%% fisher info results of coherence-combined
dataName = fullfile(dataFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData');
load(dataName);
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);
%%
nSession = numel(results_all);
CI_level = 68;
prc = [(100-CI_level)/2, 100 - (100-CI_level)/2];
[fisher_cardinal_real, fisher_oblique_real, fisher_cardinal_cross, fisher_oblique_cross] = deal(zeros(nSession, 1));
[fisher_cardinal_real_CI, fisher_oblique_real_CI, fisher_cardinal_cross_CI, fisher_oblique_cross_CI] = deal(zeros(nSession, 2));
[b_PF, cardinal_prior, oblique_prior, cardinal_delta, oblique_delta] = deal(zeros(nSession, 1));
for n = 1:numel(results_all)
    fisher_cardinal_real_sample = results_all(n).combine_fisher_cardinal_cardinal_sample;
    fisher_oblique_real_sample =   results_all(n).combine_fisher_oblique_oblique_sample;
    fisher_cardinal_cross_sample = results_all(n).combine_fisher_cardinal_oblique_sample;
    fisher_oblique_cross_sample  = results_all(n).combine_fisher_oblique_cardinal_sample;



    fisher_cardinal_real(n)     = median(fisher_cardinal_real_sample);
    fisher_cardinal_cross(n)    = median(fisher_cardinal_cross_sample);
    fisher_oblique_real(n)      = median(fisher_oblique_real_sample);
    fisher_oblique_cross(n)     = median(fisher_oblique_cross_sample);

    fisher_cardinal_real_CI(n,:) = prctile(fisher_cardinal_real_sample, prc);
    fisher_oblique_real_CI(n,:)  = prctile(fisher_oblique_real_sample, prc);

    fisher_cardinal_cross_CI(n,:) = prctile(fisher_cardinal_cross_sample, prc);
    fisher_oblique_cross_CI(n,:)  = prctile(fisher_oblique_cross_sample, prc);

    tokens = regexp(results_all(n).sessionStr, 'Model_bPF_([-]?[\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params            = tokens{1}; % Extract matched tokens
    b_PF_str                    = strrep(extracted_params{1}, '_', '.'); % Replace _ with .
    cardinal_delta_str          = strrep(extracted_params{2}, '_', '.');
    cardinal_prior_str          = strrep(extracted_params{3}, '_', '.');
    oblique_delta_str           = strrep(extracted_params{4}, '_', '.');
    oblique_prior_str           = strrep(extracted_params{5}, '_', '.');
    
    b_PF(n)                     = str2double(b_PF_str);
    cardinal_delta(n)           = str2double(cardinal_delta_str);
    cardinal_prior(n)           = str2double(cardinal_prior_str);
    oblique_delta(n)            = str2double(oblique_delta_str);
    oblique_prior(n)            = str2double(oblique_prior_str);       
end
[~,i_sort] = sort(b_PF,'ascend');
b_PF = b_PF(i_sort);
cardinal_delta = cardinal_delta(i_sort);
cardinal_prior = cardinal_prior(i_sort);
oblique_delta = oblique_delta(i_sort);
oblique_prior = oblique_prior(i_sort);
fisher_cardinal_cross = fisher_cardinal_cross(i_sort);
fisher_oblique_cross   = fisher_oblique_cross(i_sort);
fisher_cardinal_cross_CI = fisher_cardinal_cross_CI(i_sort,:);
fisher_oblique_cross_CI   = fisher_oblique_cross_CI(i_sort,:);

b_PF_list = unique(b_PF); 
prior_list = unique(cardinal_prior);
delta_list = unique(cardinal_delta);
%%%% show effect of delta
idx(1) = find(cardinal_prior == 1 & oblique_prior == 1 & cardinal_delta == 0.08 & oblique_delta == 0.05 & b_PF == 0);
%%%% show effect of prior
idx(2) = find(cardinal_prior == 1 & oblique_prior == 0.5 & cardinal_delta == 0.08 & oblique_delta == 0.08 & b_PF == 0);
%%%% 
idx(3) = find(cardinal_prior == 1 & oblique_prior == 1 & cardinal_delta == 0.08 & oblique_delta == 0.08 & b_PF == 0.8);
plotOptions = struct();
plotOptions.errorbar = 'SEM_session';
plotOptions.dottest = false;
plotOptions.plotShuffle = false;


figure
set(gcf,'unit','inches','position',[0,0,15,4])
ftsize = 20;
save_folder = '../../figures/figures_informal/cross_fisher';
for i = 1:3
    subplot(1,3,i)
    %fig.plot_bar_Info(results_all, session_list_rolo_good, plotOptions); title('MonkeyR Good sessions');
    fig.plot_bar_Info(results_all, {results_all(idx(i)).sessionStr}, plotOptions);
    if i > 1
        ylabel('')
    end
end
save_name = fullfile(save_folder,'example_newsimulation_fisher.svg');
print(save_name, '-dsvg');
% plotOptions = struct();
% plotOptions.errorbar = 'SEM_session';
% plotOptions.dottest = false;
% plotOptions.plotShuffle = false;
% figure
% ftsize = 20;
% set(gcf,'Units','inches','position',[0,0,10,4])
% 
% %fig.plot_bar_Info(results_all, session_list_rolo_good, plotOptions); title('MonkeyR Good sessions');
% fig.plot_bar_Info(results_all, {results_all(idx).sessionStr}, plotOptions);
% %% show effect of b_PF
% 
% idx = 
% plotOptions = struct();
% plotOptions.errorbar = 'SEM_session';
% plotOptions.dottest = false;
% plotOptions.plotShuffle = false;
% figure
% ftsize = 20;
% set(gcf,'Units','inches','position',[0,0,10,4])
% 
% %fig.plot_bar_Info(results_all, session_list_rolo_good, plotOptions); title('MonkeyR Good sessions');
% fig.plot_bar_Info(results_all, {results_all(idx).sessionStr}, plotOptions);