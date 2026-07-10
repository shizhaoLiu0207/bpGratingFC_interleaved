clear all
clc
close all

global   bpGlobal  
bpGratingFCGlobal();
%% 1. load simulation results
plot_moreRepeats = true;
filter_folder   = '../../results/filtered_neuron_synthetic';

b_PF            = 0.8;
subject_code    = 'gremlin';


b_PF_str        = strrep(sprintf('%.2f',b_PF),'.','_');

if plot_moreRepeats
    versionName = sprintf('sampled_subset_empirical_%s_nTotal_128_nSample_64_b_PF_%s_random_250_moreRepeats',subject_code, b_PF_str);
else
    versionName = sprintf('sampled_subset_empirical_%s_nTotal_128_nSample_64_b_PF_%s_random_250',subject_code, b_PF_str);
end

save_folder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_modelInterleaved_versionControl/%s',versionName);

load(fullfile(save_folder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData'));

results_all_simulation = get_sample_CI_cross(results_cross_sizeControl);

nSession = numel(results_all_simulation);

[cardinal_delta, cardinal_prior, oblique_delta, oblique_prior] = deal(zeros(nSession,1));
% if plot_moreRepeats
%     repeat = zeros(nSession, 1);
% end
%%%% key measure
%%% (I_cross - I_real) / I_cross; cardinal and oblique
%%%% I_{redundancy, within, pct} - I_{redundancy, cross, pct}; cardinal and
%%%% oblique
[diff_fisher_cardinal, diff_fisher_oblique, diff_redundancy_cardinal, diff_redundancy_oblique] = deal(zeros(nSession,1));
[redundancy_cardinal_within, redundancy_cardinal_cross, redundancy_oblique_within, redundancy_oblique_cross] =  deal(zeros(nSession,1));
for n = 1:nSession
    sessionStr = results_all_simulation(n).sessionStr;
    if plot_moreRepeats
         tokens = regexp(sessionStr, 'Model_bPF_([\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)_rep_([\d_]+)','tokens');
    else
        tokens = regexp(sessionStr, 'Model_bPF_([\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)','tokens');

    end
    extracted_params        = tokens{1}; % Extract matched tokens
    cardinal_delta_str      = strrep(extracted_params{2}, '_', '.'); 
    cardinal_prior_str      = strrep(extracted_params{3}, '_', '.'); 
    oblique_delta_str       = strrep(extracted_params{4}, '_', '.'); 
    oblique_prior_str       = strrep(extracted_params{5}, '_', '.'); 

    cardinal_delta(n)       = str2double(cardinal_delta_str);
    cardinal_prior(n)       = str2double(cardinal_prior_str);
    oblique_delta(n)        = str2double(oblique_delta_str);
    oblique_prior(n)        = str2double(oblique_prior_str);

    % if plot_moreRepeats
    %     repeat(n) = str2double(extracted_params{6});
    % end
    

    

    redundancy_cardinal_within(n)   = results_all_simulation(n).delta_percent_cardinal_cardinal_median;
    redundancy_cardinal_cross(n)    = results_all_simulation(n).delta_percent_cardinal_oblique_median;
    redundancy_oblique_within(n)    = results_all_simulation(n).delta_percent_oblique_oblique_median;
    redundancy_oblique_cross(n)     = results_all_simulation(n).delta_percent_oblique_cardinal_median;

    diff_fisher_cardinal(n) = results_all_simulation(n).diff_info_percent_cardinal_median;
    diff_fisher_oblique(n)  = results_all_simulation(n).diff_info_percent_oblique_median;

    diff_redundancy_cardinal(n)     = results_all_simulation(n).diff_delta_percent_cardinal_median;
    diff_redundancy_oblique(n)      = results_all_simulation(n).diff_delta_percent_oblique_median;
end

%diff_cardinal_oblique_fisher = diff_fisher_cardinal - diff_fisher_oblique;
assymmetry_redundancy = diff_redundancy_cardinal - diff_redundancy_oblique;

assymmetry_delta = cardinal_delta - oblique_delta;
assymmetry_prior = cardinal_prior - oblique_prior;

%%
if plot_moreRepeats
    params = [cardinal_delta, cardinal_prior, oblique_delta, oblique_prior];
    [unique_params, ~, group_id] = unique(params, 'rows');
    
    diff_redundancy_cardinal_reorganize     = nan(size(unique_params,1), 5);
    diff_redundancy_oblique_reorganize      = nan(size(unique_params,1), 5);

    assymmetry_redundancy_reorganize        = nan(size(unique_params,1), 5);
    
    for i = 1:size(unique_params,1)
        idx = find(group_id == i);
        diff_redundancy_cardinal_reorganize(i, :) = diff_redundancy_cardinal(idx);
    end
    
    for i = 1:size(unique_params,1)
        idx = find(group_id == i);
        diff_redundancy_oblique_reorganize(i, :) = diff_redundancy_oblique(idx);
    end
  
    for i = 1:size(unique_params,1)
        idx = find(group_id == i);
        assymmetry_redundancy_reorganize(i, :) = assymmetry_redundancy(idx);
    end


    cardinal_delta_reorganize   = unique_params(:,1);
    cardinal_prior_reorganize   = unique_params(:,2);
    oblique_delta_reorganize    = unique_params(:,3);
    oblique_prior_reorganize    = unique_params(:,4);

    assymmetry_redundancy = diff_redundancy_cardinal_reorganize - diff_redundancy_oblique_reorganize;
    
    assymmetry_delta = cardinal_delta_reorganize - oblique_delta_reorganize;
    assymmetry_prior = cardinal_prior_reorganize - oblique_prior_reorganize;

end
%% 2. load empirical data
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
save_folder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);

load(fullfile(save_folder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all_empirical = get_sample_CI_cross(results_cross_sizeControl);

switch subject_code
    case 'rolo'
        idx_session = contains({results_all_empirical(:).sessionStr},'Ro');
    case 'gremlin'
        idx_session = contains({results_all_empirical(:).sessionStr},'Gr');
end

diff_redundancy_cardinal_empirical  = cell2mat({results_all_empirical(idx_session).diff_delta_percent_cardinal_median});
diff_redundancy_oblique_empirical   = cell2mat({results_all_empirical(idx_session).diff_delta_percent_oblique_median});

diff_fisher_cardinal_empirical      = cell2mat({results_all_empirical(idx_session).diff_info_percent_cardinal_median});
diff_fisher_oblique_empirical       = cell2mat({results_all_empirical(idx_session).diff_info_percent_oblique_median});

diff_cardinal_oblique_redundancy_empirical = diff_redundancy_cardinal_empirical - diff_redundancy_oblique_empirical;

%% match between simulation and empirical data - average

if plot_moreRepeats
        dist_feature = sqrt((mean(diff_redundancy_cardinal_reorganize, 2) - mean(diff_redundancy_cardinal_empirical)).^ 2 + ...
                (mean(diff_redundancy_oblique_reorganize, 2) - mean(diff_redundancy_oblique_empirical)).^ 2 );
        
        dist_feature_cardinal   = abs(mean(diff_redundancy_cardinal_reorganize, 2) - mean(diff_redundancy_cardinal_empirical));
        dist_feature_oblique    = abs(mean(diff_redundancy_oblique_reorganize, 2) - mean(diff_redundancy_oblique_empirical));

        cardinal_delta_plot = cardinal_delta_reorganize;
        cardinal_prior_plot = cardinal_prior_reorganize;
        oblique_delta_plot = oblique_delta_reorganize;
        oblique_prior_plot = oblique_prior_reorganize;
else
        dist_feature = sqrt((diff_redundancy_cardinal - mean(diff_redundancy_cardinal_empirical)).^ 2 + ...
                    (diff_redundancy_oblique - mean(diff_redundancy_oblique_empirical)).^ 2 );

        cardinal_delta_plot = cardinal_delta;
        cardinal_prior_plot = cardinal_prior;
        oblique_delta_plot = oblique_delta;
        oblique_prior_plot = oblique_prior;

end



plot_option = 'oblique'; % 'combined', 'cardinal', 'oblique'
switch plot_option
    case 'combined'
        dist_plot = dist_feature;
        mapName = 'turbo';
    case 'cardinal'
        dist_plot = dist_feature_cardinal;
        mapName = 'hot';
    case 'oblique'
        dist_plot = dist_feature_oblique;
        mapName = 'cool';
end
rgb_colors = fig_it.vals2colormap(dist_plot, mapName);
figure
plotOptions.mapName = mapName;
plotOptions.max_min = 'min';
set(gcf, 'Units','inches','Position',[0,0,12,5])
save_folder = '../../figures/figures_final/model_fisher_largescale';
save_name = fullfile(save_folder ,sprintf('scatter_params_match_neural_monkey_%s_%s.svg',subject_code, plot_option));

ax_1 = subplot(1,4,1);
set(ax_1,'position',get(ax_1,'position') + [-0.065, 0.1, 0.02, -0.15])
p1 = cardinal_delta_plot;
p2 = cardinal_prior_plot;
plotOptions.p1_name = '\color{red}\delta_{cardinal}';
plotOptions.p2_name = '\color{red}prior_{cardinal}';

fig_it.make_scatter_match(p1, p2, dist_plot, rgb_colors, plotOptions);
xlim([0.03, 0.09]); ylim([0.4,1.1]);


ax_2 = subplot(1,4,2);
set(ax_2,'position',get(ax_2,'position') + [-0.025, 0.1, 0.02, -0.15])
p1 = oblique_delta_plot;
p2 = oblique_prior_plot;
plotOptions.p1_name = '\color{blue}\delta_{oblique}';
plotOptions.p2_name = '\color{blue}prior_{oblique}';
fig_it.make_scatter_match(p1, p2, dist_plot, rgb_colors, plotOptions);
xlim([0.03, 0.09]); ylim([0.4,1.1]);

ax_3 = subplot(1,4,3);
set(ax_3, 'position',get(ax_3, 'position')+[0.025, 0.1, 0.02, -0.15])
p1 = cardinal_delta_plot;
p2 = oblique_delta_plot;
plotOptions.p1_name = '\color{red}\delta_{cardinal}';
plotOptions.p2_name = '\color{blue}\delta_{oblique}';
fig_it.make_scatter_match(p1, p2, dist_plot, rgb_colors, plotOptions);
xlim([0.03, 0.09]); ylim([0.03, 0.09]);


ax_4 = subplot(1,4,4);
set(ax_4, 'position',get(ax_4, 'position')+[0.07, 0.1, 0.02, -0.15])
p1 = cardinal_prior_plot;
p2 = oblique_prior_plot;
plotOptions.p1_name = '\color{red}prior_{cardinal}';
plotOptions.p2_name = '\color{blue}prior_{oblique}';
fig_it.make_scatter_match(p1, p2, dist_plot, rgb_colors, plotOptions);
xlim([0.4, 1.1]); ylim([0.4,1.1]);

switch subject_code
    case 'rolo'
        sgtitle('Monkey R','fontsize',25,'fontweight','bold')
    case 'gremlin'
        sgtitle('Monkey G','fontsize',25,'fontweight','bold')
end

print(save_name,'-dsvg');

%% plot asymmetry between two tasks as a function of assymetry 
p1 = assymmetry_delta;
p2 = assymmetry_prior;

y_plot = mean(assymmetry_redundancy_reorganize, 2);

plotOptions.max_min = 'max_abs';
plotOptions.p1_name = '\color{red}\delta_{cardinal}  \color{black}- \color{blue}\delta_{oblique}';
plotOptions.p2_name = '\color{red}prior_{cardinal}  \color{black}- \color{blue}prior_{oblique}';
plotOptions.mapName =   'turbo';  
rgb_colors = fig_it.vals2colormap(y_plot, plotOptions.mapName);
plotOptions.doContour = true;
plotOptions.empirical_value = mean(diff_cardinal_oblique_redundancy_empirical);


save_name = fullfile(save_folder ,sprintf('scatter_params_match_neural_assymetry_monkey_%s.svg',subject_code));

figure;
set(gcf, 'unit','inches','position',[0,0,6,5])
fig_it.make_scatter_match(p1, p2, y_plot, rgb_colors, plotOptions);
title('Asymmetry in $I_\textrm{redundancy}$ (\%)','Interpreter','latex')

xlim([-0.05, 0.05]); ylim([-0.6,0.6]);
print(save_name,'-dsvg');
%% Also save the match result as a struct for future use
if plot_moreRepeats
    neural_match_result.dist_feature    = dist_feature;
    neural_match_result.delta_cardinal  = cardinal_delta_plot;
    neural_match_result.prior_cardinal  = cardinal_prior_plot;
    neural_match_result.delta_oblique   = oblique_delta_plot;
    neural_match_result.prior_oblique   = oblique_prior_plot; 
    switch subject_code
        case 'rolo'
            data_save_name = fullfile('../../results','model_match_neural_rolo_moreRepeats');
            neural_match_result_rolo = neural_match_result;
            save(data_save_name, 'neural_match_result_rolo');
        case 'gremlin'
            data_save_name = fullfile('../../results','model_match_neural_gremlin_moreRepeats');
            neural_match_result_gremlin = neural_match_result;
            save(data_save_name, 'neural_match_result_gremlin');
    end
end
%%

%function make_scatter_match(p1, p2, p1_name, p2_name, dist_feature, rgb_colors, mapName)



%% visualize scatter plot


% %% visualize the heat map
% 
% 
% feature_target_str_list = {'diff_redundancy_cardinal';'diff_redundancy_oblique'};
% 
% feature_title_list = {'\color{red}Diff. Redundancy_{cardinal}';'\color{blue}Diff. Redundancy_{oblique}'};
% 
% for k = 1:2
%     plotOptions.ftsize = 12;
%     plotOptions.feature_title = feature_title_list{k};
%     plotOptions.doPeak = true;
%     eval(sprintf('plotOptions.f_obs  =  mean(%s_empirical);', feature_target_str_list{k}));
%     eval(sprintf('feature_target = %s;',feature_target_str_list{k}));
%     figure
%     set(gcf,'Units','inches','position',[0,0,12,5])
%     fig_it.make_heatmap_4D(feature_target, cardinal_delta,cardinal_prior, oblique_delta,oblique_prior, plotOptions)
%     save_folder = '../../figures/figures_final/model_fisher_largescale';
%     save_name = fullfile(save_folder,sprintf('heatmap_4D_%s_%s_matched.svg',feature_target_str_list{k}, subject_code));
%     print(save_name,'-dsvg')
% end
% 
% %% visualize the heat map between difference in variables
% f = diff_cardinal_oblique_redundancy;
% p1 = diff_cardinal_oblique_delta;
% p2 = diff_cardinal_oblique_prior;
% 
% p1_name = '\color{red}\delta_{cardinal} \color{black}- \color{blue}\delta_{oblique}';
% p2_name =  '\color{red}prior_{cardinal} \color{black}- \color{blue}prior_{oblique}';
% figure
% plotOptions = struct();
% plotOptions.ftsize = 12;
% 
% plotOptions.doPeak = false;
% plotOptions.f_obs = mean(diff_cardinal_oblique_redundancy_empirical);
% set(gcf,'Units','inches','position',[0,0,5,4])
% fig_it.make_heatmap(p1,p2,f,p1_name,p2_name, plotOptions);
% 
% set(gca,'fontsize',14)
% title('\color{red}Diff. Redundancy_{cardinal} - \color{blue}Diff. Redundancy_{oblique}','Interpreter','tex',...
%     'FontWeight','bold')
% 
% save_name = fullfile(save_folder,sprintf('heatmap_diff_diff_redundancy_%s_matched.svg',subject_code));
% print(save_name,'-dsvg','-vector')