clear all
clc
close all
%%
global   bpGlobal  ftsize
bpGratingFCGlobal();
filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_%s', filter_name);


load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
results_all = results_cross_sizeControl;
results_all = get_sample_CI_cross(results_all);

idx_use = strcmp({results_all(:).sessionType},'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0;
results_all = results_all(idx_use);
%%  correlation between neural and behavior
%%%% add learning index to results
behav_results_folder = '/Users/liushizhao/Documents/projects/bpGratingEx/results/behavior';
psyKernel_table_rolo         = load(fullfile(behav_results_folder, 'Rolo_psyKernel_table_final.mat'));
psyKernel_table_gremlin      = load(fullfile(behav_results_folder, 'Gremlin_psyKernel_table'));
%%%% all interleaved session
session_list_rolo = bpGlobal.rolo.session_list.switching;
session_list_gremlin = bpGlobal.gremlin.session_list.interleaved_real;

%%%

for n = 1:numel(results_all)
    sessionStr = results_all(n).sessionStr;
     %%% add learning index values
    if contains(sessionStr,'Ro')
        predllPath = fullfile(behav_results_folder, 'psyKernel_interleaved_crossvalidate/Ro');

        idx = strcmp({psyKernel_table_rolo.psyKernel_table(:).sessionName},sessionStr);
        results_all(n).cardinal_match_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_match_amplitude;
        results_all(n).oblique_match_amplitude   =  psyKernel_table_rolo.psyKernel_table(idx).oblique_match_amplitude;
        results_all(n).cardinal_template_match  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_template_match;
        results_all(n).oblique_template_match  = psyKernel_table_rolo.psyKernel_table(idx).oblique_template_match;
        results_all(n).cardinal_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).cardinal_amplitude;
        results_all(n).oblique_amplitude  = psyKernel_table_rolo.psyKernel_table(idx).oblique_amplitude;
        
        [fitCpredC_ll,fitOpredO_ll,fitOpredC_ll, fitCpredO_ll, fitOpredC_ll_final, fitCpredO_ll_final] = ...
            read_predll_result({sessionStr}, predllPath); 
        results_all(n).fitCpredC_ll = mean(fitCpredC_ll);
        results_all(n).fitOpredO_ll = mean(fitOpredO_ll);
        results_all(n).fitOpredC_ll = mean(fitOpredC_ll_final);
        results_all(n).fitCpredO_ll = mean(fitCpredO_ll_final);

     
    elseif contains(sessionStr, 'Gr')
        predllPath = fullfile(behav_results_folder, 'psyKernel_interleaved_crossvalidate/Gr');
       
        idx = strcmp({psyKernel_table_gremlin.psyKernel_table(:).sessionName},sessionStr);
        results_all(n).cardinal_match_amplitude  = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_match_amplitude;
        results_all(n).oblique_match_amplitude   =  psyKernel_table_gremlin.psyKernel_table(idx).oblique_match_amplitude;
        results_all(n).cardinal_template_match   = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_template_match;
        results_all(n).oblique_template_match    = psyKernel_table_gremlin.psyKernel_table(idx).oblique_template_match;
        results_all(n).cardinal_amplitude        = psyKernel_table_gremlin.psyKernel_table(idx).cardinal_amplitude;
        results_all(n).oblique_amplitude         = psyKernel_table_gremlin.psyKernel_table(idx).oblique_amplitude;

        [fitCpredC_ll,fitOpredO_ll,fitOpredC_ll, fitCpredO_ll, fitOpredC_ll_final,fitCpredO_ll_final] = ...
            read_predll_result({sessionStr}, predllPath); 
        results_all(n).fitCpredC_ll = mean(fitCpredC_ll);
        results_all(n).fitOpredO_ll = mean(fitOpredO_ll);
        results_all(n).fitOpredC_ll = mean(fitOpredC_ll_final);
        results_all(n).fitCpredO_ll = mean(fitCpredO_ll_final);


    end
end
%%% corr(I_redudancy(within), learning index (within))

%%% corr(I_cross - I_real, learning index (within))
%%% corr(I_redundacy(within) - I_redundacy(cross), learning index (within))

%%% corr(I_cross - I_real, learning index (within - cross))
%%% corr(I_redundacy(within) - I_redundacy(cross), learning index (within - cross)
%%
clear  results_neural results_behav
idx_rolo    = contains({results_all(:).sessionStr},'Ro');
idx_gremlin = contains({results_all(:).sessionStr},'Gr');
sessionStr_all = {results_all(:).sessionStr};

%%%%% all neural measure

results_neural.I_cardinal          = cell2mat({results_all(:).fisher_cardinal_cardinal_median})';
I_cardinal_cross    = cell2mat({results_all(:).fisher_cardinal_oblique_median})';
results_neural.I_oblique           = cell2mat({results_all(:).fisher_oblique_oblique_median})';
I_oblique_cross     = cell2mat({results_all(:).fisher_oblique_cardinal_median})';

results_neural.I_redundancy_cardinal = cell2mat({results_all(:).delta_cardinal_cardinal_median})';
I_redundancy_cardinal_cross = cell2mat({results_all(:).delta_cardinal_oblique_median})';
results_neural.I_redundancy_oblique = cell2mat({results_all(:).delta_oblique_oblique_median})';
I_redundancy_oblique_cross  = cell2mat({results_all(:).delta_oblique_cardinal_median})';

results_neural.I_redundancy_percent_cardinal    = cell2mat({results_all(:).delta_percent_cardinal_cardinal_median})';
I_redundancy_percent_cardinal_cross     = cell2mat({results_all(:).delta_percent_cardinal_oblique_median})';
results_neural.I_redundancy_percent_oblique     = cell2mat({results_all(:).delta_percent_oblique_oblique_median})';
I_redundancy_percent_oblique_cross      = cell2mat({results_all(:).delta_percent_oblique_cardinal_median})';

results_neural.I_cardinal_diff     = I_cardinal_cross - results_neural.I_cardinal;
results_neural.I_oblique_diff      = I_oblique_cross - I_oblique;

results_neural.I_redundacy_cardinal_diff   = results_neural.I_redundancy_cardinal - I_redundancy_cardinal_cross;
results_neural.I_redundacy_oblique_diff    = results_neural.I_redundancy_oblique - I_redundancy_oblique_cross;

results_neural.I_redundacy_percent_cardinal_diff   = results_neural.I_redundancy_percent_cardinal - I_redundancy_percent_cardinal_cross;
results_neural.I_redundacy_percent_oblique_diff    = results_neural.I_redundancy_percent_oblique - I_redundancy_percent_oblique_cross;

%%% all behav measure
results_behav.learning_cardinal     = cell2mat({results_all(:).cardinal_match_amplitude})';
results_behav.learning_oblique      = cell2mat({results_all(:).oblique_match_amplitude})';

%results_behav.learning_diff         = results_behav.learning_cardinal - results_behav.learning_oblique;

% results_behav.match_cardinal        = cell2mat({results_all(:).cardinal_template_match})';
% results_behav.match_oblique         = cell2mat({results_all(:).oblique_template_match})';
% 
% results_behav.match_diff            = results_behav.match_cardinal - results_behav.match_oblique;
% 
% results_behav.amplitude_cardinal    = cell2mat({results_all(:).cardinal_amplitude})';
% results_behav.amplitude_oblique     = cell2mat({results_all(:).oblique_amplitude})';
% results_behav.amplitude_diff        = results_behav.amplitude_cardinal - results_behav.amplitude_oblique;
% 
% results_behav.fitCpredC             = cell2mat({results_all(:).fitCpredC_ll})'; 
% results_behav.fitOpredC             = cell2mat({results_all(:).fitOpredC_ll})'; 
% results_behav.predC_diff            = results_behav.fitCpredC - results_behav.fitOpredC;
% results_behav.fitOpredO             = cell2mat({results_all(:).fitOpredO_ll})'; 
% results_behav.fitCpredO             = cell2mat({results_all(:).fitCpredO_ll})';
% results_behav.predO_diff            = results_behav.fitOpredO - results_behav.fitCpredO;

results_behav.predC_diff  =  cell2mat({results_all(:).fitCpredC_ll})'  - cell2mat({results_all(:).fitOpredC_ll})'; 
results_behav.predO_diff  =  cell2mat({results_all(:).fitOpredO_ll})'  - cell2mat({results_all(:).fitCpredO_ll})'; 


table_neural    = struct2table(results_neural);
tabel_behav     = struct2table(results_behav);

varNames = table_neural.Properties.VariableNames;

table_neural_norm = table_neural;  % create a copy to store normalized data
subjectID_list = {'Ro';'Gr'};
for i = 1:numel(varNames)
    col = table_neural.(varNames{i});
    norm_col = zeros(size(col));
    
    for s = subjectID_list  % for each subject
        idx = contains(sessionStr_all, s);  % logical indices for subject s
        data = col(idx);
        
        % Normalize: (x - mean) / std
        norm_col(idx) = (data - mean(data)) / std(data);
    end
    
    table_neural_norm.(varNames{i}) = norm_col;
end

varNames = tabel_behav.Properties.VariableNames;
table_behav_norm = tabel_behav;  % create a copy to store normalized data
subjectID_list = {'Ro';'Gr'};
for i = 1:numel(varNames)
    col = tabel_behav.(varNames{i});
    norm_col = zeros(size(col));
    
    for s = subjectID_list  % for each subject
        idx = contains(sessionStr_all, s);  % logical indices for subject s
        data = col(idx);
        
        % Normalize: (x - mean) / std
        norm_col(idx) = (data - mean(data)) / std(data);
    end
    
    table_behav_norm.(varNames{i}) = norm_col;
end
%% 


% Convert tables to matrices
X = table2array(table_neural_norm);  % size: [n × p]
Y = table2array(table_behav_norm);   % size: [n × q]

% Compute correlation coefficients and p-values
[R, P] = corr(X, Y);  % Pearson correlation



% Mask non-significant correlations
R_masked = R;
R_masked(P >= 0.05) = NaN;

% Get variable names for labeling axes
neuralVars = table_neural_norm.Properties.VariableNames;
behavVars = table_behav_norm.Properties.VariableNames;

% Plot the heatmap
figure;
set(gcf,'Units','inches','Position',[0,0,12,10])
h = imagesc(R_masked);  % display correlation matrix
colorbar;    % show color scale
caxis([-1 1]);  % set color limits

% Use AlphaData to make NaNs transparent (shows background = white)
set(h, 'AlphaData', ~isnan(R_masked));

% Set figure background to white (in case transparency shows through)
set(gca, 'Color', 'w');

% Set axis ticks and labels
xticks(1:length(behavVars));
xticklabels(behavVars);
xtickangle(45);

yticks(1:length(neuralVars));
yticklabels(neuralVars);
set(gca,'TickLabelInterpreter','none')
% Add title
title('Pairwise Correlations: Neural vs Behavioral Variables');

% Annotate significant values only
for i = 1:size(R,1)
    for j = 1:size(R,2)
        if ~isnan(R_masked(i,j))
            text(j, i, sprintf('r = %.2f, p=%.2f', R_masked(i,j), P(i,j)), ...
                'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 14);
        end
    end
end
set(gca,'fontsize',16)


save_folder = '../../figures/figures_informal';
save_name = fullfile(save_folder,'pairwise_correlation_neural_behav.svg');
print(save_name,'-dsvg','-vector')
%%
