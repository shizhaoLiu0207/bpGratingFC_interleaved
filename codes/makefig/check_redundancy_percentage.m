clear all
clc
close all
%%
global bpGlobal
bpGratingFCGlobal()

%% redundancy percentage, interleaved epoch
versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
[session_list, results_neural]  = data_prep(versionName);


plotOptions.subjectCode = 'Ro';
plotOptions.task_plot = 'cardinal';
session_list_plot = session_list.rolo_switching;
Y_rolo_cardinal_interleaved = print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Ro';
plotOptions.task_plot = 'oblique';
session_list_plot = session_list.rolo_switching;
Y_rolo_oblique_interleaved = print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Gr';
plotOptions.task_plot = 'cardinal';
session_list_plot = session_list.gremlin_switching;
Y_gremlin_cardinal_interleaved = print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Gr';
plotOptions.task_plot = 'oblique';
session_list_plot = session_list.gremlin_switching;
Y_gremlin_oblique_interleaved = print_pct_mean(results_neural, session_list_plot, plotOptions);

%% redundancy percentage, learning epoch
versionName = 'all_trials_coef1_hVis2_FR1_sizeControl';
[session_list, results_neural]  = data_prep(versionName);


plotOptions.subjectCode = 'Ro';
plotOptions.task_plot = 'cardinal';
session_list_plot = bpGlobal.rolo.session_list.cardinal_late;
Y_rolo_cardinal_learning = print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Ro';
plotOptions.task_plot = 'oblique';
session_list_plot = bpGlobal.rolo.session_list.oblique_late;
Y_rolo_oblique_learning= print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Gr';
plotOptions.task_plot = 'cardinal';
session_list_plot = bpGlobal.gremlin.session_list.cardinal_late;
Y_gremlin_cardinal_learning= print_pct_mean(results_neural, session_list_plot, plotOptions);


plotOptions.subjectCode = 'Gr';
plotOptions.task_plot = 'oblique';
session_list_plot = bpGlobal.gremlin.session_list.oblique_late;
Y_gremlin_oblique_learning = print_pct_mean(results_neural, session_list_plot, plotOptions);

% %% redundancy percentage zero, interleaved epoch
% filter_name = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
% saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);
% 
% 
% load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
% results_all = results_cross_sizeControl;
% results_all = get_sample_CI_cross(results_all);
% 
% plotOptions.subjectCode = 'Ro';
% plotOptions.task_plot = 'cardinal';
% session_list_plot = session_list.rolo_switching;
% Y = print_pct_mean(results_all, session_list_plot, plotOptions);
% 
% 
% plotOptions.subjectCode = 'Ro';
% plotOptions.task_plot = 'oblique';
% session_list_plot = session_list.rolo_switching;
% Y = print_pct_mean(results_all, session_list_plot, plotOptions);
% 
% 
% plotOptions.subjectCode = 'Gr';
% plotOptions.task_plot = 'cardinal';
% session_list_plot = session_list.gremlin_switching;
% Y = print_pct_mean(results_all, session_list_plot, plotOptions);
% 
% 
% plotOptions.subjectCode = 'Gr';
% plotOptions.task_plot = 'oblique';
% session_list_plot = session_list.gremlin_switching;
% Y = print_pct_mean(results_all, session_list_plot, plotOptions);

%%

function Y = print_pct_mean(results_neural, session_list_plot, plotOptions)
idx = ismember({results_neural(:).sessionStr}, session_list_plot) & ...
         strcmp({results_neural(:).sessionType},'mainTask') &  ...
            strcmp({results_neural(:).taskType}, plotOptions.task_plot) & ...
            cell2mat({results_neural(:).timeWinIndex}) == 0;

switch plotOptions.subjectCode
    case 'Ro'
        animal_str = 'monkey R';
    case 'Gr'
        animal_str = 'monkey G';
end


Y = cell2mat({results_neural(idx).combine_delta_percent});
fprintf('%s, %s, redundancy percent: mean = %.2f, std = %.2f \n',animal_str, plotOptions.task_plot, mean(Y), std(Y))

end