clear all
clc
close all
%%
filter_name_highC_highO = 'all_trials_highC_highO_passiveViewing_sizeControl';
filter_name_lowC_lowO = 'all_trials_lowC_lowO_passiveViewing_sizeControl';
filter_name_highC_lowO = 'all_trials_highC_lowO_passiveViewing_sizeControl';
filter_name_lowC_highO = 'all_trials_lowC_highO_passiveViewing_sizeControl';



results_all_highC_highO = load_data(filter_name_highC_highO);
results_all_lowC_lowO = load_data(filter_name_lowC_lowO);

results_all_highC_lowO = load_data(filter_name_highC_lowO);
results_all_lowC_highO = load_data(filter_name_lowC_highO);


%% find sessions that appear in all conditions
session_list_highC_highO = {results_all_highC_highO(:).sessionStr};
session_list_lowC_lowO = {results_all_lowC_lowO(:).sessionStr};
session_list_highC_lowO = {results_all_highC_lowO(:).sessionStr};
session_list_lowC_highO = {results_all_lowC_highO(:).sessionStr};


session_list_intersect = session_list_highC_highO;
session_list_intersect = intersect(session_list_intersect, session_list_lowC_lowO);
session_list_intersect = intersect(session_list_intersect, session_list_highC_lowO);
session_list_intersect = intersect(session_list_intersect, session_list_lowC_highO);


session_list_Ro = session_list_intersect(contains(session_list_intersect,'Ro'));
session_list_Gr = session_list_intersect(contains(session_list_intersect,'Gr'));
%%

[redundancy_change_cardinal.highC_highO_Ro, redundancy_change_oblique.highC_highO_Ro, redundancy_change.highC_highO_Ro] = ...
    extract_data(results_all_highC_highO, session_list_Ro);

[redundancy_change_cardinal.highC_highO_Gr, redundancy_change_oblique.highC_highO_Gr, redundancy_change.highC_highO_Gr] = ...
    extract_data(results_all_highC_highO, session_list_Gr);


[redundancy_change_cardinal.lowC_lowO_Ro, redundancy_change_oblique.lowC_lowO_Ro, redundancy_change.lowC_lowO_Ro] = ...
    extract_data(results_all_lowC_lowO, session_list_Ro);

[redundancy_change_cardinal.lowC_lowO_Gr, redundancy_change_oblique.lowC_lowO_Gr, redundancy_change.lowC_lowO_Gr] = ...
    extract_data(results_all_lowC_lowO, session_list_Gr);

[redundancy_change_cardinal.highC_lowO_Ro, redundancy_change_oblique.highC_lowO_Ro, redundancy_change.highC_lowO_Ro] = ...
    extract_data(results_all_highC_lowO, session_list_Ro);
[redundancy_change_cardinal.highC_lowO_Gr, redundancy_change_oblique.highC_lowO_Gr, redundancy_change.highC_lowO_Gr] = ...
    extract_data(results_all_highC_lowO, session_list_Gr);


[redundancy_change_cardinal.lowC_highO_Ro, redundancy_change_oblique.lowC_highO_Ro, redundancy_change.lowC_highO_Ro] = ...
    extract_data(results_all_lowC_lowO, session_list_Ro);
[redundancy_change_cardinal.lowC_highO_Gr, redundancy_change_oblique.lowC_highO_Gr, redundancy_change.lowC_highO_Gr] = ...
    extract_data(results_all_lowC_lowO, session_list_Gr);

%% paired test between four groups
data_list   = {redundancy_change.highC_highO_Gr; redundancy_change.lowC_lowO_Gr};
data_mean   = cellfun(@mean, data_list);
data_sem    = cellfun(@std, data_list) ./ sqrt(cellfun(@numel, data_list));

for i = 1:numel(data_list)
    bar(i, data_mean(i)); hold on;
    errorbar(i,data_mean(i),data_sem(i),'color','black');
end


for i = 1:numel(data_list)
    for j = i+1:numel(data_list)

        ttest(data_list{i}, data_list{j})
    end
end

%% function
function results_all = load_data(filter_name)
    saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', filter_name);
    load(fullfile(saveFolder, 'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions'));
    results_all = results_cross_sizeControl;
    results_all = get_sample_CI_cross(results_all);
end



function [redundancy_change_cardinal, redundancy_change_oblique, redundancy_change] = extract_data(results_all, session_list)
    idx = strcmp({results_all(:).sessionType}, 'mainTask') & cell2mat({results_all(:).timeWinIndex}) == 0 &...
        ismember({results_all(:).sessionStr}, session_list);

    deltaI_cardinal =  cell2mat({results_all(idx).delta_percent_cardinal_cardinal_median});
    deltaI_cardinal_cross =  cell2mat({results_all(idx).delta_percent_cardinal_oblique_median});
    deltaI_oblique = cell2mat({results_all(idx).delta_percent_oblique_oblique_median});
    deltaI_oblique_cross = cell2mat({results_all(idx).delta_percent_oblique_cardinal_median});
    
    redundancy_change_cardinal  = deltaI_cardinal - deltaI_cardinal_cross;
    redundancy_change_oblique   = deltaI_oblique - deltaI_oblique_cross;

    deltaI_within     =  cell2mat({results_all(idx).delta_percent_within_median});
    deltaI_cross    = cell2mat({results_all(idx).delta_percent_cross_median});

    redundancy_change = deltaI_within - deltaI_cross;
end