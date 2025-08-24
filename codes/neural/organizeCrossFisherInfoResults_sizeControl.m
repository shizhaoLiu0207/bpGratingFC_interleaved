clear all
clc
close all

%%
versionName_list = {'all_trials_coef1_hVis2_FR1_interleaved_sizeControl'};
%...'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl'};
for n = 1:numel(versionName_list)
    versionName = versionName_list{n};
    saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', versionName);
    
    % dataName = fullfile(saveFolder,  sprintf('fisherInfo_cross_direct_all_sessions_%s',versionName));
    % load(dataName);

    data_list = dir(fullfile(saveFolder,'individual_sessions_cross','*.mat'));
    saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');
    
    results_cross_sizeControl = util_it.run_organize_cross_fisherinfo_sizeControl(data_list);
    
    save(saveName,'results_cross_sizeControl');
end
