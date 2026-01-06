clear all
clc
close all

%%
%versionName_list = {'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl'};
versionName_list = {'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'distanceEyePosition_center_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'distanceEyePosition_center_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'distanceEyePosition_surround_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'distanceEyePosition_surround_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'eyeVelocity_high_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'eyeVelocity_high_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'eyeVelocity_low_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'eyeVelocity_low_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'pupilSize_high_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'pupilSize_high_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'pupilSize_low_prctile50_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
                    'pupilSize_low_prctile75_zero_coef1_hVis2_FR1_interleaved_forEye_sizeControl'};

%...'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl'};
doTimebin = false;
for n = 1:numel(versionName_list)
    versionName = versionName_list{n};
  
    saveFolder = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', versionName);
  
    
    % dataName = fullfile(saveFolder,  sprintf('fisherInfo_cross_direct_all_sessions_%s',versionName));
    % load(dataName);
    if doTimebin
        data_list = dir(fullfile(saveFolder,'individual_sessions_cross_timebin','*.mat'));
        %saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_timebin');
        organized_folder = fullfile(saveFolder,'results_cohr_combined_session_timebin');
        mkdir(organized_folder)
    else
        data_list = dir(fullfile(saveFolder,'individual_sessions_cross','*.mat'));
        %saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');
        organized_folder = fullfile(saveFolder,'results_cohr_combined_session');
        mkdir(organized_folder)
    
    end
    
    %results_cross_sizeControl = util_it.run_organize_cross_fisherinfo_sizeControl(data_list, organized_folder);
    util_it.run_organize_cross_fisherinfo_sizeControl(data_list, organized_folder);
    %save(saveName,'results_cross_sizeControl');

     %%%  combine per-session organize result
     if doTimebin
        saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_timebin');
        individual_folder = fullfile(saveFolder,'results_cohr_combined_session_timebin');
     else
        saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');
        individual_folder = fullfile(saveFolder,'results_cohr_combined_session');
     end
    file_list = dir(fullfile(individual_folder,'*.mat'));
    
    %results_sizeControl_combined = struct();
    for t = 1:numel(file_list)
        results_sizeControl_individual = load(fullfile(individual_folder, file_list(t).name));
    
        if t == 1
            results_cross_sizeControl = results_sizeControl_individual.results_cross_sizeControl;
        else
        results_cross_sizeControl = [results_cross_sizeControl,...
            results_sizeControl_individual.results_cross_sizeControl];
        end
    
    end
    
    save(saveName,'results_cross_sizeControl');

end
% %% if need, combine per-session organize result
% doThis = 0;
% if doThis
%     versionName = 'all_trials_coef1_hVis2_FR1_interleaved_sizeControl';
%     saveFolder =  sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s', versionName);
%     if doTimebin
%         saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_timebin');
%         individual_folder = fullfile(saveFolder,'results_cohr_combined_session_timebin');
%     else
%         saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions');
%         individual_folder = fullfile(saveFolder,'results_cohr_combined_session');
%     end
% 
% 
% 
% 
%     file_list = dir(fullfile(individual_folder,'*.mat'));
% 
%     %results_sizeControl_combined = struct();
%     for n = 1:numel(file_list)
%         results_sizeControl_individual = load(fullfile(individual_folder, file_list(n).name));
% 
%         if n == 1
%             results_cross_sizeControl = results_sizeControl_individual.results_cross_sizeControl;
%         else
%         results_cross_sizeControl = [results_cross_sizeControl,...
%             results_sizeControl_individual.results_cross_sizeControl];
%         end
% 
%     end
% 
%     save(saveName,'results_cross_sizeControl');
% end