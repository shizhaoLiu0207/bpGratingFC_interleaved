clear all
clc
close all

%%

filter_name_list = {'all_trials_coef1_hVis2_FR1_interleaved_forEye_sizeControl';...
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

%runOptions.doRun = 1;
runOptions.doRuncross = 1;

%runOptions.doReplace = 1;
runOptions.doReplace_cross = 1;

runOptions.doSmallSample = 0;
%runOptions.doOrganize = 0;
runOptions.doOrganizeCross = 0;

%runOptions.doCombineAcrossCoherence = 0;


runOptions.doTask = 1;
runOptions.doOrientation = 0;


runOptions.doMultiple_timebin = 0;

runOptions.timebinsize = 200;

combineOptions = struct();

for n = 1:numel(filter_name_list)
    filter_name = filter_name_list{n};
    run_directFisherInfoEstimate_Cross_versionControl(filter_name,runOptions);
end