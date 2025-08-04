clear all
clc
close all

%%

filter_name_list = {'all_trials_coef1_hVis2_FR1_hVisOri2_FROri2_interleaved_sizeControl'};

%runOptions.doRun = 1;
runOptions.doRuncross = 1;

%runOptions.doReplace = 1;
runOptions.doReplace_cross = 1;

%runOptions.doOrganize = 0;
runOptions.doOrganizeCross = 0;

%runOptions.doCombineAcrossCoherence = 0;


runOptions.doTask = 1;
runOptions.doOrientation = 1;


runOptions.doMultiple_timebin = 0;

runOptions.timebinsize = 200;

combineOptions = struct();

for n = 1:numel(filter_name_list)
    filter_name = filter_name_list{n};
    run_directFisherInfoEstimate_Cross_versionControl(filter_name,runOptions);
end