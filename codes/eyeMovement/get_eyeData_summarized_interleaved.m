clear all
clc
close all
%%
global bpGlobal
bpGratingFCGlobal;
session_list_ro = bpGlobal.rolo.session_list.switching;


%%%% session list of gremlin
session_list_gr = bpGlobal.gremlin.session_list.interleaved_real;


session_list_all    = [session_list_ro; session_list_gr];
nSession            = numel(session_list_all);

rolo_y_bias         =  6.4284;
%%
eyeData_summarize_all = struct();
eyePath = '../../results/eyeMovement/eyeData_avgTrial';
bootstrap_path = '../../data/Bootstrap_trialIndex';
load(fullfile(bootstrap_path, 'bootstrap_trialIndex_rolo_220623-1415.mat'))
CI_level = 68;
CI_pct   = [50 - CI_level / 2, 50 + CI_level / 2];
note_str = ['The first column of each variable correponds to the cardinal, the second row oblique. \n Confidence interval is 68'];
for n = 1:nSession
    fprintf('summarizing %d/%d \n',n,nSession)
    sessionStr      = session_list_all{n};
    subjectCode     = sessionStr(1:2);
    %%%%%% load eye data
    load(fullfile(eyePath, subjectCode,sprintf('%s_eyedata_avgTrial', sessionStr)));
    %%%%% load Bootstrap index
    switch subjectCode
        case 'Ro'
            idx  = contains({BootstrapIdx(:).sessionName},sessionStr);
            trialIdx = BootstrapIdx(idx).trialIdx;
        case 'Gr'
            bootstrap_path_gr = fullfile(bootstrap_path,'Gremlin');
            flist = dir(fullfile(bootstrap_path,'Gremlin',sprintf('%s_*.mat',sessionStr)));
            bootstrap_gr = load(fullfile(flist(1).folder,flist(1).name));
            trialIdx = bootstrap_gr.BootstrapIdx.trialIdx;
    end
    %%%% load stimulus data to get task information
    data_folder = sprintf('/Users/liushizhao/projectData_local/probinf_data/extractedData/%s/bpGratingFC', subjectCode);
    dlist = dir(fullfile(data_folder, [sessionStr,'*.mat'] ));
    load(fullfile(data_folder, dlist(1).name),'stimulus');
    
    %%%%%% bootstrap trials and get eye statstics
    nBootstrap = numel(trialIdx);
    %%%% the first column is cardinal and the second is oblique
    [x_pos_avg, y_pos_avg, x_vel_avg, y_vel_avg, pupil_avg]    = deal(zeros(nBootstrap, 2));
    [x_pos_var, y_pos_var, x_vel_var, y_vel_var, pupil_var]    = deal(zeros(nBootstrap, 2));

    [x_pos_zero_avg, y_pos_zero_avg, x_vel_zero_avg, y_vel_zero_avg, pupil_zero_avg]    = deal(zeros(nBootstrap, 2));
    [x_pos_zero_var, y_pos_zero_var, x_vel_zero_var, y_vel_zero_var, pupil_zero_var]    = deal(zeros(nBootstrap, 2));
    
    trialIdx_cardinal_all   = find(strcmp(stimulus.taskType,'C'));
    trialIdx_oblique_all    = find(strcmp(stimulus.taskType,'O'));

    trialIdx_cardinal_zero_all  = find(strcmp(stimulus.taskType,'C') & stimulus.signal' == 0);
    trialIdx_oblique_zero_all   = find(strcmp(stimulus.taskType,'O') & stimulus.signal' == 0);

    for k = 1:nBootstrap


        trialIdx_cardinal   = intersect(trialIdx{k}, trialIdx_cardinal_all);
        trialIdx_oblique    = intersect(trialIdx{k}, trialIdx_oblique_all);

        trialIdx_cardinal_zero  = intersect(trialIdx{k}, trialIdx_cardinal_zero_all);
        trialIdx_oblique_zero   = intersect(trialIdx{k}, trialIdx_oblique_zero_all);



        [~, idx_trial_cardinal] = ismember(trialIdx_cardinal, eyeData_avg.trialInd);
        [~, idx_trial_oblique]  = ismember(trialIdx_oblique,  eyeData_avg.trialInd);

        [~, idx_trial_cardinal_zero] = ismember(trialIdx_cardinal_zero, eyeData_avg.trialInd);
        [~, idx_trial_oblique_zero]  = ismember(trialIdx_oblique_zero,  eyeData_avg.trialInd);

        if sum(idx_trial_cardinal == 0) > 0
            idx_trial_cardinal(idx_trial_cardinal == 0) = [];
        end
        if sum(idx_trial_oblique == 0) > 0
            idx_trial_oblique(idx_trial_oblique == 0) = [];
        end

         if sum(idx_trial_cardinal_zero == 0) > 0
            idx_trial_cardinal_zero(idx_trial_cardinal_zero == 0) = [];
        end
        if sum(idx_trial_oblique_zero == 0) > 0
            idx_trial_oblique_zero(idx_trial_oblique_zero == 0) = [];
        end
       
        % average x-position
        x_pos_avg(k,1)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_cardinal,1), 'omitnan');
        x_pos_avg(k,2)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_oblique,1), 'omitnan');

        x_pos_zero_avg(k,1) = mean(eyeData_avg.eyePosition_corrected(idx_trial_cardinal_zero,1), 'omitnan');
        x_pos_zero_avg(k,2) = mean(eyeData_avg.eyePosition_corrected(idx_trial_oblique_zero,1), 'omitnan');
        % average y-position
        y_pos_avg(k,1)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_cardinal,2), 'omitnan');
        y_pos_avg(k,2)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_oblique,2), 'omitnan');

        y_pos_zero_avg(k,1)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_cardinal_zero,2), 'omitnan');
        y_pos_zero_avg(k,2)   = mean(eyeData_avg.eyePosition_corrected(idx_trial_oblique_zero,2), 'omitnan');
        
        % average x-speed
        x_vel_avg(k,1)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_cardinal,1), 'omitnan');
        x_vel_avg(k,2)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_oblique,1), 'omitnan');

        x_vel_zero_avg(k,1)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_cardinal_zero,1), 'omitnan');
        x_vel_zero_avg(k,2)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_oblique_zero,1), 'omitnan');

        % average y-speed
        y_vel_avg(k,1)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_cardinal,2),'omitnan');
        y_vel_avg(k,2)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_oblique,2),'omitnan');

        y_vel_zero_avg(k,1)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_cardinal_zero,2),'omitnan');
        y_vel_zero_avg(k,2)   = mean(eyeData_avg.meanEyeVelocity(idx_trial_oblique_zero,2),'omitnan');

        % average pupil size
        pupil_avg(k,1)   = mean(eyeData_avg.meanPupilDiam(idx_trial_cardinal), 'omitnan');
        pupil_avg(k,2)   = mean(eyeData_avg.meanPupilDiam(idx_trial_oblique), 'omitnan');

        pupil_zero_avg(k,1)   = mean(eyeData_avg.meanPupilDiam(idx_trial_cardinal_zero), 'omitnan');
        pupil_zero_avg(k,2)   = mean(eyeData_avg.meanPupilDiam(idx_trial_oblique_zero), 'omitnan');

        % variance of x-position
        x_pos_var(k,1)    = var(eyeData_avg.eyePosition_corrected(idx_trial_cardinal,1), 'omitnan');
        x_pos_var(k,2)    = var(eyeData_avg.eyePosition_corrected(idx_trial_oblique,1), 'omitnan');

        x_pos_zero_var(k,1)    = var(eyeData_avg.eyePosition_corrected(idx_trial_cardinal_zero,1), 'omitnan');
        x_pos_zero_var(k,2)    = var(eyeData_avg.eyePosition_corrected(idx_trial_oblique_zero,1), 'omitnan');

        % variance of y-position
        y_pos_var(k,1)    = var(eyeData_avg.eyePosition_corrected(idx_trial_cardinal,2), 'omitnan');
        y_pos_var(k,2)    = var(eyeData_avg.eyePosition_corrected(idx_trial_oblique,2), 'omitnan');

        y_pos_zero_var(k,1)    = var(eyeData_avg.eyePosition_corrected(idx_trial_cardinal_zero,2), 'omitnan');
        y_pos_zero_var(k,2)    = var(eyeData_avg.eyePosition_corrected(idx_trial_oblique_zero,2), 'omitnan');

        % variance of x-speed
        x_vel_var(k,1)    = var(eyeData_avg.meanEyeVelocity(idx_trial_cardinal,1), 'omitnan');
        x_vel_var(k,2)    = var(eyeData_avg.meanEyeVelocity(idx_trial_oblique,1), 'omitnan');

        x_vel_zero_var(k,1)    = var(eyeData_avg.meanEyeVelocity(idx_trial_cardinal_zero,1), 'omitnan');
        x_vel_zero_var(k,2)    = var(eyeData_avg.meanEyeVelocity(idx_trial_oblique_zero,1), 'omitnan');

        % variance of y-speed
        y_vel_var(k,1)    = var(eyeData_avg.meanEyeVelocity(idx_trial_cardinal,2),'omitnan');
        y_vel_var(k,2)    = var(eyeData_avg.meanEyeVelocity(idx_trial_oblique,2),'omitnan');

        y_vel_zero_var(k,1)    = var(eyeData_avg.meanEyeVelocity(idx_trial_cardinal_zero,2),'omitnan');
        y_vel_zero_var(k,2)    = var(eyeData_avg.meanEyeVelocity(idx_trial_oblique_zero,2),'omitnan');


        % variance of pupil size
        pupil_var(k,1)    = var(eyeData_avg.meanPupilDiam(idx_trial_cardinal), 'omitnan');
        pupil_var(k,2)    = var(eyeData_avg.meanPupilDiam(idx_trial_oblique), 'omitnan');

        pupil_zero_var(k,1)    = var(eyeData_avg.meanPupilDiam(idx_trial_cardinal_zero), 'omitnan');
        pupil_zero_var(k,2)    = var(eyeData_avg.meanPupilDiam(idx_trial_oblique_zero), 'omitnan');
    end
    %%%%%% summarize data into one struct
    eyeData_summarize_all(n).sessionStr       = sessionStr;
    eyeData_summarize_all(n).CI_level         = CI_level; 

    eyeData_summarize_all(n).x_pos_avg_median = median(x_pos_avg, 1);
    eyeData_summarize_all(n).x_pos_avg_CI     = prctile(x_pos_avg, CI_pct, 1);
    eyeData_summarize_all(n).y_pos_avg_median = median(y_pos_avg, 1);
    eyeData_summarize_all(n).y_pos_avg_CI     = prctile(y_pos_avg, CI_pct, 1);

    eyeData_summarize_all(n).x_vel_avg_median = median(x_vel_avg, 1);
    eyeData_summarize_all(n).x_vel_avg_CI     = prctile(x_vel_avg, CI_pct, 1);
    eyeData_summarize_all(n).y_vel_avg_median = median(y_vel_avg, 1);
    eyeData_summarize_all(n).y_vel_avg_CI     = prctile(y_vel_avg, CI_pct, 1);

    eyeData_summarize_all(n).pupil_avg_median = median(pupil_avg, 1);
    eyeData_summarize_all(n).pupil_avg_CI     = prctile(pupil_avg, CI_pct, 1);

    eyeData_summarize_all(n).x_pos_var_median = median(x_pos_var, 1);
    eyeData_summarize_all(n).x_pos_var_CI     = prctile(x_pos_var, CI_pct, 1);
    eyeData_summarize_all(n).y_pos_var_median = median(y_pos_var, 1);
    eyeData_summarize_all(n).y_pos_var_CI     = prctile(y_pos_var, CI_pct, 1);

    eyeData_summarize_all(n).x_vel_var_median = median(x_vel_var, 1);
    eyeData_summarize_all(n).x_vel_var_CI     = prctile(x_vel_var, CI_pct, 1);
    eyeData_summarize_all(n).y_vel_var_median = median(y_vel_var, 1);
    eyeData_summarize_all(n).y_vel_var_CI     = prctile(y_vel_var, CI_pct, 1);

    eyeData_summarize_all(n).pupil_var_median = median(pupil_var, 1);
    eyeData_summarize_all(n).pupil_var_CI     = prctile(pupil_var, CI_pct, 1);




    eyeData_summarize_all(n).x_pos_zero_avg_median = median(x_pos_zero_avg, 1);
    eyeData_summarize_all(n).x_pos_zero_avg_CI     = prctile(x_pos_zero_avg, CI_pct, 1);
    eyeData_summarize_all(n).y_pos_zero_avg_median = median(y_pos_zero_avg, 1);
    eyeData_summarize_all(n).y_pos_zero_avg_CI     = prctile(y_pos_zero_avg, CI_pct, 1);

    eyeData_summarize_all(n).x_vel_zero_avg_median = median(x_vel_zero_avg, 1);
    eyeData_summarize_all(n).x_vel_zero_avg_CI     = prctile(x_vel_zero_avg, CI_pct, 1);
    eyeData_summarize_all(n).y_vel_zero_avg_median = median(y_vel_zero_avg, 1);
    eyeData_summarize_all(n).y_vel_zero_avg_CI     = prctile(y_vel_zero_avg, CI_pct, 1);

    eyeData_summarize_all(n).pupil_zero_avg_median = median(pupil_zero_avg, 1);
    eyeData_summarize_all(n).pupil_zero_avg_CI     = prctile(pupil_zero_avg, CI_pct, 1);

    eyeData_summarize_all(n).x_pos_zero_var_median = median(x_pos_zero_var, 1);
    eyeData_summarize_all(n).x_pos_zero_var_CI     = prctile(x_pos_zero_var, CI_pct, 1);
    eyeData_summarize_all(n).y_pos_zero_var_median = median(y_pos_zero_var, 1);
    eyeData_summarize_all(n).y_pos_zero_var_CI     = prctile(y_pos_zero_var, CI_pct, 1);

    eyeData_summarize_all(n).x_vel_zero_var_median = median(x_vel_zero_var, 1);
    eyeData_summarize_all(n).x_vel_zero_var_CI     = prctile(x_vel_zero_var, CI_pct, 1);
    eyeData_summarize_all(n).y_vel_zero_var_median = median(y_vel_zero_var, 1);
    eyeData_summarize_all(n).y_vel_zero_var_CI     = prctile(y_vel_zero_var, CI_pct, 1);

    eyeData_summarize_all(n).pupil_zero_var_median = median(pupil_zero_var, 1);
    eyeData_summarize_all(n).pupil_zero_var_CI     = prctile(pupil_zero_var, CI_pct, 1);

end
save('../../results/eyeMovement/eyedata_summarized_interleaved','eyeData_summarize_all','note_str')