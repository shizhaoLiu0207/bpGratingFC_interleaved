function run_directFisherInfoEstimate_Cross_versionControl(filter_name,runOptions)
% clear all
% clc
% close all

%%% This script adapts from run_decodeStimulus. The main purpose is to make
%%% it easy to run the analysis with different subset of trial/neurons, and
%%% keep track of the trial/neuron filter applied.
%%% By Shizhao Liu 02/20/2023

%%%% 10/09/2024: Corrected the cross decoding part
%%
%doRun = runOptions.doRun;
doRunCross = runOptions.doRuncross;
%doOrganize = runOptions.doOrganize;
doOrganizeCross = runOptions.doOrganizeCross;
%doCombineAcrossCoherence = runOptions.doCombineAcrossCoherence;
%doReplace = runOptions.doReplace;
doReplace_cross = runOptions.doReplace_cross;
doTask        = runOptions.doTask; 
doOrientation = runOptions.doOrientation;
doMultiple_timebin  = runOptions.doMultiple_timebin;

timebinsize = runOptions.timebinsize; % in ms



global bpGlobal
bpGratingFCGlobal;
%%% data path and session list
%%%%%%%%%%%%%%%%%%%%%%%%%% specify session list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% session list of rolo
session_list_ro = bpGlobal.rolo.session_list;
%session_list_ro_all = [session_list_ro.cardinal;session_list_ro.oblique;session_list_ro.switching];
session_list_ro_all = session_list_ro.switching;
%%%% session list of gremlin
session_list_gr = bpGlobal.gremlin.session_list;
session_list_gr_all    = session_list_gr.interleaved;
bad_session_list_gr    = [session_list_gr.cardinal_hugeBias;session_list_gr.cardinal_wrongtargetEarlyoff;session_list_gr.recording_error_list];
session_list_gr_all    = setdiff(session_list_gr_all,bad_session_list_gr);



session_list_all    = [session_list_ro_all; session_list_gr_all];

interleaved_session_all = [session_list_ro.switching; setdiff(session_list_gr.interleaved,bad_session_list_gr )];
%%%%%%%%%%%%%%%%%%%%%%%%%% specify data path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfolder('/Users/liushizhao/projectData_local/probinf_data/extractedData')
    % on macbook
    extractedData_folder = '/Users/liushizhao/projectData_local/probinf_data/extractedData';
elseif isfolder('/home/shizhao/Documents/projectData/probinf_data/extractedData')
    extractedData_folder = '/home/shizhao/Documents/projectData/probinf_data/extractedData';
else
    error('No extracted data folder')
end

taskDataPath_Ro         = fullfile(extractedData_folder,'Ro/bpGratingFC');
oriDataPath_Ro          = fullfile(extractedData_folder,'Ro/rfMapping_bpGrating');
unitPropertiesPath_Ro   = fullfile(extractedData_folder,'Ro/unitProperties');

taskDataPath_Gr         = fullfile(extractedData_folder,'Gr/bpGratingFC');
oriDataPath_Gr          = fullfile(extractedData_folder,'Gr/rfMapping_bpGrating');
unitPropertiesPath_Gr   = fullfile(extractedData_folder,'Gr/unitProperties');



saveFolder          = sprintf('../../results/neural/fisherInfo_cross_direct/fisherInfo_cross_direct_%s',filter_name); 


%saveFolder_session = fullfile(saveFolder,'individual_sessions');
saveFolder_session_cross = fullfile(saveFolder,'individual_sessions_cross');
%mkdir(saveFolder_session);
mkdir(saveFolder_session_cross);

%% run the analysis
if  doRunCross
    %%% load the file that specifies which trials/neurons to use
    filter_folder   = fullfile('../../results/filter_trial_neuron',filter_name);
    filter_file     = fullfile(filter_folder,sprintf('filtered_trials_neurons_%s.mat',filter_name));
    %txt_file        = fullfile(filter_folder,sprintf('criterion_%s.txt',filter_name));
    load(filter_file);
    mkdir(saveFolder)
    copyfile(filter_file,saveFolder)
   % copyfile(txt_file, saveFolder);
    
    session_list_all = {filtered_trials_neurons(:).sessionStr};
    nSession            = numel(session_list_all); 



    for i_session = 1:nSession
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load needed data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sessionStr = session_list_all{i_session};
        sessionStr_real = sessionStr(1:10);
        sessionStr_save = sessionStr;
        subjectCode = sessionStr_real(1:2);
        dateStr     = sessionStr_real(3:end);
        eval(sprintf('taskDataPath = taskDataPath_%s;', subjectCode));
        eval(sprintf('oriDataPath  = oriDataPath_%s;', subjectCode));
        eval(sprintf('unitPropertiesPath = unitPropertiesPath_%s;', subjectCode));
       % saveName = fullfile(saveFolder_session, sprintf('%s_%s.mat',sessionStr_save,filter_name));
        saveName_cross = fullfile(saveFolder_session_cross, sprintf('%s_%s.mat',sessionStr,filter_name));

        % run_fisher = 0;
        % if doRun & (~isfile(saveName) | doReplace)
        % 
        %     run_fisher = 1;
        % 
        % end
        run_fisher_cross = 0;
        if doRunCross & (~isfile(saveName_cross) | doReplace_cross)
          
             run_fisher_cross = 1;

        end

        % if ~run_fisher & ~run_fisher_cross
        %     % do not run neither, skip
        %     continue
        % end
        
        if ~run_fisher_cross 
            continue
        end

        % if ~run_fisher & run_fisher_cross & ~ismember(sessionStr, interleaved_session_all) 
        %     % only run cross but do not run original
        %     % skip if session is not in interleaved epoch
        %     continue
        % end
        % 
       


        %dat_fisher   = [];
        dat_fisher_cross = [];
        fprintf('Running session %d/%d\n',i_session,numel(session_list_all))

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find kept neuron and trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_filter = strcmp({filtered_trials_neurons(:).sessionStr}, sessionStr_save);
        trialInd_kept = filtered_trials_neurons(idx_filter).trialInd_kept;
        neuronIdx_kept = filtered_trials_neurons(idx_filter).neuronIdx_kept;
        
         
        if doTask
            %%%%%%%%%%%%%%%%%%%%% Estimate fisher information for the main task session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% load data %%%%%%%%
    
            filelist = dir(fullfile(taskDataPath,sprintf('%s_*_%s_*extractedData.mat',sessionStr_real,'bpGratingFC')));
            load(fullfile(fullfile(taskDataPath,filelist(1).name)),'stimulus','neuro','behavior')
    
            %%%%%%% extract firing rate and stimulus information %%%%%%
           
            trialInd     = trialInd_kept; % used kept trials specified in the filter
            data_run    = struct();
            data_run.signal      = stimulus.signal(trialInd);
            data_run.orientation = stimulus.meanThetaDEG(trialInd);
            %signal_abs_list = nonzeros(unique(abs(signal)));
            data_run.behavCorrect    = behavior.correct(trialInd);
    
    
            info_run                    = struct();
            info_run.sessionStr         = sessionStr_real;
            info_run.sessionType        = 'mainTask';
    
            %%%%%%%%%%%%%%%%%% prepare neuro data %%%%%%%%%%%%%%%%%%
            % do decoding for several time window: the whole trial and every 400ms
            % use 400 ms to match the orientation session
            %stimTime    = stimulus.nImages * 100; % total length of stimulus
           
            stimTime = 1600; % 2024/12/23, always use 1600 ms even for some sessions the stimulus length is 2000

            % the first time window is the whole trial
            % take into account 50 ms latency
            latency = 50; 
            if doMultiple_timebin
                nWin        = floor(stimTime / timebinsize) + 1;
                spkWindow_list = cell(nWin,1);
                spkWindow_list{1} = [latency + 1,stimTime + latency];
                for k = 2:nWin
                    spkWindow_list{k} = [latency + 1,timebinsize + latency] + (k-2) * timebinsize;
                end
            else
                nWin = 1;
                spkWindow_list{1} = [latency + 1,stimTime + latency];
            end
            
            nTrial      = numel(trialInd);
            nNeuron     = numel(neuro.unitId);
            for i_sample = 1:numel(neuronIdx_kept)
                %%%%%%%%%% kept neurons %%%%%%%%
            

                is_to_keep    = zeros(nNeuron,1);
                nNeuron_keep  = numel(neuronIdx_kept(i_sample).unitId);
                for n_k = 1:nNeuron_keep
                     idx_keep = neuro.electrode == neuronIdx_kept(i_sample).electrode(n_k)  & ...
                    neuro.idOnChannel == neuronIdx_kept(i_sample).idOnChannel(n_k);
                     is_to_keep(idx_keep) = 1;
                end
               
                for i_win = 1:nWin
                    data_run.spikeCount  = zeros(nTrial,nNeuron);
                    spkWindow   = spkWindow_list{i_win};
               
                    for n = 1:nNeuron
                        data_run.spikeCount(:,n) = arrayfun(@(t)sum(neuro.spikeTimeMS{n,t} >= neuro.stimOnMS(t) + spkWindow(1) & ...
                            neuro.spikeTimeMS{n,t} <= neuro.stimOnMS(t) + spkWindow(2)),trialInd);
                    end
                    
                    % throw away not kept neurons
                    data_run.spikeCount(:,~boolean(is_to_keep)) = [];
        
                    info_run.timeWin                      = spkWindow;
                    info_run.i_win                        = i_win;
                    % if run_fisher
                    %     nBefore = numel(dat_fisher);
                    %     dat_fisher    = run_fisher_estimate_one_session(dat_fisher,data_run,info_run);
                    %     nCurrent = numel(dat_fisher);
                    %     idx_add = [nBefore + 1 :nCurrent];
                    % 
                    %     for i_a = 1:numel(idx_add)
                    %         dat_fisher(idx_add(i_a)).i_subSample = i_sample;
                    %     end
                    % 
                    % end
        
                    if run_fisher_cross
                        nBefore = numel(dat_fisher_cross);
                        dat_fisher_cross    = run_fisher_cross_estimate_one_session(dat_fisher_cross,data_run,info_run);
                     
                        nCurrent = numel(dat_fisher_cross);
                        idx_add = [nBefore + 1 :nCurrent];
                        
                        for i_a = 1:numel(idx_add)
                            dat_fisher_cross(idx_add(i_a)).i_subSample = i_sample;
                        end
    
                    end
                   
                end
            end


            clear stimulus neuro
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% decoding analysis for orientation mapping session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% use all trials 
        %%%%% use the same subset of neurons as the main task session
  
        filelist = dir(fullfile(oriDataPath,sprintf('%s_*_%s_*extractedTuningData.mat',sessionStr_real,'rfMapping_bpGrating')));
        if ~isempty(filelist) & doOrientation
            load(fullfile(filelist(1).folder,filelist(1).name),'tuningFiringRate','dat');
            spkWindow                   = tuningFiringRate.stimTime([1,end]);

            i_win                       = 1;
            % firingRate_tuning           = tuningFiringRate.firingrateAll';
            signal_tuning               = tuningFiringRate.signalLevel;
            orientation_tuning          = tuningFiringRate.orientation;
            idx_flip                    = ismember(orientation_tuning,[90,135]);
            signal_tuning(idx_flip)     = -signal_tuning(idx_flip);
            signal_abs_list_tuning      = nonzeros(unique(abs(signal_tuning)));

            behavCorrect                = nan * ones(size(signal_tuning));

            [nNeuron,nTrial]            = size(tuningFiringRate.firingrateAll);
            spikeCount                   = zeros(nTrial,nNeuron);

            data_run    = struct();
            data_run.signal      = signal_tuning;
            data_run.orientation = orientation_tuning;
            %signal_abs_list = nonzeros(unique(abs(signal)));
            data_run.behavCorrect    = behavCorrect;

            info_run                    = struct();
            info_run.sessionStr         = sessionStr_real;
            info_run.sessionType        = 'orientation';


            neuro.unitId = sort(nonzeros(tuningFiringRate.sortcodeLookupTable(:)));
            nUnit = numel(neuro.unitId);
            [neuro.electrode,neuro.idOnChannel] = deal(zeros(nUnit,1));
            for n = 1:numel(neuro.unitId)
                [neuro.electrode(n),neuro.idOnChannel(n)] = find(tuningFiringRate.sortcodeLookupTable == neuro.unitId(n));
            end
            
            for i_sample = 1:numel(neuronIdx_kept)
                is_to_keep    = zeros(nNeuron,1);
                nNeuron_keep  = numel(neuronIdx_kept(i_sample).unitId);
                for n_k = 1:nNeuron_keep
                     idx_keep = neuro.electrode == neuronIdx_kept(i_sample).electrode(n_k)  & ...
                    neuro.idOnChannel == neuronIdx_kept(i_sample).idOnChannel(n_k);
                     is_to_keep(idx_keep) = 1;
                end
               
    
                for t = 1:nTrial
                    i_t = dat(t).rasterTimeline >= spkWindow(1) & dat(t).rasterTimeline <= spkWindow(2);
                    spikeCount(t,:) = sum(dat(t).spikes(:,i_t),2)';
                end
    
    
                data_run.spikeCount         = spikeCount;
                data_run.spikeCount(:,~boolean(is_to_keep)) = [];
    
                info_run.timeWin                      = spkWindow;
                info_run.i_win                        = i_win;
    
                % if run_fisher    
                %     nBefore = numel(dat_fisher);
                %     dat_fisher    = run_fisher_estimate_one_session(dat_fisher,data_run,info_run);
                %     nCurrent = numel(dat_fisher);
                %     idx_add = [nBefore + 1 :nCurrent];
                % 
                %     for i_a = 1:numel(idx_add)
                %         dat_fisher(idx_add(i_a)).i_subSample = i_sample;
                %     end
                % 
                % end
                if run_fisher_cross
                    nBefore = numel(dat_fisher_cross);
                    dat_fisher_cross    = run_fisher_cross_estimate_one_session(dat_fisher_cross,data_run,info_run);
                    nCurrent = numel(dat_fisher_cross);
                    idx_add = [nBefore + 1 :nCurrent];

                    for i_a = 1:numel(idx_add)
                        dat_fisher_cross(idx_add(i_a)).i_subSample = i_sample;
                    end

                end
    
               
    
            end
             clear tuningFiringRate dat
        end

        % %%%%%%% save file %%%%%%%
        % if run_fisher
        %     save(saveName,'dat_fisher')
        % end

        if run_fisher_cross 
            save(saveName_cross,"dat_fisher_cross");
        end


    end
end
% %% organize decoding results of all individual sessions into one big struct
% if doOrganize
%     filelist = dir(fullfile(saveFolder_session,'*.mat'));
%     dat_fisher = [];
%     for n = 1:numel(filelist)
%         dat_session = load(fullfile(saveFolder_session,filelist(n).name));
%         dat_fisher = [dat_fisher,dat_session.dat_fisher];
%     end
%     saveName = fullfile(saveFolder, sprintf('fisherInfo_direct_all_sessions_%s',filter_name));
%     save( saveName,'dat_fisher');
% end
% %% organize cross-fisher results of all individual sessions into one big struct
% if doOrganizeCross
%     filelist = dir(fullfile(saveFolder_session_cross, '*.mat'));
%     dat_fisher_cross = [];
%     for n = 1:numel(filelist)
%         dat_session = load(fullfile(saveFolder_session_cross,filelist(n).name));
%         dat_fisher_cross = [dat_fisher_cross,dat_session.dat_fisher_cross];
%     end
% 
%     saveName = fullfile(saveFolder, sprintf('fisherInfo_cross_direct_all_sessions_%s',filter_name));
%     save( saveName,'dat_fisher_cross');
% end



end


