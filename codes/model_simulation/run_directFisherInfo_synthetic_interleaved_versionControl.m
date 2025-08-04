clear all
clc
close all
%%
doRun_cross = 1;
%doOrganize_cross = 0; % simply put everything into a big struct

doMultipleTimebin = 0;

doOrganize_perCohr_cross = 0;  % average across subsamples
doOrganize_combine_cross = 0;
versionName = 'subset_32_random_1000'; % version of neuron filters

%versionName = 'all_neurons';
load(sprintf('../../results/filtered_neuron_synthetic/filtered_neurons_%s.mat',versionName));

subjectCode = 'Model';


% dataFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_ec_x256_g64_nr1024_contrast_all/';
% saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_model_versionControl/%s',versionName);
%%%% on macboox
dataFolder   = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved';
%%%% on linux
if ~exist(dataFolder)
    dataFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved';
end

saveFolder = sprintf('../../results/neural/fisherInfo_direct/fisherInfo_direct_modelInterleaved_versionControl/%s',versionName);

saveFolder_session = fullfile(saveFolder, 'individual_sessions_cross');
mkdir(saveFolder_session);

% data_list = dir(fullfile(dataFolder,'*.mat')); % all sessions
%%%% specific sessions

data_list(1).name = 'synthData_use_interleaved_bPF_0_00_cardinal_delta_0_00_prior_1_00_oblique_delta_0_00_prior_1_00.mat';
nSession  = numel(data_list);
if doMultipleTimebin
    nTimeBin = 9;
else
    nTimeBin = 1;
end

if doRun_cross
    for n = 1:nSession

        dat_fisher_cross    = [];

        
        
        load(fullfile(dataFolder, data_list(n).name));
        sessionStr = synthData_interleaved(1).sessionStr;

        saveName   = fullfile(saveFolder_session,sprintf('%s.mat',sessionStr));
        if isfile(saveName)
            continue
        end
        
        timeBin_index               = cell2mat({synthData_interleaved(:).timebinIndex});
        for t = 1:nTimeBin
            data_out = get_synthetic_data_realInterleaved(synthData_interleaved, t - 1);


            info_run                    = struct();
            info_run.sessionStr         = sessionStr;
            info_run.sessionType        = 'mainTask';
            info_run.timeWin            = data_out.timebin;
            info_run.i_win              = t;

            fprintf('Running session %d timewin %d\n',n,t);

            %idx_session    = strcmp({filtered_neurons(:).sessionStr}, sessionStr);
            idx_session     = 1; % we can just use the neuron filter for the first session since neuron population is exactly the same across sessions in synthetic data
            neuronIdx_kept = filtered_neurons(idx_session).neuronIdx_kept;
            nSample         = size(neuronIdx_kept,1);
            
            data_run = data_out; % first replicate everything in data_out to data run
            for i_sample = 1:nSample
                % subsample neurons
                data_run.spikeCount = data_out.spikeCount(:, neuronIdx_kept(i_sample,:));

                nBefore = numel(dat_fisher_cross);
                dat_fisher_cross     = run_fisher_cross_estimate_one_session(dat_fisher_cross,data_run,info_run);
                nCurrent = numel(dat_fisher_cross);
                idx_add = [nBefore + 1 :nCurrent];

                for i_a = 1:numel(idx_add)
                    dat_fisher_cross(idx_add(i_a)).i_subSample = i_sample;
                end
            end



        end
        save(saveName,'dat_fisher_cross');
    end
end

% %% organize decoding results of all individual sessions into one big struct
% if doOrganize_cross
%     filelist = dir(fullfile(saveFolder_session,'*.mat'));
%     dat_fisher_cross = [];
%     for n = 1:numel(filelist)
%         dat_session = load(fullfile(saveFolder_session,filelist(n).name));
%         dat_fisher_cross = [dat_fisher_cross,dat_session.dat_fisher_cross];
%     end
%     saveName = fullfile(saveFolder, 'fisherInfo_cross_direct_all_sessions_syntheticData');
%     save( saveName,'dat_fisher_cross');
% end
%% Organize across subsamples of neurons per coherence level
if doOrganize_perCohr_cross
    %saveFolder = saveFolder_session;
    filelist = dir(fullfile(saveFolder_session ,'*.mat'));
    nSession = numel(filelist);
    organized_folder = fullfile(saveFolder, 'results_organized_perCohr_session');
    mkdir(organized_folder);
    
    
%     dataName = fullfile(saveFolder,  'fisherInfo_cross_direct_all_sessions_syntheticData');
%     load(dataName);
% 
%     session_list = unique({dat_fisher_cross(:).sessionStr});
%     nSession = numel(session_list);
%     t = 1;
    for n = 1:nSession
        t = 1;
        
        load(fullfile(saveFolder_session, filelist(n).name));
        saveName_session = fullfile(organized_folder, sprintf('results_perCohr_%s', filelist(n).name));
        fprintf('Organizing session %d/%d\n',n,nSession)
        
      
        nSample = max(cell2mat({dat_fisher_cross(:).i_subSample}));
        timebin_list = unique(cell2mat({dat_fisher_cross(:).timeWinIndex}));

        for k = 1:numel(timebin_list)
            idx_base = cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1;
            cohr_list = unique(cell2mat({dat_fisher_cross(idx_base).coherence_level}));
            

            for i = 1:numel(cohr_list)
                idx = idx_base & cell2mat({dat_fisher_cross(:).coherence_level}) == cohr_list(i);
                dat_fisher_cross_sub = dat_fisher_cross(idx);
    
                results_cross_sizeControl_perCohr(t).sessionStr     = dat_fisher_cross_sub(1).sessionStr;
                results_cross_sizeControl_perCohr(t).sessionType    = dat_fisher_cross_sub(1).sessionType;
                results_cross_sizeControl_perCohr(t).timeWin        = dat_fisher_cross_sub(1).timeWin;
                results_cross_sizeControl_perCohr(t).timeWinIndex   = dat_fisher_cross_sub(1).timeWinIndex;
                results_cross_sizeControl_perCohr(t).nUnit          = dat_fisher_cross_sub(1).N;
                results_cross_sizeControl_perCohr(t).coherence_level = dat_fisher_cross_sub(1).coherence_level;



                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_bc});
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_bc});
               
                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_bc});
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_bc});

                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_cardinal_shuffle_bc});
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_oblique_shuffle_bc});
               
                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_oblique_cardinal_shuffle_bc});
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_sub(:).fisher_cardinal_oblique_shuffle_bc});

    
                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_median = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_median   = median(results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_median  = median(results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_sample);
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_median  = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_sample);
    
                results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_median = median(results_cross_sizeControl_perCohr(t).fisher_cardinal_cardinal_shuffle_sample);
                results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_median   = median(results_cross_sizeControl_perCohr(t).fisher_oblique_oblique_shuffle_sample);

                results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_median = median(results_cross_sizeControl_perCohr(t).fisher_oblique_cardinal_shuffle_sample);
                results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_median =  median(results_cross_sizeControl_perCohr(t).fisher_cardinal_oblique_shuffle_sample);
               

                
                t = t+1;

            end
        end
        save(saveName_session, 'results_cross_sizeControl_perCohr');
        clear results_cross_sizeControl_perCohr
    end
       %%%% combine into one file
    file_list = dir(fullfile(organized_folder,'*.mat'));
    results_cross_sizeControl_perCohr = [];
    for n = 1:numel(file_list)
        data = load(fullfile(organized_folder, file_list(n).name));
        results_cross_sizeControl_perCohr = [results_cross_sizeControl_perCohr, data.results_cross_sizeControl_perCohr];
        
    end
    saveName = fullfile(saveFolder,'results_SubsampleOrganized_perCohr_fisherInfo_all_sessions_syntheticData');

    save(saveName,'results_cross_sizeControl_perCohr');
end
%% Organize across subsamples of neurons and combine across coherence level
if doOrganize_combine_cross
   % saveFolder = saveFolder_session;
    filelist = dir(fullfile(saveFolder_session ,'*.mat'));
    nSession = numel(filelist);
    organized_folder = fullfile(saveFolder, 'results_cohr_combined_session');
    mkdir(organized_folder);
    

   
    for n = 1:nSession
        t = 1;
       
        load(fullfile(saveFolder_session, filelist(n).name));
        saveName_session = fullfile(organized_folder, sprintf('results_cohrCombined_%s', filelist(n).name));
        fprintf('Organizing session %d/%d\n',n,nSession)
      
        nSample = max(cell2mat({dat_fisher_cross(:).i_subSample}));
        timebin_list = unique(cell2mat({dat_fisher_cross(:).timeWinIndex}));

        for k = 1:numel(timebin_list)
            
            idx_base = cell2mat({dat_fisher_cross(:).timeWinIndex}) == k-1;
                
       

                for i = 1:nSample
                    idx = idx_base & cell2mat({dat_fisher_cross(:).i_subSample}) == i;
                    options = struct();
                  
                    [dat_fisher_cross_combine(i), options] = combine_coherence_crossfisherInfo(dat_fisher_cross(idx),options);
                end
        
                if isempty(fieldnames(dat_fisher_cross_combine(1)))
                    clear dat_fisher_cross_combine
                    continue
                end


                results_cross_sizeControl(t).sessionStr     = dat_fisher_cross_combine(1).sessionStr;
                results_cross_sizeControl(t).sessionType    = dat_fisher_cross_combine(1).sessionType;
                results_cross_sizeControl(t).timeWin        = dat_fisher_cross_combine(1).timeWin;
                results_cross_sizeControl(t).timeWinIndex   = dat_fisher_cross_combine(1).timeWinIndex;
                results_cross_sizeControl(t).nUnit          = dat_fisher_cross_combine(1).nUnit;
   


                results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal});
                results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique});
               
                results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal});
                results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique});

                results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_cardinal_shuffle});
                results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_oblique_shuffle});
               
                results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_oblique_cardinal_shuffle});
                results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample = cell2mat({dat_fisher_cross_combine(:).combine_fisher_cardinal_oblique_shuffle});

    
                results_cross_sizeControl(t).fisher_cardinal_cardinal_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_sample);
                results_cross_sizeControl(t).fisher_oblique_oblique_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_sample);
                results_cross_sizeControl(t).fisher_oblique_cardinal_median  = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_sample);
                results_cross_sizeControl(t).fisher_cardinal_oblique_median  = median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_sample);
    
                results_cross_sizeControl(t).fisher_cardinal_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_cardinal_cardinal_shuffle_sample);
                results_cross_sizeControl(t).fisher_oblique_oblique_shuffle_median   = median(results_cross_sizeControl(t).combine_fisher_oblique_oblique_shuffle_sample);

                results_cross_sizeControl(t).fisher_oblique_cardinal_shuffle_median = median(results_cross_sizeControl(t).combine_fisher_oblique_cardinal_shuffle_sample);
                results_cross_sizeControl(t).fisher_cardinal_oblique_shuffle_median =  median(results_cross_sizeControl(t).combine_fisher_cardinal_oblique_shuffle_sample);
               

                t = t+1;
                clear dat_fisher_cross_combine

        end
        save(saveName_session, 'results_cross_sizeControl');
        clear results_cross_sizeControl
    end
    %%%% combine into one file
    file_list = dir(fullfile(organized_folder,'*.mat'));
    results_cross_sizeControl = [];
    for n = 1:numel(file_list)
        data = load(fullfile(organized_folder, file_list(n).name));
        results_cross_sizeControl = [results_cross_sizeControl, data.results_cross_sizeControl];
        
    end
   saveName = fullfile(saveFolder,'results_SubsampleCombined_combinedCohr_fisherInfo_all_sessions_syntheticData');

    save(saveName,'results_cross_sizeControl');
end