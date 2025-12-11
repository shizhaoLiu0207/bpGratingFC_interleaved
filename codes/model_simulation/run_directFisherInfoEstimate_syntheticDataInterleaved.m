clear all
clc
close all
%%
doRun = 1;
doOrganize = 1;
doCombineAcrossCoherence = 1;
combineOptions = struct();




subjectCode = 'Model';


dataFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_ec_x256_g64_nr1024_contrast_all/';
saveFolder = '../../results/neural/fisherInfo_direct/fisherInfo_direct_model_interleaved_subsetNeuron';
saveFolder_session = fullfile(saveFolder, 'individual_sessions');
saveFolder_session_cross = fullfile(saveFolder, 'individual_sessions_cross');
mkdir(saveFolder_session);
mkdir(saveFolder_session_cross);


[xx,yy] = meshgrid([1:4], [1:8]);
session_id_list = xx * 10 + yy;
session_id_list = session_id_list(:);
nSession = numel(session_id_list);
nTimeBin = 9; 
if doRun
    for n = 1:nSession
        
        dat_fisher    = [];
        dat_fisher_cross = [];
        sessionStr = sprintf('Model_session%d',session_id_list(n));
        dataName = sprintf('syntheticData_interleaved_toUse_session%d',session_id_list(n));
        saveName = fullfile(saveFolder_session,sprintf('%s_subsetNeuron',sessionStr));
        saveName_cross = fullfile(saveFolder_session_cross,sprintf('%s_subsetNeuron',sessionStr));
       
        load(fullfile(dataFolder, dataName));
        timeBin_index               = cell2mat({synthData_use(:).timebinIndex}); 
        for t = 1:nTimeBin
           
            
            
            
            data_run = get_synthetic_data_cross(synthData_use, t-1);

            info_run                    = struct();
            info_run.sessionStr         = sessionStr;
            info_run.sessionType        = 'mainTask';

            idx_data                    = find(cell2mat({synthData_use(:).timebinIndex}) == t-1);
            info_run.timeWin            = synthData_use(idx_data(1)).timebin;
            info_run.i_win              = t;
            
            fprintf('Running session %d timewin %d\n',n,t);
            dat_fisher    = run_fisher_estimate_one_session(dat_fisher,data_run,info_run);
            dat_fisher_cross = run_fisher_cross_estimate_one_session(dat_fisher_cross, data_run, info_run);
           
            
        end
        save(saveName,'dat_fisher');
        save(saveName_cross, 'dat_fisher_cross');
    end
end

%% organize decoding results of all individual sessions into one big struct

filelist = dir(fullfile(saveFolder_session,'*.mat'));
dat_fisher = [];
for n = 1:numel(filelist)
    dat_session = load(fullfile(saveFolder_session,filelist(n).name));
    dat_fisher = [dat_fisher,dat_session.dat_fisher];
end
saveName = fullfile(saveFolder, 'fisherInfo_direct_all_sessions_syntheticData');
save(saveName, 'dat_fisher');
%% Organize cross fisher
filelist = dir(fullfile(saveFolder_session_cross,'*.mat'));
dat_fisher_cross = [];
for n = 1:numel(filelist)
    dat_session = load(fullfile(saveFolder_session_cross, filelist(n).name));
    dat_fisher_cross = [dat_fisher_cross,dat_session.dat_fisher_cross];
end
saveName = fullfile(saveFolder, 'fisherInfo_cross_direct_all_sessions_syntheticData');
save(saveName, 'dat_fisher_cross');
%% Combine across coherence

dataName = fullfile(saveFolder,  'fisherInfo_direct_all_sessions_syntheticData');
load(dataName);
saveName = fullfile(saveFolder,'CoherenceCombine_fisherInfo_all_sessions_syntheticData');

[dat_fisher_combine,combineOptions]    = combine_coherence_fisherInfo(dat_fisher,combineOptions);
save(saveName, 'dat_fisher_combine','combineOptions');



