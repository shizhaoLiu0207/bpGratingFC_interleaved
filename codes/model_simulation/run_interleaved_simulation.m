clear all
clc
close all
%%% This function runs large scale simulation of interleaved data with
%%% different b_PF
%% 
%%% on macbook
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/run_interleaved';

% on linux
if ~exist(saveFolder)
    saveFolder = fullfile('/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/run_interleaved');
end
%%%% 1. March 15. First batch of simulation
% stimulus_contrast_list    = {[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12]};
% prior_list                = [1,0.75,0.5];
% b_PF_list                 = [-1.5, -0.8, 0, 0.8, 1.5];
% delta_list                = [0.05, 0.08];
% dimension_X              = 256;
% dimension_G              = 64;
% number_repetitions        = 512;

%%%%% 2. March 25. Second batch of simulation with more dense sample of prior
%%%%% and delta. More realistic list of b_PF. With smaller number of trials and neurons
% stimulus_contrast_list    = {[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12]};
% prior_list                = [0.5:0.1:1]; % prior for the correct task
% b_PF_list                 = [0.8, 0];
% delta_list                = [0.04:0.01:0.08];
% dimension_X              = 128;
% dimension_G              = 32;
% number_repetitions        = 256;

%%%% 3. generate a no-learning session
stimulus_contrast_list  = {[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12]};
prior_list              = [1]; % prior for the correct task
b_PF_list               = [0];
delta_list              = [0];
dimension_X              = 256;
dimension_G              = 64;
number_repetitions        = 512;
%%
nPrior      = numel(prior_list);
nStim       = numel(stimulus_contrast_list);
nPF         = numel(b_PF_list);
nDelta      = numel(delta_list);

for i = 1:nPF
    b_PF = b_PF_list(i);
    for j = 1:nPrior
        prior                       = prior_list(j);
  
       
        for n = 1:nStim
            stimulus_contrast = stimulus_contrast_list{n};
            for m = 1:nDelta
                delta = delta_list(m);
                prior_str                   = sprintf('%.2f',prior);
                b_PF_str                    = sprintf('%.2f', b_PF);
                delta_str                   = sprintf('%.2f', delta);
            
                prior_str                   = strrep(prior_str, '.', '_');
                b_PF_str                    = strrep(b_PF_str, '.', '_');
                delta_str                   = strrep(delta_str, '.', '_');

                save_name_cardinal = sprintf('synthetic_data_cardinal_bPF_%s_taskprior_%s_delta_%s_contrast_%d_%d.mat', ...
                                        b_PF_str, prior_str, delta_str, stimulus_contrast(1),  stimulus_contrast(2));

                save_name_oblique = sprintf('synthetic_data_oblique_bPF_%s_taskprior_%s_delta_%s_contrast_%d_%d.mat', ...
                                        b_PF_str, prior_str, delta_str, stimulus_contrast(1),  stimulus_contrast(2));

                if ~isfile(fullfile(saveFolder, save_name_cardinal))
                    prior_task = [prior, 1-prior];
                    image_task = 'cardinal';
                    P = S_Exp_Para('run-interleaved-nonuniform', 'I.stimulus_contrast',stimulus_contrast,...
                                                        'G.prior_task',prior_task,...
                                                        'I.image_task',image_task,...
                                                        'G.b_PF',b_PF,...
                                                        'G.delta',delta,...
                                                        'G.dimension_X', dimension_X,...
                                                        'G.dimension_G', dimension_G,...
                                                        'S.number_repetitions', number_repetitions);
             
  
                    dat = S_Experiment(P);
                    save(fullfile(saveFolder, save_name_cardinal),'dat');
                end
                if ~isfile(fullfile(saveFolder, save_name_oblique))
                    prior_task = [1-prior, prior];
                    image_task = 'oblique';
                    P = S_Exp_Para('run-interleaved-nonuniform', 'I.stimulus_contrast',stimulus_contrast,...
                                                        'G.prior_task',prior_task,...
                                                        'I.image_task',image_task,...
                                                        'G.b_PF',b_PF,...
                                                        'G.delta',delta,...
                                                        'G.dimension_X', dimension_X,...
                                                        'G.dimension_G', dimension_G,...
                                                        'S.number_repetitions', number_repetitions);
                    %%%%% By Shizhao 07/31/2025; Replace these parameters
                    %%%%% by user set-up at the beginning of this file
                    P.G.dimension_X = dimension_X;
                    P.G.dimension_G = dimension_G;
                    P.S.number_repetitions = number_repetitions;

                    dat = S_Experiment(P);
                    save(fullfile(saveFolder, save_name_oblique),'dat');
                end
            end
        end
        
    end
end