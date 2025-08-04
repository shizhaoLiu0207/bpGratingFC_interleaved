clear all
clc
close all
%%%%% Generate a single session of interleaved synthetic data by combining data 
%%%%% of single task simulated with certain parameter set
%%
%%% on macbook
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved';

% on linux
if ~exist(saveFolder)
    saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved';
end


b_PF            = 0;
cardinal_prior  = 1;
cardinal_delta  = 0;
oblique_prior   = 1;
oblique_delta   = 0;

bPF_str                     = strrep(sprintf('%.2f', b_PF), '.', '_');
cardinal_delta_str          = strrep(sprintf('%.2f', cardinal_delta), '.', '_');
oblique_delta_str           = strrep(sprintf('%.2f', oblique_delta), '.', '_');
cardinal_prior_str          = strrep(sprintf('%.2f', cardinal_prior), '.', '_');
oblique_prior_str           = strrep(sprintf('%.2f', oblique_prior), '.', '_');

cardinal_data_name = sprintf('synthData_use_imagetask_cardinal_bPF_%s_delta_%s_taskprior_%s.mat',...
                            bPF_str, cardinal_delta_str, cardinal_prior_str);
oblique_data_name  = sprintf('synthData_use_imagetask_oblique_bPF_%s_delta_%s_taskprior_%s.mat',...
                            bPF_str, oblique_delta_str, oblique_prior_str);

dat_cardinal = load(fullfile(saveFolder, cardinal_data_name));
dat_oblique  = load(fullfile(saveFolder, oblique_data_name)); 

synthData_interleaved = [dat_cardinal.synthData_use,dat_oblique.synthData_use];


%%%% to make it easier for future use, add a scaler version of
%             %%%% prior task and stimulus contrast (should have done this on the previous step)
%             for k = 1:numel(synthData_interleaved)
%                 synthData_interleaved(k).contrast_signed = ...
%                         synthData_interleaved(k).stimulus_contrast(2) - synthData_interleaved(k).stimulus_contrast(1);
%                 synthData_interleaved(k).prior_cardinal = synthData_interleaved(k).prior_task(1);
%                 synthData_interleaved(k).prior_oblique  = synthData_interleaved(k).prior_task(2);
%             end
save_name = sprintf('synthData_use_interleaved_bPF_%s_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s',...
    bPF_str, cardinal_delta_str, cardinal_prior_str, oblique_delta_str, oblique_prior_str);

sessionStr = sprintf('Model_bPF_%s_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s',...
    bPF_str, cardinal_delta_str, cardinal_prior_str, oblique_delta_str, oblique_prior_str);

for k = 1:numel(synthData_interleaved)
    synthData_interleaved(k).sessionStr = sessionStr;
end
save(fullfile(saveFolder, 'real_interleaved', save_name),'synthData_interleaved');
