clear all
clc
close all
%%%%% This scripts generate real synthetic interleaved sessions by
%%%%% combining data of single task simulated with various parameter set
%%
figFolder = '../../figures/figures_informal/Bayesian_model_simulation';
%%% on macbook
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved';

% on linux
if ~exist(saveFolder)
    saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved';
end

filename_list  = dir(fullfile(saveFolder,'*.mat'));
nFile = numel(filename_list);
[bPF, taskprior, delta, contrast_signed] = deal(zeros(nFile, 1));
imagetask  = cell(nFile,1);
for n = 1:nFile
% Extract numbers using regular expressions
    tokens = regexp(filename_list(n).name, 'synthData_use_imagetask_([a-zA-Z0-9]+)_bPF_([-]?[\d_]+)_delta_([\d_]+)_taskprior_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params        = tokens{1}; % Extract matched tokens
    imagetask_str           = extracted_params{1};
    bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
    taskprior_str           = strrep(extracted_params{3}, '_', '.');
    delta_str               = strrep(extracted_params{4}, '_', '.');

    
    imagetask{n}               = imagetask_str;
    bPF(n)        = str2double(bPF_str);
    taskprior(n)               = str2double(taskprior_str);
    delta(n)                   = str2double(delta_str);
end
imagetask_list         = unique(imagetask); 
bPF_list   = unique(bPF);
taskprior_list          = unique(taskprior);
delta_list              = unique(delta);   


nPrior                  = numel(taskprior_list);
nDownscale              = numel(bPF_list); 
nDelta                  = numel(delta_list); 
%%
% %%%% for each bPF, get combination of each delta_task prior
% of each task
for i = 1:numel(bPF_list)


    idx_cardinal = find(bPF == bPF_list(i) & strcmp(imagetask, 'cardinal'));
    idx_oblique =  find(bPF == bPF_list(i) & strcmp(imagetask, 'oblique'));

    
    for n1 = 1:numel(idx_cardinal)
        for n2 = 1:numel(idx_oblique)
            dat_cardinal = load(fullfile(saveFolder, filename_list(idx_cardinal(n1)).name));
            dat_oblique  = load(fullfile(saveFolder, filename_list(idx_oblique(n2)).name)); 

            synthData_interleaved = [dat_cardinal.synthData_use,dat_oblique.synthData_use];
            
            bPF_str                         = strrep(sprintf('%.2f', bPF_list(i)), '.', '_');
            delta_cardinal_str              = strrep(sprintf('%.2f', dat_cardinal.synthData_use(1).learning_strength), '.', '_');
            delta_oblique_str               = strrep(sprintf('%.2f', dat_oblique.synthData_use(1).learning_strength), '.', '_');
            taskprior_cardinal_str          = strrep(sprintf('%.2f', dat_cardinal.synthData_use(1).prior_task(1)), '.', '_');
            taskprior_oblique_str           = strrep(sprintf('%.2f', dat_oblique.synthData_use(1).prior_task(2)), '.', '_');
            %%%% to make it easier for future use, add a scaler version of
%             %%%% prior task and stimulus contrast (should have done this on the previous step)
%             for k = 1:numel(synthData_interleaved)
%                 synthData_interleaved(k).contrast_signed = ...
%                         synthData_interleaved(k).stimulus_contrast(2) - synthData_interleaved(k).stimulus_contrast(1);
%                 synthData_interleaved(k).prior_cardinal = synthData_interleaved(k).prior_task(1);
%                 synthData_interleaved(k).prior_oblique  = synthData_interleaved(k).prior_task(2);
%             end
            save_name = sprintf('synthData_use_interleaved_bPF_%s_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s',...
                bPF_str, delta_cardinal_str, taskprior_cardinal_str, delta_oblique_str, taskprior_oblique_str);
           
            sessionStr = sprintf('Model_bPF_%s_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s',...
                bPF_str, delta_cardinal_str, taskprior_cardinal_str, delta_oblique_str, taskprior_oblique_str);

            for k = 1:numel(synthData_interleaved)
                synthData_interleaved(k).sessionStr = sessionStr;
            end
            save(fullfile(saveFolder, 'real_interleaved', save_name),'synthData_interleaved');
        end
    end

end


