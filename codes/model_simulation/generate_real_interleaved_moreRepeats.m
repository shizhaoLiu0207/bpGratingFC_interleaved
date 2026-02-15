clear all
clc
close all
%%%%% This scripts generate real synthetic interleaved sessions by
%%%%% combining data of single task simulated with various parameter set
%%
% figFolder = '../../figures/figures_informal/Bayesian_model_simulation';
% %%% on macbook
% saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved';
% 
% % on linux
% if ~exist(saveFolder)
%     saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved';
% end
saveFolder = '/Volumes/T7/dataTransfer/BayesianModel_interleaved_simulation/synthData_use_interleaved_moreRepeats';
if ~isfolder(saveFolder)
    mkdir(saveFolder)
end

filename_list  = dir(fullfile(saveFolder,'synthData_use_*.mat'));
nFile = numel(filename_list);
[taskprior, delta, contrast_signed, repeats] = deal(zeros(nFile, 1));
imagetask  = cell(nFile,1);
for n = 1:nFile
% Extract numbers using regular expressions
    tokens = regexp(filename_list(n).name, 'synthData_use_imagetask_([a-zA-Z0-9]+)_bPF_0_80_delta_([\d_]+)_taskprior_([\d_]+)_rep_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params        = tokens{1}; % Extract matched tokens
    imagetask_str           = extracted_params{1};
    %bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
    delta_str               = strrep(extracted_params{2}, '_', '.');
    taskprior_str           = strrep(extracted_params{3}, '_', '.');
    repeat_str              = extracted_params{4};
    
    
    imagetask{n}               = imagetask_str;
    taskprior(n)               = str2double(taskprior_str);
    delta(n)                   = str2double(delta_str);
    repeats(n)                   = str2double(repeat_str); 
end
imagetask_list          = unique(imagetask); 
taskprior_list          = unique(taskprior);
delta_list              = unique(delta);   
repeat_list             = unique(repeats); 


nPrior                  = numel(taskprior_list);
nDelta                  = numel(delta_list); 
nRep                    = numel(repeat_list); 
%%
% %%%% for each bPF, get combination of each delta_task prior
% of each task

for i = 1:5


    idx_cardinal = find(repeats == repeat_list(i) & strcmp(imagetask, 'cardinal'));
    idx_oblique =  find(repeats == repeat_list(i) & strcmp(imagetask, 'oblique'));

    
    for n1 = 1:numel(idx_cardinal)
        for n2 = 1:numel(idx_oblique)
            dat_cardinal = load(fullfile(saveFolder, filename_list(idx_cardinal(n1)).name));
            dat_oblique  = load(fullfile(saveFolder, filename_list(idx_oblique(n2)).name)); 

            synthData_interleaved = [dat_cardinal.synthData_use,dat_oblique.synthData_use];
            
           % bPF_str                         = strrep(sprintf('%.2f', bPF_list(i)), '.', '_');
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
            save_name = sprintf('synthData_use_interleaved_bPF_0_80_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s_rep_%d',...
                 delta_cardinal_str, taskprior_cardinal_str, delta_oblique_str, taskprior_oblique_str, repeat_list(i));
           
            sessionStr = sprintf('Model_bPF_0_80_cardinal_delta_%s_prior_%s_oblique_delta_%s_prior_%s_rep_%d',...
                delta_cardinal_str, taskprior_cardinal_str, delta_oblique_str, taskprior_oblique_str, repeat_list(i));

            for k = 1:numel(synthData_interleaved)
                synthData_interleaved(k).sessionStr = sessionStr;
            end
            save(fullfile(saveFolder, 'real_interleaved', save_name),'synthData_interleaved');
        end
    end

end


