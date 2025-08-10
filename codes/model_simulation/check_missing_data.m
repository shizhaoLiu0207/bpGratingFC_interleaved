clear all
clc
close all
%%
data_folder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved/real_interleaved/batch_2';

prior_list                = [0.5:0.1:1]; 
b_PF_list                 = [0.8, 0];
delta_list                = [0.04:0.01:0.08];

data_list = dir(fullfile(data_folder,'*.mat'));
%% target combinations
k = 1;
for n1 = 1:numel(b_PF_list)
    for n2 = 1:numel(delta_list)
        for n3 = 1:numel(prior_list)
                for n4 = 1:numel(delta_list)
                    for n5 = 1:numel(prior_list)
                        
                        target.b_PF(k,1) = b_PF_list(n1);
                        target.cardinal_delta(k,1)  = delta_list(n2);
                        target.cardinal_prior(k,1)  = prior_list(n3);
                        target.oblique_delta(k,1)   = delta_list(n4);
                        target.oblique_prior(k,1)   = prior_list(n5);
                        k = k+1;
                    end
                end
        end
    end
end
%% combination of data
nFile = numel(data_list);
[file.b_PF, file.cardinal_delta, file.cardinal_prior, file.oblique_delta, file.oblique_prior] = deal(zeros(nFile,1));
for n = 1:nFile
    tokens = regexp(data_list(n).name, 'synthData_use_interleaved_bPF_([-]?[\d_]+)_cardinal_delta_([\d_]+)_prior_([\d_]+)_oblique_delta_([\d_]+)_prior_([\d_]+)', 'tokens');
    
    extracted_params        = tokens{1}; % Extract matched tokens
 
    bPF_str                 = strrep(extracted_params{1}, '_', '.'); % Replace _ with .
    cardinal_delta_str      = strrep(extracted_params{2}, '_', '.');
    cardinal_prior_str      = strrep(extracted_params{3}, '_', '.');
    oblique_delta_str       = strrep(extracted_params{4}, '_', '.');
    oblique_prior_str       = strrep(extracted_params{5}, '_', '.');

    

    file.b_PF(n)            = str2double(bPF_str);
    file.cardinal_delta(n)  = str2double(cardinal_delta_str);
    file.cardinal_prior(n)  = str2double(cardinal_prior_str);
    file.oblique_delta(n)   = str2double(oblique_delta_str);
    file.oblique_prior(n)   = str2double(oblique_prior_str);

end
%% find the combination that is missing

idx_missing = [];
for n = 1:numel(target.b_PF)
    idx = file.b_PF == target.b_PF(n) & ...
        file.cardinal_delta == target.cardinal_delta(n) & ...
        file.cardinal_prior == target.cardinal_prior(n) & ...
        file.oblique_delta  == target.oblique_delta(n) & ...
        file.oblique_prior  == target.oblique_prior(n);
    if sum(idx) == 0
        idx_missing = [idx_missing;n];
    end
end

target_params_table = struct2table(target);

target_params_table(idx_missing,:)

