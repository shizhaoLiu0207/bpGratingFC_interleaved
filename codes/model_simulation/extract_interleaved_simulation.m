clear all
clc
close all
%% extract synthetic data
%%%% list of parameters
%%% on macbook
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/run_interleaved';

% on linux
if ~exist(saveFolder)
    saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/run_interleaved';
end


filename_list  = dir(fullfile(saveFolder,'*.mat'));
nFile = numel(filename_list);
[b_PF, taskprior, delta, contrast_signed] = deal(zeros(nFile, 1));
imagetask  = cell(nFile,1);
for n = 1:nFile
% Extract numbers using regular expressions
    tokens = regexp(filename_list(n).name, 'synthetic_data_([a-zA-Z0-9]+)_bPF_([-]?[\d_]+)_taskprior_([\d_]+)_delta_([\d_]+)_contrast_([\d_]+)_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params        = tokens{1}; % Extract matched tokens
    imagetask_str           = extracted_params{1};
    bPF_str                 = strrep(extracted_params{2}, '_', '.'); % Replace _ with .
    taskprior_str           = strrep(extracted_params{3}, '_', '.');
    delta_str               = strrep(extracted_params{4}, '_', '.');
    contrast_1_str          = extracted_params{5};
    contrast_2_str          = extracted_params{6};
    
    imagetask{n}               = imagetask_str;
    b_PF(n)                    = str2double(bPF_str);
    taskprior(n)               = str2double(taskprior_str);
    delta(n)                   = str2double(delta_str);
    contrast_signed(n)         = str2double(contrast_2_str) - str2double(contrast_1_str);
end
imagetask_list          = unique(imagetask); 
bPF_list                = unique(b_PF);
taskprior_list          = unique(taskprior);
delta_list              = unique(delta); 

contrast_signed_list    = unique(contrast_signed); 

nTask                   = numel(imagetask_list); 
nPrior                  = numel(taskprior_list);
nPF                     = numel(bPF_list); 
nDelta                  = numel(delta_list); 
%%
%%%% For each combination of image_task, downscaleoblique, taskprior,
%%%% delta, find all the files (contrasts) of this parameter of
%%%% combination, organize into one data file
organizeData_saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/synthData_use_interleaved';

% on linux
if ~exist(saveFolder)
    organizeData_saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved';

end

%organizeData_saveFolder = '/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/synthData_use_interleaved';
timeBin_list{1} = [4:99];
% T_total = numel(timeBin_list{1});
% nTimebin = 8;
% for n = 1:nTimebin
%     binSize = T_total / nTimebin;
%     timeBin_list{n+1} = [4 + (n-1) * binSize: n*binSize + 3];
% end

for i = 1:nTask
    for j = 1:nPrior
        for n = 1:nPF
            for m = 1:nDelta
                %%%%% find all sessions of this parameter combination
                idx = find(strcmp(imagetask, imagetask_list{i}) & ...
                      taskprior == taskprior_list(j) & ...
                      b_PF == bPF_list(n) & ...
                      delta == delta_list(m));

                  bPF_str               = strrep(sprintf('%.2f', bPF_list(n)), '.', '_');
                  delta_str             = strrep(sprintf('%.2f', delta_list(m)), '.', '_');
                  task_prior            = strrep(sprintf('%.2f', taskprior_list(j)), '.', '_');
                  session_name          = sprintf('synthData_use_imagetask_%s_bPF_%s_delta_%s_taskprior_%s',...
                                            imagetask_list{i}, bPF_str, delta_str, task_prior);
                  save_name = fullfile(organizeData_saveFolder, session_name);
                  %%%%% load and extract useful information for each
                  %%%%% session
                  synthData_use = struct();
                  d = 1;
                  for k = 1:numel(idx)
                      for t = 1:numel(timeBin_list)
                        load(fullfile(saveFolder, filename_list(idx(k)).name));
                        synthData_use(d).sessionStr             = session_name;  
                        synthData_use(d).image_task             = dat.Projection.image_task;
                        synthData_use(d).b_PF                   = dat.Projection.b_PF;
                        synthData_use(d).learning_strength      = dat.Projection.delta;
                       
                        synthData_use(d).prior_task             = dat.Projection.prior_task;
                        synthData_use(d).cardinal_prior         = synthData_use(d).prior_task(1);
                        synthData_use(d).oblique_prior          = synthData_use(d).prior_task(2);
                        synthData_use(d).timebin                = timeBin_list{t};
                        synthData_use(d).timebinIndex           = t-1;
                        synthData_use(d).X_response             = sum(dat.X(:,:,timeBin_list{t}),3);
                        
                        synthData_use(d).decision               = (dat.O(:,3,end) > 0.5) + 1;
  
                        % stimulus contrast
                        synthData_use(d).stimulus_contrast      = dat.Projection.stimulus_contrast;
                        synthData_use(d).contrast_signed        = synthData_use(d).stimulus_contrast(2)- synthData_use(d).stimulus_contrast(1);
                        % neurons's preferred orientation
                        synthData_use(d).phi_x                  = dat.Projection.phi_x;
                        d = d+1;
                      end
                  end
                  save(save_name, 'synthData_use')
                 
                
            end
        end
    end
end
              

