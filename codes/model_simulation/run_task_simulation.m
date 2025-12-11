clear all
clc
close all

savepath = '/Users/liushizhao/projectData_local/probinf_synthetic/syntheticData_interleaved/run_dual_task';
if ~isfolder(savepath)
    mkdir(savepath)
end
stimulus_contrast_list =  {[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12]};


prior_task_dual_list{1} = [0.8,0.5];
image_task_list{1}      = 'cardinal';

prior_task_dual_list{2} = [0.2,0.5];
image_task_list{2}      = 'oblique';

%%%% Generate dual task data
nDual = numel(prior_task_dual_list);
for n  = 1:nDual
    prior_task_cardinal = prior_task_dual_list{n}(1);
    prior_task_oblique  = prior_task_dual_list{n}(2);
    image_task          = image_task_list{n};
    prior_cardinal_str  = strrep(sprintf('%.1f',prior_task_cardinal), '.', '_');
    prior_oblique_str   = strrep(sprintf('%.1f',prior_task_oblique), '.', '_');
    
    
    for i = 1:numel(stimulus_contrast_list)
        stimulus_contrast = stimulus_contrast_list{i};
        contrast_str = [sprintf('%d',stimulus_contrast(1)),'_',sprintf('%d',stimulus_contrast(2))];

        savename = fullfile(savepath,sprintf('Dual_image_%s_prior_cardinal_%s_prior_oblique_%s_contrast_%s',...
            image_task, prior_cardinal_str, prior_oblique_str, contrast_str));
        
        P   = S_Exp_Para('run-dualTask', 'I.stimulus_contrast',stimulus_contrast,'I.image_task',image_task,...
            'G.prior_task_cardinal', prior_task_cardinal, 'G.prior_task_oblique', prior_task_oblique);

        
        dat      = S_Experiment(P);
        save(savename, 'dat')
    end
    
    clear dat_out
end
