clear all
clc
close all
%%%%% This script tests model simulation for interleaved task after
%%%%% implementation of oblique task image. But the problem of imbalanced orientation
%%%%% signal had not been fixed by when this was run

%%% [+c, 0], cardinal: 0 degree
%%% [0, +c], cardinal: 90 degree
%%% [+c, 0], oblique: 45 degree
%%% [0, +c], oblique: 135 degree

%% run simulation for interleaved （non-uniform phi_x was not implemented when this was ran）
doRun = 1;
if doRun
    saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample';
    stimulus_contrast_list = {[15,0],[12,0],[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9],[0,12],[0,15]};
    prior_task_list = {[1, 0], [0.8, 0.2], [0.5, 0.5], [0.2, 0.8], [0, 1]};
    image_task_list = {'cardinal','oblique'};
    nStim       = numel(stimulus_contrast_list);
    nPrior   = numel(prior_task_list);
    nTask       = numel(image_task_list); 
    for i = 1:nTask
        for j = 1:nPrior
            for k = 1:nStim
                fprintf('Simulating task %d/%d strategy %d/%d, stimulus %d/%d \n',i,nTask,j,nPrior,k,nStim)
                
                image_task = image_task_list{i};
                prior_task = prior_task_list{j};
                stimulus_contrast = stimulus_contrast_list{k};
                

                save_name = fullfile(saveFolder, sprintf('synthetic_data_task%d_strategy_%d_stim_%d',i,j,k));
               
                P = S_Exp_Para('test-interleaved', 'I.stimulus_contrast',stimulus_contrast,...
                                'G.prior_task',prior_task,'I.image_task',image_task);
                dat = S_Experiment(P);
                save(save_name, 'dat');
            end
        end
    end
end
%% plot psychometric curves of interleaved simulation (should be uniform phi_x)
doOrganize = 0;
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample';
if doOrganize 
  
    t = 1;
    for i = 1:nTask
        for j = 1:nPrior
            for k = 1:nStim
                 save_name = fullfile(saveFolder, sprintf('synthetic_data_task%d_strategy_%d_stim_%d',i,j,k));
                 load(save_name);
                 decision = (dat.O(:,3,end)>0.5) + 1; 
    
                 behav_results(t).image_task = dat.Projection.image_task;

                 prior_task =  dat.Projection.prior_task;
                 behav_results(t).prior_cardinal = prior_task(1);

                 stimulus_contrast = dat.Projection.stimulus_contrast;
                 behav_results(t).stimulus_contrast = stimulus_contrast(2) - stimulus_contrast(1); 

                 behav_results(t).decision = decision;
                 pChoice_ori2 = sum(decision == 2) / numel(decision);
                 behav_results(t).pChoice_ori2 = pChoice_ori2;
                 behav_results(t).semChoice_ori2 = sqrt(pChoice_ori2 .* (1 - pChoice_ori2) /size(dat.O,1));

                 t = t+1;
            end
        end
    end
    semChoice_ori2 = sqrt(pChoice_ori2 .* (1 - pChoice_ori2) /size(dat.O,1));
    save(fullfile(saveFolder, 'choice_all_sessions'), 'behav_results');
else
    load(fullfile(saveFolder, 'choice_all_sessions'));
end



prior_cardinal = cell2mat({behav_results(:).prior_cardinal});
prior_oblique  = 1 - prior_cardinal;  

prior_cardinal_list = unique(prior_cardinal);
prior_oblique_list  = unique(prior_oblique);
nPrior = numel(prior_oblique_list);

figure;
for i = 1:nPrior
    subplot(1,nPrior,i)
    idx = prior_cardinal == prior_cardinal_list(i) & strcmp({behav_results(:).image_task} , 'cardinal');
    x = cell2mat({behav_results(idx).stimulus_contrast});
    pChoice_ori2 = cell2mat({behav_results(idx).pChoice_ori2});
    semChoice_ori2 = cell2mat({behav_results(idx).semChoice_ori2});
    
    errorbar(x, pChoice_ori2, semChoice_ori2,'LineWidth',2,'color','red'); hold on
    
    
    idx = prior_oblique == prior_oblique_list(i) & strcmp({behav_results(:).image_task} , 'oblique');
    x = cell2mat({behav_results(idx).stimulus_contrast});
    pChoice_ori2 = cell2mat({behav_results(idx).pChoice_ori2});
    semChoice_ori2 = cell2mat({behav_results(idx).semChoice_ori2});
    
    errorbar(x, pChoice_ori2, semChoice_ori2,'LineWidth',2,'color','blue'); 

    box off
    set(gca,'fontsize',18);
    xlabel('contrast');ylabel('Prob. Choice 2');
    legend('Cardinal','Oblique')
    title(sprintf('Prior for correct task = %.1f',prior_oblique_list(i)))
end
