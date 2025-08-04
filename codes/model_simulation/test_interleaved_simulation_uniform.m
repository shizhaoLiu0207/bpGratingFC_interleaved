clear all
clc
close all
%%%%
%%% This script was to confirm no bias exists in neuron's response and
%%% psychometric curves after the input image was fixed such that
%%% orientation signal in cardinal and oblique image is equal
%%

doThis = 1;
if doThis
    nTrial = 1000;
    nBootstrap = 1000;
    nNeuron = 64;
    stim_signal = 12;
    
    P = S_Exp_Para('test-interleaved','G.fct','nxN','G.dimension_X',nNeuron);
    projective_fields = C_Projection(P.G.fct, P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations);
    
    [response_cardinal_1,response_cardinal_2,response_oblique_1,response_oblique_2] = deal(zeros(nTrial,nNeuron));
    
    
    
    im_type     = P.I.fct;
    n_locs      = P.G.number_locations;
    im_height   = P.G.ny;
    nNeuron     = size(projective_fields.G,2);
    
    for t = 1:nTrial
        %%%%% cardinal image
        image_task  = 'cardinal';
        stimulus_contrast = [stim_signal,0];
        stim_cardinal_1 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
        image_task  = 'cardinal';
        stimulus_contrast = [0,stim_signal];
        stim_cardinal_2 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
    
        %%%%% cardinal signal
        response_cardinal_1(t,:) = arrayfun(@(n)stim_cardinal_1(:)' * projective_fields.G(:,n), [1:nNeuron]);
        response_cardinal_2(t,:) = arrayfun(@(n)stim_cardinal_2(:)' * projective_fields.G(:,n), [1:nNeuron]);
        %cardinal_signal(i,t) = mean(abs(response_cardinal_1 - response_cardinal_2));
    
    
        %%%%% oblique image
        image_task  = 'oblique';
        stimulus_contrast = [stim_signal,0];
        stim_oblique_1 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
        stimulus_contrast = [0,stim_signal];
        stim_oblique_2 = InputImage(im_type, n_locs, im_height, stimulus_contrast, image_task);
        %%%%% oblique signal
        response_oblique_1(t,:) = arrayfun(@(n)stim_oblique_1(:)' * projective_fields.G(:,n), [1:nNeuron]);
        response_oblique_2(t,:) = arrayfun(@(n)stim_oblique_2(:)' * projective_fields.G(:,n), [1:nNeuron]);
        %oblique_signal(i,t) = mean(abs(response_oblique_1 - response_oblique_2));
    
    
    
    
    end
    [dprime_cardinal, dprime_oblique] = deal(zeros(nBootstrap,nNeuron));
    for t = 1:nBootstrap
        idx = randsample(nTrial,nTrial/2);
        dprime_cardinal(t,:) = abs(mean(response_cardinal_1(idx,:),1) - mean(response_cardinal_2(idx,:),1)) ./ ...
            sqrt((var(response_cardinal_1(idx,:), [], 1) + var(response_cardinal_2(idx,:), [], 1)) / 2);
        dprime_oblique(t,:) = abs(mean(response_oblique_1(idx,:),1) - mean(response_oblique_2(idx,:),1)) ./ ...
            sqrt((var(response_oblique_1(idx,:), [], 1) + var(response_oblique_2(idx,:), [], 1)) / 2);
    end

     figure;
    subplot(2,2,1);hold on
    plot(projective_fields.phi_x,mean(response_cardinal_1,1),'-o');
    
    plot(projective_fields.phi_x,mean(response_cardinal_2,1),'-o');
    plot(projective_fields.phi_x,mean(response_oblique_1,1),'-o');
    plot(projective_fields.phi_x,mean(response_oblique_2,1),'-o');
    set(gca,'fontsize',18); box off
    xlabel('Phi'); ylabel('Response')
    legend('0','90','45','135');

    subplot(2,2,2);hold on
    plot(projective_fields.phi_x,mean(dprime_cardinal,1),'-o');
    plot(projective_fields.phi_x,mean(dprime_oblique,1),'-o');
    set(gca,'fontsize',18); box off
    xlabel('Phi'); ylabel('dprime')
    legend('cardinal task','oblique task')

    subplot(2,2,3); hold on
    cdfplot([mean(response_cardinal_1,1), mean(response_cardinal_2,1)]);
    cdfplot([mean(response_oblique_1,1), mean(response_oblique_2,1)]);
    set(gca,'fontsize',18); box off
    xlabel('Response'); ylabel('CDF. nNeuron')
    legend('cardinal task','oblique task')

    subplot(2,2,4);  hold on
    cdfplot(mean(dprime_cardinal,1));
    cdfplot(mean(dprime_oblique,1));
    set(gca,'fontsize',18); box off
    xlabel('dprime'); ylabel('CDF. nNeuron')
    legend('cardinal task','oblique task')
    sgtitle('Uniform phi_{x}', 'fontsize', 20)
   
end
%% plot projective field
figure
for n = 1:nNeuron
    subplot(8,8,n)
    imagesc(reshape(projective_fields.G(:,n),im_height,im_height))
end
%% run simulation for interleaved （non-uniform phi_x was not implemented when this was ran）
doRun = 0;
if doRun
    saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample_uniformity';
    mkdir(saveFolder);
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
               
                P = S_Exp_Para('test-interleaved', 'G.fct','nxN','I.stimulus_contrast',stimulus_contrast,...
                                'G.prior_task',prior_task,'I.image_task',image_task);
                dat = S_Experiment(P);
                save(save_name, 'dat');
            end
        end
    end
end

%%
doOrganize = 0;
saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample_uniformity';
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


% dprime_cardinal = squeeze(dprime_cardinal);
% dprime_oblique  = squeeze(dprime_oblique);
% 
% subplot(nRow, nCol, (i-1)* nCol + j)
% plot(downscale_oblique_list,  mean(dprime_cardinal,2), 'color','red','linewidth',2); hold on
% plot(downscale_oblique_list , mean(dprime_oblique,2), 'color','blue','linewidth',2)
% % errorbar(obl_list, mean(dprime_cardinal,2), std(dprime_cardinal,[],2), 'color','red','linewidth',2); hold on
% % errorbar(obl_list, mean(dprime_oblique,2),std(dprime_oblique,[],2),'color','blue','linewidth',2);
% xlabel('downscale oblique'); ylabel('Avg. task dprime')
% set(gca,'fontsize',18)
% legend('Cardinal','Oblique')
% title(sprintf('Signal level = %d, nNeuron = %d',stim_signal, nNeuron))
% box off