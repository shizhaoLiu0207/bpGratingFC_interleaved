clear all
clc
close all
%%%% This script was to test how neuron's response and psychometric curves
%%%% change with b_PF: a scale factor that controls the shape of von mises
%%%% function - the proportion of neurons against their orientation
%%%% perference.
%%

doThis = 1;
if doThis
    figfolder = '../../figures/figures_informal/Bayesian_model_simulation';
    b_PF_list = [-1.5, 0, 0.8, 1.5];

    for ib = 1:numel(b_PF_list)
        b_PF = b_PF_list(ib);
        nTrial = 1000;
        nBootstrap = 1000;
        nNeuron = 64;
        stim_signal = 12;
        
        
        P = S_Exp_Para('test-interleaved-nonuniform','G.b_PF',b_PF,'G.dimension_X',nNeuron);
        projective_fields = C_Projection(P.G.fct, P.G.nx, P.G.dimension_X, P.G.dimension_G, P.G.number_locations, P.G.b_PF);
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
        set(gcf,'Units','normalized','Position',[0,0,0.5,0.5]);
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
        legend('cardinal task','oblique task','location','northwest')

        subplot(2,2,4);  hold on
        cdfplot(mean(dprime_cardinal,1));
        cdfplot(mean(dprime_oblique,1));
        set(gca,'fontsize',18); box off
        xlabel('dprime'); ylabel('CDF. nNeuron')
        legend('cardinal task','oblique task','location','northwest')
        sgtitle(sprintf('b_{PF} = %.1f',b_PF), 'fontsize', 20)
        figname = fullfile(figfolder, sprintf('neuronResponse_nonuniform_bPF_%.1f.png',b_PF));
        print(gcf, figname ,'-dpng')
    end

    
end
%%  test how downscale_oblique affects model's psychometric curves
doThis = 0;
if doThis
    
    %%% on macbook
    saveFolder = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample_nonuniform';

    % on linux
    if ~exist(saveFolder)
        saveFolder = fullfile('/home/shizhao/Documents/projectData/probinf_data/syntheticData_interleaved/test_smallsample_nonuniform');
    end

    stimulus_contrast_list  = {[9,0],[6,0],[3,0],[0,0],[0,3],[0,6],[0,9]};
    prior_list              = [1, 0.75, 0.5]; % prior for the correct task
    b_PF_list               = [0.8, 1.5];
    
    nPrior      = numel(prior_list);
    nStim       = numel(stimulus_contrast_list);
    nPF         = numel(b_PF_list);
    
    for i = 1:nPF
        b_PF = b_PF_list(i);
        for j = 1:nPrior
            prior                       = prior_list(j);
            [dat_cardinal, dat_oblique] = deal(cell(nStim, 1));
    
            prior_str                   = sprintf('%.2f',prior);
            b_PF_str                    = sprintf('%.2f', b_PF);
            prior_str                   = strrep(prior_str, '.', '_');
            b_PF_str                    = strrep(b_PF_str, '.', '_');
            save_name = sprintf(sprintf('synthetic_data_bPF_%s_taskprior_%s.mat', b_PF_str, prior_str));
            for n = 1:nStim
                stimulus_contrast = stimulus_contrast_list{n};
               
                prior_task = [prior, 1-prior];
                image_task = 'cardinal';
                P = S_Exp_Para('test-interleaved-nonuniform','G.b_PF',b_PF, 'I.stimulus_contrast',stimulus_contrast,...
                                            'G.prior_task',prior_task,'I.image_task',image_task);
                dat_cardinal{n} = S_Experiment(P);
            
                prior_task = [1-prior, prior];
                image_task = 'oblique';
                P = S_Exp_Para('test-interleaved-nonuniform','G.b_PF',b_PF, 'I.stimulus_contrast',stimulus_contrast,...
                                            'G.prior_task',prior_task,'I.image_task',image_task);
                dat_oblique{n} = S_Experiment(P);
            end
            save(fullfile(saveFolder, save_name),'dat_cardinal','dat_oblique','stimulus_contrast_list');
        end
    end
end
%% plot downscale_oblique and model performance
saveFolder          = '/Users/liushizhao/projectData_local/probinf_data/syntheticData_interleaved/test_smallsample_nonuniform';
figureFolder        = '../../figures/figures_informal/Bayesian_model_simulation';
filename_list       = dir(fullfile(saveFolder,'*.mat'));
nFile               = numel(filename_list);
[b_PF, taskprior]   = deal(zeros(nFile, 1));
for n = 1:nFile
% Extract numbers using regular expressions
    tokens = regexp(filename_list(n).name, 'synthetic_data_bPF_([\d_]+)_taskprior_([\d_]+)', 'tokens');
    % Convert extracted strings back to numbers
    extracted_params            = tokens{1}; % Extract matched tokens
    b_PF_str                    = strrep(extracted_params{1}, '_', '.'); % Replace _ with .
    taskprior_str               = strrep(extracted_params{2}, '_', '.');
    
    b_PF(n)                     = str2double(b_PF_str);
    taskprior(n)                = str2double(taskprior_str);
end
task_prior_list         = unique(taskprior);
b_PF_list               = unique(b_PF);
nPrior                  = numel(task_prior_list);
nPF                     = numel(b_PF_list); 
%%%% psychometric curves of each b_PF
%%%% one figure for one task prior
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:nPrior

    for j = 1:nPF
        idx = taskprior == task_prior_list(i) & b_PF == b_PF_list(j);
        load(fullfile(saveFolder, filename_list(idx).name));
        nStim = numel(stimulus_contrast_list);
        [pChoice_ori2_cardinal, semChoice_ori2_cardinal,...
            pChoice_ori2_oblique, semChoice_ori2_oblique ,x] = deal(zeros(nStim,1));
        
        for n = 1:nStim
            decision_cardinal = (dat_cardinal{n}.O(:,3,end)>0.5) + 1;
            pChoice_ori2_cardinal(n) = sum(decision_cardinal == 2) / numel(decision_cardinal);

            semChoice_ori2_cardinal(n) = ...
                sqrt(pChoice_ori2_cardinal(n) .* (1 - pChoice_ori2_cardinal(n)) /size(dat_cardinal{n}.O,1));

            decision_oblique = (dat_oblique{n}.O(:,3,end)>0.5) + 1;
            pChoice_ori2_oblique(n) = sum(decision_oblique == 2) / numel(decision_oblique);

            semChoice_ori2_oblique(n) = ...
                sqrt(pChoice_ori2_oblique(n) .* (1 - pChoice_ori2_oblique(n)) /size(dat_cardinal{n}.O,1));
            stimulus_contrast = dat_oblique{n}.Projection.stimulus_contrast;
            x(n) = stimulus_contrast(2) - stimulus_contrast(1);
        end
        subplot(nPrior,nPF,(i-1) * nPF + j)
        errorbar(x, pChoice_ori2_cardinal, semChoice_ori2_cardinal,'LineWidth',2,'color','red'); hold on
        errorbar(x, pChoice_ori2_oblique, semChoice_ori2_oblique,'LineWidth',2,'color','blue'); hold on
    
        box off
        set(gca,'fontsize',18);
        xlabel('contrast');ylabel('Prob. Choice 2');
        legend('Cardinal','Oblique','Location','southeast')
        title(sprintf('Task prior = %.2f, b PF = %.1f',task_prior_list(i), b_PF_list(j)));
    end
    % sgtitle(sprintf('Task prior = %.2f',task_prior_list(i)), 'fontsize', 20);
    % prior_str                   = sprintf('%.2f',task_prior_list(i)); 
    % prior_str                   = strrep(prior_str, '.', '_');
   
end
print(fullfile(figureFolder,'Psychometric_curve_test_nonuniform.png'),'-dpng');
%%
%%%%% plot choice accuracy as a function of b_PF for each
%%%%% contrast level

[cardinal_accuracy, oblique_accuracy, sem_cardinal_accuracy, sem_oblique_accuracy] = deal(zeros(nPrior,nPF, nStim));
for i = 1:nPrior
    for j = 1:nPF

        idx_file = taskprior == task_prior_list(i) & b_PF == b_PF_list(j);
        load(fullfile(saveFolder, filename_list(idx_file).name));
        contrast_signed         = cellfun(@(x)x(2) - x(1),stimulus_contrast_list); 
        contrast_signed_list    = unique(abs(contrast_signed));
        nContrast               = numel(abs(contrast_signed_list));

        for n = 1:nContrast
            idx = find(abs(contrast_signed) == contrast_signed_list(n));
            decision_cardinal = arrayfun(@(k) (dat_cardinal{k}.O(:,3,end)>0.5) + 1, idx, 'UniformOutput', false);
            decision_oblique  = arrayfun(@(k) (dat_oblique{k}.O(:,3,end)>0.5) + 1, idx, 'UniformOutput', false);

            decision_cardinal = cat(2, decision_cardinal{:});
            decision_oblique  = cat(2, decision_oblique{:});

            if contrast_signed_list(n) == 0
                stim_truth = [1,1];
            else
                stim_truth       = 1.5 + sign(contrast_signed(idx))/2;
            end

            cardinal_correct = decision_cardinal == repmat(stim_truth,size(decision_cardinal,1),1);
            oblique_correct  = decision_oblique == repmat(stim_truth,size(decision_cardinal,1),1);

            cardinal_accuracy(i,j,n) = sum(cardinal_correct(:)) / numel(cardinal_correct(:));
            oblique_accuracy(i,j,n)  = sum(oblique_correct(:)) / numel(oblique_correct(:));

            sem_cardinal_accuracy(i,j,n) =  sqrt(cardinal_accuracy(i,j,n) .* (1 - cardinal_accuracy(i,j,n)) / numel(cardinal_correct(:)));
            sem_oblique_accuracy(i,j,n) =  sqrt(oblique_accuracy(i,j,n) .* (1 - oblique_accuracy(i,j,n)) / numel(oblique_correct(:)));
        end
    end
end

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:nPrior

    for n = 1:nContrast
        subplot(nPrior,nContrast,(i-1)*nContrast + n)
        errorbar(b_PF_list,cardinal_accuracy(i,:,n), sem_cardinal_accuracy(i,:,n),'LineWidth',2,'color','red'); hold on
        errorbar(b_PF_list,oblique_accuracy(i,:,n), sem_oblique_accuracy(i,:,n),'LineWidth',2,'color','blue'); hold on
        box off
        set(gca,'fontsize',18);
        xlabel('b PF');ylabel('Accuracy');
        legend('Cardinal','Oblique','Location','southeast')
        title(sprintf('Task prior = %.2f, Contrast = %d',task_prior_list(i), contrast_signed_list(n)));
       % ylim([0,1])
    end
end
